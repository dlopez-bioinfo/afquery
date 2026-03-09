import duckdb
from pyroaring import BitMap

from .bitmaps import deserialize
from .constants import normalize_chrom


def annotate_vcf(
    engine,
    input_vcf: str,
    output_vcf: str,
    icd10_codes: list[str],
    sex_filter: str = "both",
) -> dict:
    """Annotate a VCF with AFQUERY_AC / AFQUERY_AN / AFQUERY_AF INFO fields.

    Returns a summary dict: {"n_variants": int, "n_annotated": int, "n_uncovered": int}.
    """
    import cyvcf2

    icd10_bitmap = engine._build_icd10_bitmap(icd10_codes)
    sex_bitmap = engine._build_sex_bitmap(sex_filter)

    vcf = cyvcf2.VCF(input_vcf)
    vcf.add_info_to_header({
        "ID": "AFQUERY_AC", "Number": "A", "Type": "Integer",
        "Description": "Allele count in eligible samples",
    })
    vcf.add_info_to_header({
        "ID": "AFQUERY_AN", "Number": "1", "Type": "Integer",
        "Description": "Allele number in eligible samples",
    })
    vcf.add_info_to_header({
        "ID": "AFQUERY_AF", "Number": "A", "Type": "Float",
        "Description": "Allele frequency in eligible samples",
    })
    writer = cyvcf2.Writer(output_vcf, vcf)

    stats = {"n_variants": 0, "n_annotated": 0, "n_uncovered": 0}
    buffer: list = []
    current_chrom: str | None = None

    for variant in vcf:
        norm_chrom = normalize_chrom(variant.CHROM)
        if norm_chrom != current_chrom:
            _flush_buffer(buffer, current_chrom, engine, icd10_bitmap, sex_bitmap, writer, stats)
            buffer = []
            current_chrom = norm_chrom
        buffer.append(variant)

    _flush_buffer(buffer, current_chrom, engine, icd10_bitmap, sex_bitmap, writer, stats)

    writer.close()
    vcf.close()
    return stats


def _flush_buffer(
    buffer: list,
    chrom: str | None,
    engine,
    icd10_bitmap: BitMap,
    sex_bitmap: BitMap,
    writer,
    stats: dict,
) -> None:
    if not buffer or chrom is None:
        return

    unique_positions = list({v.POS for v in buffer})

    # Compute eligible / AN per unique position
    pos_data: dict[int, tuple[BitMap, int]] = {}
    for pos in unique_positions:
        eligible, AN = engine._compute_eligible(chrom, pos, icd10_bitmap, sex_bitmap)
        pos_data[pos] = (eligible, AN)

    valid_positions = [p for p in unique_positions if pos_data[p][1] > 0]

    # Batch-fetch variant rows from Parquet
    variant_data: dict[tuple[int, str, str], tuple[bytes, bytes]] = {}
    parquet_path = engine._db / "variants" / f"{chrom}.parquet"
    if valid_positions and parquet_path.exists():
        con = duckdb.connect()
        placeholders = ", ".join("?" * len(valid_positions))
        rows = con.execute(
            f"SELECT pos, ref, alt, het_bitmap, hom_bitmap"
            f" FROM read_parquet(?) WHERE pos IN ({placeholders})",
            [str(parquet_path)] + valid_positions,
        ).fetchall()
        con.close()
        for pos, ref, alt, het_bytes, hom_bytes in rows:
            variant_data[(pos, ref, alt)] = (bytes(het_bytes), bytes(hom_bytes))

    for variant in buffer:
        pos = variant.POS
        eligible, AN = pos_data[pos]
        stats["n_variants"] += 1

        if AN == 0:
            stats["n_uncovered"] += 1
            variant.INFO["AFQUERY_AN"] = 0
            writer.write_record(variant)
            continue

        ac_list = []
        af_list = []
        any_found = False

        for alt in variant.ALT:
            key = (pos, variant.REF, alt)
            if key in variant_data:
                het_bytes, hom_bytes = variant_data[key]
                het_bm = deserialize(het_bytes)
                hom_bm = deserialize(hom_bytes)
                AC = len(het_bm & eligible) + 2 * len(hom_bm & eligible)
                any_found = True
            else:
                AC = 0  # covered position, absent variant
            ac_list.append(AC)
            af_list.append(AC / AN)

        # cyvcf2 Number=A fields must be set as comma-separated strings
        variant.INFO["AFQUERY_AC"] = ",".join(str(v) for v in ac_list)
        variant.INFO["AFQUERY_AN"] = AN
        variant.INFO["AFQUERY_AF"] = ",".join("%.10g" % v for v in af_list)
        if any_found:
            stats["n_annotated"] += 1
        writer.write_record(variant)
