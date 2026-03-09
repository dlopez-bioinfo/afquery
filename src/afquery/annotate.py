import os
import duckdb

from .bitmaps import deserialize
from .constants import normalize_chrom
from .models import SampleFilter


def _compute_chunk_annotations(
    db_path: str,
    chrom: str,
    bucket_id: int,
    records: list[tuple[int, str, list[str]]],
    sf: SampleFilter,
) -> dict[tuple[int, str, str], tuple[int, int, bool]]:
    """Compute (AC, AN, in_parquet) for each (pos, ref, alt) in a (chrom, bucket) chunk.

    No cyvcf2 dependency — safe to run in a subprocess worker.
    """
    from pathlib import Path
    from .query import QueryEngine

    engine = QueryEngine(db_path)
    sample_bm = engine._build_sample_bitmap(sf)

    unique_positions = list({pos for pos, _ref, _alts in records})

    pos_data: dict[int, tuple] = {}
    for pos in unique_positions:
        eligible, AN = engine._compute_eligible(chrom, pos, sample_bm)
        pos_data[pos] = (eligible, AN)

    valid_positions = [p for p in unique_positions if pos_data[p][1] > 0]

    variant_data: dict[tuple[int, str, str], tuple[bytes, bytes]] = {}
    _db = Path(db_path)
    bucket_start = bucket_id * 1_000_000
    bucket_end = (bucket_id + 1) * 1_000_000 - 1

    if chrom in engine._partitioned_chroms:
        parquet_file = _db / "variants" / chrom / f"bucket_{bucket_id}.parquet"
        if valid_positions and parquet_file.exists():
            con = duckdb.connect()
            placeholders = ", ".join("?" * len(valid_positions))
            rows = con.execute(
                f"SELECT pos, ref, alt, het_bitmap, hom_bitmap"
                f" FROM read_parquet('{parquet_file}') WHERE pos IN ({placeholders})",
                valid_positions,
            ).fetchall()
            con.close()
            for pos, ref, alt, het_bytes, hom_bytes in rows:
                variant_data[(pos, ref, alt)] = (bytes(het_bytes), bytes(hom_bytes))
    else:
        parquet_file = _db / "variants" / f"{chrom}.parquet"
        if valid_positions and parquet_file.exists():
            con = duckdb.connect()
            rows = con.execute(
                f"SELECT pos, ref, alt, het_bitmap, hom_bitmap"
                f" FROM read_parquet('{parquet_file}') WHERE pos BETWEEN ? AND ?",
                [bucket_start, bucket_end],
            ).fetchall()
            con.close()
            valid_pos_set = set(valid_positions)
            for pos, ref, alt, het_bytes, hom_bytes in rows:
                if pos in valid_pos_set:
                    variant_data[(pos, ref, alt)] = (bytes(het_bytes), bytes(hom_bytes))

    result: dict[tuple[int, str, str], tuple[int, int, bool]] = {}
    for pos, ref, alts in records:
        eligible, AN = pos_data[pos]
        for alt in alts:
            key = (pos, ref, alt)
            if key in result:
                continue  # dedup
            if AN == 0:
                result[key] = (0, 0, False)
            elif key in variant_data:
                het_bytes, hom_bytes = variant_data[key]
                het_bm = deserialize(het_bytes)
                hom_bm = deserialize(hom_bytes)
                AC = len(het_bm & eligible) + 2 * len(hom_bm & eligible)
                result[key] = (AC, AN, True)
            else:
                result[key] = (0, AN, False)

    return result


def annotate_vcf(
    engine,
    input_vcf: str,
    output_vcf: str,
    sf: SampleFilter,
    n_workers: int | None = None,
) -> dict:
    """Annotate a VCF with AFQUERY_AC / AFQUERY_AN / AFQUERY_AF INFO fields.

    Returns a summary dict: {"n_variants": int, "n_annotated": int, "n_uncovered": int}.
    """
    import cyvcf2
    from concurrent.futures import ProcessPoolExecutor, as_completed

    # Phase 1: collect all variants per (chrom, bucket_id)
    work_order: list[tuple[str, int]] = []
    variant_buffers: dict[tuple[str, int], list] = {}
    query_buffers: dict[tuple[str, int], list[tuple[int, str, list[str]]]] = {}

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

    for variant in vcf:
        norm = normalize_chrom(variant.CHROM)
        bucket = variant.POS // 1_000_000
        key = (norm, bucket)
        if key not in variant_buffers:
            work_order.append(key)
            variant_buffers[key] = []
            query_buffers[key] = []
        variant_buffers[key].append(variant)
        query_buffers[key].append((variant.POS, variant.REF, list(variant.ALT)))

    # Phase 2: parallel compute
    db_path = str(engine._db)
    n_units = len(work_order)
    if n_workers is None:
        effective = max(1, min(os.cpu_count() or 1, n_units)) if n_units > 0 else 1
    else:
        effective = max(1, min(n_workers, n_units)) if n_units > 0 else 1

    unit_annotations: dict[tuple[str, int], dict] = {}
    if effective == 1:
        for chrom, bucket in work_order:
            unit_annotations[(chrom, bucket)] = _compute_chunk_annotations(
                db_path, chrom, bucket, query_buffers[(chrom, bucket)], sf
            )
    else:
        with ProcessPoolExecutor(max_workers=effective) as executor:
            futures = {
                executor.submit(
                    _compute_chunk_annotations,
                    db_path, chrom, bucket, query_buffers[(chrom, bucket)], sf,
                ): (chrom, bucket)
                for chrom, bucket in work_order
            }
            for future in as_completed(futures):
                unit_annotations[futures[future]] = future.result()

    # Phase 3: write annotated records in VCF input order
    stats = {"n_variants": 0, "n_annotated": 0, "n_uncovered": 0}
    for (chrom, bucket) in work_order:
        ann = unit_annotations[(chrom, bucket)]
        for variant in variant_buffers[(chrom, bucket)]:
            pos = variant.POS
            stats["n_variants"] += 1

            # AN is per-position; read from any alt's entry (same value for all alts)
            AN = 0
            if variant.ALT:
                AN = ann.get((pos, variant.REF, variant.ALT[0]), (0, 0, False))[1]

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
                AC, _, in_parquet = ann.get(key, (0, AN, False))
                if in_parquet:
                    any_found = True
                ac_list.append(AC)
                af_list.append(AC / AN)

            # cyvcf2 Number=A fields must be set as comma-separated strings
            variant.INFO["AFQUERY_AC"] = ",".join(str(v) for v in ac_list)
            variant.INFO["AFQUERY_AN"] = AN
            variant.INFO["AFQUERY_AF"] = ",".join("%.10g" % v for v in af_list)
            if any_found:
                stats["n_annotated"] += 1
            writer.write_record(variant)

    writer.close()
    vcf.close()
    return stats
