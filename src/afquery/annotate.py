import logging
import os
import warnings
import duckdb

from .bitmaps import deserialize
from .constants import normalize_chrom
from .models import AfqueryWarning, SampleFilter
from .ploidy import split_ploidy

logger = logging.getLogger(__name__)


def _compute_chunk_annotations(
    db_path: str,
    chrom: str,
    bucket_id: int,
    records: list[tuple[int, str, list[str]]],
    sf: SampleFilter,
) -> dict[tuple[int, str, str], tuple[int, int, bool, int, int, int, int]]:
    """Compute (AC, AN, in_parquet, N_FAIL, N_HET, N_HOM_ALT, N_HOM_REF) for each (pos, ref, alt).

    No cyvcf2 dependency — safe to run in a subprocess worker.
    """
    from pathlib import Path
    from .query import QueryEngine

    engine = QueryEngine(db_path)
    if chrom not in engine._all_known_chroms:
        warnings.warn(
            f"Chromosome {chrom!r} has no data — variants on this chromosome will get AN=0.",
            AfqueryWarning, stacklevel=2,
        )
    sample_bm = engine._build_sample_bitmap(sf)
    male_bm = engine._male_bm
    female_bm = engine._female_bm

    unique_positions = list({pos for pos, _ref, _alts in records})

    pos_data: dict[int, tuple] = {}
    for pos in unique_positions:
        eligible, AN = engine._compute_eligible(chrom, pos, sample_bm)
        pos_data[pos] = (eligible, AN)

    valid_positions = [p for p in unique_positions if pos_data[p][1] > 0]

    variant_data: dict[tuple[int, str, str], tuple[bytes, bytes, bytes]] = {}
    _db = Path(db_path)
    bucket_start = bucket_id * 1_000_000
    bucket_end = (bucket_id + 1) * 1_000_000 - 1

    if chrom in engine._partitioned_chroms:
        parquet_file = _db / "variants" / chrom / f"bucket_{bucket_id}.parquet"
        if valid_positions and parquet_file.exists():
            con = duckdb.connect()
            placeholders = ", ".join("?" * len(valid_positions))
            rows = con.execute(
                f"SELECT pos, ref, alt, het_bitmap, hom_bitmap, fail_bitmap"
                f" FROM read_parquet('{parquet_file}') WHERE pos IN ({placeholders})",
                valid_positions,
            ).fetchall()
            con.close()
            for row in rows:
                pos, ref, alt, het_bytes, hom_bytes, fail_bytes = row
                variant_data[(pos, ref, alt)] = (bytes(het_bytes), bytes(hom_bytes), bytes(fail_bytes))
    else:
        parquet_file = _db / "variants" / f"{chrom}.parquet"
        if valid_positions and parquet_file.exists():
            con = duckdb.connect()
            rows = con.execute(
                f"SELECT pos, ref, alt, het_bitmap, hom_bitmap, fail_bitmap"
                f" FROM read_parquet('{parquet_file}') WHERE pos BETWEEN ? AND ?",
                [bucket_start, bucket_end],
            ).fetchall()
            con.close()
            valid_pos_set = set(valid_positions)
            for row in rows:
                pos, ref, alt, het_bytes, hom_bytes, fail_bytes = row
                if pos in valid_pos_set:
                    variant_data[(pos, ref, alt)] = (bytes(het_bytes), bytes(hom_bytes), bytes(fail_bytes))

    result: dict[tuple[int, str, str], tuple[int, int, bool, int, int, int, int]] = {}
    for pos, ref, alts in records:
        eligible, AN = pos_data[pos]
        for alt in alts:
            key = (pos, ref, alt)
            if key in result:
                continue  # dedup
            if AN == 0:
                result[key] = (0, 0, False, 0, 0, 0, 0)
            elif key in variant_data:
                het_bytes, hom_bytes, fail_bytes = variant_data[key]
                het_bm = deserialize(het_bytes)
                hom_bm = deserialize(hom_bytes)
                N_HET = len(het_bm & eligible)
                N_HOM_ALT = len(hom_bm & eligible)
                haploid_elig, diploid_elig = split_ploidy(
                    eligible, male_bm, female_bm, chrom, pos, engine._genome_build
                )
                het_elig = het_bm & eligible
                hom_elig = hom_bm & eligible
                AC = (len((het_elig | hom_elig) & haploid_elig)
                      + len(het_elig & diploid_elig)
                      + 2 * len(hom_elig & diploid_elig))
                N_FAIL: int = len(deserialize(fail_bytes) & eligible)
                N_HOM_REF = len(eligible) - N_HET - N_HOM_ALT - N_FAIL
                result[key] = (AC, AN, True, N_FAIL, N_HET, N_HOM_ALT, N_HOM_REF)
            else:
                result[key] = (0, AN, False, 0, 0, 0, len(eligible))

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
    vcf.add_info_to_header({
        "ID": "AFQUERY_N_HET", "Number": "A", "Type": "Integer",
        "Description": "Count of heterozygous eligible samples per alt allele",
    })
    vcf.add_info_to_header({
        "ID": "AFQUERY_N_HOM_ALT", "Number": "A", "Type": "Integer",
        "Description": "Count of homozygous alt eligible samples per alt allele",
    })
    vcf.add_info_to_header({
        "ID": "AFQUERY_N_HOM_REF", "Number": "A", "Type": "Integer",
        "Description": "Count of homozygous ref eligible samples per alt allele",
    })
    vcf.add_info_to_header({
        "ID": "AFQUERY_N_FAIL", "Number": "1", "Type": "Integer",
        "Description": "Eligible samples with FILTER!=PASS for this variant",
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

    n_variants = sum(len(v) for v in variant_buffers.values())
    logger.info("[annotate] Loaded %d variant(s) in %d chunk(s) from %s",
                n_variants, len(work_order), os.path.basename(input_vcf))

    # Phase 2: parallel compute
    db_path = str(engine._db)
    n_units = len(work_order)
    if n_workers is None:
        effective = max(1, min(os.cpu_count() or 1, n_units)) if n_units > 0 else 1
    else:
        effective = max(1, min(n_workers, n_units)) if n_units > 0 else 1

    logger.debug("[annotate] Using %d worker(s) for %d chunk(s).", effective, len(work_order))

    unit_annotations: dict[tuple[str, int], dict] = {}
    if effective == 1:
        for chrom, bucket in work_order:
            unit_annotations[(chrom, bucket)] = _compute_chunk_annotations(
                db_path, chrom, bucket, query_buffers[(chrom, bucket)], sf
            )
            logger.debug("  [annotate] chunk %s/bucket_%d: done", chrom, bucket)
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
                chrom, bucket = futures[future]
                unit_annotations[(chrom, bucket)] = future.result()
                logger.debug("  [annotate] chunk %s/bucket_%d: done", chrom, bucket)

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
            n_het_list = []
            n_hom_alt_list = []
            n_hom_ref_list = []
            any_found = False
            fail_count_total = 0

            for alt in variant.ALT:
                key = (pos, variant.REF, alt)
                AC, _, in_parquet, n_fail, n_het, n_hom_alt, n_hom_ref = ann.get(
                    key, (0, AN, False, 0, 0, 0, 0)
                )
                if in_parquet:
                    any_found = True
                ac_list.append(AC)
                af_list.append(AC / AN)
                n_het_list.append(n_het)
                n_hom_alt_list.append(n_hom_alt)
                n_hom_ref_list.append(n_hom_ref)
                fail_count_total += n_fail

            # cyvcf2 Number=A fields must be set as comma-separated strings
            variant.INFO["AFQUERY_AC"] = ",".join(str(v) for v in ac_list)
            variant.INFO["AFQUERY_AN"] = AN
            variant.INFO["AFQUERY_AF"] = ",".join("%.10g" % v for v in af_list)
            variant.INFO["AFQUERY_N_HET"] = ",".join(str(v) for v in n_het_list)
            variant.INFO["AFQUERY_N_HOM_ALT"] = ",".join(str(v) for v in n_hom_alt_list)
            variant.INFO["AFQUERY_N_HOM_REF"] = ",".join(str(v) for v in n_hom_ref_list)
            variant.INFO["AFQUERY_N_FAIL"] = fail_count_total
            if any_found:
                stats["n_annotated"] += 1
            writer.write_record(variant)

    writer.close()
    vcf.close()

    logger.info("[annotate] Complete: %d annotated, %d uncovered, %d total.",
                stats["n_annotated"], stats["n_uncovered"], stats["n_variants"])

    return stats
