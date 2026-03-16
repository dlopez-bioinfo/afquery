import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import tqdm

import pyarrow as pa
import pyarrow.parquet as pq

from ..constants import normalize_chrom
from ..models import Sample

logger = logging.getLogger(__name__)

INGEST_SCHEMA = pa.schema([
    ("chrom",       pa.utf8()),
    ("pos",         pa.uint32()),
    ("ref",         pa.utf8()),
    ("alt",         pa.utf8()),
    ("gt_ac",       pa.uint8()),
    ("sample_id",   pa.uint32()),
    ("filter_pass", pa.bool_()),
])


class IngestError(RuntimeError):
    pass


def ingest_sample(sample_id: int, vcf_path: str, tmp_dir: str) -> tuple[str, float]:
    """Parse one VCF, write Parquet file. Returns (output_path, elapsed_seconds)."""
    import cyvcf2  # local import for worker-safety

    t0 = time.monotonic()
    out_path = os.path.join(tmp_dir, f"sample_{sample_id}.parquet")

    chroms: list[str] = []
    positions: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    gt_acs: list[int] = []
    sample_ids: list[int] = []
    filter_passes: list[bool] = []

    vcf = cyvcf2.VCF(vcf_path)
    for variant in vcf:
        if not (variant.is_snp or variant.is_indel):
            continue

        gt = variant.genotypes[0]
        alleles = [a for a in gt[:-1] if a >= 0]

        chrom = normalize_chrom(variant.CHROM)
        pos = variant.POS
        ref = variant.REF
        # cyvcf2: FILTER is None for PASS/missing, string for others (e.g. "LowQual")
        fp = variant.FILTER is None or variant.FILTER == "PASS"

        # Missing GT (./.) at a failed site: track as N_FAIL for all ALTs
        if not alleles:
            if not fp:
                for alt_str in variant.ALT:
                    if alt_str == "*":
                        continue
                    chroms.append(chrom)
                    positions.append(pos)
                    refs.append(ref)
                    alts.append(alt_str)
                    gt_acs.append(0)
                    sample_ids.append(sample_id)
                    filter_passes.append(False)
            continue

        for idx, alt_str in enumerate(variant.ALT):
            if alt_str == "*":
                continue
            alt_allele_num = idx + 1
            ac = alleles.count(alt_allele_num)
            if ac == 0:
                continue

            chroms.append(chrom)
            positions.append(pos)
            refs.append(ref)
            alts.append(alt_str)
            gt_acs.append(ac)
            sample_ids.append(sample_id)
            filter_passes.append(fp)

    vcf.close()

    table = pa.table(
        {
            "chrom":       pa.array(chroms,        type=pa.utf8()),
            "pos":         pa.array(positions,     type=pa.uint32()),
            "ref":         pa.array(refs,          type=pa.utf8()),
            "alt":         pa.array(alts,          type=pa.utf8()),
            "gt_ac":       pa.array(gt_acs,        type=pa.uint8()),
            "sample_id":   pa.array(sample_ids,    type=pa.uint32()),
            "filter_pass": pa.array(filter_passes, type=pa.bool_()),
        },
        schema=INGEST_SCHEMA,
    )

    pq.write_table(table, out_path)
    elapsed = time.monotonic() - t0
    return out_path, elapsed


def _ingest_worker(args: tuple) -> str:
    sample_id, vcf_path, tmp_dir = args
    return ingest_sample(sample_id, vcf_path, tmp_dir)


def ingest_all(
    samples: list[Sample],
    vcf_paths: list[str],
    tmp_dir: str,
    n_workers: int = 8,
    resume: bool = True,
) -> list[str]:
    """Parallel ingestion via ProcessPoolExecutor. Collects all errors before raising."""
    args_list = [(s.sample_id, vcf_paths[i], tmp_dir) for i, s in enumerate(samples)]

    if resume:
        pending = []
        n_skip = 0
        for a in args_list:
            parquet_path = os.path.join(tmp_dir, f"sample_{a[0]}.parquet")
            if os.path.exists(parquet_path):
                n_skip += 1
            else:
                pending.append(a)
        if n_skip:
            logger.info("[ingest] Skipping %d already-ingested sample(s); ingesting %d remaining.",
                        n_skip, len(pending))
        args_list = pending

    if not args_list:
        logger.info("[ingest] All samples already ingested, nothing to do.")
        return []

    logger.info("[ingest] Ingesting %d VCF(s) with %d worker(s)...", len(args_list), n_workers)
    t0 = time.monotonic()

    results: list[str] = []
    errors: list[str] = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_args = {
            executor.submit(_ingest_worker, args): args for args in args_list
        }
        with tqdm.tqdm(
            total=len(args_list),
            unit="VCF",
            desc="[ingest]",
            dynamic_ncols=True,
            disable=not sys.stderr.isatty(),
        ) as bar:
            for future in as_completed(future_to_args):
                args = future_to_args[future]
                try:
                    out_path, elapsed = future.result()
                    results.append(out_path)
                    bar.update(1)
                    logger.debug("  [ingest] sample_id=%d (%s) done in %.2fs",
                                 args[0], os.path.basename(args[1]), elapsed)
                except Exception as e:
                    errors.append(f"Sample {args[0]} ({args[1]}): {e}")

    logger.info("[ingest] Complete: %d ok, %d error(s) (%.1fs)",
                len(results), len(errors), time.monotonic() - t0)

    if errors:
        raise IngestError("VCF ingestion failed:\n" + "\n".join(errors))

    return results
