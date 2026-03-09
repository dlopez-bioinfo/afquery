import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import pyarrow as pa
import pyarrow.parquet as pq

from ..constants import normalize_chrom
from ..models import Sample

logger = logging.getLogger(__name__)

INGEST_SCHEMA = pa.schema([
    ("chrom",     pa.utf8()),
    ("pos",       pa.uint32()),
    ("ref",       pa.utf8()),
    ("alt",       pa.utf8()),
    ("gt_ac",     pa.uint8()),
    ("sample_id", pa.uint32()),
])


class IngestError(RuntimeError):
    pass


def ingest_sample(sample_id: int, vcf_path: str, tmp_dir: str) -> str:
    """Parse one VCF, write Parquet file. Returns output path."""
    import cyvcf2  # local import for worker-safety

    out_path = os.path.join(tmp_dir, f"sample_{sample_id}.parquet")

    chroms: list[str] = []
    positions: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    gt_acs: list[int] = []
    sample_ids: list[int] = []

    vcf = cyvcf2.VCF(vcf_path)
    for variant in vcf:
        if not (variant.is_snp or variant.is_indel):
            continue

        gt = variant.genotypes[0]
        alleles = [a for a in gt[:-1] if a >= 0]
        if not alleles:
            continue

        chrom = normalize_chrom(variant.CHROM)
        pos = variant.POS
        ref = variant.REF

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

    vcf.close()

    table = pa.table(
        {
            "chrom":     pa.array(chroms,     type=pa.utf8()),
            "pos":       pa.array(positions,  type=pa.uint32()),
            "ref":       pa.array(refs,       type=pa.utf8()),
            "alt":       pa.array(alts,       type=pa.utf8()),
            "gt_ac":     pa.array(gt_acs,     type=pa.uint8()),
            "sample_id": pa.array(sample_ids, type=pa.uint32()),
        },
        schema=INGEST_SCHEMA,
    )

    pq.write_table(table, out_path)
    return out_path


def _ingest_worker(args: tuple) -> str:
    sample_id, vcf_path, tmp_dir = args
    return ingest_sample(sample_id, vcf_path, tmp_dir)


def ingest_all(
    samples: list[Sample],
    vcf_paths: list[str],
    tmp_dir: str,
    n_workers: int = 8,
) -> list[str]:
    """Parallel ingestion via ProcessPoolExecutor. Collects all errors before raising."""
    args_list = [(s.sample_id, vcf_paths[i], tmp_dir) for i, s in enumerate(samples)]

    results: list[str] = []
    errors: list[str] = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_args = {
            executor.submit(_ingest_worker, args): args for args in args_list
        }
        for future in as_completed(future_to_args):
            args = future_to_args[future]
            try:
                results.append(future.result())
            except Exception as e:
                errors.append(f"Sample {args[0]} ({args[1]}): {e}")

    if errors:
        raise IngestError("VCF ingestion failed:\n" + "\n".join(errors))

    return results
