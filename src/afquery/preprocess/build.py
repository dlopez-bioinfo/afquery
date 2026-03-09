import glob as glob_module
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from pyroaring import BitMap

from ..bitmaps import serialize
from ..constants import ALL_CHROMS

logger = logging.getLogger(__name__)

BUCKET_SIZE = 1_000_000

PARQUET_SCHEMA = pa.schema([
    ("pos",        pa.uint32()),
    ("ref",        pa.large_utf8()),
    ("alt",        pa.large_utf8()),
    ("het_bitmap", pa.large_binary()),
    ("hom_bitmap", pa.large_binary()),
])


def get_chroms_in_temp_files(tmp_dir: str) -> list[str]:
    """SELECT DISTINCT chrom from all sample_*.parquet in tmp_dir."""
    parquet_files = glob_module.glob(os.path.join(tmp_dir, "sample_*.parquet"))
    if not parquet_files:
        return []

    glob_pattern = os.path.join(tmp_dir, "sample_*.parquet").replace("'", "''")
    con = duckdb.connect()
    try:
        rows = con.execute(
            f"SELECT DISTINCT chrom FROM read_parquet('{glob_pattern}')"
        ).fetchall()
    except Exception:
        return []
    finally:
        con.close()

    return [r[0] for r in rows if r[0] in ALL_CHROMS]


def _make_table(positions, refs, alts, het_bitmaps, hom_bitmaps) -> pa.Table:
    return pa.table(
        {
            "pos":        pa.array(positions,   type=pa.uint32()),
            "ref":        pa.array(refs,        type=pa.large_utf8()),
            "alt":        pa.array(alts,        type=pa.large_utf8()),
            "het_bitmap": pa.array(het_bitmaps, type=pa.large_binary()),
            "hom_bitmap": pa.array(hom_bitmaps, type=pa.large_binary()),
        },
        schema=PARQUET_SCHEMA,
    )


def _write_flat(
    chrom: str,
    variants_dir: str,
    positions, refs, alts, het_bitmaps, hom_bitmaps,
    row_group_size: int,
) -> None:
    table = _make_table(positions, refs, alts, het_bitmaps, hom_bitmaps)
    out_path = os.path.join(variants_dir, f"{chrom}.parquet")
    tmp_path = out_path + ".tmp"
    pq.write_table(table, tmp_path, row_group_size=row_group_size)
    os.rename(tmp_path, out_path)


def _write_partitioned(
    chrom: str,
    variants_dir: str,
    positions, refs, alts, het_bitmaps, hom_bitmaps,
    row_group_size: int,
) -> None:
    chrom_dir = os.path.join(variants_dir, chrom)
    os.makedirs(chrom_dir, exist_ok=True)

    # Group rows by bucket (pos // BUCKET_SIZE)
    buckets: dict[int, tuple[list, list, list, list, list]] = {}
    for pos, ref, alt, het_bm, hom_bm in zip(
        positions, refs, alts, het_bitmaps, hom_bitmaps
    ):
        bucket = pos // BUCKET_SIZE
        if bucket not in buckets:
            buckets[bucket] = ([], [], [], [], [])
        b = buckets[bucket]
        b[0].append(pos)
        b[1].append(ref)
        b[2].append(alt)
        b[3].append(het_bm)
        b[4].append(hom_bm)

    for bucket_id, (bpos, brefs, balts, bhets, bhoms) in buckets.items():
        out_path = os.path.join(chrom_dir, f"bucket_{bucket_id}.parquet")
        tmp_path = out_path + ".tmp"
        table = _make_table(bpos, brefs, balts, bhets, bhoms)
        pq.write_table(table, tmp_path, row_group_size=row_group_size)
        os.rename(tmp_path, out_path)


def build_chromosome_parquet(
    chrom: str,
    tmp_dir: str,
    variants_dir: str,
    row_group_size: int = 100_000,
    partitioned: bool = False,
) -> int:
    """Build Parquet for one chromosome. Returns variant count."""
    glob_pattern = os.path.join(tmp_dir, "sample_*.parquet").replace("'", "''")

    con = duckdb.connect()
    rows = con.execute(
        f"""
        SELECT
            pos, ref, alt,
            list(sample_id ORDER BY sample_id) AS sample_ids,
            list(gt_ac     ORDER BY sample_id) AS gt_acs
        FROM read_parquet('{glob_pattern}')
        WHERE chrom = ?
        GROUP BY pos, ref, alt
        ORDER BY pos, alt
        """,
        [chrom],
    ).fetchall()
    con.close()

    if not rows:
        return 0

    positions = []
    refs = []
    alts = []
    het_bitmaps = []
    hom_bitmaps = []

    for pos, ref, alt, sample_ids, gt_acs in rows:
        het_ids = [sid for sid, ac in zip(sample_ids, gt_acs) if ac == 1]
        hom_ids = [sid for sid, ac in zip(sample_ids, gt_acs) if ac == 2]
        positions.append(pos)
        refs.append(ref)
        alts.append(alt)
        het_bitmaps.append(serialize(BitMap(het_ids)))
        hom_bitmaps.append(serialize(BitMap(hom_ids)))

    if partitioned:
        _write_partitioned(
            chrom, variants_dir,
            positions, refs, alts, het_bitmaps, hom_bitmaps,
            row_group_size,
        )
    else:
        _write_flat(
            chrom, variants_dir,
            positions, refs, alts, het_bitmaps, hom_bitmaps,
            row_group_size,
        )

    return len(rows)


def _build_chrom_worker(
    chrom: str,
    tmp_dir: str,
    variants_dir: str,
    row_group_size: int,
    partitioned: bool,
) -> tuple[str, int]:
    """Top-level worker for ProcessPoolExecutor (must be picklable)."""
    count = build_chromosome_parquet(
        chrom, tmp_dir, variants_dir, row_group_size, partitioned
    )
    return chrom, count


def build_all_parquets(
    tmp_dir: str,
    variants_dir: str,
    row_group_size: int = 100_000,
    n_workers: int | None = None,
    partitioned: bool = True,
) -> dict[str, int]:
    """Build Parquet for all discovered chromosomes. Returns {chrom: count}.

    Uses ProcessPoolExecutor for parallel chromosome builds.
    n_workers=None uses os.cpu_count(). partitioned=True writes bucketed files.
    """
    chroms = get_chroms_in_temp_files(tmp_dir)
    valid_chroms = []
    for chrom in chroms:
        if chrom in ALL_CHROMS:
            valid_chroms.append(chrom)
        else:
            logger.warning("Skipping unknown contig: %s", chrom)

    if not valid_chroms:
        return {}

    result: dict[str, int] = {}
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(
                _build_chrom_worker,
                chrom, tmp_dir, variants_dir, row_group_size, partitioned,
            ): chrom
            for chrom in valid_chroms
        }
        for fut in as_completed(futures):
            chrom, count = fut.result()  # propagates exceptions
            if count > 0:
                result[chrom] = count

    return result
