#!/usr/bin/env python3
"""Generate synthetic demo data for the AFQuery tutorial.

Creates:
- 10 single-sample VCFs (~100 variants each) in demo_output/vcfs/
- A manifest.tsv in demo_output/
- BED files for 2 WES technologies in demo_output/beds/

All data is synthetic and not based on real patient genomes.
"""

import gzip
import os
import random
from pathlib import Path

SEED = 42
OUTPUT_DIR = Path(__file__).resolve().parent / "demo_output"

# --- Sample definitions -------------------------------------------------------

SAMPLES = [
    # (name, sex, tech, phenotype_codes)
    ("DEMO_001", "female", "wgs",      "E11.9,I10"),
    ("DEMO_002", "male",   "wgs",      "E11.9"),
    ("DEMO_003", "female", "wgs",      "I10"),
    ("DEMO_004", "male",   "wgs",      "control"),
    ("DEMO_005", "female", "wes_v1",   "E11.9,I10"),
    ("DEMO_006", "male",   "wes_v1",   "control"),
    ("DEMO_007", "female", "wes_v1",   "I10"),
    ("DEMO_008", "male",   "wes_v2",   "E11.9"),
    ("DEMO_009", "female", "wes_v2",   "control"),
    ("DEMO_010", "male",   "wes_v2",   "I10,control"),
]

# --- Variant pool -------------------------------------------------------------
# Shared pool of positions; each sample gets a random subset with random GTs.

CHROMS_AUTOSOME = ["chr1", "chr2", "chr7", "chr12", "chr15"]
BASES = ["A", "C", "G", "T"]

# WES capture regions (shared by both kits with slight differences)
WES_V1_REGIONS = [
    # (chrom, start, end)
    ("chr1",  900000,  950000),
    ("chr1",  1200000, 1250000),
    ("chr2",  1000000, 1080000),
    ("chr7",  5000000, 5100000),
    ("chr12", 2000000, 2050000),
    ("chr15", 4000000, 4100000),
]

WES_V2_REGIONS = WES_V1_REGIONS + [
    # V2 adds extra regions
    ("chr1",  1800000, 1850000),
    ("chr7",  7000000, 7050000),
    ("chr12", 3000000, 3050000),
]


def _alt_for(ref: str) -> str:
    """Pick a random ALT allele different from REF."""
    choices = [b for b in BASES if b != ref]
    return random.choice(choices)


def _generate_variant_pool(n_per_chrom: int = 20) -> list[tuple[str, int, str, str]]:
    """Generate a pool of (chrom, pos, ref, alt) variants."""
    pool = []

    for chrom in CHROMS_AUTOSOME:
        for _ in range(n_per_chrom):
            pos = random.randint(900000, 10000000)
            ref = random.choice(BASES)
            alt = _alt_for(ref)
            pool.append((chrom, pos, ref, alt))

    # Add chrX, chrY, chrM variants
    for _ in range(10):
        pos = random.randint(5000000, 100000000)
        ref = random.choice(BASES)
        pool.append(("chrX", pos, ref, _alt_for(ref)))

    for _ in range(5):
        pos = random.randint(2700000, 10000000)
        ref = random.choice(BASES)
        pool.append(("chrY", pos, ref, _alt_for(ref)))

    for _ in range(5):
        pos = random.randint(100, 16500)
        ref = random.choice(BASES)
        pool.append(("chrM", pos, ref, _alt_for(ref)))

    # Sort by chrom, pos
    chrom_order = {c: i for i, c in enumerate(
        CHROMS_AUTOSOME + ["chrX", "chrY", "chrM"]
    )}
    pool.sort(key=lambda v: (chrom_order.get(v[0], 99), v[1]))

    return pool


def _write_vcf(path: Path, sample_name: str, sex: str,
               variants: list[tuple[str, int, str, str]]) -> None:
    """Write a minimal single-sample VCF."""
    header = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FILTER=<ID=LowQual,Description="Low quality">\n'
    )

    # Add contig headers
    contigs = sorted({v[0] for v in variants})
    for c in contigs:
        header += f"##contig=<ID={c}>\n"

    header += f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"

    lines = []
    for chrom, pos, ref, alt in variants:
        # Random genotype
        r = random.random()
        if r < 0.60:
            gt = "0/0"  # hom ref (most common)
        elif r < 0.85:
            gt = "0/1"  # het
        elif r < 0.95:
            gt = "1/1"  # hom alt
        else:
            continue  # skip (not called)

        # Haploid on chrY/chrM for males, chrY excluded for females
        if chrom == "chrY":
            if sex == "female":
                continue
            gt = gt.replace("0/", "").replace("/0", "").replace("/1", "")
            if gt == "0":
                gt = "0"
            elif "1" in gt:
                gt = "1"

        if chrom == "chrM":
            gt = "1" if "1" in gt else "0"

        # Small chance of LowQual filter
        filt = "LowQual" if random.random() < 0.05 else "PASS"

        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filt}\t.\tGT\t{gt}\n")

    with gzip.open(path, "wt") as f:
        f.write(header)
        f.writelines(lines)


def _write_bed(path: Path, regions: list[tuple[str, int, int]]) -> None:
    """Write a BED file."""
    with open(path, "w") as f:
        for chrom, start, end in sorted(regions):
            f.write(f"{chrom}\t{start}\t{end}\n")


def main() -> None:
    random.seed(SEED)

    vcf_dir = OUTPUT_DIR / "vcfs"
    bed_dir = OUTPUT_DIR / "beds"
    vcf_dir.mkdir(parents=True, exist_ok=True)
    bed_dir.mkdir(parents=True, exist_ok=True)

    # Generate variant pool
    pool = _generate_variant_pool(n_per_chrom=20)

    # Write VCFs
    for name, sex, tech, pheno in SAMPLES:
        # Each sample gets ~80% of the pool (random subset)
        sample_variants = [v for v in pool if random.random() < 0.80]
        vcf_path = vcf_dir / f"{name}.vcf.gz"
        _write_vcf(vcf_path, name, sex, sample_variants)

    # Write BED files
    _write_bed(bed_dir / "wes_v1.bed", WES_V1_REGIONS)
    _write_bed(bed_dir / "wes_v2.bed", WES_V2_REGIONS)

    # Write manifest
    manifest_path = OUTPUT_DIR / "manifest.tsv"
    with open(manifest_path, "w") as f:
        f.write("sample_name\tvcf_path\tsex\ttech_name\tphenotype_codes\n")
        for name, sex, tech, pheno in SAMPLES:
            vcf_abs = str((vcf_dir / f"{name}.vcf.gz").resolve())
            f.write(f"{name}\t{vcf_abs}\t{sex}\t{tech}\t{pheno}\n")

    print(f"Demo data created in: {OUTPUT_DIR}")
    print(f"  VCFs:     {vcf_dir} ({len(SAMPLES)} files)")
    print(f"  BED files: {bed_dir} (wes_v1.bed, wes_v2.bed)")
    print(f"  Manifest: {manifest_path}")
    print()
    print("Next steps:")
    print(f"  afquery create-db --manifest {manifest_path} --output-dir ./demo_db/ "
          f"--genome-build GRCh38 --bed-dir {bed_dir}")


if __name__ == "__main__":
    main()
