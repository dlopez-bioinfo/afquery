from pathlib import Path

import pytest

from afquery.preprocess.synth import generate_synthetic_manifest


def test_generates_manifest_file(tmp_path):
    """generate_synthetic_manifest returns a manifest.tsv Path that exists."""
    manifest_path = generate_synthetic_manifest(
        tmp_path, n_samples=5, n_variants_per_chrom=10, chroms=("chr1",)
    )
    assert manifest_path.exists()
    assert manifest_path.name == "manifest.tsv"


def test_manifest_has_correct_sample_count(tmp_path):
    """manifest.tsv has a header + n_samples data lines."""
    manifest_path = generate_synthetic_manifest(
        tmp_path, n_samples=8, n_variants_per_chrom=5, chroms=("chr1",)
    )
    lines = [l for l in manifest_path.read_text().splitlines() if l.strip()]
    assert lines[0].startswith("sample_name")  # header
    assert len(lines) == 9  # 1 header + 8 samples


def test_vcf_files_are_created(tmp_path):
    """A VCF file is created for each sample."""
    n = 6
    manifest_path = generate_synthetic_manifest(
        tmp_path, n_samples=n, n_variants_per_chrom=10, chroms=("chr1",)
    )
    vcf_dir = tmp_path / "vcfs"
    vcf_files = list(vcf_dir.glob("*.vcf"))
    assert len(vcf_files) == n


def test_sex_alternates(tmp_path):
    """Sex alternates male/female across samples."""
    manifest_path = generate_synthetic_manifest(
        tmp_path, n_samples=4, n_variants_per_chrom=5, chroms=("chr1",)
    )
    lines = manifest_path.read_text().splitlines()[1:]  # skip header
    sexes = [l.split("\t")[1] for l in lines if l.strip()]
    assert sexes == ["male", "female", "male", "female"]


def test_seed_reproducibility(tmp_path):
    """Same seed produces identical VCF content (paths differ but genotypes match)."""
    generate_synthetic_manifest(
        tmp_path / "a", n_samples=3, n_variants_per_chrom=5, seed=42
    )
    generate_synthetic_manifest(
        tmp_path / "b", n_samples=3, n_variants_per_chrom=5, seed=42
    )
    # Compare VCF content (not manifest, which contains absolute paths)
    vcf_a = (tmp_path / "a" / "vcfs" / "synth_000000.vcf").read_text()
    vcf_b = (tmp_path / "b" / "vcfs" / "synth_000000.vcf").read_text()
    assert vcf_a == vcf_b


def test_different_seeds_differ(tmp_path):
    """Different seeds produce different VCF content."""
    generate_synthetic_manifest(
        tmp_path / "a", n_samples=2, n_variants_per_chrom=20, seed=1
    )
    generate_synthetic_manifest(
        tmp_path / "b", n_samples=2, n_variants_per_chrom=20, seed=2
    )
    vcf_a = (tmp_path / "a" / "vcfs" / "synth_000000.vcf").read_text()
    vcf_b = (tmp_path / "b" / "vcfs" / "synth_000000.vcf").read_text()
    assert vcf_a != vcf_b


def test_synth_integrates_with_preprocess(tmp_path, tmp_path_factory):
    """Synthetic manifest can be processed by run_preprocess end-to-end."""
    from afquery.preprocess import run_preprocess
    from afquery.database import Database

    manifest_path = generate_synthetic_manifest(
        tmp_path / "synth",
        n_samples=5,
        n_variants_per_chrom=50,
        chroms=("chr1",),
    )

    db_path = str(tmp_path_factory.mktemp("synth_db"))
    run_preprocess(
        manifest_path=str(manifest_path),
        output_dir=db_path,
        genome_build="GRCh37",
        n_threads=2,
    )

    db = Database(db_path)
    info = db.info()
    assert info["sample_count"] == 5
    assert info["genome_build"] == "GRCh37"
