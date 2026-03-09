import json
import sys

import click

from .database import Database


@click.group()
def cli():
    """afquery: genomic allele frequency query engine."""


@cli.command()
@click.option("--db",    required=True, help="Path to database directory.")
@click.option("--chrom", required=True, help="Chromosome (e.g. 1, chrX).")
@click.option("--pos",   required=True, type=int, help="1-based position.")
@click.option("--phenotype", multiple=True, help="Phenotype code(s). Repeatable. If omitted, all samples are used.")
@click.option("--sex",   default="both", type=click.Choice(["male", "female", "both"]))
@click.option("--ref",   default=None, help="Reference allele (optional filter).")
@click.option("--alt",   default=None, help="Alternate allele (optional filter).")
@click.option("--format", "fmt", default="text", type=click.Choice(["text", "json", "tsv"]))
def query(db, chrom, pos, phenotype, sex, ref, alt, fmt):
    """Query allele frequency at a single position."""
    database = Database(db)
    results = database.query(chrom=chrom, pos=pos, phenotype=list(phenotype), sex=sex, ref=ref, alt=alt)

    if fmt == "json":
        out = []
        for r in results:
            out.append({
                "chrom": r.variant.chrom, "pos": r.variant.pos,
                "ref": r.variant.ref, "alt": r.variant.alt,
                "AC": r.AC, "AN": r.AN,
                "AF": r.AF, "n_eligible": r.n_samples_eligible,
            })
        click.echo(json.dumps(out, indent=2))
    elif fmt == "tsv":
        click.echo("chrom\tpos\tref\talt\tAC\tAN\tAF\tn_eligible")
        for r in results:
            af = f"{r.AF:.6f}" if r.AF is not None else "NA"
            click.echo(
                f"{r.variant.chrom}\t{r.variant.pos}\t{r.variant.ref}\t"
                f"{r.variant.alt}\t{r.AC}\t{r.AN}\t{af}\t{r.n_samples_eligible}"
            )
    else:  # text
        if not results:
            click.echo("No variants found at this position for the given filters.")
            return
        for r in results:
            af = f"{r.AF:.4f}" if r.AF is not None else "NA"
            click.echo(
                f"{r.variant.chrom}:{r.variant.pos} {r.variant.ref}>{r.variant.alt}  "
                f"AC={r.AC}  AN={r.AN}  AF={af}  n_eligible={r.n_samples_eligible}"
            )


@cli.command("query-batch")
@click.option("--db",       required=True, help="Path to database directory.")
@click.option("--chrom",    required=True, help="Chromosome (e.g. 1, chrX).")
@click.option("--variants", required=True, type=click.Path(exists=True), help="TSV file with columns: pos ref alt (no header, whitespace-separated).")
@click.option("--phenotype",    multiple=True, help="Phenotype code(s). Repeatable. If omitted, all samples are used.")
@click.option("--sex",      default="both", type=click.Choice(["male", "female", "both"]))
@click.option("--format",   "fmt", default="text", type=click.Choice(["text", "json", "tsv"]))
def query_batch(db, chrom, variants, phenotype, sex, fmt):
    """Query allele frequency at multiple variants (read from TSV file: pos ref alt)."""
    variant_list = []
    with open(variants) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            variant_list.append((int(parts[0]), parts[1], parts[2]))

    database = Database(db)
    results = database.query_batch(chrom=chrom, variants=variant_list, phenotype=list(phenotype), sex=sex)

    if fmt == "json":
        out = []
        for r in results:
            out.append({
                "chrom": r.variant.chrom, "pos": r.variant.pos,
                "ref": r.variant.ref, "alt": r.variant.alt,
                "AC": r.AC, "AN": r.AN,
                "AF": r.AF, "n_eligible": r.n_samples_eligible,
            })
        click.echo(json.dumps(out, indent=2))
    elif fmt == "tsv":
        click.echo("chrom\tpos\tref\talt\tAC\tAN\tAF\tn_eligible")
        for r in results:
            af = f"{r.AF:.6f}" if r.AF is not None else "NA"
            click.echo(
                f"{r.variant.chrom}\t{r.variant.pos}\t{r.variant.ref}\t"
                f"{r.variant.alt}\t{r.AC}\t{r.AN}\t{af}\t{r.n_samples_eligible}"
            )
    else:
        if not results:
            click.echo("No variants found for the given filters.")
            return
        for r in results:
            af = f"{r.AF:.4f}" if r.AF is not None else "NA"
            click.echo(
                f"{r.variant.chrom}:{r.variant.pos} {r.variant.ref}>{r.variant.alt}  "
                f"AC={r.AC}  AN={r.AN}  AF={af}  n_eligible={r.n_samples_eligible}"
            )


@cli.command()
@click.option("--db",      required=True, help="Path to database directory.")
@click.option("--input",   "input_vcf",  required=True, help="Input VCF (plain or .gz).")
@click.option("--output",  "output_vcf", required=True, help="Output annotated VCF.")
@click.option("--phenotype",   multiple=True, help="Phenotype code(s). Repeatable. If omitted, all available phenotypes are used.")
@click.option("--sex",     default="both", type=click.Choice(["male", "female", "both"]))
@click.option("--threads", default=None, type=int,
              help="Worker processes for parallel annotation. Default: all CPU cores.")
def annotate(db, input_vcf, output_vcf, phenotype, sex, threads):
    """Annotate a VCF with AFQUERY_AC / AFQUERY_AN / AFQUERY_AF INFO fields."""
    database = Database(db)
    # If no phenotype specified, use all available phenotypes
    phenotype_codes = list(phenotype) if phenotype else database.get_all_phenotypes()
    stats = database.annotate_vcf(input_vcf, output_vcf, phenotype_codes, sex, n_workers=threads)
    click.echo(
        f"Annotated {stats['n_annotated']} variants "
        f"({stats['n_uncovered']} uncovered, {stats['n_variants']} total)."
    )


@cli.command()
@click.option("--db", required=True, help="Path to database directory.")
def info(db):
    """Show database metadata."""
    database = Database(db)
    data = database.info()
    for k, v in data.items():
        click.echo(f"{k}: {v}")


@cli.command()
@click.option("--manifest",     required=True, help="Path to TSV manifest file.")
@click.option("--output-dir",   required=True, help="Path to output database directory.")
@click.option("--genome-build", required=True, type=click.Choice(["GRCh37", "GRCh38"]))
@click.option("--threads",      default=8, show_default=True, type=int)
@click.option("--tmp-dir",      default=None, help="Temporary directory for Arrow IPC files.")
@click.option("--bed-dir",      default=None, help="Directory containing BED files for WES technologies.")
def preprocess(manifest, output_dir, genome_build, threads, tmp_dir, bed_dir):
    """Preprocess VCFs from manifest into the query database."""
    from .preprocess import run_preprocess
    from .preprocess.manifest import ManifestError
    from .preprocess.ingest import IngestError
    try:
        run_preprocess(
            manifest_path=manifest,
            output_dir=output_dir,
            genome_build=genome_build,
            bed_dir=bed_dir,
            threads=threads,
            tmp_dir=tmp_dir,
        )
        click.echo(f"Database written to {output_dir}")
    except (ManifestError, IngestError) as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)


@cli.command("add-samples")
@click.option("--db",           required=True, help="Path to database directory.")
@click.option("--manifest",     required=True, type=click.Path(exists=True), help="Path to TSV manifest of new samples.")
@click.option("--threads",      default=8, show_default=True, type=int)
@click.option("--tmp-dir",      default=None, help="Temporary directory for ingestion files.")
@click.option("--bed-dir",      default=None, help="Directory containing BED files for WES technologies.")
@click.option("--genome-build", default=None, type=click.Choice(["GRCh37", "GRCh38"]), help="Validate DB genome build matches.")
def add_samples_cmd(db, manifest, threads, tmp_dir, bed_dir, genome_build):
    """Add new samples to an existing database without full rebuild."""
    from .preprocess.update import add_samples, UpdateError
    from .preprocess.manifest import ManifestError
    from .preprocess.ingest import IngestError
    try:
        result = add_samples(
            db, manifest, threads=threads, tmp_dir=tmp_dir,
            bed_dir=bed_dir, genome_build=genome_build,
        )
        click.echo(
            f"Added {result['new_samples']} sample(s): "
            f"{result['new_variants']} new variants, "
            f"{result['updated_variants']} updated variants."
        )
    except (UpdateError, ManifestError, IngestError) as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)


@cli.command("remove-samples")
@click.option("--db",           required=True, help="Path to database directory.")
@click.option("--sample-names", required=True, help="Comma-separated list of sample names to remove.")
def remove_samples_cmd(db, sample_names):
    """Remove samples from the database by clearing their bitmap bits."""
    from .preprocess.update import remove_samples, UpdateError
    names = [n.strip() for n in sample_names.split(",") if n.strip()]
    try:
        result = remove_samples(db, names)
        click.echo(f"Removed {len(result['removed'])} sample(s): {', '.join(result['removed'])}")
    except UpdateError as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)


@cli.command("check")
@click.option("--db", required=True, help="Path to database directory.")
def check_cmd(db):
    """Validate database integrity."""
    from .preprocess.update import check_database
    results = check_database(db)
    has_error = False
    for r in results:
        if r.severity == "error":
            prefix = "[ERROR]"
            has_error = True
        elif r.severity == "warning":
            prefix = "[WARN ]"
        else:
            prefix = "[INFO ]"
        click.echo(f"{prefix} {r.message}")
    if has_error:
        raise SystemExit(1)


@cli.command()
@click.argument("db_path")
@click.option("--dry-run", is_flag=True, help="Report what would be done without modifying files.")
def compact(db_path, dry_run):
    """Remove dead bits from removed samples to reclaim disk space."""
    from .preprocess.compact import compact_database
    from pathlib import Path

    if dry_run:
        click.echo("Dry-run mode: no files will be modified.")
        # Count parquets and active samples for informational output
        import sqlite3
        db = Path(db_path)
        con = sqlite3.connect(str(db / "metadata.sqlite"))
        n_active = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        con.close()
        variants_dir = db / "variants"
        n_flat = len(list(variants_dir.glob("*.parquet")))
        n_bucket = sum(1 for d in variants_dir.iterdir() if d.is_dir()
                       for _ in d.glob("bucket_*.parquet"))
        click.echo(
            f"Would compact {n_flat + n_bucket} Parquet file(s) "
            f"against {n_active} active sample(s)."
        )
        return

    stats = compact_database(Path(db_path))
    click.echo(
        f"Compact complete: {stats['files_rewritten']} file(s) rewritten, "
        f"{stats['rows_removed']} empty row(s) removed, "
        f"{stats['rows_processed']} row(s) processed. "
        f"Size: {stats['size_before']} → {stats['size_after']} bytes."
    )


@cli.command()
@click.option("--n-samples",  default=1000, show_default=True, type=int, help="Number of synthetic samples.")
@click.option("--n-variants", default=10_000, show_default=True, type=int, help="Variants per chromosome.")
@click.option("--output",     default="benchmark_report.json", show_default=True, help="Output JSON report path.")
@click.option("--db-dir",     default=None, help="Use an existing DB instead of generating synthetic data.")
def benchmark(n_samples, n_variants, output, db_dir):
    """Run performance benchmark suite and output a timing report."""
    from .benchmark import run_benchmark, run_benchmark_with_synth
    from pathlib import Path

    if db_dir is not None:
        click.echo(f"Benchmarking existing DB at {db_dir} ...")
        results = run_benchmark(Path(db_dir))
        import json
        with open(output, "w") as f:
            json.dump(results, f, indent=2)
    else:
        click.echo(
            f"Generating synthetic DB ({n_samples} samples, {n_variants} variants/chrom) ..."
        )
        results = run_benchmark_with_synth(
            n_samples=n_samples,
            n_variants=n_variants,
            output_report=output,
        )

    if "error" in results:
        click.echo(f"Benchmark error: {results['error']}", err=True)
        raise SystemExit(1)

    click.echo(f"Report written to {output}")
    click.echo(f"  point cold:  {results['point_query_cold_ms']:.1f} ms")
    click.echo(f"  point warm:  {results['point_query_warm_ms']:.1f} ms")
    click.echo(f"  batch-100:   {results['batch_100_ms']:.1f} ms")
    click.echo(f"  batch-1000:  {results['batch_1000_ms']:.1f} ms")
    targets = results.get("targets", {})
    ok = all(targets.values())
    click.echo("Targets: " + ("ALL MET" if ok else "SOME MISSED"))
