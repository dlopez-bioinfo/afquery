import json
import logging
import sys

import click

from .database import Database


def _expand_tokens(values: tuple[str, ...]) -> list[str]:
    """Expand comma-separated tokens into a flat list."""
    result = []
    for v in values:
        result.extend(t.strip() for t in v.split(",") if t.strip())
    return result


def _configure_logging(verbose: bool) -> None:
    """Configure logging to stderr with optional DEBUG verbosity."""
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter("%(message)s"))
    root = logging.getLogger("afquery")
    root.setLevel(level)
    root.handlers.clear()
    root.addHandler(handler)
    root.propagate = False


@click.group()
def cli():
    """afquery: genomic allele frequency query engine."""


@cli.command()
@click.option("--db",    required=True, help="Path to database directory.")
@click.option("--chrom", required=True, help="Chromosome (e.g., chr1, chrX).")
@click.option("--pos",   required=True, type=int, help="1-based genomic position.")
@click.option("--phenotype", multiple=True, help="Phenotype to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--sex",   default="both", type=click.Choice(["male", "female", "both"]), help="Restrict to specified sex. Options: male, female, both. (default: both)")
@click.option("--ref",   default=None, help="Filter to specific reference allele.")
@click.option("--alt",   default=None, help="Filter to specific alternate allele.")
@click.option("--tech", multiple=True, help="Technology filter to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--format", "fmt", default="text", type=click.Choice(["text", "json", "tsv"]), help="Output format. Options: text, json, tsv. (default: text)")
def query(db, chrom, pos, phenotype, sex, ref, alt, tech, fmt):
    """Query allele frequency at a single position."""
    database = Database(db)
    results = database.query(
        chrom=chrom, pos=pos,
        phenotype=_expand_tokens(phenotype), sex=sex,
        ref=ref, alt=alt, tech=_expand_tokens(tech),
    )

    if fmt == "json":
        out = []
        for r in results:
            entry = {
                "chrom": r.variant.chrom, "pos": r.variant.pos,
                "ref": r.variant.ref, "alt": r.variant.alt,
                "AC": r.AC, "AN": r.AN, "AF": r.AF, "n_eligible": r.n_samples_eligible,
                "N_HET": r.N_HET, "N_HOM_ALT": r.N_HOM_ALT, "N_HOM_REF": r.N_HOM_REF,
            }
            if r.N_FAIL is not None:
                entry["N_FAIL"] = r.N_FAIL
            out.append(entry)
        click.echo(json.dumps(out, indent=2))
    elif fmt == "tsv":
        has_fail = any(r.N_FAIL is not None for r in results)
        header = "chrom\tpos\tref\talt\tAC\tAN\tAF\tn_eligible\tN_HET\tN_HOM_ALT\tN_HOM_REF"
        if has_fail:
            header += "\tN_FAIL"
        click.echo(header)
        for r in results:
            af = f"{r.AF:.6f}" if r.AF is not None else "NA"
            line = (
                f"{r.variant.chrom}\t{r.variant.pos}\t{r.variant.ref}\t"
                f"{r.variant.alt}\t{r.AC}\t{r.AN}\t{af}\t{r.n_samples_eligible}\t"
                f"{r.N_HET}\t{r.N_HOM_ALT}\t{r.N_HOM_REF}"
            )
            if has_fail:
                line += f"\t{r.N_FAIL if r.N_FAIL is not None else '?'}"
            click.echo(line)
    else:  # text
        if not results:
            click.echo("No variants found at this position for the given filters.")
            return
        for r in results:
            af = f"{r.AF:.4f}" if r.AF is not None else "NA"
            fail_str = f"  N_FAIL={r.N_FAIL}" if r.N_FAIL is not None else ""
            click.echo(
                f"{r.variant.chrom}:{r.variant.pos} {r.variant.ref}>{r.variant.alt}  "
                f"AC={r.AC}  AN={r.AN}  AF={af}  n_eligible={r.n_samples_eligible}  "
                f"N_HET={r.N_HET}  N_HOM_ALT={r.N_HOM_ALT}  N_HOM_REF={r.N_HOM_REF}{fail_str}"
            )


@cli.command("query-batch")
@click.option("--db",       required=True, help="Path to database directory.")
@click.option("--chrom",    required=True, help="Chromosome (e.g., chr1, chrX).")
@click.option("--variants", required=True, type=click.Path(exists=True), help="TSV file with columns: pos ref alt (no header, whitespace-separated).")
@click.option("--phenotype",    multiple=True, help="Phenotype code to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--sex",      default="both", type=click.Choice(["male", "female", "both"]), help="Restrict to specified sex. Options: male, female, both. (default: both)")
@click.option("--tech", multiple=True, help="Technology filter to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--format",   "fmt", default="text", type=click.Choice(["text", "json", "tsv"]), help="Output format. Options: text, json, tsv. (default: text)")
def query_batch(db, chrom, variants, phenotype, sex, tech, fmt):
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
    results = database.query_batch(
        chrom=chrom, variants=variant_list,
        phenotype=_expand_tokens(phenotype), sex=sex, tech=_expand_tokens(tech),
    )

    if fmt == "json":
        out = []
        for r in results:
            entry = {
                "chrom": r.variant.chrom, "pos": r.variant.pos,
                "ref": r.variant.ref, "alt": r.variant.alt,
                "AC": r.AC, "AN": r.AN, "AF": r.AF, "n_eligible": r.n_samples_eligible,
                "N_HET": r.N_HET, "N_HOM_ALT": r.N_HOM_ALT, "N_HOM_REF": r.N_HOM_REF,
            }
            if r.N_FAIL is not None:
                entry["N_FAIL"] = r.N_FAIL
            out.append(entry)
        click.echo(json.dumps(out, indent=2))
    elif fmt == "tsv":
        has_fail = any(r.N_FAIL is not None for r in results)
        header = "chrom\tpos\tref\talt\tAC\tAN\tAF\tn_eligible\tN_HET\tN_HOM_ALT\tN_HOM_REF"
        if has_fail:
            header += "\tN_FAIL"
        click.echo(header)
        for r in results:
            af = f"{r.AF:.6f}" if r.AF is not None else "NA"
            line = (
                f"{r.variant.chrom}\t{r.variant.pos}\t{r.variant.ref}\t"
                f"{r.variant.alt}\t{r.AC}\t{r.AN}\t{af}\t{r.n_samples_eligible}\t"
                f"{r.N_HET}\t{r.N_HOM_ALT}\t{r.N_HOM_REF}"
            )
            if has_fail:
                line += f"\t{r.N_FAIL if r.N_FAIL is not None else '?'}"
            click.echo(line)
    else:
        if not results:
            click.echo("No variants found for the given filters.")
            return
        for r in results:
            af = f"{r.AF:.4f}" if r.AF is not None else "NA"
            fail_str = f"  N_FAIL={r.N_FAIL}" if r.N_FAIL is not None else ""
            click.echo(
                f"{r.variant.chrom}:{r.variant.pos} {r.variant.ref}>{r.variant.alt}  "
                f"AC={r.AC}  AN={r.AN}  AF={af}  n_eligible={r.n_samples_eligible}  "
                f"N_HET={r.N_HET}  N_HOM_ALT={r.N_HOM_ALT}  N_HOM_REF={r.N_HOM_REF}{fail_str}"
            )


@cli.command("query-region")
@click.option("--db",    required=True, help="Path to database directory.")
@click.option("--chrom", required=True, help="Chromosome (e.g., chr1, chrX).")
@click.option("--start", required=True, type=int, help="1-based start position (inclusive).")
@click.option("--end",   required=True, type=int, help="1-based end position (inclusive).")
@click.option("--phenotype", multiple=True, help="Phenotype code to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--sex",   default="both", type=click.Choice(["male", "female", "both"]), help="Restrict to specified sex. Options: male, female, both. (default: both)")
@click.option("--tech",  multiple=True, help="Technology filter to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--format", "fmt", default="text", type=click.Choice(["text", "json", "tsv"]), help="Output format. Options: text, json, tsv. (default: text)")
def query_region(db, chrom, start, end, phenotype, sex, tech, fmt):
    """Query allele frequencies for all variants in a chromosomal region [start, end]."""
    database = Database(db)
    results = database.query_region(
        chrom=chrom, start=start, end=end,
        phenotype=_expand_tokens(phenotype), sex=sex, tech=_expand_tokens(tech),
    )

    if fmt == "json":
        out = []
        for r in results:
            entry = {
                "chrom": r.variant.chrom, "pos": r.variant.pos,
                "ref": r.variant.ref, "alt": r.variant.alt,
                "AC": r.AC, "AN": r.AN, "AF": r.AF, "n_eligible": r.n_samples_eligible,
                "N_HET": r.N_HET, "N_HOM_ALT": r.N_HOM_ALT, "N_HOM_REF": r.N_HOM_REF,
            }
            if r.N_FAIL is not None:
                entry["N_FAIL"] = r.N_FAIL
            out.append(entry)
        click.echo(json.dumps(out, indent=2))
    elif fmt == "tsv":
        has_fail = any(r.N_FAIL is not None for r in results)
        header = "chrom\tpos\tref\talt\tAC\tAN\tAF\tn_eligible\tN_HET\tN_HOM_ALT\tN_HOM_REF"
        if has_fail:
            header += "\tN_FAIL"
        click.echo(header)
        for r in results:
            af = f"{r.AF:.6f}" if r.AF is not None else "NA"
            line = (
                f"{r.variant.chrom}\t{r.variant.pos}\t{r.variant.ref}\t"
                f"{r.variant.alt}\t{r.AC}\t{r.AN}\t{af}\t{r.n_samples_eligible}\t"
                f"{r.N_HET}\t{r.N_HOM_ALT}\t{r.N_HOM_REF}"
            )
            if has_fail:
                line += f"\t{r.N_FAIL if r.N_FAIL is not None else '?'}"
            click.echo(line)
    else:
        if not results:
            click.echo("No variants found in this region for the given filters.")
            return
        for r in results:
            af = f"{r.AF:.4f}" if r.AF is not None else "NA"
            fail_str = f"  N_FAIL={r.N_FAIL}" if r.N_FAIL is not None else ""
            click.echo(
                f"{r.variant.chrom}:{r.variant.pos} {r.variant.ref}>{r.variant.alt}  "
                f"AC={r.AC}  AN={r.AN}  AF={af}  n_eligible={r.n_samples_eligible}  "
                f"N_HET={r.N_HET}  N_HOM_ALT={r.N_HOM_ALT}  N_HOM_REF={r.N_HOM_REF}{fail_str}"
            )


@cli.command()
@click.option("--db",      required=True, help="Path to database directory.")
@click.option("--input",   "input_vcf",  required=True, help="Input VCF file (plain or .gz).")
@click.option("--output",  "output_vcf", required=True, help="Output annotated VCF file.")
@click.option("--phenotype",   multiple=True, help="Phenotype code to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--sex",     default="both", type=click.Choice(["male", "female", "both"]), help="Restrict to specified sex. Options: male, female, both. (default: both)")
@click.option("--tech", multiple=True, help="Technology filter to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--threads", default=None, type=int,
              help="Number of worker threads for parallel annotation. (default: all available CPU cores)")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output with per-item progress. (default: false)")
def annotate(db, input_vcf, output_vcf, phenotype, sex, tech, threads, verbose):
    """Annotate a VCF with AFQUERY_AC / AFQUERY_AN / AFQUERY_AF INFO fields."""
    _configure_logging(verbose)
    database = Database(db)
    stats = database.annotate_vcf(
        input_vcf, output_vcf,
        phenotype=_expand_tokens(phenotype), sex=sex,
        tech=_expand_tokens(tech), n_workers=threads,
    )
    click.echo(
        f"Annotated {stats['n_annotated']} variants "
        f"({stats['n_uncovered']} uncovered, {stats['n_variants']} total)."
    )


@cli.command()
@click.option("--db",       required=True, help="Path to database directory.")
@click.option("--samples",  is_flag=True,  help="List all samples with metadata. (default: false)")
@click.option("--changelog", "show_changelog", is_flag=True,
              help="Show full changelog history. (default: false)")
@click.option("--format", "fmt", default="table",
              type=click.Choice(["table", "tsv", "json"]),
              help="Output format. Options: table, tsv, json. (default: table)")
def info(db, samples, show_changelog, fmt):
    """Show database metadata and introspection."""
    database = Database(db)

    if samples:
        sample_list = database.list_samples()
        if fmt == "json":
            click.echo(json.dumps(sample_list, indent=2))
        elif fmt == "tsv":
            click.echo("sample_id\tsample_name\tsex\ttech\tphenotypes\tvcf_path\tingested_at")
            for s in sample_list:
                click.echo(
                    f"{s['sample_id']}\t{s['sample_name']}\t{s['sex']}\t{s['tech']}\t"
                    f"{','.join(s['phenotypes'])}\t{s['vcf_path'] or ''}\t{s['ingested_at'] or ''}"
                )
        else:
            click.echo(f"{'ID':<6} {'Name':<30} {'Sex':<8} {'Tech':<15} Phenotypes")
            click.echo("-" * 80)
            for s in sample_list:
                click.echo(
                    f"{s['sample_id']:<6} {s['sample_name']:<30} {s['sex']:<8} "
                    f"{s['tech']:<15} {','.join(s['phenotypes'])}"
                )
        return

    if show_changelog:
        cl = database.changelog()
        if fmt == "json":
            click.echo(json.dumps(cl, indent=2))
        elif fmt == "tsv":
            click.echo("event_id\tevent_type\tevent_time\tsample_names\tnotes")
            for e in cl:
                names = ','.join(e['sample_names']) if e['sample_names'] else ''
                click.echo(
                    f"{e['event_id']}\t{e['event_type']}\t{e['event_time']}\t"
                    f"{names}\t{e['notes'] or ''}"
                )
        else:
            for e in cl:
                n_str = f" ({len(e['sample_names'])} samples)" if e['sample_names'] else ""
                click.echo(f"  [{e['event_time']}] {e['event_type']:<20} — {e['notes'] or ''}{n_str}")
        return

    data = database.info()
    if fmt == "json":
        click.echo(json.dumps(data, indent=2))
    elif fmt == "tsv":
        for k, v in data.items():
            click.echo(f"{k}\t{v}")
    else:
        click.echo(f"Database:     {data['db_path']}")
        click.echo(f"Version:      {data['db_version']}")
        click.echo(f"Genome build: {data['genome_build']}   Schema: {data['schema_version']}")
        click.echo(f"Created: {data.get('created_at') or 'N/A'}    Updated: {data.get('updated_at') or 'N/A'}")
        click.echo()
        click.echo(f"Samples: {data['sample_count']} total")
        if data.get("by_sex"):
            sex_str = "  ".join(f"{k}={v}" for k, v in sorted(data["by_sex"].items()))
            click.echo(f"  By sex:       {sex_str}")
        if data.get("by_tech"):
            tech_str = "  ".join(f"{k}={v}" for k, v in sorted(data["by_tech"].items()))
            click.echo(f"  By tech:      {tech_str}")
        if data.get("by_phenotype"):
            pheno_str = "  ".join(f"{k}={v}" for k, v in sorted(data["by_phenotype"].items()))
            click.echo(f"  By phenotype: {pheno_str}")
        recent = data.get("changelog_recent", [])
        if recent:
            click.echo()
            click.echo("Recent changes:")
            for e in recent:
                click.echo(f"  [{e['event_time']}] {e['event_type']:<20} — {e['notes'] or ''}")


@cli.group()
def version():
    """Manage the database version label."""


@version.command("show")
@click.option("--db", required=True, help="Path to database directory.")
def version_show(db):
    """Show the current database version."""
    database = Database(db)
    click.echo(database.info().get("db_version", "unknown"))


@version.command("set")
@click.option("--db", required=True, help="Path to database directory.")
@click.argument("new_version")
def version_set(db, new_version):
    """Set the database version label."""
    database = Database(db)
    database.set_db_version(new_version)
    click.echo(f"Database version set to: {new_version}")


@cli.command()
@click.option("--manifest",     required=True, help="Path to TSV manifest file.")
@click.option("--output-dir",   required=True, help="Path to output database directory.")
@click.option("--genome-build", required=True, type=click.Choice(["GRCh37", "GRCh38"]), help="Reference genome build. Options: GRCh37, GRCh38.")
@click.option("--threads",      default=None, type=int, help="Number of worker threads for parallel processing. (default: all available CPU cores)")
@click.option("--build-threads", default=None, type=int, help="Max parallel workers for the build phase. (default: min(cpu_count, n_buckets))")
@click.option("--build-memory", default="2GB", help="DuckDB memory limit per build worker. Increase for WGS or dense regions. (default: 2GB)")
@click.option("--tmp-dir",      default=None, help="Temporary directory for intermediate files. (default: {output_dir}/.tmp_preprocess)")
@click.option("--bed-dir",      default=None, help="Directory containing BED files for WES technologies.")
@click.option("--force", is_flag=True, default=False, help="Delete any partial results and restart from scratch. (default: False)")
@click.option("--db-version", "db_version", default="1.0", help="Version label for this database. (default: 1.0)")
@click.option("--include-all-filters", "include_all_filters", is_flag=True, default=False,
              help="Count all variants regardless of FILTER field (default: PASS-only).")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output with per-item progress. (default: false)")
def preprocess(manifest, output_dir, genome_build, threads: int | None, build_threads: int | None, build_memory: str, tmp_dir, bed_dir, force, db_version, include_all_filters, verbose):
    """Preprocess VCFs from manifest into the query database."""
    _configure_logging(verbose)
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
            build_threads=build_threads,
            build_memory=build_memory,
            tmp_dir=tmp_dir,
            force=force,
            db_version=db_version,
            include_all_filters=include_all_filters,
        )
        click.echo(f"Database written to {output_dir}")
    except (ManifestError, IngestError) as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)


@cli.command("add-samples")
@click.option("--db",           required=True, help="Path to database directory.")
@click.option("--manifest",     required=True, type=click.Path(exists=True), help="Path to TSV manifest of new samples.")
@click.option("--threads",      default=None, type=int, help="Number of worker threads for parallel processing. (default: all available CPU cores)")
@click.option("--tmp-dir",      default=None, help="Temporary directory for intermediate files. (default: system temp)")
@click.option("--bed-dir",      default=None, help="Directory containing BED files for WES technologies.")
@click.option("--genome-build", default=None, type=click.Choice(["GRCh37", "GRCh38"]), help="Validate database genome build matches. Options: GRCh37, GRCh38.")
@click.option("--db-version", "db_version", default=None, help="New version label after adding samples. (default: auto-increment current version)")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output with per-item progress. (default: false)")
def add_samples_cmd(db, manifest, threads: int | None, tmp_dir, bed_dir, genome_build, db_version, verbose):
    """Add new samples to an existing database without full rebuild."""
    _configure_logging(verbose)
    from .preprocess.update import add_samples, UpdateError
    from .preprocess.manifest import ManifestError
    from .preprocess.ingest import IngestError
    try:
        result = add_samples(
            db, manifest, threads=threads, tmp_dir=tmp_dir,
            bed_dir=bed_dir, genome_build=genome_build, db_version=db_version,
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
@click.option("--sample-names", required=True, multiple=True,
              help="Sample name to remove. Repeatable; comma-separated or multiple flags.")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output with per-item progress. (default: false)")
def remove_samples_cmd(db, sample_names, verbose):
    """Remove samples from the database by clearing their bitmap bits."""
    _configure_logging(verbose)
    from .preprocess.update import remove_samples, UpdateError
    names = _expand_tokens(sample_names)
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
@click.option("--db", "db_path", required=True, type=click.Path(exists=True),
              help="Path to database directory.")
@click.option("--dry-run", is_flag=True, help="Dry-run mode: report changes without modifying files. (default: false)")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output with per-item progress. (default: false)")
def compact(db_path, dry_run, verbose):
    """Remove dead bits from removed samples to reclaim disk space."""
    _configure_logging(verbose)
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
@click.option("--n-samples",  default=1000, type=int, help="Number of synthetic samples to generate. (default: 1000)")
@click.option("--n-variants", default=10_000, type=int, help="Number of variants per chromosome. (default: 10000)")
@click.option("--output",     default="benchmark_report.json", help="Output path for JSON benchmark report. (default: benchmark_report.json)")
@click.option("--db-dir",     default=None, help="Use existing database instead of generating synthetic data.")
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
