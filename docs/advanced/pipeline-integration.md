# Pipeline Integration

AFQuery integrates into automated pipelines as a command-line annotation step. The database is a read-only file directory, safe for concurrent access from multiple processes.

---

## Nextflow

### Basic Annotation Process

```groovy
process AFQUERY_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(afquery_db)

    output:
    tuple val(meta), path("*.annotated.vcf.gz"), path("*.annotated.vcf.gz.tbi"), emit: vcf

    script:
    """
    afquery annotate \
        --db ${afquery_db} \
        --input ${vcf} \
        --output ${meta.id}.annotated.vcf.gz \
        --threads ${task.cpus}

    tabix -p vcf ${meta.id}.annotated.vcf.gz
    """
}
```

### Usage in a Workflow

```groovy
workflow VARIANT_ANNOTATION {
    take:
    ch_vcfs       // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_afquery_db // channel: path(db_dir)

    main:
    AFQUERY_ANNOTATE(ch_vcfs, ch_afquery_db.collect())

    emit:
    vcf = AFQUERY_ANNOTATE.out.vcf
}
```

### nf-core Integration Pattern

For Sarek-like pipelines, add AFQuery as a post-annotation step after VEP/SnpEff:

```groovy
// In your main workflow, after VEP annotation:
AFQUERY_ANNOTATE(
    VEP.out.vcf,
    Channel.fromPath(params.afquery_db)
)
```

Add to `nextflow.config`:

```groovy
params {
    afquery_db = null  // Path to AFQuery database directory
}
```

---

## Snakemake

### Basic Rule

```python
rule afquery_annotate:
    input:
        vcf="results/variants/{sample}.vcf.gz",
        db=config["afquery_db"],
    output:
        vcf="results/annotated/{sample}.annotated.vcf.gz",
    threads: 8
    shell:
        """
        afquery annotate \
            --db {input.db} \
            --input {input.vcf} \
            --output {output.vcf} \
            --threads {threads}
        """
```

### Per-Chromosome Parallelism

For large VCFs, split by chromosome and annotate in parallel:

```python
rule afquery_annotate_chrom:
    input:
        vcf="results/split/{sample}.{chrom}.vcf.gz",
        db=config["afquery_db"],
    output:
        vcf="results/annotated/{sample}.{chrom}.annotated.vcf.gz",
    threads: 4
    shell:
        """
        afquery annotate \
            --db {input.db} \
            --input {input.vcf} \
            --output {output.vcf} \
            --threads {threads}
        """

rule merge_annotated:
    input:
        vcfs=expand(
            "results/annotated/{{sample}}.{chrom}.annotated.vcf.gz",
            chrom=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"],
        ),
    output:
        vcf="results/annotated/{sample}.annotated.vcf.gz",
    shell:
        "bcftools concat {input.vcfs} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}"
```

---

## Parallelism Notes

| Level | Mechanism | Notes |
|-------|-----------|-------|
| Per-file | `--threads N` in `afquery annotate` | Controls internal read-ahead and query parallelism |
| Per-sample | Pipeline-level (Nextflow channels, Snakemake wildcards) | Each sample runs as an independent process |
| Per-chromosome | Split VCF → annotate → merge | Useful for very large single-sample VCFs |

### Concurrent Read Access

The AFQuery database directory is **read-only during queries and annotation**. Multiple processes can read the same database simultaneously without locking or corruption. This means:

- Multiple Nextflow tasks can share one database via `collect()`
- Snakemake rules can reference the same `config["afquery_db"]` path
- No need for database copies or per-worker instances

!!! warning "Do not run `create-db`, `update`, or `compact` while queries are in progress"
    Write operations (database creation and updates) modify files in the database directory. Run these as dedicated pipeline steps with no concurrent readers.

---

## Shared Database as a Pipeline Resource

A common pattern is to maintain one AFQuery database per institution and reference it from multiple pipelines:

```
/shared/databases/afquery/
  cohort_v2/              ← current production database
    manifest.json
    metadata.sqlite
    variants/
    capture/
  cohort_v1/              ← previous version (kept for reproducibility)
```

Point all pipelines to the current version:

```groovy
// nextflow.config
params.afquery_db = "/shared/databases/afquery/cohort_v2"
```

```yaml
# Snakemake config.yaml
afquery_db: /shared/databases/afquery/cohort_v2
```

---

## Next Steps

- [Annotate a VCF](../guides/annotate-vcf.md) — full annotation CLI reference
- [Performance Tuning](performance.md) — thread and memory optimization
- [Multi-cohort Strategies](multi-cohort-strategies.md) — managing multiple databases
