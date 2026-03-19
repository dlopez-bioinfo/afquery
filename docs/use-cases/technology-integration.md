NOTA_PARA_LLM: Reescribir completamente esta sección para destacar que AFQuery permite crear bases de datos locales a partir de muestras secuenciaciadas con diferentes tecnologías, calculando correctamente la AN según la zona secuenciada y cubierta por el bed de captura de cada muestra. Esto, hacerlo a manualmente es muy complicado, especialmente si hay un gran número de paneles. De hecho, es muy común aún que los hospitales hagan paneles de genes para los distintos grupos de enfermedad en lugar de hacer WES o WGS. Este problema también se tiene cuando se usan distintas versiones de un mismo kit de captura, que tienen pequeñas diferencias. Si no se tiene en cuenta la versión del kit de captura para el cálculo de AN, estos valores estarán sesgados. Calcular las AF teniendo en cuenta las regiones secuenciadas de cada kit de captura es un proceso muy laborioso, y no existen herramientas bioinformáticas que lo hagan de forma automática. Esta es una de las features más importantes de esta herramienta. Reescribe completamente esta sección para plasmar esta idea como una guía más.

# Technology Integration (WES + WGS)

## Scenario

Your cohort contains samples sequenced with whole-genome sequencing (WGS), two different whole-exome kits (WES_v1 and WES_v2), and a targeted gene panel. Each technology covers a different set of genomic positions. You want to compute AF at a position covered only by WGS and WES_v1 — without contaminating the AN estimate with panel samples that do not cover that position.

## Why Standard Databases Fall Short

When mixing sequencing technologies, a variant in a gene not covered by a WES kit produces GT=missing for those samples. Without technology-aware filtering, these missing genotypes are either:
- Excluded from the VCF entirely (most callers) → AN appears correct but WES samples are implicitly treated as non-carriers
- Included as no-call → inflates AN if counted

Neither approach is principled. AFQuery solves this by associating each sample with its technology's BED file and excluding WES samples from AN computation at positions outside their capture region.

## How AFQuery Handles Technology Coverage

When you build the database:

```bash
afquery create-db \
  --manifest manifest.tsv \
  --bed-dir ./beds/ \          # Directory with <tech_name>.bed files
  --output-dir ./db/ \
  --genome-build GRCh38
```

For each non-WGS technology, AFQuery builds an interval tree from the BED file. At query time, a WES sample is eligible at a position only if that position falls within its technology's capture regions.

## Step-by-Step Example

NOTA_PARA_LLM: La sección Step-by-Step Example está repetida en múltiples sitios en la documentación, incluyendo todas las guías de clinical workflows. Piensa una manera de optimizar la documentación para que no sea tan redundante y facilite la lectura de toda la documentación.

### 1. Manifest with mixed technologies

```tsv
sample_name	vcf_path	sex	tech_name	phenotype_codes
SAMP_001	vcfs/SAMP_001.vcf.gz	female	wgs	control
SAMP_002	vcfs/SAMP_002.vcf.gz	male	wgs	control
SAMP_003	vcfs/SAMP_003.vcf.gz	female	wes_v1	rare_disease
SAMP_004	vcfs/SAMP_004.vcf.gz	male	wes_v1	rare_disease
SAMP_005	vcfs/SAMP_005.vcf.gz	female	wes_v2	rare_disease
SAMP_006	vcfs/SAMP_006.vcf.gz	male	panel_cardio	control
```

### 2. BED files in `./beds/`

```
beds/
├── wes_v1.bed          # Coverage for WES kit v1 samples
├── wes_v2.bed          # Coverage for WES kit v2 samples
└── panel_cardio.bed    # Coverage for cardiac gene panel
```

No BED file is needed for `wgs` — WGS samples are always eligible.

### 3. Query at a position covered by WES_v1 but not by panel

```bash
afquery query --db ./db/ --locus chr1:925952
```

```
chr1:925952
  REF=G  ALT=A  AC=8  AN=200  AF=0.0400
```

AN=200 reflects only samples whose technology covers this position (WGS + WES_v1 — not WES_v2 if not in capture, and not panel_cardio if not in capture). The panel samples contribute AN=0 at this position automatically.

### 4. Query restricted to one technology

```bash
# Only WES_v1 samples
afquery query --db ./db/ --locus chr1:925952 --tech wes_v1

# Only WGS samples
afquery query --db ./db/ --locus chr1:925952 --tech wgs
```

### 5. Detect technology artifacts

Compare AF across technologies for the same variant:

```bash
for tech in wgs wes_v1 wes_v2; do
  echo -n "$tech: "
  afquery query --db ./db/ --locus chr1:925952 --tech $tech --format tsv | \
    awk 'NR>1 {print "AF=" $6 " AN=" $5}'
done
```

If AF differs substantially between WGS and WES, this may indicate:
- A capture artifact (systematic over- or under-calling in WES)
- A WES coverage boundary effect (positions near capture edges)
- True biological difference (different disease distributions between WGS/WES groups)

### 6. Position not covered by any technology

```bash
afquery query --db ./db/ --locus chr1:999999999
```

```
chr1:999999999 — no results (AN=0 for all variants)
```

AN=0 means no eligible samples at this position. For WES-only cohorts, most positions return AN=0 (only captured exonic regions contribute).

## Biological Interpretation

Technology-aware AN ensures that reported AF reflects only samples for which genotype calling at that position was attempted and expected to succeed. This avoids the misleading situation where a variant at an uncaptured position appears at AF=0 (rather than AF=unknown/not-covered).

**Best practices:**
- Always provide BED files for all non-WGS technologies
- Verify with `afquery check` that all BED files were registered correctly
- Use `--tech wgs` for gene-agnostic cohort-wide analyses
- Report `AN` alongside `AF` — AN < 100 indicates insufficient coverage to interpret AF

## Related Features

- [Manifest Format](../guides/manifest-format.md) — tech_name and BED file setup
- [Create a Database](../guides/create-database.md) — `--bed-dir` option
- [Sample Filtering](../guides/sample-filtering.md) — `--tech` filter syntax
- [Cohort Stratification](cohort-stratification.md) — compare WGS vs. WES AF
