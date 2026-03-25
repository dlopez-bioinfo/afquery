#!/usr/bin/env bash
# Download 1000 Genomes Phase 3 chr22 and split into per-sample VCFs.
#
# Prerequisites: bcftools, wget/curl, tabix
# Output: per-sample VCFs + manifest.tsv in ONEKG_DIR (see config.py)
#
# Runtime: ~30-60 min depending on network and CPU (splitting 2504 samples).

set -euo pipefail

# Load required modules from HPC
module load BCFtools/1.18-GCC-12.3.0
module load parallel

# ---------------------------------------------------------------------------
# Configuration — must match config.py
# ---------------------------------------------------------------------------
DATA_DIR="/mnt/lustre/home/dlopez/projects/afquery_bench_data"
ONEKG_DIR="${DATA_DIR}/1kg"
VCF_DIR="${ONEKG_DIR}/vcfs"
THREADS=16

FTP_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
CHR22_VCF="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL_FILE="integrated_call_samples_v3.20130502.ALL.panel"

# ---------------------------------------------------------------------------
# Create directories
# ---------------------------------------------------------------------------
mkdir -p "${ONEKG_DIR}" "${VCF_DIR}"

# ---------------------------------------------------------------------------
# Step 1: Download chr22 multi-sample VCF + index + panel
# ---------------------------------------------------------------------------
echo "[1/4] Downloading 1KG chr22 VCF..."
if [ ! -f "${ONEKG_DIR}/${CHR22_VCF}" ]; then
    wget -q --show-progress -O "${ONEKG_DIR}/${CHR22_VCF}" \
        "${FTP_BASE}/${CHR22_VCF}"
    wget -q --show-progress -O "${ONEKG_DIR}/${CHR22_VCF}.tbi" \
        "${FTP_BASE}/${CHR22_VCF}.tbi"
else
    echo "  Already downloaded: ${CHR22_VCF}"
fi

echo "[2/4] Downloading panel (sample metadata)..."
if [ ! -f "${ONEKG_DIR}/${PANEL_FILE}" ]; then
    wget -q --show-progress -O "${ONEKG_DIR}/${PANEL_FILE}" \
        "${FTP_BASE}/${PANEL_FILE}"
else
    echo "  Already downloaded: ${PANEL_FILE}"
fi

# ---------------------------------------------------------------------------
# Step 2: Extract sample list
# ---------------------------------------------------------------------------
echo "[3/4] Splitting into per-sample VCFs (${THREADS} parallel jobs)..."

SAMPLE_LIST="${ONEKG_DIR}/sample_list.txt"
bcftools query -l "${ONEKG_DIR}/${CHR22_VCF}" > "${SAMPLE_LIST}"
N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo "  Found ${N_SAMPLES} samples"

# ---------------------------------------------------------------------------
# Step 3: Split multi-sample VCF into per-sample VCFs (parallel)
# ---------------------------------------------------------------------------
split_one_sample() {
    local sample="$1"
    local vcf_dir="$2"
    local merged_vcf="$3"
    local out_vcf="${vcf_dir}/${sample}.vcf.gz"

    if [ -f "${out_vcf}" ]; then
        return 0
    fi

    # Extract sample, remove hom-ref (0/0 and ./.) to keep files sparse,
    # then compress with bgzip
    bcftools view -s "${sample}" --min-ac 1:nref "${merged_vcf}" \
        -Oz -o "${out_vcf}" 2>/dev/null
}
export -f split_one_sample

# Use GNU parallel if available, otherwise xargs
if command -v parallel &>/dev/null; then
    parallel -j "${THREADS}" --bar \
        split_one_sample {} "${VCF_DIR}" "${ONEKG_DIR}/${CHR22_VCF}" \
        :::: "${SAMPLE_LIST}"
else
    echo "  (GNU parallel not found, falling back to xargs -P)"
    cat "${SAMPLE_LIST}" | xargs -P "${THREADS}" -I{} bash -c \
        'split_one_sample "$@"' _ {} "${VCF_DIR}" "${ONEKG_DIR}/${CHR22_VCF}"
fi

# ---------------------------------------------------------------------------
# Step 4: Generate manifest.tsv
# ---------------------------------------------------------------------------
echo "[4/4] Generating manifest.tsv..."

MANIFEST="${ONEKG_DIR}/manifest.tsv"
echo -e "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes" > "${MANIFEST}"

# Panel columns: sample pop super_pop gender
# gender: male/female
while IFS=$'\t' read -r sample pop super_pop gender; do
    # Skip header
    [ "${sample}" = "sample" ] && continue

    vcf_path="${VCF_DIR}/${sample}.vcf.gz"
    [ -f "${vcf_path}" ] || continue

    # Map 1KG gender to afquery sex
    sex="${gender}"

    # Use super_population as phenotype code
    echo -e "${sample}\t${sex}\tWGS\t${vcf_path}\t${super_pop}"
done < "${ONEKG_DIR}/${PANEL_FILE}" >> "${MANIFEST}"

MANIFEST_LINES=$(tail -n +2 "${MANIFEST}" | wc -l)
echo "Done. Manifest has ${MANIFEST_LINES} samples: ${MANIFEST}"
