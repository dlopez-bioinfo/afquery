#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Load HPC modules
module load BCFtools/1.18-GCC-12.3.0
module load BEDTools/2.31.0-GCC-12.3.0
module load htslib/1.18-GCC-12.3.0
module load parallel

PYTHON="/mnt/lustre/home/dlopez/micromamba/envs/afquery/bin/python"

echo "=========================================="
echo "Capture Kit Mixing Benchmark Pipeline"
echo "=========================================="

echo ""
echo "[Step 1] Prepare samples (split, mask, manifests)"
${PYTHON} 01_prepare_samples.py

echo ""
echo "[Step 2] Build AFQuery databases (1 WGS + 3 WES)"
${PYTHON} 02_build_databases.py

echo ""
echo "[Step 3] Compute AF metrics (NADI, AF error)"
${PYTHON} 03_compute_metrics.py

echo ""
echo "[Step 4] ACMG classification analysis"
${PYTHON} 04_classify_acmg.py

echo ""
echo "[Step 5] Generate figures"
${PYTHON} 05_plot_figures.py

echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "Figures:  ${SCRIPT_DIR}/../figures/"
echo "Results:  /mnt/lustre/home/dlopez/projects/afquery_bench_data/capture_kit/results/"
echo "=========================================="
