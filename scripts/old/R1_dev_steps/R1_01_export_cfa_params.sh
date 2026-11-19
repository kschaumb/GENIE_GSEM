#!/usr/bin/env bash
# R&R step 01: export CFA parameter tables from existing RDS (no refitting).
# Usage (from anywhere):  bash R1_01_export_cfa_params.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM"
Rscript scripts/6.CFA_export_params.R1.R 2>&1 | tee gsem/gsem_output/log_files/R1_01_export_cfa_params.log
echo "Done. Outputs in gsem/gsem_output/{disorder_cfa,symptoms_cfa}/params_R1/"
