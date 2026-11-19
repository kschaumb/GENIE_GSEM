#!/usr/bin/env bash
# R&R step 06: rerun disorder- and symptom-level EFAs with geomin rotation (preregistered) and
# communality-based variance explained (scripts/5.EFA.R1.r).
# Outputs: gsem/gsem_output/{disorder_efa_R1,symptoms_efa_R1}/   (original EFA outputs untouched)
# Requires: GPArotation  (install.packages("GPArotation") if missing)
# Usage:   bash R1_06_efa_geomin.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM"
LOG=gsem/gsem_output/log_files
Rscript -e 'if (!requireNamespace("GPArotation", quietly = TRUE)) install.packages("GPArotation", repos = "https://cloud.r-project.org")'
Rscript scripts/5.EFA.R1.r 2>&1 | tee "$LOG/R1_06_efa_geomin.log"
echo "Done. Check gsem/gsem_output/{disorder_efa_R1,symptoms_efa_R1}/"
