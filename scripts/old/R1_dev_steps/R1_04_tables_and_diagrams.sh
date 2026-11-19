#!/usr/bin/env bash
# R&R step 04: rebuild CFA path diagrams and supplemental tables from the R1 CFA outputs,
# with the Child Physical Development factor reflected for reporting.
# Outputs: gsem/gsem_output/{disorder_cfa_R1,symptoms_cfa_R1}/diagrams/*.svg|png
#          gsem/gsem_output/Supplemental_Tables_R1.xlsx   (original workbook untouched)
# Usage:   bash R1_04_tables_and_diagrams.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM"
LOG=gsem/gsem_output/log_files

Rscript scripts/7.CFA_Diagrams.R1.r 2>&1 | tee "$LOG/R1_04_diagrams.log"
Rscript gsem/gsem_output/Supplemental_Tables.R1.R 2>&1 | tee "$LOG/R1_04_tables.log"

echo "Done."
