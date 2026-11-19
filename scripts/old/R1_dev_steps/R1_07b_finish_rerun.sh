#!/usr/bin/env bash
# R&R step 07b: finish the corrected-ANX-N rerun after Aim 3 (diagrams, supplemental tables,
# before/after comparison). Use if R1_07 stopped after step 11.
# Usage:   bash R1_07b_finish_rerun.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM"
LOG=gsem/gsem_output/log_files
run() { local tag=$1; shift; echo; echo "=== $tag ==="; Rscript -e "suppressPackageStartupMessages(library(dplyr)); $*" 2>&1 | tee "$LOG/R1_07_${tag}.log"; }
run 13_diagrams     'source("scripts/7.CFA_Diagrams.R1.r")'
run 14_supp_tables  'source("gsem/gsem_output/Supplemental_Tables.R1.R")'
run 15_compare      'source("scripts/R1_compare_rerun.R")'
echo; echo "Done."
