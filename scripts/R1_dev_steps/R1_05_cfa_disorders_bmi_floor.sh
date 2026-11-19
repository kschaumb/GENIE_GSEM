#!/usr/bin/env bash
# R&R step 05: refit the disorder-level CFA with the residual floor on BMI only
# (scripts/6.CFA_disorders.R1.R), export parameters, then rebuild diagrams + supplemental tables.
# Usage:   bash R1_05_cfa_disorders_bmi_floor.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM"
LOG=gsem/gsem_output/log_files

Rscript -e 'suppressPackageStartupMessages(library(dplyr)); source("scripts/6.CFA_disorders.R1.R")' 2>&1 | tee "$LOG/R1_05_cfa_disorders.log"

Rscript -e '
suppressPackageStartupMessages({library(dplyr); library(readr)})
d <- "gsem/gsem_output/disorder_cfa_R1"
obj <- readRDS(file.path(d, "cfa_DWLS_results.RDS")); od <- file.path(d, "params_R1")
dir.create(od, showWarnings = FALSE)
for (f in names(obj$fits)) write_csv(obj$fits[[f]]$results, file.path(od, paste0(f, "_params.csv")))
cat("GenomicSEM version:", as.character(packageVersion("GenomicSEM")), "\n")
' 2>&1 | tee "$LOG/R1_05_export_params.log"

Rscript scripts/7.CFA_Diagrams.R1.r 2>&1 | tee "$LOG/R1_05_diagrams.log"
Rscript gsem/gsem_output/Supplemental_Tables.R1.R 2>&1 | tee "$LOG/R1_05_tables.log"
echo "Done."
