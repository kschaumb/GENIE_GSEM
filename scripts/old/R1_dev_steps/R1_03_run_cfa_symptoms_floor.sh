#!/usr/bin/env bash
# R&R step 03: refit the symptom-level CFA with the BMI residual floor applied
# (guard_indicators = "BMI" in scripts/6.CFA_symtpoms.R1.R), then export parameter tables.
# Outputs overwrite gsem/gsem_output/symptoms_cfa_R1/ (R1 outputs only; originals untouched).
# Usage:   bash R1_03_run_cfa_symptoms_floor.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM"
LOG=gsem/gsem_output/log_files

Rscript -e 'suppressPackageStartupMessages(library(dplyr)); source("scripts/6.CFA_symtpoms.R1.R")' 2>&1 | tee "$LOG/R1_03_cfa_symptoms.log"

Rscript -e '
suppressPackageStartupMessages({library(dplyr); library(readr)})
d <- "gsem/gsem_output/symptoms_cfa_R1"
obj <- readRDS(file.path(d, "cfa_DWLS_results.RDS")); od <- file.path(d, "params_R1")
dir.create(od, showWarnings = FALSE)
for (f in names(obj$fits)) write_csv(obj$fits[[f]]$results, file.path(od, paste0(f, "_params.csv")))
' 2>&1 | tee "$LOG/R1_03_export_params.log"

echo "Done. Check gsem/gsem_output/symptoms_cfa_R1/"
