#!/usr/bin/env bash
# R&R step 02: refit disorder- and symptom-level CFAs with the Internalizing ~~ Development
# correlation freely estimated (R1 scripts). Original scripts/outputs are untouched.
# Outputs: gsem/gsem_output/disorder_cfa_R1/ and gsem/gsem_output/symptoms_cfa_R1/
# Usage:   bash R1_02_run_cfa.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM"
LOG=gsem/gsem_output/log_files

# dplyr is loaded first because the scripts use %>% before their own library() calls
Rscript -e 'suppressPackageStartupMessages(library(dplyr)); source("scripts/6.CFA_disorders.R1.R")' 2>&1 | tee "$LOG/R1_02_cfa_disorders.log"
Rscript -e 'suppressPackageStartupMessages(library(dplyr)); source("scripts/6.CFA_symtpoms.R1.R")'  2>&1 | tee "$LOG/R1_02_cfa_symptoms.log"

# export parameter tables for the new fits (same helper as step 01, pointed at the _R1 folders)
Rscript -e '
suppressPackageStartupMessages({library(dplyr); library(readr)})
for (d in c("gsem/gsem_output/disorder_cfa_R1","gsem/gsem_output/symptoms_cfa_R1")) {
  obj <- readRDS(file.path(d, "cfa_DWLS_results.RDS")); od <- file.path(d, "params_R1")
  dir.create(od, showWarnings = FALSE)
  for (f in names(obj$fits)) write_csv(obj$fits[[f]]$results, file.path(od, paste0(f, "_params.csv")))
  writeLines(c("## original", obj$models$original_clean, "", "## constrained", obj$models$constrained_clean),
             file.path(od, "model_syntax.txt"))
}' 2>&1 | tee "$LOG/R1_02_export_params.log"

echo "Done. Check gsem/gsem_output/{disorder_cfa_R1,symptoms_cfa_R1}/"
