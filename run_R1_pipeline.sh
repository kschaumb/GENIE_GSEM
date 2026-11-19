#!/usr/bin/env bash
# =====================================================================================
# GENIE Genomic SEM — full R1 (JCPP revision) pipeline
#
# Runs every step behind the revised manuscript, from GWAS settings to supplemental
# tables and figures. All R1 steps read/write *_R1 locations; the original (submission)
# scripts and outputs are left in place for comparison.
#
#   01 settings corrections -> data/gwas_sumstats/GWAS_sumstats_settings_R1.xlsx
#   02 settings             -> gsem/gsem_settings_R1/
#   03 munge (ANX only)     -> data/gwas_sumstats/munged/hm3_munged_R1/
#   04 LDSC                 -> gsem/mvLD_R1/
#   05 genetic correlations -> gsem/gsem_output/RGs_R1/            (Figure 1)
#   06 EFA (geomin)         -> gsem/gsem_output/*_efa_R1/          (Tables S3, S4, S6, S7)
#   07-08 CFA               -> gsem/gsem_output/*_cfa_R1/          (Tables S2, S5, S8)
#   09 export CFA params    -> gsem/gsem_output/*_cfa_R1/params_R1/
#   10-11 Aim 3             -> gsem/gsem_output/aim3_conditional_R1/ (Figure 4, Table S9)
#   12 path diagrams        -> gsem/gsem_output/*_cfa_R1/diagrams/ (drafts for Figures 2-3)
#   13 supplemental tables  -> gsem/gsem_output/Supplemental_Tables_R1.xlsx
#   14 supplemental figures -> gsem/gsem_output/Supplemental_Figures.R1.pdf (needs Quarto)
#
# Requires the GWAS summary statistics and LD scores under data/ (see README.md).
# Usage (from anywhere):  bash run_R1_pipeline.sh
# =====================================================================================
set -euo pipefail
cd "$(dirname "$0")"
LOG=gsem/gsem_output/log_files; mkdir -p "$LOG"
run() { local tag=$1; shift; echo; echo "=== $tag ==="; Rscript -e "suppressPackageStartupMessages(library(dplyr)); $*" 2>&1 | tee "$LOG/R1_${tag}.log"; }

run 01_settings_corrections 'source("scripts/0.Settings_corrections.R1.R")'
run 02_settings             'source("scripts/0.Settings.R1.R")'
run 03_munge_ANX            'source("scripts/1.Munge.R1.R")'
run 04_ldsc                 'source("scripts/2.LDSC.R1.R")'
run 05_rgs                  'source("scripts/4.RGs.R1.R")'
run 06_efa                  'source("scripts/5.EFA.R1.r")'
run 07_cfa_disorders        'source("scripts/6.CFA_disorders.R1.R")'
run 08_cfa_symptoms         'source("scripts/6.CFA_symtpoms.R1.R")'
run 09_export_params '
for (d in c("gsem/gsem_output/disorder_cfa_R1","gsem/gsem_output/symptoms_cfa_R1")) {
  obj <- readRDS(file.path(d, "cfa_DWLS_results.RDS")); od <- file.path(d, "params_R1")
  dir.create(od, showWarnings = FALSE)
  for (f in names(obj$fits)) readr::write_csv(obj$fits[[f]]$results, file.path(od, paste0(f, "_params.csv")))
}'
run 10_aim3_trios           'source("scripts/8.Conditional_Overlap.R1.r")'
run 11_aim3_gmr             'source("scripts/9.GeneticMR.R1.r")'
run 12_diagrams             'source("scripts/7.CFA_Diagrams.R1.r")'
run 13_supp_tables          'source("gsem/gsem_output/Supplemental_Tables.R1.R")'

echo; echo "=== 14_supp_figures ==="
if command -v quarto >/dev/null 2>&1; then
  (cd gsem/gsem_output && quarto render Supplemental_Figures.R1.qmd --to pdf) 2>&1 | tee "$LOG/R1_14_supp_figures.log"
else
  echo "Quarto not found; skipping Supplemental_Figures.R1.pdf"
fi
echo; echo "Done."
