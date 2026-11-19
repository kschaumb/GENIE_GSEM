#!/usr/bin/env bash
# R&R step 07: full rerun with corrected settings (ANX effective N ≈ 11,625; BE_noAN sample
# prevalence 0.5; ST1 label fixes). All originals untouched; everything writes to *_R1 locations.
#   settings  -> data/gwas_sumstats/GWAS_sumstats_settings_R1.xlsx, gsem/gsem_settings_R1/
#   munge     -> data/gwas_sumstats/munged/hm3_munged_R1/ANX.sumstats.gz
#   LDSC      -> gsem/mvLD_R1/
#   downstream-> RGs_R1, *_efa_R1, *_cfa_R1, aim3_conditional_R1, Supplemental_Tables_R1.xlsx
#   check     -> gsem/gsem_output/R1_rerun_comparison/
# LDSC on the 22-trait set is the slow step. Raw ANX sumstats and LD scores are pulled from Box
# on first access, so keep Box Drive online.
# Usage:   bash R1_07_full_rerun_ANX_N.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM"
LOG=gsem/gsem_output/log_files
run() { local tag=$1; shift; echo; echo "=== $tag ==="; Rscript -e "suppressPackageStartupMessages(library(dplyr)); $*" 2>&1 | tee "$LOG/R1_07_${tag}.log"; }

# snapshot current R1 CFA estimates so the comparison has a "before"
SNAP=gsem/gsem_output/R1_rerun_comparison/pre_rerun
mkdir -p "$SNAP"
for k in disorder_cfa_R1 symptoms_cfa_R1; do
  [ -d "$SNAP/$k" ] || cp -R "gsem/gsem_output/$k/params_R1" "$SNAP/$k"
done

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
# (8.Sensitivity.R1.R removed: exploratory, never run in the original analysis, not reported)
run 13_diagrams             'source("scripts/7.CFA_Diagrams.R1.r")'
run 14_supp_tables          'source("gsem/gsem_output/Supplemental_Tables.R1.R")'
run 15_compare              'source("scripts/R1_compare_rerun.R")'
echo; echo "Done. Send me the R1_07_15_compare.log (or just say done)."
