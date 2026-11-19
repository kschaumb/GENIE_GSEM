# =====================================================================
# 6.CFA_export_params.R1.R  (R&R revision, Nov 2026)
# Read-only helper: exports the full parameter tables for every fitted
# CFA (original/constrained x residual floor) from the saved RDS files,
# so loadings, residual variances and factor correlations can be
# reported/checked without refitting. Does NOT refit any model.
# Run from the project root (GENIE_GSEM/).
# =====================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr) })

panels <- c(disorders = "gsem/gsem_output/disorder_cfa",
            symptoms  = "gsem/gsem_output/symptoms_cfa")

for (p in names(panels)) {
  in_rds  <- file.path(panels[[p]], "cfa_DWLS_results.RDS")
  out_dir <- file.path(panels[[p]], "params_R1")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  obj <- readRDS(in_rds)
  for (fit_name in names(obj$fits)) {
    res <- obj$fits[[fit_name]]$results
    if (!is.data.frame(res)) { message("No $results for ", p, "/", fit_name); next }
    res <- res %>% mutate(panel = p, fit = fit_name, .before = 1)
    write_csv(res, file.path(out_dir, paste0(fit_name, "_params.csv")))
    message("Wrote ", file.path(out_dir, paste0(fit_name, "_params.csv")))
  }
  # also save the model syntax actually fitted
  writeLines(c("## original", obj$models$original_clean, "",
               "## constrained", obj$models$constrained_clean),
             file.path(out_dir, "model_syntax.txt"))
}
