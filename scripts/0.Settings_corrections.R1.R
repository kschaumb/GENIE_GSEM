# ============================================================
# R&R (R1): documented corrections to GWAS_sumstats_settings.xlsx
# Reads the original settings workbook, applies corrections IN CODE, and writes
#   data/gwas_sumstats/GWAS_sumstats_settings_R1.xlsx
# The original workbook is not modified.
#
# Corrections
#  ANX (Otowa et al., 2016, case-control analysis): counts had been taken from the
#     Table 1 total row (7,016 / 14,745) rather than the case-control meta-analysis that
#     produced the summary statistics used here (3,695 cases / 13,615 controls; N = 17,310).
#     gsem_n is recomputed as full effective N, 4 / (1/Ncase + 1/Ncontrol) ≈ 11,625,
#     matching the effective-N convention used for every other binary trait
#     (the previous value, 9,507.9, was half-Neff on the incorrect counts).
#  BE_noAN: gsem_samplePrevalence was missing; set to 0.5 because the trait is munged
#     with full effective N (Neff_half column, doubled by munge), as for AN and BEBROAD.
#  ANX_FS: type "binary" -> "Continuous" (factor-score GWAS; label only).
#  Dep_Sx: Year 2025 -> 2016 (Okbay et al., 2016; label only).
#  Neurotic: name "Neuroticisim" -> "Neuroticism" (label only).
# ============================================================
suppressPackageStartupMessages({ library(readxl); library(openxlsx); library(dplyr) })

in_xlsx  <- file.path("data", "gwas_sumstats", "GWAS_sumstats_settings.xlsx")
out_xlsx <- file.path("data", "gwas_sumstats", "GWAS_sumstats_settings_R1.xlsx")

s <- readxl::read_excel(in_xlsx)
set_val <- function(df, code, col, value) {
  i <- which(df$code == code); stopifnot(length(i) == 1, col %in% names(df))
  old <- df[[col]][i]
  if (is.character(df[[col]])) value <- as.character(value)
  df[[col]][i] <- value
  message(sprintf("  %-9s %-26s %s -> %s", code, col, format(old), format(value)))
  df
}

anx_case <- 3695; anx_ctrl <- 13615
anx_neff <- 4 / (1 / anx_case + 1 / anx_ctrl)

message("Applying R1 corrections:")
s <- set_val(s, "ANX", "n_case", anx_case)
s <- set_val(s, "ANX", "n_control", anx_ctrl)
s <- set_val(s, "ANX", "n_total", anx_case + anx_ctrl)
s <- set_val(s, "ANX", "gsem_n", round(anx_neff, 1))
s <- set_val(s, "ANX", "sumstats_samplePrevalence", round(anx_case / (anx_case + anx_ctrl), 4))
s <- set_val(s, "BE_noAN", "gsem_samplePrevalence", 0.5)
s <- set_val(s, "ANX_FS", "type", "Continuous")
s <- set_val(s, "Dep_Sx", "Year", 2016)
s <- set_val(s, "Neurotic", "name", "Neuroticism")

openxlsx::write.xlsx(s, out_xlsx, overwrite = TRUE)
message("Wrote ", out_xlsx)
