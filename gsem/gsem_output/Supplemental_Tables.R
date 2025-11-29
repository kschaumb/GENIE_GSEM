# ===========================
# Build Supplemental_Tables.xlsx (simple version)
# ===========================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(openxlsx); library(tibble); library(lavaan)
})

source('gsem/gsem_output/Supplemental_Table_Functions.R')
# ---- Paths ----
out_dir   <- file.path("gsem", "gsem_output")
out_xlsx  <- file.path(out_dir, "Supplemental_Tables.xlsx")
gwas_file <- file.path("data", "gwas_sumstats", "GWAS_sumstats_settings.xlsx")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



# ===========================================
# Build workbook
# ===========================================
wb <- createWorkbook()

# ST1
st1 <- openxlsx::read.xlsx(gwas_file) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>% arrange(across(1)) 
add_sheet_simple(wb, "ST1_GWAS_Sumstats_Settings", st1)

# ST2: Disorders EFA fit
fit_disorders <- read_if("disorder_efa/efa_fit_psych_simplified.csv") %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  rename(Factors = Factors, `Chi-Square` = ChiSq, `Degrees_of_Freedom` = df,
         TLI = TLI, BIC = BIC, `Total_Variance` = CumVar, `Delta_Variance_pp` = DeltaVar)
add_sheet_simple(wb, "ST2_EFA_FitStats_Disorders", fit_disorders)

# ST3: Disorders loadings
L1_d <- load_efa_csv("disorder_efa/efa_loadings_1fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L2_d <- load_efa_csv("disorder_efa/efa_loadings_2fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L3_d <- load_efa_csv("disorder_efa/efa_loadings_3fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L4_d <- load_efa_csv("disorder_efa/efa_loadings_4fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
addWorksheet(wb, "ST3_Loadings_Disorders_1to4"); setColWidths(wb, "ST3_Loadings_Disorders_1to4", 1:50, "auto")
r <- 1; r <- write_block(wb, "ST3_Loadings_Disorders_1to4", "1-Factor Loadings (promax)", L1_d, r)
r <- write_block(wb, "ST3_Loadings_Disorders_1to4", "2-Factor Loadings (promax)", L2_d, r)
r <- write_block(wb, "ST3_Loadings_Disorders_1to4", "3-Factor Loadings (promax)", L3_d, r)
r <- write_block(wb, "ST3_Loadings_Disorders_1to4", "4-Factor Loadings (promax)", L4_d, r)

# ST4: EFA Fit Stats — Symptoms (no TLI)
fit_symptoms <- read_if("symptoms_efa/efa_fit_psych_simplified.csv") %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  rename(
    Factors = Factors,
    `Chi-Square` = ChiSq,
    `Degrees_of_Freedom` = df,
    BIC = BIC,
    `Total_Variance` = CumVar,
    `Delta_Variance_pp` = DeltaVar
  ) %>%
  select(-any_of(c("TLI")))
add_sheet_simple(wb, "ST4_EFA_FitStats_Symptoms", fit_symptoms)

# ST5: Symptoms loadings
L1_s <- load_efa_csv("symptoms_efa/efa_loadings_1fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L2_s <- load_efa_csv("symptoms_efa/efa_loadings_2fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L3_s <- load_efa_csv("symptoms_efa/efa_loadings_3fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L4_s <- load_efa_csv("symptoms_efa/efa_loadings_4fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L5_s <- load_efa_csv("symptoms_efa/efa_loadings_5fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
addWorksheet(wb, "ST5_Loadings_Symptoms_1to5"); setColWidths(wb, "ST5_Loadings_Symptoms_1to5", 1:50, "auto")
r <- 1; r <- write_block(wb, "ST5_Loadings_Symptoms_1to5", "1-Factor Loadings (promax)", L1_s, r)
r <- write_block(wb, "ST5_Loadings_Symptoms_1to5", "2-Factor Loadings (promax)", L2_s, r)
r <- write_block(wb, "ST5_Loadings_Symptoms_1to5", "3-Factor Loadings (promax)", L3_s, r)
r <- write_block(wb, "ST5_Loadings_Symptoms_1to5", "4-Factor Loadings (promax)", L4_s, r)
r <- write_block(wb, "ST5_Loadings_Symptoms_1to5", "5-Factor Loadings (promax)", L5_s, r)

# ST7: Combined CFA fits (Symptoms + Disorders), floor001 omitted, NO SHADING
load_fit_summary <- function(xlsx_path, model_label) {
  sh <- openxlsx::getSheetNames(xlsx_path); sheet <- if ("fit_summary" %in% sh) "fit_summary" else sh[1]
  raw <- openxlsx::read.xlsx(xlsx_path, sheet = sheet)
  nm <- names(raw)
  tab <- raw[, c(pick_col(nm, c("Model","model")), pick_col(nm, c("ChiSq","Chisq","chisq")),
                 pick_col(nm, c("df","DF")), pick_col(nm, c("CFI","cfi"))), drop = FALSE]
  names(tab) <- c("Model","ChiSq","df","CFI")
  tab %>%
    filter(!str_detect(Model, "floor001")) %>%
    mutate(ChiSq = as.numeric(ChiSq), df = as.integer(round(as.numeric(df))),
           CFI = as.numeric(CFI), CFA_Model = model_label) %>%
    select(CFA_Model, everything())
}
sym_xlsx <- file.path(out_dir, "symptoms_cfa", "cfa_DWLS_results.xlsx")
dis_xlsx <- file.path(out_dir, "disorder_cfa",  "cfa_DWLS_results.xlsx")
st7_both <- bind_rows(load_fit_summary(sym_xlsx, "Symptoms"),
                      load_fit_summary(dis_xlsx, "Disorders")) %>%
  arrange(CFA_Model, desc(str_detect(Model, "^orig_")), Model)

addWorksheet(wb, "ST7_CFA_DWLS_Fits")
setColWidths(wb, "ST7_CFA_DWLS_Fits", 1:ncol(st7_both), "auto")
writeData(wb, "ST7_CFA_DWLS_Fits",
          "A. CFA Fit Indices (DWLS): Unconstrained vs. Constrained — Symptoms & Disorders (floor001 omitted)",
          startRow = 1, startCol = 1)
mergeCells(wb, "ST7_CFA_DWLS_Fits", rows = 1, cols = 1:ncol(st7_both))
addStyle(wb, "ST7_CFA_DWLS_Fits", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, "ST7_CFA_DWLS_Fits", st7_both, startRow = 2, startCol = 1, withFilter = TRUE)
freezePane(wb, "ST7_CFA_DWLS_Fits", firstRow = TRUE)

# ST8: Disorders — Loadings (with SE/p if present), Residuals, Φ
dis_rds <- file.path(out_dir, "disorder_cfa", "cfa_DWLS_results.RDS")
if (file.exists(dis_rds)) {
  L_dis_wide_raw <- final_constrained_loadings(dis_rds)
  L_dis_tidy     <- normalize_cfa_output(L_dis_wide_raw)
  
  Phi_dis <- extract_phi_matrix(dis_rds)
  
  L_dis_for_resid <- L_dis_tidy %>%
    dplyr::select(Trait, Factor, Loading) %>%
    tidyr::pivot_wider(names_from = Factor, values_from = Loading)
  
  Resid_dis <- compute_residuals(L_dis_for_resid, Phi_dis)
  
  addWorksheet(wb, "ST8_CFA_Disorders")
  r <- 1
  r <- write_long_center0(
    wb, "ST8_CFA_Disorders",
    "A. Standardized Loadings — Final Constrained CFA (Disorders; DWLS)",
    L_dis_tidy, r
  )
  r <- write_block_plain(
    wb, "ST8_CFA_Disorders",
    paste0("B. Residual Variances (1 − h²)", if (is.null(Phi_dis)) " — computed assuming Φ = I" else ""),
    Resid_dis, r
  )
  r <- write_corr_center0(
    wb, "ST8_CFA_Disorders",
    "C. Latent Factor Correlations (Φ)", Phi_dis, r
  )
}

# ST9: Symptoms — Loadings (with SE/p if present), Residuals, Φ
sym_rds <- file.path(out_dir, "symptoms_cfa", "cfa_DWLS_results.RDS")
if (file.exists(sym_rds)) {
  L_sym_wide_raw <- final_constrained_loadings(sym_rds)
  L_sym_tidy     <- normalize_cfa_output(L_sym_wide_raw)
  
  Phi_sym <- extract_phi_matrix(sym_rds)
  
  L_sym_for_resid <- L_sym_tidy %>%
    dplyr::select(Trait, Factor, Loading) %>%
    tidyr::pivot_wider(names_from = Factor, values_from = Loading)
  
  Resid_sym <- compute_residuals(L_sym_for_resid, Phi_sym)
  
  addWorksheet(wb, "ST9_CFA_Symptoms")
  r <- 1
  r <- write_long_center0(
    wb, "ST9_CFA_Symptoms",
    "A. Standardized Loadings — Final Constrained CFA (Symptoms; DWLS)",
    L_sym_tidy, r
  )
  r <- write_block_plain(
    wb, "ST9_CFA_Symptoms",
    paste0("B. Residual Variances (1 − h²)", if (is.null(Phi_sym)) " — computed assuming Φ = I" else ""),
    Resid_sym, r
  )
  r <- write_corr_center0(
    wb, "ST9_CFA_Symptoms",
    "C. Latent Factor Correlations (Φ)", Phi_sym, r
  )
}

# README
addWorksheet(wb, "README")
readme <- c(
  "Supplemental Tables (Genetic EFA & CFA on LDSC-derived models)",
  "",
  "ST1_GWAS_Sumstats_Settings: Settings/metadata for all GWAS summary statistics.",
  "ST2_EFA_FitStats_Disorders: EFA fit statistics (disorders).",
  "ST3_Loadings_Disorders_1to4: EFA loadings (disorders; blue–white–red).",
  "ST4_EFA_FitStats_Symptoms: EFA fit statistics (symptoms).",
  "ST5_Loadings_Symptoms_1to5: EFA loadings (symptoms; blue–white–red).",
  "ST7_CFA_DWLS_Fits: CFA fit indices for BOTH models (Symptoms & Disorders); floor001 omitted; includes CFA_Model; no shading.",
  "ST8_CFA_Disorders: (A) standardized loadings (SE/p in scientific if present), (B) residual variances, (C) factor correlations Φ.",
  "ST9_CFA_Symptoms: (A) standardized loadings (SE/p in scientific if present), (B) residual variances, (C) factor correlations Φ."
)
writeData(wb, "README", readme); setColWidths(wb, "README", 1, 120)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
