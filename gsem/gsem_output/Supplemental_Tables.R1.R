# ===========================
# Build Supplemental_Tables.xlsx (simple version)
# ===========================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(openxlsx); library(tibble); library(lavaan)
})

source('gsem/gsem_output/Supplemental_Table_Functions.R')

# ---- R1 (JCPP R&R): reporting conventions --------------------------------
# CFA results are read from the *_cfa_R1 folders (final model: free Internalizing~~Development
# correlation, BMI residual floor). The Child Physical Development factor is REFLECTED for
# reporting (loadings and its factor correlations x -1) so that higher scores = earlier/greater
# physical development (higher childhood BMI, earlier menarche). Model fit is unaffected.
REFLECT_FACTOR <- "Child Physical Development"
reflect_loadings <- function(L, factor = REFLECT_FACTOR) {
  L$Loading[L$Factor == factor] <- -L$Loading[L$Factor == factor]; L
}
# R1: standardized residual variances taken directly from the fitted model (not recomputed as
# max(0, 1 - h2), which hid the negative BMI residual in the original tables). Parameters fixed
# by the residual floor have no SE and are flagged.
fitted_residuals <- function(res_long) {
  latents <- unique(res_long$lhs[res_long$op == "=~"])
  th <- res_long[res_long$op == "~~" & res_long$lhs == res_long$rhs & !(res_long$lhs %in% latents), ]
  se <- suppressWarnings(as.numeric(th$STD_Genotype_SE))
  tibble::tibble(Trait = th$lhs,
                 Residual = round(as.numeric(th$STD_Genotype), 3),
                 Note = ifelse(is.na(se), "fixed at .01 (residual floor)", ""))
}
reflect_phi <- function(Phi, factor = REFLECT_FACTOR) {
  if (is.null(Phi) || !factor %in% rownames(Phi)) return(Phi)
  i <- which(rownames(Phi) == factor); j <- which(colnames(Phi) == factor)
  Phi[i, -j] <- -Phi[i, -j]; Phi[-i, j] <- -Phi[-i, j]; Phi
}
# ---- Paths ----
out_dir   <- file.path("gsem", "gsem_output")
out_xlsx  <- file.path(out_dir, "Supplemental_Tables_R1.xlsx")  # R1: separate output; original workbook untouched
gwas_file <- file.path("data", "gwas_sumstats", "GWAS_sumstats_settings_R1.xlsx")  # R1: corrected settings (scripts/0.Settings_corrections.R1.R)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



# ===========================================
# Build workbook
# ===========================================
wb <- createWorkbook()

# ST1
st1 <- openxlsx::read.xlsx(gwas_file) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>% arrange(across(1)) %>%
  select(
    "Trait Name" = name,
    "Trait Code" = code,
    Type = type,
    Ancestry = ancestry,
    Year,
    "Population Prevalence" = populationPrevalence,
    "Sumstats Sample Prevalence" = sumstats_samplePrevalence,
    "Case N" = n_case,
    "Control N" = n_control,
    "Total N" = n_total,
    "GSEM N" = gsem_n,
    Reference = reference_doi
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange("Trait Code") %>%
  mutate(Type = str_to_title(Type))  %>%
  mutate(
    Reference = str_squish(Reference),  # remove stray spaces
    Reference = case_when(
      
      # If it's already a DOI URL
      str_detect(Reference, "^https://doi\\.org/") ~ Reference,
      
      # medRxiv URL (leave untouched)
      str_detect(Reference, "medrxiv") ~ Reference,
      
      # Extract and convert raw DOI or DOI with prefix
      str_detect(str_to_lower(Reference), "doi") ~ 
        paste0("https://doi.org/", str_remove_all(str_to_lower(Reference), "doi:|doi|https://doi\\.org/")),
      
      # PMIDs → standard "PMID: ######"
      str_detect(Reference, "PMID") ~ Reference,
      
      TRUE ~ Reference
    )
  )
add_sheet_simple(wb, "ST1_GWAS_Sumstats_Settings", st1)

# ST3: Disorders EFA — variance explained (R1: geomin EFA; chi-square/TLI/BIC dropped)
efa_var_table <- function(f) {
  read_if(f) %>%
    transmute(Factors,
              `Cumulative common variance explained` = round(CumVar, 3),
              `Increment (percentage points)` = round(DeltaVar, 1))
}
EFA_NOTE <- paste("Note. Exploratory factor analysis (maximum likelihood, oblique geomin rotation) of the LDSC genetic",
                  "correlation matrix. Common variance explained is the mean communality across indicators.",
                  "Chi-square-based fit indices are not reported because they depend on an arbitrary sample size for",
                  "a genetic correlation matrix; factor retention considered scree plots, parallel analysis, eigenvalues and interpretability of loadings.")
fit_disorders <- efa_var_table("disorder_efa_R1/efa_fit_psych_simplified.csv")
add_sheet_simple(wb, "ST3_EFA_Variance_Disorders", fit_disorders)
writeData(wb, "ST3_EFA_Variance_Disorders", EFA_NOTE, startRow = nrow(fit_disorders) + 3, startCol = 1)

# ST4: Disorders loadings
L1_d <- load_efa_csv("disorder_efa_R1/efa_loadings_1fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L2_d <- load_efa_csv("disorder_efa_R1/efa_loadings_2fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L3_d <- load_efa_csv("disorder_efa_R1/efa_loadings_3fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L4_d <- load_efa_csv("disorder_efa_R1/efa_loadings_4fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
addWorksheet(wb, "ST4_Loadings_Disorders_1to4"); setColWidths(wb, "ST4_Loadings_Disorders_1to4", 1:50, "auto")
r <- 1; r <- write_block(wb, "ST4_Loadings_Disorders_1to4", "1-Factor Loadings (geomin)", L1_d, r)
r <- write_block(wb, "ST4_Loadings_Disorders_1to4", "2-Factor Loadings (geomin)", L2_d, r)
r <- write_block(wb, "ST4_Loadings_Disorders_1to4", "3-Factor Loadings (geomin)", L3_d, r)
r <- write_block(wb, "ST4_Loadings_Disorders_1to4", "4-Factor Loadings (geomin)", L4_d, r)

# ST6: Symptoms EFA — variance explained (R1)
fit_symptoms <- efa_var_table("symptoms_efa_R1/efa_fit_psych_simplified.csv")
add_sheet_simple(wb, "ST6_EFA_Variance_Symptoms", fit_symptoms)
writeData(wb, "ST6_EFA_Variance_Symptoms", EFA_NOTE, startRow = nrow(fit_symptoms) + 3, startCol = 1)

# ST7: Symptoms loadings
L1_s <- load_efa_csv("symptoms_efa_R1/efa_loadings_1fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L2_s <- load_efa_csv("symptoms_efa_R1/efa_loadings_2fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L3_s <- load_efa_csv("symptoms_efa_R1/efa_loadings_3fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L4_s <- load_efa_csv("symptoms_efa_R1/efa_loadings_4fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
L5_s <- load_efa_csv("symptoms_efa_R1/efa_loadings_5fac_dx.csv") %>% round_df(2) %>% order_by_top_factor()
addWorksheet(wb, "ST7_Loadings_Symptoms_1to5"); setColWidths(wb, "ST7_Loadings_Symptoms_1to5", 1:50, "auto")
r <- 1; r <- write_block(wb, "ST7_Loadings_Symptoms_1to5", "1-Factor Loadings (geomin)", L1_s, r)
r <- write_block(wb, "ST7_Loadings_Symptoms_1to5", "2-Factor Loadings (geomin)", L2_s, r)
r <- write_block(wb, "ST7_Loadings_Symptoms_1to5", "3-Factor Loadings (geomin)", L3_s, r)
r <- write_block(wb, "ST7_Loadings_Symptoms_1to5", "4-Factor Loadings (geomin)", L4_s, r)
r <- write_block(wb, "ST7_Loadings_Symptoms_1to5", "5-Factor Loadings (geomin)", L5_s, r)

# ST2: Combined CFA fits (Symptoms + Disorders), floor001 omitted, NO SHADING
load_fit_summary <- function(xlsx_path, model_label) {
  sh <- openxlsx::getSheetNames(xlsx_path)
  sheet <- if ("fit_summary" %in% sh) "fit_summary" else sh[1]
  raw <- openxlsx::read.xlsx(xlsx_path, sheet = sheet)
  nm <- names(raw)
  tab <- raw[, c(
    pick_col(nm, c("Model","model")),
    pick_col(nm, c("ChiSq","Chisq","chisq")),
    pick_col(nm, c("df","DF")),
    pick_col(nm, c("CFI","cfi"))
  ), drop = FALSE]
  names(tab) <- c("Model","ChiSq","df","CFI")
  tab %>%
    filter(str_detect(Model, "floor001")) %>%   # R1: report the floored (final) models, matching ST5/ST8
    mutate(
      ChiSq = as.numeric(ChiSq),
      df = as.integer(round(as.numeric(df))),
      CFI = as.numeric(CFI),
      CFA_Model = model_label
    ) %>%
    select(CFA_Model, everything())
}

sym_xlsx <- file.path(out_dir, "symptoms_cfa_R1", "cfa_DWLS_results.xlsx")
dis_xlsx <- file.path(out_dir, "disorder_cfa_R1",  "cfa_DWLS_results.xlsx")

st2_both <- bind_rows(
  load_fit_summary(sym_xlsx, "Symptoms"),
  load_fit_summary(dis_xlsx, "Disorders")
) %>%
  arrange(CFA_Model, desc(str_detect(Model, "^orig_")), Model) %>%
  # reader-facing row labels (base = no residual-variance floor; floor sensitivity rows are dropped above)
  mutate(Model = case_when(
    str_detect(Model, "^orig_")        ~ "Full CFA (EFA-based specification)",
    str_detect(Model, "^constrained_") ~ "Final CFA",
    str_detect(Model, "^Delta")        ~ "Difference (Final \u2212 Full)",
    TRUE ~ Model
  )) %>%
  mutate(order = match(Model, c("Full CFA (EFA-based specification)", "Final CFA", "Difference (Final \u2212 Full)"))) %>%
  arrange(CFA_Model, order) %>% select(-order)

addWorksheet(wb, "ST2_CFA_DWLS_Fits")
setColWidths(wb, "ST2_CFA_DWLS_Fits", 1:ncol(st2_both), "auto")
writeData(
  wb, "ST2_CFA_DWLS_Fits",
  "CFA Fit Indices (DWLS) for the Full and Final Models — Disorders and Symptoms (residual variance of BMI fixed at .01 to resolve a Heywood case)",
  startRow = 1, startCol = 1
)
mergeCells(wb, "ST2_CFA_DWLS_Fits", rows = 1, cols = 1:ncol(st2_both))
addStyle(wb, "ST2_CFA_DWLS_Fits", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, "ST2_CFA_DWLS_Fits", st2_both, startRow = 2, startCol = 1, withFilter = TRUE)
freezePane(wb, "ST2_CFA_DWLS_Fits", firstRow = TRUE)
# R1: table note
writeData(wb, "ST2_CFA_DWLS_Fits",
  paste("Note. The Final CFA removes cross-loadings of childhood BMI and binge eating on the Adult Anthropometric",
        "factor that were below |.20| and non-significant in the Full CFA; all factor correlations are freely estimated.",
        "In Genomic SEM, model chi-square is computed with the full sampling covariance matrix rather than the DWLS",
        "fit function, so chi-square differences between nested models are not likelihood-ratio tests and may be negative."),
  startRow = nrow(st2_both) + 4, startCol = 1)

# ST5: Disorders — Loadings (with SE/p if present), Residuals, Φ
dis_rds <- file.path(out_dir, "disorder_cfa_R1", "cfa_DWLS_results.RDS")
if (file.exists(dis_rds)) {
  L_dis_wide_raw <- final_constrained_loadings(dis_rds)
  L_dis_tidy     <- normalize_cfa_output(L_dis_wide_raw)
  
  L_dis_tidy <- L_dis_tidy %>%
    dplyr::mutate(
      Factor = dplyr::recode(
        Factor,
        "Development"     = "Child Physical Development",
        "Anthropometric"  = "Adult Anthropometric",
        "Internalizing"   = "Internalizing"
      )
    )
  
  Phi_dis <- extract_phi_matrix(dis_rds)
  old_to_new_dis <- c(
    "Development"    = "Child Physical Development",
    "Anthropometric" = "Adult Anthropometric",
    "Internalizing"  = "Internalizing"
  )
  rownames(Phi_dis) <- old_to_new_dis[rownames(Phi_dis)]
  colnames(Phi_dis) <- old_to_new_dis[colnames(Phi_dis)]
  # R1: reflect Child Physical Development for reporting
  L_dis_tidy <- reflect_loadings(L_dis_tidy)
  Phi_dis    <- reflect_phi(Phi_dis)
  
  L_dis_for_resid <- L_dis_tidy %>%
    dplyr::select(Trait, Factor, Loading) %>%
    tidyr::pivot_wider(names_from = Factor, values_from = Loading)
  
  Resid_dis <- fitted_residuals(L_dis_wide_raw)  # R1: from fitted model
  
  addWorksheet(wb, "ST5_CFA_Disorders")
  r <- 1
  r <- write_long_center0(
    wb, "ST5_CFA_Disorders",
    "A. Standardized Loadings — Final Constrained CFA (Disorders; DWLS)",
    L_dis_tidy, r
  )
  r <- write_block_plain(
    wb, "ST5_CFA_Disorders",
    "B. Standardized Residual Variances (fitted model)",
    Resid_dis, r
  )
  r <- write_corr_center0(
    wb, "ST5_CFA_Disorders",
    "C. Latent Factor Correlations (Φ)", Phi_dis, r
  )
}

# ST8: Symptoms — Loadings (with SE/p if present), Residuals, Φ
sym_rds <- file.path(out_dir, "symptoms_cfa_R1", "cfa_DWLS_results.RDS")
if (file.exists(sym_rds)) {
  L_sym_wide_raw <- final_constrained_loadings(sym_rds)
  L_sym_tidy     <- normalize_cfa_output(L_sym_wide_raw)
  
  L_sym_tidy <- L_sym_tidy %>%
    dplyr::mutate(
      Factor = dplyr::recode(
        Factor,
        "Development"     = "Child Physical Development",
        "Anthropometric"  = "Adult Anthropometric",
        "Internalizing"   = "Internalizing"
      )
    )
  
  Phi_sym <- extract_phi_matrix(sym_rds)
  old_to_new_sym <- c(
    "Development"    = "Child Physical Development",
    "Anthropometric" = "Adult Anthropometric",
    "Internalizing"  = "Internalizing"
  )
  rownames(Phi_sym) <- old_to_new_sym[rownames(Phi_sym)]
  colnames(Phi_sym) <- old_to_new_sym[colnames(Phi_sym)]
  # R1: reflect Child Physical Development for reporting
  L_sym_tidy <- reflect_loadings(L_sym_tidy)
  Phi_sym    <- reflect_phi(Phi_sym)
  
  L_sym_for_resid <- L_sym_tidy %>%
    dplyr::select(Trait, Factor, Loading) %>%
    tidyr::pivot_wider(names_from = Factor, values_from = Loading)
  
  Resid_sym <- fitted_residuals(L_sym_wide_raw)  # R1: from fitted model
  
  addWorksheet(wb, "ST8_CFA_Symptoms")
  r <- 1
  r <- write_long_center0(
    wb, "ST8_CFA_Symptoms",
    "A. Standardized Loadings — Final Constrained CFA (Symptoms; DWLS)",
    L_sym_tidy, r
  )
  r <- write_block_plain(
    wb, "ST8_CFA_Symptoms",
    "B. Standardized Residual Variances (fitted model)",
    Resid_sym, r
  )
  r <- write_corr_center0(
    wb, "ST8_CFA_Symptoms",
    "C. Latent Factor Correlations (Φ)", Phi_sym, r
  )
}

# ST9 (R1, new): Aim 3 conditional genomic regression results (R2 #10)
trait_labels <- c(
  AN = "Anorexia nervosa", BEBROAD = "Binge eating (broad)", BE_noAN = "Binge eating (cohorts not ascertained for AN)",
  OCD = "Obsessive-compulsive disorder", MDD_WtGain_Case = "MDD with weight gain", MDD_WtLoss_Case = "MDD with weight loss",
  Tanner_StageMidAdol = "Tanner stage, mid-adolescence", Age_Menarche = "Age at menarche",
  Insulin = "Fasting insulin (BMI-adjusted)", BMI = "Adult BMI", BMI_Childhood = "Childhood BMI")
lab <- function(x) ifelse(x %in% names(trait_labels), trait_labels[x], x)
fmt_p <- function(x) ifelse(is.na(x), NA_character_, formatC(x, format = "e", digits = 2))
st9 <- read_if("aim3_conditional_R1/gmr_trios/gmr_trios_BMIfocused.csv") %>%
  transmute(
    `ED outcome`      = lab(Outcome),
    `Focal trait`     = lab(Focal),
    `Conditioned on`  = lab(Controlled_for),
    `b (unadjusted)`  = round(b_raw, 3),
    `SE (unadjusted)` = round(se_raw, 3),
    `Std. β (unadjusted)` = round(beta_raw_std, 3),
    `p (unadjusted)`  = fmt_p(p_raw),
    `FDR p (unadjusted)` = fmt_p(p_raw_FDR),
    `b (BMI-adjusted)`  = round(b_adj, 3),
    `SE (BMI-adjusted)` = round(se_adj, 3),
    `Std. β (BMI-adjusted)` = round(beta_adj_std, 3),
    `p (BMI-adjusted)`  = fmt_p(p_adj),
    `FDR p (BMI-adjusted)` = fmt_p(p_adj_FDR)
  ) %>%
  arrange(`ED outcome`, `Focal trait`, `Conditioned on`)
add_sheet_simple(wb, "ST9_Aim3_Conditional", st9)
writeData(wb, "ST9_Aim3_Conditional",
  paste("Note. Genomic multiple regression of each ED outcome on a focal trait, before (unadjusted) and after",
        "conditioning on adult or childhood BMI. Trait trios were selected where the focal trait showed |rg| > .25",
        "with both an ED phenotype and a BMI trait. FDR = Benjamini-Hochberg false discovery rate correction within Aim 3.",
        "Values correspond to Figure 4c."),
  startRow = nrow(st9) + 3, startCol = 1)

# README
if ("README" %in% names(wb)) removeWorksheet(wb, "README")
addWorksheet(wb, "README")
readme <- c(
  "Supplementary Tables S1-S9: Mapping Shared and Unique Genetic Liability to Eating and Internalizing Psychopathology",
  "All tables are generated by gsem/gsem_output/Supplemental_Tables.R1.R (https://github.com/kschaumb/GENIE_GSEM).",
  "Trait codes used throughout are defined in Table S1 (columns 'Trait Name' and 'Trait Code').",
  "",
  "ST1_GWAS_Sumstats_Settings - Table S1. GWAS summary statistics included: trait name and code, phenotype type, ancestry, publication year, population and sample prevalence (binary traits), case/control/total N, N used in Genomic SEM, and reference.",
  "ST2_CFA_DWLS_Fits - Table S2. Fit of the full (EFA-based) and final confirmatory factor models (DWLS) for the disorder-level and symptom-level panels, and the difference between them. See table note on interpreting chi-square differences in Genomic SEM.",
  "ST3_EFA_Variance_Disorders - Table S3. Disorder-level exploratory factor analysis: cumulative common variance explained (mean communality) and increment for 1-5 factor solutions.",
  "ST4_Loadings_Disorders_1to4 - Table S4. Disorder-level EFA standardized loadings (maximum likelihood, oblique geomin rotation) for 1-4 factor solutions.",
  "ST5_CFA_Disorders - Table S5. Final disorder-level CFA: (A) standardized loadings with SE and p, (B) standardized residual variances from the fitted model, (C) latent factor correlations.",
  "ST6_EFA_Variance_Symptoms - Table S6. Symptom-level exploratory factor analysis: cumulative common variance explained and increment for 1-5 factor solutions.",
  "ST7_Loadings_Symptoms_1to5 - Table S7. Symptom-level EFA standardized loadings (maximum likelihood, oblique geomin rotation) for 1-5 factor solutions.",
  "ST8_CFA_Symptoms - Table S8. Final symptom-level CFA: (A) standardized loadings with SE and p, (B) standardized residual variances from the fitted model, (C) latent factor correlations.",
  "ST9_Aim3_Conditional - Table S9. Aim 3 conditional genomic regression results: associations of focal traits with ED phenotypes before and after conditioning on adult or childhood BMI (values shown in Figure 4c).",
  "",
  "Conventions: The Child Physical Development factor is scored so that higher values reflect higher childhood BMI and earlier menarche (loadings and correlations involving this factor are reflected relative to the raw model output; model fit is unaffected). In final CFA models, the residual variance of adult BMI was fixed at .01 to resolve a Heywood case."
)
writeData(wb, "README", readme)
setColWidths(wb, "README", 1, 120)

# R1: order sheets README, ST1..ST9
sh <- names(wb)
st <- sort(sh[grepl("^ST[0-9]", sh)])
worksheetOrder(wb) <- match(c("README", st), sh)

# save workbook
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

