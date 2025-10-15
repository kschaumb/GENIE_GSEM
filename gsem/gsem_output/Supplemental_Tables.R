# ===========================
# Build Supplemental_Tables.xlsx (simple version)
# ===========================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(openxlsx); library(tibble); library(lavaan)
})

# ---- Paths ----
out_dir   <- file.path("gsem", "gsem_output")
out_xlsx  <- file.path(out_dir, "Supplemental_Tables.xlsx")
gwas_file <- file.path("data", "gwas_sumstats", "GWAS_sumstats_settings.xlsx")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Small helpers ----
read_if <- function(fname) readr::read_csv(file.path(out_dir, fname), show_col_types = FALSE)

load_efa_csv <- function(file) {
  df <- read_if(file)
  if (is.na(names(df)[1]) || names(df)[1] %in% c("", "...1")) names(df)[1] <- "Trait"
  df %>% column_to_rownames("Trait")
}

round_df <- function(df, digits = 2) {
  num <- sapply(df, is.numeric); df[num] <- lapply(df[num], round, digits); df
}

order_by_top_factor <- function(L) {
  A <- as.matrix(L)
  ord <- order(apply(abs(A), 1, which.max), -apply(abs(A), 1, max), rownames(L))
  L[ord, , drop = FALSE]
}

fmt_sci2 <- function(x) if (is.numeric(x)) sprintf("%.2e", x) else x

add_sheet_simple <- function(wb, sheet, df) {
  addWorksheet(wb, sheet); writeData(wb, sheet, df, withFilter = TRUE)
  freezePane(wb, sheet, firstRow = TRUE); setColWidths(wb, sheet, 1:ncol(df), "auto")
}

write_block <- function(wb, sheet, title, df, start_row) {
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:(ncol(df) + 1))
  tbl <- data.frame(Trait = rownames(df), df, row.names = NULL, check.names = FALSE)
  writeData(wb, sheet, tbl, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  num_cols <- which(sapply(tbl, is.numeric))
  if (length(num_cols)) conditionalFormatting(
    wb, sheet, cols = num_cols,
    rows = (start_row + 2):(start_row + 1 + nrow(tbl)),
    style = c("#4575b4", "#ffffff", "#d73027"), type = "colorScale"
  )
  start_row + 1 + nrow(tbl) + 2
}

write_block_center0 <- function(wb, sheet, title, df, start_row) {
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:(ncol(df) + 1))
  tbl <- data.frame(Trait = rownames(df), df, row.names = NULL, check.names = FALSE)
  writeData(wb, sheet, tbl, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  
  num_cols <- which(sapply(tbl, is.numeric))
  if (length(num_cols)) {
    vals <- unlist(tbl[num_cols], use.names = FALSE)
    M <- max(abs(vals), na.rm = TRUE)
    if (is.finite(M) && M > 0) {
      conditionalFormatting(
        wb, sheet, cols = num_cols,
        rows = (start_row + 2):(start_row + 1 + nrow(tbl)),
        type  = "colorScale",
        style = c("#4575b4", "#ffffff", "#d73027"),
        rule  = c(-M, 0, M)
      )
    }
  }
  setColWidths(wb, sheet, 1:ncol(tbl), "auto")
  start_row + 1 + nrow(tbl) + 2
}

write_block_plain <- function(wb, sheet, title, df, start_row) {
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:ncol(df))
  writeData(wb, sheet, df, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  setColWidths(wb, sheet, 1:ncol(df), "auto")
  start_row + 1 + nrow(df) + 2
}

write_long_center0 <- function(wb, sheet, title, df, start_row) {
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:5)
  writeData(wb, sheet, df, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  if (nrow(df)) {
    M <- max(abs(df$Loading), na.rm = TRUE)
    if (is.finite(M) && M > 0) conditionalFormatting(
      wb, sheet, cols = 3, rows = (start_row + 2):(start_row + 1 + nrow(df)),
      type = "colorScale", style = c("#4575b4", "#ffffff", "#d73027"), rule = c(-M, 0, M)
    )
  }
  setColWidths(wb, sheet, 1:5, "auto")
  start_row + 1 + nrow(df) + 2
}

write_corr_center0 <- function(wb, sheet, title, Phi, start_row) {
  if (is.null(Phi)) {
    writeData(wb, sheet, paste0(title, " (not available)"), startCol = 1, startRow = start_row)
    return(start_row + 2)
  }
  Phi <- round(Phi, 3)
  if (is.null(colnames(Phi))) colnames(Phi) <- paste0("F", seq_len(ncol(Phi)))
  if (is.null(rownames(Phi))) rownames(Phi) <- colnames(Phi)
  df <- tibble::rownames_to_column(as.data.frame(Phi, check.names = FALSE), "Factor")
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:ncol(df))
  writeData(wb, sheet, df, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  conditionalFormatting(
    wb, sheet, cols = 2:ncol(df),
    rows = (start_row + 2):(start_row + 1 + nrow(df)),
    type = "colorScale", style = c("#4575b4", "#ffffff", "#d73027"), rule = c(-1, 0, 1)
  )
  setColWidths(wb, sheet, 1:ncol(df), "auto")
  start_row + 1 + nrow(df) + 2
}

# Minimal picker / null-coalescer
pick_col <- function(nm, choices) { hit <- intersect(choices, nm); if (length(hit)) hit[1] else NA_character_ }
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Choose best fit from $fits ----
.pick_fit <- function(nms) {
  for (pat in c("^constrained_DWLS_floor001$",
                "^constrained_DWLS_base$",
                "^constrained",
                "^orig_DWLS_floor001$",
                "^orig_DWLS_base$",
                "^orig")) {
    hit <- grep(pat, nms, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  nms[1]
}

# ---- Final constrained loader (reads from $fits; returns long 'results' if available) ----
final_constrained_loadings <- function(rds_path) {
  obj  <- readRDS(rds_path)
  fits <- if (is.list(obj$fits)) obj$fits else obj
  if (!length(fits)) return(data.frame())
  
  fit_name <- .pick_fit(names(fits))
  fit <- fits[[fit_name]]
  
  if (is.list(fit) && is.data.frame(fit$results) && all(c("lhs","op","rhs") %in% names(fit$results))) {
    return(fit$results)  # long parameter table; normalize later
  }
  
  L <- fit$loadings %||% fit$L %||% fit$std_loadings %||% fit$load
  if (!is.null(L)) return(as.data.frame(L, check.names = FALSE))
  
  L <- obj$loadings %||% obj$L %||% obj$std_loadings %||% obj$load
  as.data.frame(L, check.names = FALSE)
}

# ---- Φ extractor (reads from $fits; builds from latent '~~' in results) ----
extract_phi_matrix <- function(rds_path) {
  obj  <- readRDS(rds_path)
  fits <- if (is.list(obj$fits)) obj$fits else obj
  if (!length(fits)) return(NULL)
  
  fit_name <- .pick_fit(names(fits))
  fit <- fits[[fit_name]]
  
  if (is.matrix(fit$Phi)) return(fit$Phi)
  if (!is.null(fit$lavaan_output) && inherits(fit$lavaan_output, "lavaan")) {
    ph <- try(lavaan::inspect(fit$lavaan_output, "cor.lv"), silent = TRUE)
    if (!inherits(ph, "try-error") && is.matrix(ph)) return(ph)
  }
  
  if (!is.data.frame(fit$results) || !all(c("lhs","op","rhs") %in% names(fit$results))) return(NULL)
  res <- fit$results
  
  latents <- unique(res$lhs[res$op == "=~"])
  latents <- latents[!is.na(latents) & nzchar(latents)]
  if (!length(latents)) return(NULL)
  
  std_col <- pick_col(names(res), c("STD_Genotype","Std.all","STD_All","std.all"))
  if (is.na(std_col)) return(NULL)
  
  phi_long <- subset(res, op == "~~" & lhs %in% latents & rhs %in% latents,
                     select = c("lhs","rhs", std_col))
  names(phi_long)[3] <- "val"
  if (!nrow(phi_long)) return(NULL)
  
  Phi <- matrix(0, nrow = length(latents), ncol = length(latents),
                dimnames = list(latents, latents))
  diag(Phi) <- 1
  for (i in seq_len(nrow(phi_long))) {
    a <- phi_long$lhs[i]; b <- phi_long$rhs[i]
    v <- suppressWarnings(as.numeric(phi_long$val[i]))
    if (is.finite(v)) { Phi[a,b] <- v; Phi[b,a] <- v }
  }
  Phi
}

# ---- Normalize to 5 cols (robust SE/p; keeps strings like "< 5e-300") ----
# ---- Normalize to 5 cols (ONLY loadings; robust SE/p; keeps strings like "< 5e-300") ----
normalize_cfa_output <- function(df) {
  nm <- names(df)
  
  fmt_numeric_or_keep <- function(x, digits = 2) {
    if (is.null(x)) return(rep(NA_character_, nrow(df)))
    x_num <- suppressWarnings(as.numeric(x))
    if (all(is.na(x_num)) && is.character(x)) {
      x_num <- suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", x)))
    }
    out <- as.character(x)
    ok <- is.finite(x_num)
    out[ok] <- sprintf(paste0("%.", digits, "e"), x_num[ok])
    out
  }
  
  if (all(c("lhs","rhs") %in% nm)) {
    # LONG FORM: keep ONLY loadings rows
    df2 <- if ("op" %in% nm) df[df$op == "=~", , drop = FALSE] else df
    
    est <- pick_col(names(df2), c("std_loading","STD_Genotype","STD_All","Std.all","std.all","est","Estimate"))
    se  <- pick_col(names(df2), c("STD_Genotype_SE","SE","se","Unstand_SE"))
    pv  <- pick_col(names(df2), c("p_value","p","pvalue"))
    
    out <- df2[, c("rhs","lhs", est, se, pv)]
    names(out) <- c("Trait","Factor","Loading","SE","p")
    
    out$Loading <- suppressWarnings(as.numeric(out$Loading))
    out$SE <- fmt_numeric_or_keep(out$SE)
    out$p  <- fmt_numeric_or_keep(out$p)
    
  } else {
    # WIDE FORM
    if (!"Trait" %in% nm && !is.null(rownames(df))) df <- tibble::rownames_to_column(df, "Trait")
    num_cols <- names(df)[sapply(df, is.numeric)]
    num_cols <- setdiff(num_cols, c("SE","se","p","p_value","pvalue"))
    out <- tidyr::pivot_longer(df, cols = dplyr::all_of(num_cols),
                               names_to = "Factor", values_to = "Loading") %>%
      mutate(SE = NA_character_, p = NA_character_) %>%
      select(Trait, Factor, Loading, SE, p)
    out$Loading <- as.numeric(out$Loading)
  }
  
  # Order nicely within factor
  dplyr::arrange(out, Factor, dplyr::desc(abs(Loading)), Trait)
}


# ---- Residuals with Φ aligned to loading columns ----
compute_residuals <- function(loads_wide, phi = NULL) {
  if (!"Trait" %in% names(loads_wide)) stop("compute_residuals() needs a Trait column")
  
  factors <- setdiff(names(loads_wide), "Trait")
  A <- as.matrix(loads_wide[, factors, drop = FALSE])
  A[is.na(A)] <- 0
  
  # align/build Phi to match 'factors'
  if (is.null(phi)) {
    Phi <- diag(length(factors))
    dimnames(Phi) <- list(factors, factors)
  } else {
    Phi0 <- as.matrix(phi)
    Phi <- matrix(0, nrow = length(factors), ncol = length(factors),
                  dimnames = list(factors, factors))
    diag(Phi) <- 1
    rn <- rownames(Phi0); cn <- colnames(Phi0)
    if (!is.null(rn) && !is.null(cn)) {
      common <- intersect(factors, intersect(rn, cn))
      if (length(common)) {
        Phi[common, common] <- Phi0[common, common, drop = FALSE]
      }
    } else if (nrow(Phi0) == ncol(Phi0) && ncol(Phi0) == length(factors)) {
      Phi <- Phi0
      dimnames(Phi) <- list(factors, factors)
    }
  }
  
  h2 <- rowSums((A %*% Phi) * A)
  tibble::tibble(Trait = loads_wide$Trait, Residual = round(pmax(0, 1 - h2), 3))
}

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
