# ============================================================
# CFA (GenomicSEM, DWLS): Original vs Constrained + Sensitivity
# ============================================================
# Outputs -> gsem/gsem_output/disorder_cfa/
# Usage   -> source("scripts/cfa_dwls_pipeline.R")
# Notes   -> - DWLS only
#           - Constrained model: Development ~~ 0*Internalizing (via label+constraint)
#           - Optional sensitivity: residual floors on BMI/PTSD
#           - Writes ONE consolidated .RDS (everything) and optionally ONE .xlsx
# ============================================================
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

# --------------------
# 0) Parameters/Paths
# --------------------
base_path   <- "~/Library/CloudStorage/Box-Box/GENIE_GSEM"
out_dir     <- file.path(base_path, "gsem/gsem_output/disorder_cfa")
# Filter for European-ancestry traits
disorders_fa <- sumstats %>% 
  filter(Disorder_GSEM == 1)
disorders_fa <- disorders_fa$code

# Dynamically construct the expected mvLD filename
# (Assumes consistent naming scheme: "mvLD_<traits joined by _>.Rds")
mvld_path <- file.path(folderpaths$gsem.mvLD,
  paste0("mvLD_", paste(disorders_fa, collapse = "_"), ".Rds")
)

# Confirm the path and read
if (!file.exists(mvld_path)) {
  stop("❌ mvLD file not found: ", mvld_path)
}


# Diagnostics thresholds
THRESH_LOAD <- 1.00         # flag standardized loading > 1
THRESH_CORR <- 0.95         # flag |factor correlation| >= 0.95

# Sensitivity floors to apply to BMI/PTSD residual variances
residual_floors <- c(0.00, 0.01)

# Indicators to guard (from diagnostics insight)
guard_indicators <- c("BMI", "PTSD")

# Output mode
write_xlsx <- TRUE  # set FALSE for RDS-only
xlsx_name  <- "cfa_DWLS_results.xlsx"
rds_name   <- "cfa_DWLS_results.RDS"

# --------------------
# 1) Setup & Libraries
# --------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(GenomicSEM)
  library(lavaan)
})

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(base_path)

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
log_msg("Starting DWLS CFA pipeline (consolidated outputs)")

# --------------------
# 2) Load Inputs
# --------------------
load(file.path(base_path, "gsem/gsem_settings/folderpaths.RData"))
load(file.path(base_path, "gsem/gsem_settings/sumstats_metadata.RData"))
load(file.path(base_path, "gsem/gsem_settings/variant_paths.RData"))

ldsc_full <- readRDS(mvld_path)
ldsc <- if (!is.null(ldsc_full$hm3)) ldsc_full$hm3 else ldsc_full
stopifnot(!is.null(ldsc$S), !is.null(ldsc$V))
obs_traits <- colnames(ldsc$S)

# --------------------
# 3) Helpers
# --------------------
clean_model <- function(x){
  x <- gsub("\r", "", x)
  x <- gsub("[\u00A0]", " ", x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  lines <- strsplit(x, "\n", fixed = TRUE)[[1]]
  lines <- sub("\\s*#.*$", "", lines)
  lines <- trimws(lines)
  paste(lines[nzchar(lines)], collapse = "\n")
}
num_safe   <- function(x) suppressWarnings(as.numeric(x))
pick_col <- function(df, candidates) {
  hits <- intersect(candidates, names(df))
  if (length(hits)) hits[1] else NA_character_
}

fit_dwls <- function(model_string, std.lv = TRUE){
  out <- usermodel(covstruc = ldsc, estimation = "DWLS", model = model_string, std.lv = std.lv)
  out$model <- model_string
  out$estimation <- "DWLS"
  out
}

get_tables <- function(fit){
  if ("lavaan_output" %in% names(fit)) {
    pe <- lavaan::parameterEstimates(fit$lavaan_output, standardized = TRUE)
    if (!"STD_Genotype" %in% names(pe) && "std.all" %in% names(pe)) {
      pe$STD_Genotype <- pe$std.all
    }
    return(pe)
  }
  if ("results" %in% names(fit) && is.data.frame(fit$results)) return(fit$results)
  stop("Could not find parameter tables in fit object")
}

row_from_fit <- function(nm, f) {
  mf <- if ("modelfit" %in% names(f)) f$modelfit else NULL
  if (is.null(mf)) return(data.frame(Model = nm, Chisq = NA, df = NA, p_chisq = NA,
                                     AIC = NA, CFI = NA, SRMR = NA, Wald_p = NA))
  data.frame(Model = nm,
             Chisq = mf$chisq, df = mf$df, p_chisq = mf$p_chisq,
             AIC = mf$AIC, CFI = mf$CFI, SRMR = mf$SRMR,
             Wald_p = NA_real_)
}

diag_one_model <- function(model_name, fit_obj, obs_names,
                           THRESH_LOAD = 1.00, THRESH_CORR = 0.95) {
  pe <- get_tables(fit_obj)
  
  # Loadings
  L <- subset(pe, op == "=~")
  stdL_col <- pick_col(L, c("STD_Genotype","std.all","Std.all","STD.ALL","std"))
  if (is.na(stdL_col)) stdL_col <- pick_col(L, c("est","EST","Estimate","estimate"))
  L$std_loading <- num_safe(L[[stdL_col]])
  L$flag_loading_gt1 <- abs(L$std_loading) > THRESH_LOAD
  
  # Residuals (observed)
  Th <- subset(pe, op == "~~" & lhs == rhs & lhs %in% obs_names)
  est_col_Th <- pick_col(Th, c("est","EST","Estimate","estimate","Est"))
  if (is.na(est_col_Th)) Th$est <- NA_real_ else Th$est <- num_safe(Th[[est_col_Th]])
  stdT <- pick_col(Th, c("std.all","STD_Genotype","Std.all","STD.ALL","std"))
  Th$std_resid <- if (!is.na(stdT)) num_safe(Th[[stdT]]) else NA_real_
  Th$flag_resid_neg_unstd <- !is.na(Th$est) & Th$est < 0
  Th$flag_resid_neg_std   <- !is.na(Th$std_resid) & Th$std_resid < 0
  
  # Factor correlations (among latents only)
  latents <- unique(L$lhs)
  R <- subset(pe, op == "~~" & lhs %in% latents & rhs %in% latents & lhs != rhs)
  stdR_col <- pick_col(R, c("STD_Genotype","std.all","Std.all","STD.ALL","std"))
  if (is.na(stdR_col)) stdR_col <- pick_col(R, c("est","EST","Estimate","estimate"))
  R$std_corr <- num_safe(R[[stdR_col]])
  R$flag_hi_corr <- !is.na(R$std_corr) & abs(R$std_corr) >= THRESH_CORR
  
  # Worst offenders
  worst_unstd <- if (nrow(Th) && any(!is.na(Th$est))) Th[order(Th$est), ][1, c("lhs","est"), drop=FALSE] else data.frame(lhs=NA, est=NA)
  worst_std   <- if (nrow(Th) && any(!is.na(Th$std_resid))) Th[order(Th$std_resid), ][1, c("lhs","std_resid"), drop=FALSE] else data.frame(lhs=NA, std_resid=NA)
  worst_load  <- if (nrow(L)  && any(!is.na(L$std_loading))) L[order(-abs(L$std_loading)), ][1, c("lhs","rhs","std_loading"), drop=FALSE] else data.frame(lhs=NA, rhs=NA, std_loading=NA)
  worst_corr  <- if (nrow(R)  && any(!is.na(R$std_corr)))    R[order(-abs(R$std_corr)),   ][1, c("lhs","rhs","std_corr"),     drop=FALSE] else data.frame(lhs=NA, rhs=NA, std_corr=NA)
  
  summary_row <- data.frame(
    Model = model_name,
    n_neg_resid_unstd = sum(Th$flag_resid_neg_unstd, na.rm=TRUE),
    n_neg_resid_std   = sum(Th$flag_resid_neg_std,   na.rm=TRUE),
    min_resid_unstd   = round(worst_unstd$est, 6),
    min_resid_unstd_var = worst_unstd$lhs,
    min_resid_std     = round(worst_std$std_resid, 6),
    min_resid_std_var = worst_std$lhs,
    n_load_gt1        = sum(L$flag_loading_gt1, na.rm=TRUE),
    max_abs_loading   = round(max(abs(L$std_loading), na.rm=TRUE), 6),
    max_loading_item  = ifelse(is.na(worst_load$lhs), NA, paste0(worst_load$lhs," -> ",worst_load$rhs)),
    n_high_fac_corr   = sum(R$flag_hi_corr, na.rm=TRUE),
    max_abs_fac_corr  = round(max(abs(R$std_corr), na.rm=TRUE), 6),
    max_fac_corr_pair = ifelse(is.na(worst_corr$lhs), NA, paste0(worst_corr$lhs," ~~ ",worst_corr$rhs)),
    stringsAsFactors = FALSE
  )
  
  neg_unstd <- subset(Th, flag_resid_neg_unstd)[, c("lhs","est"), drop=FALSE]    |> dplyr::rename(value = est)
  if (nrow(neg_unstd)) neg_unstd$issue <- "neg_resid_unstd"
  
  neg_std   <- subset(Th, flag_resid_neg_std)[, c("lhs","std_resid"), drop=FALSE] |> dplyr::rename(value = std_resid)
  if (nrow(neg_std))   neg_std$issue   <- "neg_resid_std"
  
  gt1 <- subset(L, flag_loading_gt1)[, c("lhs","rhs","std_loading"), drop=FALSE] |> dplyr::rename(value = std_loading)
  if (nrow(gt1)) gt1$issue <- "loading_gt_1"
  
  hic <- subset(R, flag_hi_corr)[, c("lhs","rhs","std_corr"), drop=FALSE] |> dplyr::rename(value = std_corr)
  if (nrow(hic)) hic$issue <- "high_factor_corr"
  
  issues <- dplyr::bind_rows(neg_unstd, neg_std, gt1, hic)
  if (nrow(issues)) issues$Model <- model_name
  
  list(summary = summary_row, issues = issues,
       tables = list(loadings=L, theta=Th, fac_corr=R))
}

# --------------------
# 4) Model Templates
# --------------------
model_original <- "
Anthropometric  =~ BMI + BMI_Childhood + BFP + AN + BE_noAN + OCD + FFM
Development =~ Age_Menarche + BE_noAN + BMI_Childhood 
Internalizing     =~ MDD + PTSD + OCD + AN + BE_noAN + ANX
Development ~~ Internalizing
Anthropometric  ~~ Development
Anthropometric  ~~ Internalizing
"

model_constrained <- "
Anthropometric  =~ BMI + BFP + AN + OCD + FFM
Development =~ Age_Menarche + BE_noAN + BMI_Childhood 
Internalizing     =~ MDD + PTSD + OCD + AN + BE_noAN + ANX
Development ~~ cPI*Internalizing
Anthropometric  ~~ Development
Anthropometric  ~~ Internalizing
cPI == 0
"


# Parse sanity
invisible(lavaan::lavaanify(model_original))
invisible(lavaan::lavaanify(model_constrained))

# -------------------------------------------------
# 5) Run fits across residual floors (sensitivity)
# -------------------------------------------------
append_resid_floors <- function(model_txt, floor_value, indicators) {
  if (isTRUE(all.equal(floor_value, 0))) return(model_txt)  # no change
  lines <- paste(sprintf("%s ~~ %g*%s", indicators, floor_value, indicators), collapse = "\n")
  paste(model_txt, lines, sep = "\n")
}

fits <- list()
fit_rows <- list()
issues_all <- list()
diag_summary_rows <- list()
per_fit_tables <- list()  # kept in-memory; placed in RDS

for (f in residual_floors) {
  suffix <- if (isTRUE(all.equal(f, 0))) "base" else paste0("floor", gsub("\\.", "", as.character(f)))
  m_orig_f <- append_resid_floors(model_original,    f, guard_indicators)
  m_con_f  <- append_resid_floors(model_constrained, f, guard_indicators)
  
  log_msg("Fitting original (DWLS), floor=", f)
  fit_o <- fit_dwls(m_orig_f)
  log_msg("Fitting constrained (DWLS), floor=", f)
  fit_c <- fit_dwls(m_con_f)
  
  name_o <- paste0("orig_DWLS_", suffix)
  name_c <- paste0("constrained_DWLS_", suffix)
  fits[[name_o]] <- fit_o
  fits[[name_c]] <- fit_c
  
  fit_rows[[name_o]] <- row_from_fit(name_o, fit_o)
  fit_rows[[name_c]] <- row_from_fit(name_c, fit_c)
  
  d_o <- diag_one_model(name_o, fit_o, obs_traits, THRESH_LOAD, THRESH_CORR)
  d_c <- diag_one_model(name_c, fit_c, obs_traits, THRESH_LOAD, THRESH_CORR)
  diag_summary_rows[[name_o]] <- d_o$summary
  diag_summary_rows[[name_c]] <- d_c$summary
  
  if (nrow(d_o$issues)) issues_all[[name_o]] <- d_o$issues
  if (nrow(d_c$issues)) issues_all[[name_c]] <- d_c$issues
  
  # keep detailed tables per fit only in the RDS to avoid file clutter
  per_fit_tables[[name_o]] <- d_o$tables
  per_fit_tables[[name_c]] <- d_c$tables
}

# -----------------------------------------
# 6) Wald test (unconstrained, base floor)
# -----------------------------------------
wald_p <- NA_real_
nm_base <- "orig_DWLS_base"
wald_obj_print <- NULL
if (nm_base %in% names(fits) && "lavaan_output" %in% names(fits[[nm_base]])) {
  wobj <- try(lavaan::lavTestWald(fits[[nm_base]]$lavaan_output, "Development~~Internalizing == 0"), silent = TRUE)
  if (!inherits(wobj, "try-error")) {
    wald_p <- tryCatch({
      if (is.data.frame(wobj) && "p.value" %in% names(wobj)) as.numeric(wobj$p.value[1])
      else if (!is.null(wobj[["p.value"]])) as.numeric(wobj[["p.value"]][1])
      else if (!is.null(attr(wobj, "p.value"))) as.numeric(attr(wobj, "p.value"))
      else NA_real_
    }, error = function(e) NA_real_)
    wald_obj_print <- capture.output(wobj)
  }
}

# -------------------------------------------------------
# 7) Build one tidy summary (with Δ rows per floor)
# -------------------------------------------------------
fit_rows_df <- dplyr::bind_rows(fit_rows)
fit_rows_df$Wald_p[fit_rows_df$Model == nm_base] <- wald_p

mk_delta <- function(df, suffix) {
  r_o <- df %>% dplyr::filter(Model == paste0("orig_DWLS_", suffix))
  r_c <- df %>% dplyr::filter(Model == paste0("constrained_DWLS_", suffix))
  if (!nrow(r_o) || !nrow(r_c)) return(NULL)
  out <- r_c
  out$Model   <- paste0("Delta(constrained - orig)_", suffix)
  out$Chisq   <- r_c$Chisq   - r_o$Chisq
  out$df      <- r_c$df      - r_o$df
  out$p_chisq <- r_c$p_chisq - r_o$p_chisq
  out$AIC     <- r_c$AIC     - r_o$AIC
  out$CFI     <- r_c$CFI     - r_o$CFI
  out$SRMR    <- r_c$SRMR    - r_o$SRMR
  out$Wald_p  <- NA_real_
  out
}
suffixes <- ifelse(residual_floors == 0, "base", paste0("floor", gsub("\\.", "", residual_floors)))
delta_rows <- dplyr::bind_rows(lapply(suffixes, function(suf) mk_delta(fit_rows_df, suf)))

fits_one <- dplyr::bind_rows(fit_rows_df, delta_rows)
num_cols <- setdiff(names(fits_one), "Model")
fits_one[num_cols] <- lapply(fits_one[num_cols], function(x) round(as.numeric(x), 4))

# -------------------------------------------------
# 8) Diagnostics (compact + long-form issues)
# -------------------------------------------------
diag_summary_df <- dplyr::bind_rows(diag_summary_rows)
issues_df <- if (length(issues_all)) dplyr::bind_rows(issues_all) else
  data.frame(message = "No flagged issues under current thresholds.", stringsAsFactors = FALSE)

# ----------------------
# 9) Provenance/Session
# ----------------------
provenance <- list(
  session_info = capture.output(sessionInfo()),
  obs_traits   = obs_traits,
  ldsc_shapes  = list(S = dim(ldsc$S), V_diag_len = length(diag(ldsc$V))),
  thresholds   = list(THRESH_LOAD = THRESH_LOAD, THRESH_CORR = THRESH_CORR),
  residual_floors = residual_floors,
  guard_indicators = guard_indicators
)

# ----------------------
# 10) Consolidated SAVE
# ----------------------
consolidated <- list(
  models = list(
    original_clean     = model_original,
    constrained_clean  = model_constrained
  ),
  fits = fits,                       # full GenomicSEM fit objects
  tables = list(
    summary = fits_one,              # fit indices + deltas (+ Wald p in base)
    diagnostics_summary = diag_summary_df,
    issues = issues_df
  ),
  per_fit_tables = per_fit_tables,   # loadings/theta/factor_corr per fit
  wald = list(
    contrast = "Development~~Internalizing == 0 (unconstrained, base floor)",
    p_value  = wald_p,
    printout = wald_obj_print
  ),
  provenance = provenance
)

saveRDS(consolidated, file.path(out_dir, rds_name))
log_msg("Saved consolidated RDS -> ", file.path(out_dir, rds_name))

# Optionally write a single xlsx with key tables
if (isTRUE(write_xlsx)) {
  ok <- requireNamespace("openxlsx", quietly = TRUE)
  if (ok) {
    wb <- openxlsx::createWorkbook()
    add_sheet <- function(name, df) {
      openxlsx::addWorksheet(wb, name)
      openxlsx::writeData(wb, name, df)
    }
    add_sheet("fit_summary", consolidated$tables$summary)
    add_sheet("diagnostics_summary", consolidated$tables$diagnostics_summary)
    # Only write issues if it's a data.frame of issues; otherwise write the message
    if (is.data.frame(consolidated$tables$issues) && !"message" %in% names(consolidated$tables$issues)) {
      add_sheet("issues", consolidated$tables$issues)
    } else {
      openxlsx::addWorksheet(wb, "issues")
      openxlsx::writeData(wb, "issues",
                          if (is.data.frame(consolidated$tables$issues)) consolidated$tables$issues else
                            data.frame(message = "No flagged issues under current thresholds."))
    }
    # Provenance sheet
    openxlsx::addWorksheet(wb, "provenance")
    prov_txt <- paste(
      "==== SESSION INFO ====",
      paste(consolidated$provenance$session_info, collapse = "\n"),
      "",
      "==== OBS TRAITS ====",
      paste(consolidated$provenance$obs_traits, collapse = ", "),
      "",
      sprintf("==== LDSC SHAPES ====\nS: %s x %s\nV diag length: %s",
              consolidated$provenance$ldsc_shapes$S[1],
              consolidated$provenance$ldsc_shapes$S[2],
              consolidated$provenance$ldsc_shapes$V_diag_len),
      "",
      sprintf("Thresholds -> THRESH_LOAD: %s | THRESH_CORR: %s",
              THRESH_LOAD, THRESH_CORR),
      sprintf("Residual floors: %s", paste(residual_floors, collapse = ", ")),
      sprintf("Guard indicators: %s", paste(guard_indicators, collapse = ", ")),
      sep = "\n"
    )
    openxlsx::writeData(wb, "provenance", prov_txt)
    out_xlsx <- file.path(out_dir, xlsx_name)
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    log_msg("Saved Excel workbook -> ", out_xlsx)
  } else {
    log_msg("Package 'openxlsx' not available; skipping .xlsx export")
  }
}

log_msg("Done. Consolidated outputs written to: ", out_dir)
