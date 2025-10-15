# ============================================================
# Genomic multiple regression: BMI-focused table (2dp, no bold)
# - Keeps trios where a predictor is BMI or BMI_Childhood
# - Shows bivariate focal effect and BMI/BMI_Childhood-partialled effect
# - Adds Model R², semi-partial R², and BH-FDR (within AN & BEBROAD)
# - Outputs: numeric CSV + wide HTML (no bold; all numeric 2dp, p-values 3dp/<.001)
# ============================================================

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(dplyr)
  library(readr)
  library(purrr)
  library(gt)
})

# ---------- paths ----------
base_dir  <- "~/Library/CloudStorage/Box-Box/GENIE_GSEM"
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")
# Filter for European-ancestry traits
all <- sumstats %>% 
  filter(Symptom_RGs == 1 | Disorder_RGs == 1)
all <- all$code

# Dynamically construct the expected mvLD filename
# (Assumes consistent naming scheme: "mvLD_<traits joined by _>.Rds")
ldsc_path <- file.path(folderpaths$gsem.mvLD,
                       paste0("mvLD_", paste(all, collapse = "_"), ".Rds")
)
# Confirm the path and read
if (!file.exists(ldsc_path)) {
  stop("❌ mvLD file not found: ", ldsc_path)
}

ldsc_rds <- readRDS(ldsc_path)


trio_path <- file.path(base_dir, "gsem/gsem_output/aim3_conditional/aim3_trait_trios.csv")
out_dir   <- file.path(base_dir, "gsem/gsem_output/aim3_conditional/gmr_trios")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_out  <- file.path(out_dir, "gmr_trios_BMIfocused.csv")
html_out <- file.path(out_dir, "gmr_trios_BMIfocused.html")

# ---------- load LDSC cov structure ----------
ldsc <- readRDS(ldsc_path)
S <- ldsc$hm3$S
rownames(S) <- colnames(S)

# ---------- read trios (ED=Y, Anthro=X1, Q=X2) ----------
trios <- read_csv(trio_path, show_col_types = FALSE) |>
  select(ED, Anthro, Q) |>
  rename(Y = ED, X1 = Anthro, X2 = Q) |>
  filter(Y %in% colnames(S), X1 %in% colnames(S), X2 %in% colnames(S))
stopifnot(nrow(trios) > 0)

# keep only those where a predictor is BMI or BMI_Childhood
BMI_SET <- c("BMI", "BMI_Childhood")
trios_bmi <- trios |>
  filter(X1 %in% BMI_SET | X2 %in% BMI_SET) |>
  mutate(
    BMI_var   = ifelse(X1 %in% BMI_SET, X1, X2),
    Focal_var = ifelse(X1 %in% BMI_SET, X2, X1)
  )
stopifnot(nrow(trios_bmi) > 0)

# ---------- helpers (formatting: 2dp for numbers, 3dp/<.001 for p) ----------
fmt_num2 <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.2f", x))
fmt_p    <- function(p) ifelse(is.na(p), NA_character_, ifelse(p < .001, "<.001", sprintf("%.3f", p)))

# ---------- run one trio ----------
run_one <- function(Y, Focal, BMIv, ldsc, S) {
  # Bivariate: Y ~ Focal
  fit_biv <- GenomicSEM::usermodel(
    covstruc = ldsc$hm3,
    model    = sprintf("%s ~ %s", Y, Focal),
    std.lv   = FALSE
  )
  pe_biv <- fit_biv$results
  biv <- pe_biv |>
    filter(lhs == Y, op == "~", rhs == Focal) |>
    transmute(b_raw = Unstand_Est, se_raw = Unstand_SE, p_raw = p_value,
              beta_raw_std = STD_Genotype)
  
  # Trio (BMI-controlled): Y ~ Focal + BMIv
  fit_tri <- GenomicSEM::usermodel(
    covstruc = ldsc$hm3,
    model    = sprintf("%s ~ %s + %s", Y, Focal, BMIv),
    std.lv   = FALSE
  )
  pe_tri <- fit_tri$results
  
  foc <- pe_tri |>
    filter(lhs == Y, op == "~", rhs == Focal) |>
    transmute(b_adj = Unstand_Est, se_adj = Unstand_SE, p_adj = p_value,
              beta_adj_std = STD_Genotype)
  
  # Model R² from residual variance; semi-partial R² for focal
  resid_var <- pe_tri |> filter(lhs == Y, op == "~~", rhs == Y) |> pull(Unstand_Est)
  R2_tri <- as.numeric(1 - resid_var / S[Y, Y])
  semi_R2 <- as.numeric(foc$b_adj * S[Focal, Y] / S[Y, Y])
  
  tibble(
    Outcome = Y,
    Focal   = Focal,
    Controlled_for = BMIv,
    b_raw = biv$b_raw,  se_raw = biv$se_raw,  p_raw = biv$p_raw,
    b_adj = foc$b_adj,  se_adj = foc$se_adj,  p_adj = foc$p_adj,
    beta_raw_std = biv$beta_raw_std,
    beta_adj_std = foc$beta_adj_std,
    semi_R2 = semi_R2,
    Model_R2 = R2_tri
  )
}

# ---------- run and collect ----------
res <- pmap_dfr(
  list(trios_bmi$Y, trios_bmi$Focal_var, trios_bmi$BMI_var),
  ~run_one(..1, ..2, ..3, ldsc, S)
)

# ---------- add FDR p-values (within outcomes AN and BEBROAD) ----------
res <- res |>
  group_by(Outcome) |>
  mutate(
    p_raw_FDR = if_else(Outcome %in% c("AN","ANR", "BEBROAD", "BE_noAN"),
                        p.adjust(p_raw, method = "BH"),
                        as.numeric(NA)),
    p_adj_FDR = if_else(Outcome %in% c("AN","ANR", "BEBROAD", "BE_noAN"),
                        p.adjust(p_adj, method = "BH"),
                        as.numeric(NA))
  ) |>
  ungroup()

# ---------- save numeric CSV ----------
write_csv(res, csv_out)

# ---------- build formatted WIDE HTML table ----------
# Pre-compute strings with 2dp and 95% CIs
res_disp <- res |>
  mutate(
    b_raw_ci = paste0(fmt_num2(b_raw), " [", fmt_num2(b_raw - 1.96*se_raw), ", ", fmt_num2(b_raw + 1.96*se_raw), "]"),
    b_adj_ci = paste0(fmt_num2(b_adj), " [", fmt_num2(b_adj - 1.96*se_adj), ", ", fmt_num2(b_adj + 1.96*se_adj), "]"),
    p_raw_FDR_txt = fmt_p(p_raw_FDR),
    p_adj_FDR_txt = fmt_p(p_adj_FDR)
  )

# One bivariate line per (Outcome, Focal)
biv <- res_disp |>
  group_by(Outcome, Focal) |>
  slice(1) |>
  ungroup() |>
  transmute(
    Outcome, Focal,
    `Bivariate b [95% CI]` = b_raw_ci,
    `Bivariate p (FDR)`    = p_raw_FDR_txt,
    `β_raw (std-geno)`     = fmt_num2(beta_raw_std),
    `semi-R²`              = fmt_num2(semi_R2),
    `Model R²`             = fmt_num2(Model_R2)
  )

# Adjusted effects wide under spanners
adj_wide <- res_disp |>
  select(Outcome, Focal, Controlled_for, b_adj_ci, p_adj_FDR_txt, beta_adj_std) |>
  mutate(
    b_col  = ifelse(Controlled_for == "BMI", "BMI: b [95% CI]", "BMI_Childhood: b [95% CI]"),
    p_col  = ifelse(Controlled_for == "BMI", "BMI: p (FDR)",   "BMI_Childhood: p (FDR)"),
    bt_col = ifelse(Controlled_for == "BMI", "BMI: β (std-geno)", "BMI_Childhood: β (std-geno)"),
    beta_fmt = fmt_num2(beta_adj_std)
  ) |>
  select(Outcome, Focal, b_col, p_col, bt_col, b_adj_ci, p_adj_FDR_txt, beta_fmt)

adj_b <- adj_wide |>
  select(Outcome, Focal, b_col, b_adj_ci) |>
  tidyr::pivot_wider(names_from = b_col, values_from = b_adj_ci)

adj_p <- adj_wide |>
  select(Outcome, Focal, p_col, p_adj_FDR_txt) |>
  tidyr::pivot_wider(names_from = p_col, values_from = p_adj_FDR_txt)

adj_beta <- adj_wide |>
  select(Outcome, Focal, bt_col, beta_fmt) |>
  tidyr::pivot_wider(names_from = bt_col, values_from = beta_fmt)

tbl_wide <- biv |>
  left_join(adj_b,    by = c("Outcome","Focal")) |>
  left_join(adj_p,    by = c("Outcome","Focal")) |>
  left_join(adj_beta, by = c("Outcome","Focal"))

# Order rows by strongest BMI-controlled absolute effect within outcome
order_helper <- res |>
  group_by(Outcome, Focal) |>
  summarize(max_abs_b_adj = max(abs(b_adj), na.rm = TRUE), .groups = "drop")
tbl_wide <- tbl_wide |>
  left_join(order_helper, by = c("Outcome","Focal")) |>
  arrange(Outcome, desc(max_abs_b_adj)) |>
  select(-max_abs_b_adj)

# Build gt table (no bolding)
gt_tbl <- tbl_wide |>
  gt(groupname_col = "Outcome", rowname_col = "Focal") |>
  tab_header(
    title = md("**Genomic multiple regression: focal effects (BMI-controlled)**"),
    subtitle = "Bivariate vs. BMI/BMI_Childhood-partialled; BH-FDR within AN, ANR, BEBROAD, and BE_noAN"
  ) |>
  tab_spanner(label = "Controlled for BMI",
              columns = c(`BMI: b [95% CI]`, `BMI: p (FDR)`, `BMI: β (std-geno)`)) |>
  tab_spanner(label = "Controlled for BMI_Childhood",
              columns = c(`BMI_Childhood: b [95% CI]`, `BMI_Childhood: p (FDR)`, `BMI_Childhood: β (std-geno)`)) |>
  cols_label(
    `Bivariate b [95% CI]` = md("*b* [95% CI]"),
    `Bivariate p (FDR)`    = md("*p* (FDR)"),
    `β_raw (std-geno)`     = md("&beta;<sub>raw</sub>"),
    `BMI: b [95% CI]`      = md("*b* [95% CI]"),
    `BMI: p (FDR)`         = md("*p* (FDR)"),
    `BMI: β (std-geno)`    = md("&beta;"),
    `BMI_Childhood: b [95% CI]` = md("*b* [95% CI]"),
    `BMI_Childhood: p (FDR)`    = md("*p* (FDR)"),
    `BMI_Childhood: β (std-geno)`= md("&beta;"),
    `semi-R²` = md("sr&sup2;"),
    `Model R²` = md("R&sup2;")
  ) |>
  cols_align(align = "center") |>
  opt_table_outline()

gtsave(gt_tbl, html_out)

message("✅ Wrote:\n- ", csv_out, "\n- ", html_out)
