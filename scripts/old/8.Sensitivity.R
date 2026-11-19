# ============================================================
# Aim 3 – Sensitivity analysis
# Genomic multiple regression: BMI ~ AN + OCD
# ------------------------------------------------------------
# - Uses the same mvLD (hm3) structure as other Aim 3 scripts
# - Outcome: BMI
# - Predictors: AN and OCD
# - Outputs: one CSV + RDS with unstandardized + std-geno betas & model R²
# ============================================================

suppressPackageStartupMessages({
  source("scripts/0.Packages.R")
  library(GenomicSEM)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ---------- paths ----------
base_dir  <- "~/Library/CloudStorage/Box-Box/GENIE_GSEM"
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

# Filter for European-ancestry traits (same as other scripts)
all <- sumstats %>%
  filter(Symptom_RGs == 1 | Disorder_RGs == 1) %>%
  pull(code)

ldsc_path <- file.path(
  folderpaths$gsem.mvLD,
  paste0("mvLD_", paste(all, collapse = "_"), ".Rds")
)

if (!file.exists(ldsc_path)) {
  stop("❌ mvLD file not found: ", ldsc_path)
}

ldsc <- readRDS(ldsc_path)
covstruc <- ldsc$hm3
S <- covstruc$S
rownames(S) <- colnames(S)

message("✅ Loaded covstruc with ", ncol(S), " traits")

# ---------- check required traits ----------
needed <- c("BMI", "AN", "OCD")
if (!all(needed %in% colnames(S))) {
  missing <- needed[!needed %in% colnames(S)]
  stop("❌ Missing required traits in S: ", paste(missing, collapse = ", "))
}

# ---------- fit model: BMI ~ AN + OCD ----------
model_BMI <- "BMI ~ AN + OCD"

fit_BMI <- GenomicSEM::usermodel(
  covstruc = covstruc,
  model    = model_BMI,
  std.lv   = FALSE
)

pe <- fit_BMI$results

# Extract regression paths into BMI
betas <- pe %>%
  filter(lhs == "BMI", op == "~", rhs %in% c("AN", "OCD")) %>%
  transmute(
    outcome      = lhs,
    predictor    = rhs,
    b            = Unstand_Est,
    se           = Unstand_SE,
    z            = Z,
    p            = p_value,
    beta_std_gen = STD_Genotype
  )

# Residual variance & model R² for BMI
resid_var <- pe %>%
  filter(lhs == "BMI", op == "~~", rhs == "BMI") %>%
  pull(Unstand_Est)

R2_BMI <- as.numeric(1 - resid_var / S["BMI", "BMI"])

sens_tbl <- betas %>%
  mutate(model_R2 = R2_BMI)

# ---------- save outputs ----------
out_dir  <- file.path(base_dir, "gsem/gsem_output/aim3_conditional")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_out <- file.path(out_dir, "aim3_sensitivity_BMI_AN_OCD.csv")
rds_out <- file.path(out_dir, "aim3_sensitivity_BMI_AN_OCD.Rds")

write_csv(sens_tbl, csv_out)

saveRDS(
  list(
    fit     = fit_BMI,
    results = pe,
    summary = sens_tbl,
    model   = model_BMI
  ),
  rds_out
)

cat("\n✅ Sensitivity model: BMI ~ AN + OCD\n")
print(sens_tbl)
cat("\nResults written to:\n  ", csv_out, "\n  ", rds_out, "\n")
