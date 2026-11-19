#!/usr/bin/env Rscript
# ============================================================
# master_efa.R
# EFA on LDSC Genetic Covariance — Symptoms + Disorders
# - Preserves all original output paths/filenames
# - Adds CLI: --which symptoms,disorders  (default: both)
#   e.g., Rscript scripts/master_efa.R --which disorders
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(psych)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# -----------------------
# Setup (unchanged)
# -----------------------
dir.create("gsem/gsem_output/symptoms_efa",  recursive = TRUE, showWarnings = FALSE)
dir.create("gsem/gsem_output/disorder_efa",  recursive = TRUE, showWarnings = FALSE)

source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")
load("gsem/gsem_settings/variant_paths.RData")

stamp <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

# -----------------------
# CLI: --which symptoms,disorders
# -----------------------
args <- commandArgs(trailingOnly = TRUE)
which_to_run <- c("symptoms","disorders")
if (length(args)) {
  for (i in seq_along(args)) {
    if (identical(args[[i]], "--which") && i < length(args)) {
      which_to_run <- strsplit(args[[i+1]], ",", fixed = TRUE)[[1]]
    }
  }
}
which_to_run <- intersect(which_to_run, c("symptoms","disorders"))
if (length(which_to_run) == 0) which_to_run <- c("symptoms","disorders")

# -----------------------
# Helpers (shared)
# -----------------------

# Read LDSC hm3 (if present), stabilize S, return list(S_pd, R, trait_names, ldsc_raw)
read_and_prepare <- function(codes) {
  fname <- paste0("mvLD_", paste(codes, collapse = "_"), ".Rds")
  candidate1 <- file.path(folderpaths$gsem.mvLD, fname)
  candidate2 <- file.path("gsem", "mvld", fname)   # fallback if needed
  
  chosen <- if (file.exists(candidate1)) candidate1 else candidate2
  if (!file.exists(chosen)) stop("❌ mvLD file not found: ", candidate1, " (or fallback ", candidate2, ")")
  
  ldsc <- readRDS(chosen)
  if (!is.null(ldsc$hm3)) ldsc <- ldsc$hm3   # keep your original behavior
  trait_names <- colnames(ldsc$S)
  
  S_raw <- as.matrix(ldsc$S)
  S_pd  <- as.matrix(Matrix::nearPD(S_raw, corr = FALSE)$mat)
  S_pd  <- as.matrix(Matrix::nearPD(S_pd,  corr = FALSE)$mat)
  R     <- cov2cor(S_pd)
  
  list(S_pd = S_pd, R = R, trait_names = trait_names, ldsc_raw = ldsc)
}

# Effective N via harmonic mean, using ldsc$N if present, else sumstats_metadata
effective_N <- function(ldsc_raw, trait_names) {
  HM <- function(x) length(x) / sum(1 / x)
  if (!is.null(ldsc_raw$N)) {
    Ns <- as.numeric(ldsc_raw$N); names(Ns) <- trait_names
  } else if (exists("sumstats_metadata")) {
    Ns <- sumstats_metadata$N[match(trait_names, sumstats_metadata$trait)]
  } else stop("No per-trait Ns found (need ldsc$N or sumstats_metadata).")
  Ns <- Ns[is.finite(Ns) & Ns > 0]
  N_eff <- floor(HM(Ns))
  stopifnot(is.numeric(N_eff), length(N_eff) == 1, is.finite(N_eff))
  N_eff
}

# Scree plot (paths provided by caller)
make_and_save_scree <- function(R, title_text, out_path) {
  eig <- eigen(R, symmetric = TRUE)$values
  p_scree <- ggplot(data.frame(Component = seq_along(eig), Eigenvalue = eig),
                    aes(Component, Eigenvalue)) +
    geom_point() + geom_line() + geom_hline(yintercept = 1, color = "red") +
    labs(x = "Component", y = "Eigenvalue", title = title_text) +
    theme(text = element_text(size = 18)) +
    embarktools::embark_theme_a
  ggsave(out_path, p_scree, width = 10, height = 10, units = "in")
  p_scree
}

# Parallel analysis (paths provided by caller)
make_and_save_parallel <- function(R, N_eff, title_text, out_path) {
  pa <- psych::fa.parallel(R, n.obs = N_eff, fa = "fa", fm = "ml", plot = FALSE)
  df_pa <- data.frame(
    Component = seq_along(pa$fa.values),
    Observed  = pa$fa.values,
    Random    = pa$fa.sim
  ) |>
    pivot_longer(cols = c(Observed, Random), names_to = "Series", values_to = "Eigenvalue")
  
  p_pa <- ggplot(df_pa, aes(Component, Eigenvalue, group = Series, color = Series)) +
    geom_point() +
    geom_line(size = 1) +
    labs(title = title_text, x = "Component", y = "Eigenvalue") +
    scale_color_manual(values = c("Observed" = "#1A4F66", "Random" = "#C2B824")) +
    theme(text = element_text(size = 18)) +
    embarktools::embark_theme_a +
    coord_cartesian(ylim = c(0, NA))
  ggsave(out_path, p_pa, width = 10, height = 10, units = "in")
  p_pa
}

# Fit EFAs (psych::fa 1–5) + simplified table and variance summary
efa_fit_and_table <- function(R, N_eff, out_csv_path) {
  fit_psych <- lapply(1:5, function(f)
    psych::fa(R, nfactors = f, fm = "ml", rotate = "oblimin", n.obs = N_eff)
  )
  fit_tab_psych <- data.frame(
    Factors = 1:5,
    ChiSq   = sapply(fit_psych, \(x) unname(x$STATISTIC)),
    df      = sapply(fit_psych,  \(x) unname(x$dof)),
    TLI     = sapply(fit_psych,  \(x) x$TLI),
    BIC     = sapply(fit_psych,  \(x) x$BIC)
  )
  
  var_from_loadings <- function(fit) {
    L <- as.matrix(fit$loadings)
    p <- nrow(L)
    ss <- colSums(L^2, na.rm = TRUE)
    cum_var <- sum(ss) / p
    data.frame(CumVar = cum_var)
  }
  variance_list <- lapply(fit_psych, var_from_loadings)
  variance_summary <- do.call(rbind, variance_list)
  variance_summary$Factors <- seq_len(nrow(variance_summary))
  variance_summary <- variance_summary |>
    mutate(
      DeltaVar = c(NA, diff(CumVar) * 100),
      CumVar  = round(CumVar, 3),
      DeltaVar = round(DeltaVar, 1)
    ) |>
    select(Factors, CumVar, DeltaVar)
  
  fit_tab_full <- left_join(fit_tab_psych, variance_summary, by = "Factors") |>
    mutate(
      ChiSq = round(ChiSq, 3),
      TLI   = round(TLI, 3),
      BIC   = round(BIC, 3)
    )
  
  write.csv(fit_tab_full, out_csv_path, row.names = FALSE)
  fit_tab_full
}

# Export loadings via factanal for 1–5 factors (exact file paths provided by caller)
export_factanal_loadings <- function(S_pd, N_eff,
                                     out_f1, out_f2, out_f3, out_f4, out_f5,
                                     out_phi3 = NULL, out_phi4 = NULL,
                                     out_all  = NULL) {
  mk_fit <- function(f) {
    fit <- factanal(covmat = S_pd, factors = f, rotation = "promax", n.obs = N_eff)
    L <- as.matrix(fit$loadings)
    keep <- setdiff(rownames(L), c("SS loadings", "Proportion Var", "Cumulative Var"))
    L <- round(L[keep, , drop = FALSE], 2)
    Phi <- tryCatch(fit$Phi, error = function(e) NULL)
    list(L = L, Phi = Phi, fit = fit)
  }
  
  f1 <- mk_fit(1); L1 <- f1$L; colnames(L1) <- "1 Factor"
  f2 <- mk_fit(2); L2 <- f2$L; colnames(L2) <- c("2 Factor 1", "2 Factor 2")
  f3 <- mk_fit(3); L3 <- f3$L; colnames(L3) <- c("3 Factor 1", "3 Factor 2", "3 Factor 3")
  f4 <- mk_fit(4); L4 <- f4$L; colnames(L4) <- paste0("4 Factor ", 1:4)
  f5 <- mk_fit(5); L5 <- f5$L; colnames(L5) <- paste0("5 Factor ", 1:5)
  
  write.csv(L1, out_f1)
  write.csv(L2, out_f2)
  write.csv(L3, out_f3)
  write.csv(L4, out_f4)
  write.csv(L5, out_f5)
  
  if (!is.null(out_phi3) && !is.null(f3$Phi)) write.csv(round(f3$Phi, 3), out_phi3)
  if (!is.null(out_phi4) && !is.null(f4$Phi)) write.csv(round(f4$Phi, 3), out_phi4)
  
  if (!is.null(out_all)) {
    # symptoms: cbind L1..L5 ; disorders: cbind L1..L4 (match your originals)
    if (nchar(out_all) && grepl("efa_loadings_all_sx\\.csv$", out_all)) {
      efa_all <- cbind(L1, L2, L3, L4, L5)
    } else {
      efa_all <- cbind(L1, L2, L3, L4)
    }
    write.csv(efa_all, out_all)
  }
}

# -----------------------
# SYMPTOMS pipeline (paths preserved)
# -----------------------
run_symptoms <- function() {
  stamp("Running: SYMPTOMS EFA")
  symptoms_fa <- (sumstats %>% filter(Symptoms_GSEM == 1))$code
  prep <- read_and_prepare(symptoms_fa)
  N_eff <- effective_N(prep$ldsc_raw, prep$trait_names)
  message(sprintf("Using effective N (harmonic mean): %s", format(N_eff, big.mark=",")))
  
  # Scree + Parallel (keep filenames the same)
  p_scree <- make_and_save_scree(
    R = prep$R,
    title_text = "Scree plot of eigenvalues\nsymptoms dataset",
    out_path = "gsem/gsem_output/symptoms_efa/Scree_plot_eigenvalues_symptoms.png"
  )
  print(p_scree)
  
  message("== Running parallel analysis (symptoms) ==")
  p_pa <- make_and_save_parallel(
    R = prep$R, N_eff = N_eff,
    title_text = "Parallel analysis: Symptoms dataset",
    out_path = "gsem/gsem_output/symptoms_efa/Parallel_analysis_plot_symptoms.png"
  )
  print(p_pa)
  
  # Fit table (psych::fa), keep filename
  fit_tab_full <- efa_fit_and_table(
    R = prep$R, N_eff = N_eff,
    out_csv_path = "gsem/gsem_output/symptoms_efa/efa_fit_psych_simplified.csv"
  )
  message("== Simplified fit + variance summary (psych::fa) ==")
  print(fit_tab_full)
  
  # Loadings exports (keep filenames)
  export_factanal_loadings(
    S_pd = prep$S_pd, N_eff = N_eff,
    out_f1  = "gsem/gsem_output/symptoms_efa/efa_loadings_1fac_dx.csv",
    out_f2  = "gsem/gsem_output/symptoms_efa/efa_loadings_2fac_dx.csv",
    out_f3  = "gsem/gsem_output/symptoms_efa/efa_loadings_3fac_dx.csv",
    out_f4  = "gsem/gsem_output/symptoms_efa/efa_loadings_4fac_dx.csv",
    out_f5  = "gsem/gsem_output/symptoms_efa/efa_loadings_5fac_dx.csv",
    out_phi3 = "gsem/gsem_output/symptoms_efa/phi_3factor.csv",
    out_phi4 = "gsem/gsem_output/symptoms_efa/phi_4factor.csv",
    out_all  = "gsem/gsem_output/symptoms_efa/efa_loadings_all_sx.csv"
  )
  
  message("== Done (symptoms) ==")
}

# -----------------------
# DISORDERS pipeline (paths preserved)
# -----------------------
run_disorders <- function() {
  stamp("Running: DISORDERS EFA")
  disorders_fa <- (sumstats %>% filter(Disorder_GSEM == 1))$code
  prep <- read_and_prepare(disorders_fa)
  N_eff <- effective_N(prep$ldsc_raw, prep$trait_names)
  message(sprintf("Using effective N (harmonic mean): %s", format(N_eff, big.mark=",")))
  
  # Scree + Parallel (keep filenames the same)
  p_scree <- make_and_save_scree(
    R = prep$R,
    title_text = "Scree plot of eigenvalues\ndisorders dataset",
    out_path = "gsem/gsem_output/disorder_efa/Scree_plot_eigenvalues_disorders.png"
  )
  
  message("== Running parallel analysis (disorders) ==")
  p_pa <- make_and_save_parallel(
    R = prep$R, N_eff = N_eff,
    title_text = "Parallel analysis: Disorders dataset",
    out_path = "gsem/gsem_output/disorder_efa/Parallel_analysis_plot.png"
  )
  
  # Fit table (psych::fa), keep filename
  fit_tab_full <- efa_fit_and_table(
    R = prep$R, N_eff = N_eff,
    out_csv_path = "gsem/gsem_output/disorder_efa/efa_fit_psych_simplified.csv"
  )
  message("== Simplified fit + variance summary (psych::fa) ==")
  print(fit_tab_full)
  
  # Loadings exports (keep filenames)
  export_factanal_loadings(
    S_pd = prep$S_pd, N_eff = N_eff,
    out_f1  = "gsem/gsem_output/disorder_efa/efa_loadings_1fac_dx.csv",
    out_f2  = "gsem/gsem_output/disorder_efa/efa_loadings_2fac_dx.csv",
    out_f3  = "gsem/gsem_output/disorder_efa/efa_loadings_3fac_dx.csv",
    out_f4  = "gsem/gsem_output/disorder_efa/efa_loadings_4fac_dx.csv",
    out_f5  = "gsem/gsem_output/disorder_efa/efa_loadings_5fac_dx.csv",
    out_phi3 = "gsem/gsem_output/disorder_efa/phi_3factor.csv",
    out_phi4 = "gsem/gsem_output/disorder_efa/phi_4factor.csv",
    out_all  = "gsem/gsem_output/disorder_efa/efa_loadings_all_dx.csv"
  )
  
  message("== Done (disorders) ==")
}

# -----------------------
# Execute selection
# -----------------------
if ("symptoms" %in% which_to_run)  run_symptoms()
if ("disorders" %in% which_to_run) run_disorders()
stamp("All done.")
