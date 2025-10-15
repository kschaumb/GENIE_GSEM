# ===========================
# EFA on LDSC Genetic Covariance + Parallel Analysis + Simplified Fit Table
# ===========================

# --- Setup & packages ---
options(stringsAsFactors = FALSE)

dir.create("gsem/gsem_output/symptoms_efa", recursive = TRUE, showWarnings = FALSE)

source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")
load("gsem/gsem_settings/variant_paths.RData")

library(Matrix)
library(psych)
library(dplyr)
library(ggplot2)
library(tidyr)

# --- Load LDSC object (no $hm3 assumption) ---
# Filter for European-ancestry traits
symptoms_fa <- sumstats %>% 
  filter(Symptoms_GSEM == 1)
symptoms_fa <- symptoms_fa$code

# Dynamically construct the expected mvLD filename
# (Assumes consistent naming scheme: "mvLD_<traits joined by _>.Rds")
ldsc_path <- file.path(folderpaths$gsem.mvLD,
  paste0("mvLD_", paste(symptoms_fa, collapse = "_"), ".Rds")
)

# Confirm the path and read
if (!file.exists(ldsc_path)) {
  stop("❌ mvLD file not found: ", ldsc_path)
}

ldsc <- readRDS(ldsc_path)

ldsc <- ldsc$hm3
trait_names <- colnames(ldsc$S)

# --- Stabilize S (positive-definite) ---
S_raw <- as.matrix(ldsc$S)
S_pd  <- as.matrix(Matrix::nearPD(S_raw, corr = FALSE)$mat)
S_pd  <- as.matrix(Matrix::nearPD(S_pd,  corr = FALSE)$mat)
R     <- cov2cor(S_pd)

# --- Effective N (harmonic mean across per-trait Ns) ---
HM <- function(x) length(x) / sum(1 / x)
if (!is.null(ldsc$N)) {
  Ns <- as.numeric(ldsc$N); names(Ns) <- trait_names
} else if (exists("sumstats_metadata")) {
  Ns <- sumstats_metadata$N[match(trait_names, sumstats_metadata$trait)]
} else stop("No per-trait Ns found (need ldsc$N or sumstats_metadata).")
Ns <- Ns[is.finite(Ns) & Ns > 0]
N_eff <- floor(HM(Ns))
stopifnot(is.numeric(N_eff), length(N_eff) == 1, is.finite(N_eff))
message(sprintf("Using effective N (harmonic mean): %s", format(N_eff, big.mark=",")))

# --- Scree plot ---
eig <- eigen(R, symmetric = TRUE)$values
p_scree <- ggplot(data.frame(Component = seq_along(eig), Eigenvalue = eig),
                  aes(Component, Eigenvalue)) +
  geom_point() + geom_line() + geom_hline(yintercept = 1, color = "red") +
  labs(x = "Component", y = "Eigenvalue",
       title = "Scree plot of eigenvalues\ndisorders dataset") +
  theme(text = element_text(size = 18)) +
  embarktools::embark_theme_a
ggsave("gsem/gsem_output/symptoms_efa/Scree_plot_eigenvalues_symptoms.png",
       p_scree, width = 10, height = 10, units = "in")
print(p_scree)
# --- Parallel analysis (no Kaiser line) ---
message("== Running parallel analysis ==")
pa <- psych::fa.parallel(
  R, n.obs = N_eff, fa = "fa", fm = "ml",
  plot = FALSE
)
df_pa <- data.frame(
  Component = seq_along(pa$fa.values),
  Observed  = pa$fa.values,
  Random    = pa$fa.sim
) |>
  pivot_longer(cols = c(Observed, Random), names_to = "Series", values_to = "Eigenvalue")

p_pa <- ggplot(df_pa, aes(Component, Eigenvalue, group = Series, color = Series)) +
  geom_point() +
  geom_line(size = 1) +
  labs(title = "Parallel analysis: Disorders dataset",
       x = "Component", y = "Eigenvalue") +
  scale_color_manual(values = c("Observed" = "#1A4F66", "Random" = "#C2B824")) +
  theme(text = element_text(size = 18)) +
  embarktools::embark_theme_a +
  coord_cartesian(ylim = c(0, NA))
ggsave("gsem/gsem_output/symptoms_efa/Parallel_analysis_plot_symptoms.png",
       p_pa, width = 10, height = 10, units = "in")
print(p_pa)
# --- Fit 1–5 factor EFAs (TLI, BIC; no RMSEA/p) ---
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

# --- Variance explained from loadings ---
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

# --- Merge and finalize table ---
fit_tab_full <- left_join(fit_tab_psych, variance_summary, by = "Factors") |>
  mutate(
    ChiSq = round(ChiSq, 3),
    TLI = round(TLI, 3),
    BIC = round(BIC, 3)
  )

# Save clean table
write.csv(fit_tab_full,
          "gsem/gsem_output/symptoms_efa/efa_fit_psych_simplified.csv",
          row.names = FALSE)

message("== Simplified fit + variance summary (psych::fa) ==")
print(fit_tab_full)

# --- Optional: loadings export (1–4 factors) ---
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

write.csv(L1, "gsem/gsem_output/symptoms_efa/efa_loadings_1fac_dx.csv")
write.csv(L2, "gsem/gsem_output/symptoms_efa/efa_loadings_2fac_dx.csv")
write.csv(L3, "gsem/gsem_output/symptoms_efa/efa_loadings_3fac_dx.csv")
write.csv(L4, "gsem/gsem_output/symptoms_efa/efa_loadings_4fac_dx.csv")
write.csv(L5, "gsem/gsem_output/symptoms_efa/efa_loadings_5fac_dx.csv")

if (!is.null(f3$Phi)) write.csv(round(f3$Phi, 3), "gsem/gsem_output/symptoms_efa/phi_3factor.csv")
if (!is.null(f4$Phi)) write.csv(round(f4$Phi, 3), "gsem/gsem_output/symptoms_efa/phi_4factor.csv")

efa_all <- cbind(L1, L2, L3, L4, L5)
write.csv(efa_all, "gsem/gsem_output/symptoms_efa/efa_loadings_all_sx.csv")

message("== Done ==")
