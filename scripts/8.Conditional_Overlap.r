# ============================================================
# Aim 3 (Steps 1–3): Identify candidate trait trios
# ============================================================
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")
load("gsem/gsem_settings/variant_paths.RData")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(readr); library(glue)
})

# ------------------------------
# 0) Load LDSC matrices
# ------------------------------
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

ldsc <- readRDS(ldsc_path)



S <- ldsc$hm3$S
V <- ldsc$hm3$V
colnames(S) <- rownames(S) <- make.names(colnames(S))
stopifnot(identical(rownames(S), colnames(S)))

message("✅ Loaded S (", nrow(S), "×", ncol(S), ") and V (", nrow(V), "×", ncol(V), ")")

# ------------------------------
# 1) Parameters
# ------------------------------
out_dir <- "gsem/gsem_output/aim3_conditional"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Add ANR and BE_noAN here
ED_traits <- c(
  AN      = "AN",
  BEBROAD      = "BEBROAD",
  BE_noAN = "BE_noAN"
)

anthro_adult <- "BMI"
anthro_child <- "BMI_Childhood"
rg_cut       <- 0.25

# Keep only ED traits that are present in the matrix; warn if any are missing
present_ED <- ED_traits[names(ED_traits) %in% names(ED_traits) & ED_traits %in% colnames(S)]
missing_ED <- setdiff(ED_traits, present_ED)
if (length(missing_ED)) {
  message("⚠️ Missing ED traits (not in S): ", paste(missing_ED, collapse = ", "))
}
stopifnot(length(present_ED) >= 1)

# ------------------------------
# 2) Compute correlation matrix
# ------------------------------
cov2cor_named <- function(S) {
  D <- diag(1 / sqrt(diag(S)))
  R <- D %*% S %*% D
  rownames(R) <- rownames(S); colnames(R) <- colnames(S)
  R
}
R <- cov2cor_named(S)

# ------------------------------
# 3) Candidate finder
# ------------------------------


find_candidates <- function(ed_name, R, rg_cut, anthro_adult, anthro_child,
                            ed_blocklist,
                            extra_block = c("BFP", "FFM")) {
  rg_ed   <- sort(R[, ed_name], decreasing = TRUE)
  keep_ed <- names(rg_ed)[abs(rg_ed) >= rg_cut & names(rg_ed) != ed_name]
  
  # Full blocklist: anthro + EDs + extras (+ extra ED phenos to exclude)
  blocklist <- unique(c(anthro_adult, anthro_child, ed_blocklist, extra_block))
  
  tibble(
    trait        = keep_ed,
    rg_with_ED   = R[trait, ed_name],
    rg_with_BMI  = R[trait, anthro_adult],
    rg_with_cBMI = R[trait, anthro_child],
    meets_BMI    = abs(rg_with_BMI)  >= rg_cut,
    meets_cBMI   = abs(rg_with_cBMI) >= rg_cut
  ) %>%
    mutate(meets_any_anthro = meets_BMI | meets_cBMI) %>%
    filter(meets_any_anthro) %>%
    filter(!(trait %in% blocklist)) %>%
    arrange(desc(abs(rg_with_ED)))
}

# ------------------------------
# 4) Run for all ED traits (no per-ED files; one combined CSV)
# ------------------------------
cand_list <- imap(present_ED, ~{
  ed_label <- .y    # e.g., "AN","ANR","BE","BE_noAN"
  ed_code  <- .x    # column name in S/R
  if (!ed_code %in% colnames(R)) {
    message("⚠️ Skipping ", ed_label, " (", ed_code, ") — not in R")
    return(NULL)
  }
  cand <- find_candidates(
    ed_name          = ed_code,
    R                = R,
    rg_cut           = rg_cut,
    anthro_adult     = anthro_adult,
    anthro_child     = anthro_child,
    ed_blocklist     = unname(present_ED),
    extra_block      = c("BFP", "FFM")  )
  if (is.null(cand) || nrow(cand) == 0) {
    message(glue("ℹ️ No candidates for {ed_label}"))
    return(NULL)
  }
  cand %>% mutate(ED = ed_label, ED_code = ed_code, AnthroHit = case_when(
    meets_BMI & meets_cBMI ~ "BMI & cBMI",
    meets_BMI             ~ "BMI",
    meets_cBMI            ~ "cBMI",
    TRUE                  ~ NA_character_
  ))
})

# Combine and write one candidates file
candidates_all <- cand_list %>%
  discard(is.null) %>%
  bind_rows()

cand_file <- file.path(out_dir, "aim3_candidates_all.csv")
if (nrow(candidates_all) > 0) {
  readr::write_csv(candidates_all, cand_file)
  message(glue("✅ Combined candidates saved: {basename(cand_file)} (n={nrow(candidates_all)})"))
} else {
  message("⚠️ No candidates found across ED traits.")
}

# ------------------------------
# 5) Build trio table for next step (also excludes ANBP/BENARROW)
# ------------------------------
make_trios <- function(cand_df_for_one_ed, ed_name, anthro_adult, anthro_child) {
  adult <- cand_df_for_one_ed %>% filter(meets_BMI)  %>% mutate(Anthro = anthro_adult)
  child <- cand_df_for_one_ed %>% filter(meets_cBMI) %>% mutate(Anthro = anthro_child)
  bind_rows(adult, child) %>%
    distinct(trait, Anthro) %>%
    mutate(ED = ed_name, Q = trait) %>%
    # ensure Q never includes excluded ED phenos
    select(ED, Anthro, Q)
}

if (nrow(candidates_all) > 0) {
  trio_list <- candidates_all %>%
    group_split(ED, .keep = TRUE) %>%
    map(~ make_trios(.x, unique(.x$ED), anthro_adult, anthro_child))
  
  trait_trios <- bind_rows(trio_list)
} else {
  trait_trios <- tibble(ED = character(), Anthro = character(), Q = character())
}

out_trio_file <- file.path(out_dir, "aim3_trait_trios.csv")
readr::write_csv(trait_trios, out_trio_file)

# ------------------------------
# 6) Console preview
# ------------------------------
cat("\n✅ Saved combined candidates to:", cand_file, "\n")
cat("✅ Saved trio table to:", out_trio_file, "\n\n")
print(trait_trios, n = 20)

