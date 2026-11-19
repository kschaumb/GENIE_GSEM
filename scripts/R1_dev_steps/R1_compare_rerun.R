# ============================================================
# R&R (R1): compare results before vs after the corrected-ANX-N rerun.
#  (1) LDSC: h2 and genetic correlations, gsem/mvLD vs gsem/mvLD_R1
#  (2) CFA standardized estimates: snapshot taken before rerun vs new params_R1
# Output: gsem/gsem_output/R1_dev_steps/ANX_N_rerun_comparison/*.csv (+ console summary)
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(readr) })
out <- "gsem/gsem_output/R1_dev_steps/ANX_N_rerun_comparison"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

rg_of <- function(S) { D <- diag(1 / sqrt(diag(S))); R <- D %*% S %*% D; dimnames(R) <- dimnames(S); R }
h2_out <- list()
for (f in list.files("gsem/mvLD_R1", pattern = "^mvLD_.*\\.Rds$")) {
  old_f <- file.path("gsem/mvLD", f)
  if (!file.exists(old_f)) { message("No original for ", f); next }
  o <- readRDS(old_f)$hm3; n <- readRDS(file.path("gsem/mvLD_R1", f))$hm3
  rownames(o$S) <- colnames(o$S); rownames(n$S) <- colnames(n$S)
  dR <- abs(rg_of(n$S) - rg_of(o$S))
  idx <- which(dR == max(dR, na.rm = TRUE), arr.ind = TRUE)[1, ]
  message(sprintf("%s\n   traits=%d  max |delta rg| = %.4f (%s, %s)", substr(f, 1, 70), ncol(n$S),
                  max(dR, na.rm = TRUE), colnames(n$S)[idx[1]], colnames(n$S)[idx[2]]))
  h2_out[[f]] <- data.frame(file = f, code = colnames(n$S), h2_old = diag(o$S), h2_new = diag(n$S),
                            int_old = diag(o$I), int_new = diag(n$I))
  write.csv(round(rg_of(n$S) - rg_of(o$S), 4), file.path(out, paste0("delta_rg_", ncol(n$S), "traits.csv")))
}
h2 <- bind_rows(h2_out) %>% mutate(ratio = h2_new / h2_old)
write_csv(h2, file.path(out, "h2_old_vs_new.csv"))
message("\nh2 changes (|ratio-1| > .001):"); print(h2 %>% filter(abs(ratio - 1) > .001) %>% distinct(code, .keep_all = TRUE), digits = 4)

for (kind in c("disorder_cfa_R1", "symptoms_cfa_R1")) {
  pre <- file.path("gsem/gsem_output/R1_dev_steps/ANX_N_rerun_comparison/pre_rerun", kind)
  post <- file.path("gsem/gsem_output", kind, "params_R1")
  for (p in list.files(pre, pattern = "_params\\.csv$")) {
    a <- read_csv(file.path(pre, p), show_col_types = FALSE); b <- read_csv(file.path(post, p), show_col_types = FALSE)
    est <- intersect(c("STD_Genotype", "STD_All", "std_all"), names(a))[1]
    m <- inner_join(a, b, by = c("lhs", "op", "rhs"), suffix = c("_pre", "_post")) %>%
      transmute(lhs, op, rhs, pre = .data[[paste0(est, "_pre")]], post = .data[[paste0(est, "_post")]]) %>%
      mutate(across(c(pre, post), ~ suppressWarnings(as.numeric(.))), diff = post - pre)
    write_csv(m, file.path(out, paste0(kind, "__", p)))
    message(sprintf("%s / %s: max |delta std est| = %.4f", kind, p, max(abs(m$diff), na.rm = TRUE)))
  }
}
message("\nSaved comparison files to ", out)
