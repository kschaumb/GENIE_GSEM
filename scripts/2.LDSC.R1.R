# ============================================================
# R&R (R1) version of 2.LDSC.R
#  - Uses R1 settings (corrected ANX N and BE_noAN sample prevalence).
#  - Same trait lists (all / symptoms / disorders) and LD scores as the original.
#  - Output: gsem/mvLD_R1/mvLD_<codes>.Rds (same file names as the originals in gsem/mvLD).
#  - Also writes gsem/mvLD_R1/h2_R1.csv (liability/observed-scale h2, SE, intercept).
# ============================================================
source("scripts/0.Packages.R")
library(GenomicSEM)
load("gsem/gsem_settings_R1/folderpaths.RData")
load("gsem/gsem_settings_R1/sumstats_metadata.RData")

trait_lists <- list(
  all       = sumstats$code[sumstats$Symptom_RGs %in% 1 | sumstats$Disorder_RGs %in% 1],
  symptoms  = sumstats$code[sumstats$Symptom_RGs %in% 1],
  disorders = sumstats$code[sumstats$Disorder_RGs %in% 1]
)

h2_rows <- list()
for (nm in names(trait_lists)) {
  codes <- trait_lists[[nm]]
  sel   <- sumstats[sumstats$code %in% codes, ]
  rds   <- file.path(folderpaths$gsem.mvLD, paste0("mvLD_", paste(sel$code, collapse = "_"), ".Rds"))
  message("== LDSC: ", nm, " (", nrow(sel), " traits) -> ", rds)
  res <- ldsc(traits = sel$mungedpath.hm3,
              sample.prev = suppressWarnings(as.numeric(sel$gsem_samplePrevalence)),
              population.prev = suppressWarnings(as.numeric(sel$populationPrevalence)),
              trait.names = sel$code,
              ld = folderpaths$data.ld_scores.eur_w_ld_chr)
  saveRDS(list(hm3 = res), rds)
  k <- ncol(res$S); se <- sqrt(diag(res$V))[cumsum(c(1, k:2))]   # diagonal elements of vech(S)
  h2_rows[[nm]] <- data.frame(list = nm, code = colnames(res$S), h2 = diag(res$S),
                              h2_se = se, intercept = diag(res$I))
}
write.csv(do.call(rbind, h2_rows), file.path(folderpaths$gsem.mvLD, "h2_R1.csv"), row.names = FALSE)

logs <- list.files(pattern = "_ldsc\\.log$")
if (length(logs)) file.rename(logs, file.path(folderpaths$gsem.mvLD, "log_files", logs))
message("Done.")
