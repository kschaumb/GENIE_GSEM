# ============================================================
# R&R (R1) version of 1.Munge.R — re-munge ANX only, with corrected N (gsem_n ≈ 11,625).
# Same filters as the original (INFO > 0.8, MAF > 0.01, HapMap3).
# Output: data/gwas_sumstats/munged/hm3_munged_R1/ANX.sumstats.gz (original munged file untouched)
# ============================================================
source("scripts/0.Packages.R")
library(GenomicSEM)
load("gsem/gsem_settings_R1/folderpaths.RData")
load("gsem/gsem_settings_R1/sumstats_metadata.RData")
load("gsem/gsem_settings_R1/variant_paths.RData")

to_munge <- sumstats[sumstats$code %in% c("ANX"), ]
to_munge$cleanedpath <- file.path("data/gwas_sumstats", to_munge$cleaned_file)
print(to_munge[, c("code", "cleanedpath", "gsem_n")])

munge(files = to_munge$cleanedpath, trait.name = to_munge$code,
      N = as.numeric(to_munge$gsem_n), info.filter = 0.8, maf.filter = 0.01,
      hm3 = variant_lists$w_hm3.snplist)

out_dir <- folderpaths$data.gwas_sumstats.munged.hm3_munged_R1
for (f in c(paste0(to_munge$code, ".sumstats.gz"), list.files(pattern = "_munge\\.log$")))
  if (file.exists(f)) file.rename(f, file.path(out_dir, basename(f)))
message("Munged files moved to ", out_dir)
