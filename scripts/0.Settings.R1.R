# ============================================================
# R&R (R1) version of 0.Settings.R
#  - Run from the GENIE_GSEM project root (paths are RELATIVE, so the old
#    absolute Box path no longer matters).
#  - Reads GWAS_sumstats_settings_R1.xlsx (see 0.Settings_corrections.R1.R).
#  - Saves settings to gsem/gsem_settings_R1/ (original gsem_settings untouched).
#  - LDSC output goes to gsem/mvLD_R1/.
#  - ANX is re-munged with the corrected N into data/gwas_sumstats/munged/hm3_munged_R1/;
#    all other traits use their existing munged files (their N did not change).
# ============================================================
source("scripts/0.Packages.R")

sub_paths <- c("data", "data/gwas_sumstats",
               "data/gwas_sumstats/munged/hm3_munged_info_0_8_maf_0_1",
               "data/gwas_sumstats/munged/hm3_munged_R1",
               "data/ld_scores/eur_w_ld_chr", "data/variant_lists",
               "scripts", "gsem", "gsem/mvLD_R1", "gsem/mvLD_R1/log_files",
               "gsem/gsem_output", "gsem/gsem_output/log_files", "gsem/gsem_settings_R1")
for (p in sub_paths) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
folderpaths <- setNames(as.list(sub_paths), gsub("/", ".", sub_paths))
folderpaths$gsem.mvLD <- "gsem/mvLD_R1"   # downstream scripts read folderpaths$gsem.mvLD

vl <- list.files(folderpaths$data.variant_lists)
variant_lists <- as.list(setNames(file.path(folderpaths$data.variant_lists, vl), vl))

sumstats <- readxl::read_excel(file.path(folderpaths$data.gwas_sumstats, "GWAS_sumstats_settings_R1.xlsx"))
sumstats$mungedpath.hm3 <- file.path(folderpaths$data.gwas_sumstats.munged.hm3_munged_info_0_8_maf_0_1,
                                     paste0(sumstats$code, ".sumstats.gz"))
REMUNGE <- c("ANX")
sumstats$mungedpath.hm3[sumstats$code %in% REMUNGE] <-
  file.path(folderpaths$data.gwas_sumstats.munged.hm3_munged_R1, paste0(REMUNGE, ".sumstats.gz"))
sumstats$rawpath <- file.path(folderpaths$data.gwas_sumstats, sumstats$raw_file)

save(sumstats,      file = "gsem/gsem_settings_R1/sumstats_metadata.RData")
save(folderpaths,   file = "gsem/gsem_settings_R1/folderpaths.RData")
save(variant_lists, file = "gsem/gsem_settings_R1/variant_paths.RData")
message("Saved R1 settings to gsem/gsem_settings_R1/")
