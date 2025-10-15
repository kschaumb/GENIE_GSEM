setwd("~/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")
library(stringr)
# If fread isn't already available via your packages file:
# library(data.table)
# library(dplyr)

# CLEAN MDD FILES
sumstats <- sumstats %>% filter(cleaning_script == 'Insulin_Clean')

# make a list of traits by pulling a list of all the codes in the data frame
traits <- as.list(sumstats$code)
traits <- lapply(traits, unlist)
traits <- unlist(traits)  # <- ensure it's a simple character vector

sumstats_cleaned <- list()

# make a new column in sumstats called readpath that is the relative path to the sumstats file
sumstats <- sumstats %>%
  mutate(readpath = stringr::str_remove(rawpath, ".*GENIE_GSEM/"))

# make a savepath column that is the relative version of the cleaned path
sumstats <- sumstats %>%
  mutate(cleanedpath = paste0("data/gwas_sumstats/", cleaned_file))

# read raw files
for (trait in traits){
  sumstats_cleaned[[trait]] <- fread(sumstats %>% filter(code == !!trait) %>% pull(readpath))
}


# ---------- Trait-specific cleaning ----------
cols_to_munge <- c('CHROM', 'ID', 'A1', 'A2', 'BETA', 'SE', 'P', 'MAF')

# ---------- Trait-specific cleaning ----------

# 1) Insulin (multi-ancestry)
if ("Insulin" %in% names(sumstats_cleaned)) {
  sumstats_cleaned$Insulin <- sumstats_cleaned$Insulin %>%
    dplyr::filter(stringr::str_detect(rsid, "^rs")) %>%
    dplyr::rename(
      ID    = rsid,
      CHROM = chromosome,
      A1    = effect_allele,
      A2    = other_allele,
      BETA  = beta,
      SE    = standard_error,
      P     = p_value,
      MAF   = effect_allele_frequency
    ) %>%
    dplyr::mutate(
      MAF = pmin(MAF, 1 - MAF)
    ) %>%
    dplyr::select(dplyr::all_of(cols_to_munge))
}

# 2) Insulin_Eur (EUR-only GWAS file)
if ("Insulin_Eur" %in% names(sumstats_cleaned)) {
  sumstats_cleaned$Insulin_Eur <- sumstats_cleaned$Insulin_Eur %>%
    dplyr::filter(stringr::str_detect(variant, "^rs")) %>%   # <- use 'variant'
    dplyr::rename(
      ID    = variant,
      CHROM = chromosome,
      A1    = effect_allele,
      A2    = other_allele,
      BETA  = beta,
      SE    = standard_error,
      P     = p_value,
      MAF   = effect_allele_frequency
    ) %>%
    dplyr::mutate(
      MAF = pmin(as.numeric(MAF), 1 - as.numeric(MAF))       # keep as minor allele freq
    ) %>%
    dplyr::select(dplyr::all_of(cols_to_munge))
}


# ---------- Save cleaned files (gzip) ----------
for (trait in traits){
  cleaned_path <- sumstats %>% filter(code == trait) %>% pull(cleanedpath)
  
  if (length(cleaned_path) != 1 || is.na(cleaned_path)) {
    stop(paste("Error: Invalid cleaned file path for trait:", trait))
  }
  
  write.table(sumstats_cleaned[[trait]], gzfile(cleaned_path), 
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
}
