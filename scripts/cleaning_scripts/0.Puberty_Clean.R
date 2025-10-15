
setwd("~/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

# CLEAN MDD FILES
sumstats <- sumstats %>% filter(cleaning_script == 'Puberty_Clean')
# make a list of traits by pulling a list of all the codes in the data frame
traits <- as.list(sumstats$code)
traits <- lapply(traits, unlist)

sumstats_cleaned <- list()

# make a new colum in sumstats_daner called readpath that is the relative path to the sumstats file (rawpath but starting at data/ -- use stringr to drop everything from rawpath before data/)
sumstats <- sumstats %>%
  mutate(readpath = stringr::str_remove(rawpath, ".*GENIE_GSEM/"))

# make a savepath column that is the relative version of the cleaned path
sumstats<- sumstats %>%
  mutate(cleanedpath = paste0("data/gwas_sumstats/", cleaned_file))

for (trait in traits){
  sumstats_cleaned[[trait]] <- fread(sumstats %>% filter(code == !!trait) %>% pull(readpath))
}
cols_to_munge <- c('CHROM', 'ID', 'A1', 'A2', 'BETA', 'SE', 'P', 'MAF')
sumstats_cleaned$Age_Menarche <- sumstats_cleaned$Age_Menarche %>%
  filter(str_detect(rsid, "^rs")) %>%  # Keep only rows with valid rsIDs
  rename(
    ID = rsid,
    CHROM = chromosome,
    A1 = effect_allele,
    A2 = other_allele,
    BETA = beta,
    SE = standard_error,
    P = p_value,
    MAF = effect_allele_frequency
  ) %>%
  mutate(
    MAF = pmin(MAF, 1 - MAF)  # Ensure MAF is ≤ 0.5
  ) %>%
  select(all_of(cols_to_munge))

cols_to_munge <- c('CHROM', 'ID', 'A1', 'A2', 'BETA', 'SE', 'P')

sumstats_cleaned$Tanner_StageMidAdol <- sumstats_cleaned$Tanner_StageMidAdol %>%
  filter(str_detect(RSID, "^rs")) %>%  # Keep valid rsIDs only
  rename(
    ID = RSID,
    CHROM = Chromosome,
    A1 = EA,
    A2 = NEA,
    BETA = BETA,
    SE = SE,
    P = P
  ) %>%
  select(all_of(cols_to_munge))

cols_to_munge <- c('CHROM', 'ID', 'A1', 'A2', 'BETA', 'SE', 'P', 'MAF')

sumstats_cleaned$Age_Menarche_Eur <- sumstats_cleaned$Age_Menarche_Eur %>%
    dplyr::filter(!is.na(variant_id), stringr::str_detect(variant_id, "^rs")) %>%
    dplyr::rename(
      ID    = variant_id,
      CHROM = chromosome,
      A1    = effect_allele,
      A2    = other_allele,
      BETA  = beta,
      SE    = standard_error,
      P     = p_value,
      MAF   = effect_allele_frequency
    ) %>%
    dplyr::mutate(
      MAF = pmin(as.numeric(MAF), 1 - as.numeric(MAF))
    ) %>%
    dplyr::select(dplyr::all_of(cols_to_munge))


# save the cleaned sumstats using the sumstats_sub cleanedpath, compressing each to gzip
for (trait in traits){
  cleaned_path <- sumstats %>% filter(code == trait) %>% pull(cleanedpath)
  
  if (length(cleaned_path) != 1 || is.na(cleaned_path)) {
    stop(paste("Error: Invalid cleaned file path for trait:", trait))
  }
  
  write.table(sumstats_cleaned[[trait]], gzfile(cleaned_path), 
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
}
rm(list = ls())

