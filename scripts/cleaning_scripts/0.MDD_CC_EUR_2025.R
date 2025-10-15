
#setwd("/Volumes/kschaumberg/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

# CLEAN MDD FILES
sumstats <- sumstats %>% filter(cleaning_script == 'MDD_CC_EUR_2025')
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

cols_to_munge <- c('CHR', #Chromosome number
                   'SNP', # SNP
                   'A1', # Alelle 1 (OR on this alelle)
                   'A2', # Allele 2
                   'MAF', # MAF frequency 
                   'INFO', # imputation quality
                   'OR', # odds ratio ()
                   'SE', 
                   'P', # P-value
                   'Neff_half' # effective sample size/2 -- LDSC will double it
)

# 3) MDD_Eur
if ("MDD_Eur" %in% names(sumstats_cleaned)) {
  sumstats_cleaned$MDD_Eur <- sumstats_cleaned$MDD_Eur %>%
    dplyr::filter(stringr::str_detect(ID, "^rs")) %>%
    dplyr::rename(
      CHR  = `#CHROM`,
      SNP  = ID,
      A1   = EA,
      A2   = NEA,
      P    = PVAL,
      INFO = IMPINFO
    ) %>%
    dplyr::mutate(
      # derive Ncontrol from Neff and Ncase: 4/Neff = 1/Ncase + 1/Ncontrol
      Ncontrol = 1 / ((4/NEFF) - (1/NCAS)),
      EAF = (FCAS * NCAS + FCON * Ncontrol) / (NCAS + Ncontrol),
      MAF = pmin(EAF, 1 - EAF),
      OR  = exp(BETA),       # PGC BETA is log(OR)
      Neff_half = NEFF / 2
    ) %>%
    dplyr::select(dplyr::all_of(cols_to_munge))
}



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

