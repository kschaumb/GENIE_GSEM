
setwd("/Volumes/kschaumberg/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

# CLEAN MDD FILES
sumstats_mddsx <- sumstats %>% filter(cleaning_script == 'MDD_Symptoms')
# make a list of traits by pulling a list of all the codes in the data frame
mddsx_traits <- as.list(sumstats_mddsx$code)
mddsx_traits <- lapply(mddsx_traits, unlist)

sumstats_cleaned <- list()

# make a new colum in sumstats_daner called readpath that is the relative path to the sumstats file (rawpath but starting at data/ -- use stringr to drop everything from rawpath before data/)
sumstats_mddsx <- sumstats_mddsx %>%
  mutate(readpath = stringr::str_remove(rawpath, ".*GENIE_GSEM/"))

# make a savepath column that is the relative version of the cleaned path
sumstats_mddsx<- sumstats_mddsx %>%
  mutate(cleanedpath = paste0("data/gwas_sumstats/", cleaned_file))

for (trait in mddsx_traits){
  sumstats_cleaned[[trait]] <- fread(sumstats_mddsx %>% filter(code == !!trait) %>% pull(readpath))
}

cols_to_munge <- c('CHROM', #Chromosome number
                   'POS', # Position
                   'ID', # SNP ID
                   'A1', # Alelle 1 (OR on this alelle)
                   'A2', # Allele 2
                   'MAF', # MAF Frequency
                   'INFO', # imputation quality
                   'BETA', # effect size
                   'SE', 
                   'P', # P-value
                   'NCAS', # number of cases
                   'NCON', # number of controls
                   'NEFF' # effective sample size
)

for (trait in mddsx_traits) {
  sumstats_cleaned[[trait]] <- sumstats_cleaned[[trait]] %>%
    mutate(
      MAF_raw = (FCAS * NCAS + FCON * NCON) / (NCAS + NCON),  # Raw calculation
      MAF = pmin(MAF_raw, 1 - MAF_raw)  # Ensures MAF is ≤ 0.5
    ) %>% 
    # rename #CHROM to CHROM, PVAL to P and IMPINFO to INFO
    rename(CHROM = "#CHROM", P = PVAL, INFO = IMPINFO) %>%
    select(all_of(cols_to_munge))
}


# save the cleaned sumstats using the sumstats_sub cleanedpath, compressing each to gzip
for (trait in mddsx_traits){
  cleaned_path <- sumstats_mddsx %>% filter(code == trait) %>% pull(cleanedpath)
  
  if (length(cleaned_path) != 1 || is.na(cleaned_path)) {
    stop(paste("Error: Invalid cleaned file path for trait:", trait))
  }
  
  write.table(sumstats_cleaned[[trait]], gzfile(cleaned_path), 
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
}


