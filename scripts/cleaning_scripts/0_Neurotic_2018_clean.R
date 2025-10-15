
setwd("~/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

# CLEAN MDD FILES
sumstats <- sumstats %>% filter(cleaning_script == 'Neurotic_2018_Clean')
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

cols_to_munge <- c('CHROM',  # Chromosome number
                   'ID',     # RS ID
                   'A1',     # Effect allele
                   'A2',     # Other allele
                   'BETA',   # Effect size
                   'SE',     # Standard error
                   'P',      # P-value
                   'MAF',    # Minor allele frequency
                   'INFO',    # INFO score (optional – NA if not present)
                   'N'       # Sample size
)

for (trait in traits) {
  sumstats_cleaned[[trait]] <- sumstats_cleaned[[trait]] %>%
    mutate(
      CHROM = CHR,
      ID = SNP,
      N = N_analyzed
    ) %>%
    select(all_of(cols_to_munge))
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


