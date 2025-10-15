
#setwd("/Volumes/kschaumberg/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

# CLEAN MDD FILES
sumstats_anx <- sumstats %>% filter(cleaning_script == 'ANX_2016')
# make a list of traits by pulling a list of all the codes in the data frame
anx_traits <- as.list(sumstats_anx$code)
anx_traits <- lapply(anx_traits, unlist)

sumstats_cleaned <- list()

# make a new colum in sumstats_daner called readpath that is the relative path to the sumstats file (rawpath but starting at data/ -- use stringr to drop everything from rawpath before data/)
sumstats_anx <- sumstats_anx %>%
  mutate(readpath = stringr::str_remove(rawpath, ".*GENIE_GSEM/"))

# make a savepath column that is the relative version of the cleaned path
sumstats_anx<- sumstats_anx %>%
  mutate(cleanedpath = paste0("data/gwas_sumstats/", cleaned_file))

for (trait in anx_traits){
  sumstats_cleaned[[trait]] <- fread(sumstats_anx %>% filter(code == !!trait) %>% pull(readpath))
}

cols_to_munge <- c('CHROM', #Chromosome number
                   'ID', # SNP ID
                   'A1', # Alelle 1 (OR on this alelle)
                   'A2', # Allele 2
                   'MAF', # MAF Frequency
                   'Effect', # effect size
                   'SE', 
                   'P', # P-value
                   'N' # Total N
)

for (trait in anx_traits) {
  sumstats_cleaned[[trait]] <- sumstats_cleaned[[trait]] %>%
    mutate(
      MAF = pmin(Freq1, 1 - Freq1)  # Ensures MAF is ≤ 0.5
    ) %>% 
    # rename #CHROM to CHROM, PVAL to P and IMPINFO to INFO
    rename(CHROM = CHR, P = P.value, N = TotalN, 
           SE = StdErr, A1= Allele1, A2 = Allele2, ID = SNPID) %>%
    
    select(all_of(cols_to_munge))
}


# save the cleaned sumstats using the sumstats_sub cleanedpath, compressing each to gzip
for (trait in anx_traits){
  cleaned_path <- sumstats_anx %>% filter(code == trait) %>% pull(cleanedpath)
  
  if (length(cleaned_path) != 1 || is.na(cleaned_path)) {
    stop(paste("Error: Invalid cleaned file path for trait:", trait))
  }
  
  write.table(sumstats_cleaned[[trait]], gzfile(cleaned_path), 
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
}


