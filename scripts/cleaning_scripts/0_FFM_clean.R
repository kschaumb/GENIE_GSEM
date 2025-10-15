
setwd("~/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")
library(stringr)
# CLEAN MDD FILES
sumstats <- sumstats %>% filter(cleaning_script == 'FFM_Clean')
# make a list of traits by pulling a list of all the codes in the data frame
traits <- as.list(sumstats$code)
traits <- lapply(traits, unlist)

sumstats_cleaned <- list()

sumstats <- sumstats %>%
  mutate(readpath = stringr::str_remove(rawpath, ".*GENIE_GSEM/"))

# make a savepath column that is the relative version of the cleaned path
sumstats<- sumstats %>%
  mutate(cleanedpath = paste0("data/gwas_sumstats/", cleaned_file))

for (trait in traits){
  sumstats_cleaned[[trait]] <- fread(sumstats %>% filter(code == !!trait) %>% pull(readpath))
}
library(data.table)
library(dplyr)


# 2. Keep only what LDSC/GenomicSEM need and rename
ffm_clean <- sumstats_cleaned$FFM %>%
  transmute(
    SNP   = rsid,                          # rsIDs are present – use them
    A1    = toupper(effect_allele),        # effect allele
    A2    = toupper(other_allele),         # other allele
    BETA  = beta,                          # effect estimate (already linear scale)
    SE    = standard_error,                # standard error
    P     = p_value                        # p-value
  ) %>%
  filter(nchar(A1)==1, nchar(A2)==1)       # drop indels / multi-base alleles

# 3. Save the cleaned file for munging
fwrite(ffm_clean, "data/gwas_sumstats/FFM/FFM.gz", sep = "\t", quote = FALSE)
