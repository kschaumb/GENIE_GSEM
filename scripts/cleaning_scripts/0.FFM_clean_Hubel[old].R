
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

# make a new colum in sumstats_daner called readpath that is the relative path to the sumstats file (rawpath but starting at data/ -- use stringr to drop everything from rawpath before data/)
sumstats <- sumstats %>%
  mutate(readpath = stringr::str_remove(rawpath, ".*GENIE_GSEM/"))

# make a savepath column that is the relative version of the cleaned path
sumstats<- sumstats %>%
  mutate(cleanedpath = paste0("data/gwas_sumstats/", cleaned_file))

for (trait in traits){
  sumstats_cleaned[[trait]] <- fread(sumstats %>% filter(code == !!trait) %>% pull(readpath))
}

library(data.table)
library(bigsnpr)

# Specify the directory where you want to download the data
download_1000G(dir = "data/variant_lists/1000G", overwrite = FALSE, delete_zip = TRUE)

# Read the .bim file
bim <- fread("data/variant_lists/1000G/1000G_phase3_common_norel.bim", col.names = c("CHR", "rsID", "CM", "BP", "A1", "A2"))

# Create a reference mapping
snp_ref <- bim[, .(CHR, BP, rsID)]
snp_ref <- snp_ref %>%
  rename(CHROM = CHR) 
library(dplyr)
library(tidyr)

cols_to_munge <- c('CHROM', 'ID', 'A1', 'A2', 'BETA', 'SE', 'P', 'MAF', 'INFO')

library(stringr)

ffm_raw <- fread(sumstats %>% filter(code == "FFM") %>% pull(readpath))

ffm_cleaned <- ffm_raw %>%
  filter(str_detect(MarkerName, "^[0-9XY]+:[0-9]+_[ACGTacgt]+_[ACGTacgt]+")) %>%
  separate(MarkerName, into = c("chr_pos", "A1_raw", "A2_raw"), sep = "_") %>%
  separate(chr_pos, into = c("CHROM", "BP"), sep = ":", convert = TRUE)



sumstats_cleaned$FFM <- ffm_raw %>%
  filter(str_detect(MarkerName, "^rs")) %>%  # Keep only valid rsIDs
  rename(
    ID = MarkerName,
    A1 = Allele1,
    A2 = Allele2,
    BETA = Effect,
    SE = StdErr,
    P = Pvalue,
    MAF = Freq1,
    INFO = info,
    CHROM = chr
  ) %>%
  mutate(MAF = pmin(MAF, 1 - MAF)) %>%
  semi_join(hm3_snps, by = "ID") %>%  # Optional but recommended for LDSC
  select(all_of(cols_to_munge))



# save the cleaned sumstats using the sumstats_sub cleanedpath, compressing each to gzip
for (trait in traits){
  cleaned_path <- sumstats %>% filter(code == trait) %>% pull(cleanedpath)
  
  if (length(cleaned_path) != 1 || is.na(cleaned_path)) {
    stop(paste("Error: Invalid cleaned file path for trait:", trait))
  }
  
  write.table(sumstats_cleaned[[trait]], gzfile(cleaned_path), 
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
}


