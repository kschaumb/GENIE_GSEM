
#setwd("/Volumes/kschaumberg/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

sumstats_BMI <- sumstats %>% filter(cleaning_script == 'BMI_2019')

# load the raw file (already cleaned) and change it to .gz , then save it

sumstats_BMI <- sumstats_BMI %>%
  mutate(readpath = stringr::str_remove(rawpath, ".*GENIE_GSEM/"))

sumstats_BMI<- sumstats_BMI %>%
  mutate(cleanedpath = paste0("data/gwas_sumstats/", cleaned_file))

sumstats_cleaned <- list()

BMI_traits <- as.list(sumstats_BMI$code)

BMI_traits <- lapply(BMI_traits, unlist)

for (trait in BMI_traits){
  sumstats_cleaned[[trait]] <- fread(sumstats_BMI %>% filter(code == !!trait) %>% pull(readpath))
}


for (trait in BMI_traits){
  cleaned_path <- sumstats_BMI %>% filter(code == trait) %>% pull(cleanedpath)
  
  if (length(cleaned_path) != 1 || is.na(cleaned_path)) {
    stop(paste("Error: Invalid cleaned file path for trait:", trait))
  }
  
  write.table(sumstats_cleaned[[trait]], gzfile(cleaned_path), 
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
}