base_path <- ("/Volumes/kschaumberg/GENIE_GSEM")
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")

# CLEAN DANER FILES


# traits that need N column done
daner_files <- c("AN", "ANR", "ANBP", "BENARROW", "BEBROAD", "OCSymptoms", "BE_noAN")


# subset only sumstats that need N column done
sumstats_daner <- sumstats %>% filter(code %in% daner_files)

# make empty list to store sumstats with N column
sumstats_cleaned <- list()

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

# make a new colum in sumstats_daner called readpath that is the relative path to the sumstats file (rawpath but starting at data/ -- use stringr to drop everything from rawpath before data/)
sumstats_daner <- sumstats_daner %>%
  mutate(readpath = stringr::str_remove(rawpath, ".*GENIE_GSEM/"))

# make a savepath column that is the relative version of the cleaned path
sumstats_daner <- sumstats_daner %>%
  mutate(cleanedpath = paste0("data/gwas_sumstats/", cleaned_file))

for (trait in daner_files){
  sumstats_cleaned[[trait]] <- fread(sumstats_daner %>% filter(code == !!trait) %>% pull(readpath))
}




for (trait in daner_files){
  # identify the column names of sumstats_cleaned[[trait]] that start with FRQ 
  frq_cols <- (grep("^FRQ", names(sumstats_cleaned[[trait]]), value = TRUE))
  # identify the n_cases and n_controls from the frq_cols list n_cases is the number after FRQ_A, n_controls is the number after FRQ_U
  n_cases <- as.numeric(stringr::str_extract(frq_cols, "(?<=FRQ_A_)\\d+")[1])
  n_controls <- as.numeric(stringr::str_extract(frq_cols, "(?<=FRQ_U_)\\d+")[2])
   # calculate MAF column by multiplying the column that starts with FRQ_A * n_cases + the column that starts with FRQ_U * n_controls and dividing by n_cases + n_controls
  # Calculate Minor Allele Frequency (MAF)
  sumstats_cleaned[[trait]] <- sumstats_cleaned[[trait]] %>%
    mutate(
      MAF_raw = (get(frq_cols[1]) * n_cases + get(frq_cols[2]) * n_controls) / (n_cases + n_controls),
      MAF = pmin(MAF_raw, 1 - MAF_raw)  # Ensure MAF is ≤ 0.5
    ) %>%
    select(all_of(cols_to_munge))
}

# save the cleaned sumstats using the sumstats_sub cleanedpath, compressing each to gzip
for (trait in daner_files){
  cleaned_path <- sumstats_daner %>% filter(code == trait) %>% pull(cleanedpath)
  
  if (length(cleaned_path) != 1 || is.na(cleaned_path)) {
    stop(paste("Error: Invalid cleaned file path for trait:", trait))
  }
  
  write.table(sumstats_cleaned[[trait]], gzfile(cleaned_path), 
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
}


# 