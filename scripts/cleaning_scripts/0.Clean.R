base_path <- ("/Volumes/kschaumberg/GENIE_GSEM")
source(paste0(base_path, "/scripts/0.Packages.R"))
load(paste0(base_path,"/gsem/gsem_settings/folderpaths.RData"))
load("gsem/gsem_settings/sumstats_metadata.RData")

# CLEAN DANER FILES

# traits that need N column done
daner_files <- c("AN", "ANR", "ANBP", "BENARROW", "BEBROAD", "OCSymptoms")

gz_files <- sumstats$rawpath[grepl("\\.gz$", sumstats$rawpath)]
# add ~ to the beginning of the file path

for (file in gz_files) {
  cat("\n### Preview of:", file, "###\n")
  system(paste("gunzip -c", shQuote(file), "| head -100"))
}

# subset only sumstats that need N column done
sumstats_daner <- sumstats %>% filter(code %in% daner_files)

# make empty list to store sumstats with N column
sumstats_cleaned <- list()

cols_to_munge <- c('CHR', #Chromosome number
                   'SNP', # SNP
                   'A1', # Alelle 1 (OR on this alelle)
                   'A2', # Allele 2
                   'MAF', # A1 frequency 
                   'INFO', # imputation quality
                   'OR', # odds ratio ()
                   'SE', 
                   'P', # P-value
                   'Neff_half' # effective sample size/2 -- LDSC will double it
                   )

# for each trait, read in the sumstats, add N column that mutliplies Neff_half by 2
for (trait in daner_files){
  sumstats_cleaned[[trait]] <- fread(sumstats_sub %>% filter(trait == !!trait) %>% pull(rawpath))
}

for (trait in pgc_traits){
  # identify the column names of sumstats_cleaned[[trait]] that start with FRQ 
  frq_cols <- (grep("^FRQ", names(sumstats_cleaned[[trait]]), value = TRUE))
  # identify the n_cases and n_controls from the frq_cols list n_cases is the number after FRQ_A, n_controls is the number after FRQ_U
  n_cases <- as.numeric(stringr::str_extract(frq_cols, "(?<=FRQ_A_)\\d+")[1])
  n_controls <- as.numeric(stringr::str_extract(frq_cols, "(?<=FRQ_U_)\\d+")[2])
   # calculate MAF column by multiplying the column that starts with FRQ_A * n_cases + the column that starts with FRQ_U * n_controls and dividing by n_cases + n_controls
  sumstats_cleaned[[trait]] <- sumstats_cleaned[[trait]] %>%
    mutate(MAF = (get(frq_cols[1]) * n_cases + get(frq_cols[2]) * n_controls) / (n_cases + n_controls)) %>%
    select(all_of(cols_to_munge))
}

# save the cleaned sumstats using the sumstats_sub cleanedpath, compressing each to gzip
for (trait in daner_files){
  write.table(sumstats_cleaned[[trait]], gzfile(sumstats_sub[trait, 'cleanedpath']), 
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
}


# 