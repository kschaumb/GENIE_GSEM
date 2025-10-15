
library(GenomicSEM)
source(file.path("scripts/0.Packages.R"))
load(file.path("gsem/gsem_settings/folderpaths.RData"))
load(file.path("gsem/gsem_settings/sumstats_metadata.RData"))


all <- unlist(as.list(sumstats$code))
# remove 'Height' from codes
all <- all[!all %in% c("Height")]

trait_lists <- list(
  all,
  c('AN', 'ANX', 'BEBROAD', 'BMI', 'MDD', 'OCD', 'PTSD', 'BMI_Childhood')
)
folderpaths <- lapply(folderpaths, function(x) gsub("/Volumes/kschaumberg/GENIE_GSEM/", "", x))
# remove /Volumes/kschaumberg/GENIE_GSEM from the file paths
sumstats$mungedpath.hm3 <- gsub("/Volumes/kschaumberg/GENIE_GSEM/", "", sumstats$mungedpath.hm3)

# remove from variant lists
variant_lists$w_hm3.snplist <- gsub("/Volumes/kschaumberg/GENIE_GSEM/", "", variant_lists$w_hm3.snplist)

# Create separate lists to store heritabilities and intercepts for each list and each LD dataset
heritabilities_list <- list()
intercepts_list <- list()

# Loop through each trait list
for (list in trait_lists) {
  selectedCodes <- list
  
  # Check if the RDS file for this list already exists
  rds_file_name <- paste("mvLD_", paste(list, collapse = "_"), ".Rds", sep = "")
  rds_file_path <- file.path(folderpaths$gsem.mvLD, rds_file_name)
  

    # Perform LD Score regression for the first LD dataset
    ld_result_hm3 <- ldsc(
      traits = sumstats[sumstats$code %in% selectedCodes, ]$mungedpath.hm3,
      sample.prev = sumstats[sumstats$code %in% selectedCodes, ]$samplePrevalence,
      population.prev = sumstats[sumstats$code %in% selectedCodes, ]$populationPrevalence,
      trait.names = sumstats[sumstats$code %in% selectedCodes, ]$code,
      ld = folderpaths$data.ld_scores.eur_w_ld_chr
    )
  
    # Save the results as an RDS file
    saveRDS(object = list(hm3 = ld_result_hm3), file = rds_file_path)
  
  # Extract heritabilities and intercepts for the LD datasets
  ld_result <- readRDS(file = rds_file_path)
  heritabilities_hm3 <- diag(ld_result$hm3$S)
  intercepts_hm3 <- diag(ld_result$hm3$I)

  # Assign names to heritabilities and intercepts
  names(heritabilities_hm3) <- colnames(ld_result$hm3$S)
  names(intercepts_hm3) <- colnames(ld_result$hm3$I)

  # Store the heritabilities and intercepts in separate lists
  heritabilities_list[[rds_file_name]] <- list(hm3 = heritabilities_hm3)
  intercepts_list[[rds_file_name]] <- list(hm3 = intercepts_hm3)
  
  # Save the lists as RDS files
  saveRDS(object = heritabilities_list, file = file.path(folderpaths$gsem.mvLD, "heritabilities_list.Rds"))
  saveRDS(object = intercepts_list, file = file.path(folderpaths$gsem.mvLD, "intercepts_list.Rds"))

}

rm('ld_result', 'ld_result_hm3', 'heritabilities_hm3', 'intercepts_hm3')

# # move log files
# log_files <- list.files(pattern = "log", full.names = TRUE)
# file.copy(log_files, file.path(folderpaths$gxsub.gxsub.mvLD, "log_files"))
# # delete log files
# file.remove(log_files)

