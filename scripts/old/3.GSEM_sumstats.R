source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")


all <- unlist(as.list(sumstats$code))
# remove 'Height' from codes
all <- all[!all %in% c("Height", "BE_noAN")]

symptoms <- sumstats %>%
  filter(Symptoms_GSEM == 1)
symptoms <- symptoms$code

disorders <- sumstats %>% 
  filter(Disorder_GSEM == 1)
disorders <- disorders$code

trait_lists <- list(
 # all,
  symptoms
 # disorders
)

info.filter <- 0.6
maf.filter <- 0.01
folderpaths <- lapply(folderpaths, function(x) gsub("~/GENIE_GSEM/", "", x))
  # Define the subdirectory path
  subdirectory <- file.path(folderpaths$gsem.gsem_sumstats)
  
  # Create the subdirectory if it doesn't exist
  if (!dir.exists(subdirectory)) {
    dir.create(subdirectory, recursive = TRUE)
  }
  
  for (list in trait_lists) {
    rdata_file_name <- paste(paste(list, collapse = "_"), "_sumstats.RData", sep = "")
    fp <- file.path(subdirectory, rdata_file_name)
    
    # Filter the sumstats for the current trait list
    sumstats_filtered <- sumstats %>%
      filter(code %in% list)
    
# Only run if the file doesn't already exist
      # Prepare the sumstats for the current list
      col_name <- paste("mungedpath.", sep = "")
      
      gsem_sumstats <- sumstats(
        files = as.vector(sumstats_filtered$cleanedpath),
        ref = variant_lists$reference.1000G.maf.0.005.txt.gz,
        trait.names = as.vector(sumstats_filtered$code),
        se.logit = as.vector(sumstats_filtered$gsem_se.logit),
        info.filter = info.filter,
        maf.filter = maf.filter,
        OLS = as.vector(sumstats_filtered$gsem_OLS),
        linprob = as.vector(sumstats_filtered$gsem_linprob),
        betas = NULL,
        N = as.vector(sumstats_filtered$gsem_n) #  previously used N = sumstatsToMunge$n_total - from tutorial: " If the summary statistics file includes a sample size column, then the user can also list NA if they wish to use the SNP-specific sample sizes to perform the rescaling. However ,we note again that this should be the sum of effective sample sizes for dichotomous traits. "
        
      )
      
      gsem_sumstats$MAF <- as.numeric(gsem_sumstats$MAF)
      
      # Save the prepared sumstats as an RData file
      save(gsem_sumstats, file = fp)
      print(paste("Saved sumstats for list", paste(list, collapse = "_")))
     
  
    
  }

  # # move log files
  # log_files <- list.files(pattern = "log", full.names = TRUE)
  # file.copy(log_files, file.path(folderpaths$gxsub.gsem_sumstats, "log_files"))
  # # delete log files
  # file.remove(log_files)
