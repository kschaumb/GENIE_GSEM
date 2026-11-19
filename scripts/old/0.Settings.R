
# Make list of folderpaths
folderpaths <- list()
# Define a base path
base_path <- ("~/Library/CloudStorage/Box-Box/GENIE_GSEM")
source(paste0(base_path, "/scripts/0.Packages.R"))

# Define sub-paths
sub_paths <- c("data", 
               "data/alignment_chains", #h19toHg38 chain
               "data/gwas_sumstats", # contains below subfolders
               "data/gwas_sumstats/munged", # contains below subfolders
               "data/gwas_sumstats/munged/hm3_munged_info_0_8_maf_0_1",
               "data/gwas_sumstats/raw", # contains raw sumstats 
               "data/ld_scores", #contains below subfolders
               "data/ld_scores/hc1kgp3.b38.eur.l2.jz2023.chr", #1kg build 38 ld scores european
               "data/ld_scores/eur_w_ld_chr",# european hapmap 3 ld scores
               "data/variant_lists", # hm3 and 1kg variant lists
               "scripts", #scripts
               "gsem", # directory for gwas-by-subtraction output
               "gsem/mvLD", #empty folder - mvLD ouput goes here
               "gsem/mvLD/log_files", #empty folder - mvLD ouput goes here
               "gsem/gsem_sumstats", #empty folder - gsem sumstats output goes here
               "gsem/gsem_sumstats/log_files", #empty folder - gsem sumstats output goes here
               "gsem/gsem_output", # empty folder - other gsem output goes here
               "gsem/gsem_output/log_files", # empty folder - other gsem output goes here
               "gsem/gsem_settings") # empty folder - gsem settings go here


# Detect if windows or linux
# if (Sys.info()[['sysname']] == 'Windows') {
# If windows, replace forward slashes with backslashes
# sub_paths <- gsub("/", "\\\\", sub_paths)
#}

# Create a list of folder paths by concatenating the base path with sub-paths
folder_paths <- lapply(sub_paths, function(sub_path) file.path(base_path, sub_path))

# make folders if they don't exist
lapply(folder_paths, function(path) {
  if (!dir.exists(path)) {
    suppressWarnings(dir.create(path, recursive = TRUE))
  }
})


# Assign the folder paths to the corresponding names in folderpath
folder_path_names <- gsub("/", ".", sub_paths)
names(folder_paths) <- folder_path_names
folderpaths <- folder_paths
# folderpaths$base_path <- base_path

rm(list = 'sub_paths', 'wd', 'folder_paths')

variant_lists <- data.frame(
  FileName = c(
    list.files(folderpaths$data.variant_lists, full.names = FALSE)
  ),
  FilePath = c(
    list.files(folderpaths$data.variant_lists, full.names = TRUE))
)

# Create a list with file names as names and file paths as list items
variant_lists <- as.list(setNames(variant_lists$FilePath, variant_lists$FileName))
ref_hm3 <- variant_lists$reference.1000G.maf.0.005.txt 
refpanels <- c('hm3')

# Load sumstats metadata
sumstats <- readxl::read_excel(file.path(folderpaths$data.gwas_sumstats, "GWAS_sumstats_settings.xlsx"))

#Add munged pxaths to dataset
sumstats$mungedpath.hm3 <- file.path(folderpaths$data.gwas_sumstats.munged.hm3_munged_info_0_8_maf_0_1, paste0(sumstats$code, '.sumstats.gz'))

sumstats$rawpath <- file.path(folderpaths$data.gwas_sumstats, paste0(sumstats$raw_file))

setwd(base_path)

save(sumstats, file = file.path('gsem/gsem_settings', "sumstats_metadata.RData"))
save(folderpaths, file = file.path('gsem/gsem_settings', "folderpaths.RData"))
save(variant_lists, file = file.path('gsem/gsem_settings', "variant_paths.RData"))
rm(list = ls())
