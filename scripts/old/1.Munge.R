# Munges files using the Genomic SEM package

# setwd("/Volumes/kschaumberg/GENIE_GSEM")
source("scripts/0.Packages.R")

load(file.path("gsem/gsem_settings/folderpaths.RData"))
load(file.path("gsem/gsem_settings/sumstats_metadata.RData"))
load(file.path("gsem/gsem_settings/variant_paths.RData"))

# remove /Volumes/kschaumberg/GENIE_GSEM from the file paths
sumstats$mungedpath.hm3 <- gsub("/Volumes/kschaumberg/GENIE_GSEM/", "", sumstats$mungedpath.hm3)

# remove from variant lists
variant_lists<- lapply(variant_lists, function(x) gsub("/Volumes/kschaumberg/GENIE_GSEM/", "", x))

# remove from all folderpaths
folderpaths <- lapply(folderpaths, function(x) gsub("/Volumes/kschaumberg/GENIE_GSEM/", "", x))

# setwd(folderpaths$base_path) #Set Working Directory

codes <- unlist(as.list(sumstats$code))
sumstatsToMunge<-sumstats |> 
  filter(code %in% codes)

# add data/gwas_sumstats/ to cleaned_file
sumstatsToMunge$cleanedpath <- paste0("data/gwas_sumstats/", sumstatsToMunge$cleaned_file)

munge(
  files = sumstatsToMunge$cleanedpath,
  trait.name = sumstatsToMunge$code,
  N = as.numeric(sumstatsToMunge$gsem_n),
  info.filter = 0.6,
  maf.filter = 0.01, 
  hm3 = variant_lists$w_hm3.snplist
)



# Move munged files to the munged folder
file.rename(
  from = list.files(pattern = ".gz", full.names = TRUE),
  to = file.path(folderpaths$data.gwas_sumstats.munged.hm3_munged_info_0_6, basename(list.files(pattern = ".gz", full.names = TRUE)))
)

# Move log files to the munged folder
file.rename(
  from = list.files(pattern = ".log", full.names = TRUE),
  to = file.path(folderpaths$data.gwas_sumstats.munged.hm3_munged_info_0_6, basename(list.files(pattern = ".log", full.names = TRUE)))
)
