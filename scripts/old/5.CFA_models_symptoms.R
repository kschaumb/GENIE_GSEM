base_path <- ('~/Library/CloudStorage/Box-Box/GENIE_GSEM/')
# change directory to base path
setwd(base_path)
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")
load("gsem/gsem_settings/variant_paths.RData")

library(reshape2)
library(ggplot2)
library(GenomicSEM)


ldsc <- readRDS("gsem/mvLD/mvLD_ANR_ANX_FS_BENARROW_BMI_MDD_WtGain_Case_MDD_WtLoss_Case_OCSymptoms_PTSD_PCS_BMI_Childhood_BFP_Neurotic_PA_FFM_Insulin_Age_Menarche_Tanner_StageMidAdol_Depressive_Sx.Rds")
ldsc <- ldsc$hm3
symptoms <- sumstats %>%
  filter(Symptoms_GSEM == 1)
trait_list <- symptoms$code



cfa_model_twofac <- '
WeightRelated =~ BMI + BMI_Childhood + BFP + ANR + BENARROW + MDD_WtGain_Case + MDD_WtLoss_Case + Age_Menarche + Tanner_StageMidAdol + FFM + PA
INT =~ ANX_FS + BENARROW + ANR + OCSymptoms + PTSD_PCS + Neurotic + Depressive_Sx
WeightRelated ~~ INT
'

cfa_model_threefac <- '
WeightRelated =~ BMI + BMI_Childhood + BFP + ANR + BENARROW + MDD_WtGain_Case + MDD_WtLoss_Case + PA + FFM
INT =~ ANX_FS + BENARROW + ANR + OCSymptoms + PTSD_PCS + Neurotic + Depressive_Sx
Puberty =~ BMI + Tanner_StageMidAdol + Age_Menarche + ANR + BENARROW + BMI_Childhood + Insulin
WeightRelated ~~ INT
WeightRelated ~~ Puberty
INT ~~ Puberty
'


cfa_model_threefac_pared <- '
WeightRelated =~ BMI + BFP + ANR + MDD_WtGain_Case + MDD_WtLoss_Case + PA + FFM
INT =~ ANX_FS + BENARROW + ANR + OCSymptoms + PTSD_PCS + Neurotic + Depressive_Sx
Puberty =~ BMI + Tanner_StageMidAdol + Age_Menarche + BENARROW + BMI_Childhood 
WeightRelated ~~ INT
WeightRelated ~~ Puberty
Puberty ~~ 0*INT
'

cfa_models <- list(cfa_model_twofac, cfa_model_threefac, cfa_model_threefac_pared)
cfa_results <- list()

for (cfa_model in cfa_models) {
  cfa_result <- usermodel(
    covstruc = ldsc,
   estimation = "DWLS",
    model = cfa_model,
    std.lv = TRUE
   # toler = 1e-40
  )
  # save the results
  cfa_results[[cfa_model]] <- cfa_result
}


# rename cfa_resutls to one_fac, two_fac, three_fac
names(cfa_results) <- c("two_fac_symptoms", "three_fac_symptoms", "three_fac_symptoms_pared")
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Loop over each CFA model in cfa_results
for(model_name in names(cfa_results)) {
  model_result <- cfa_results[[model_name]]
  
  # Subset only the loadings (rows where op == "=~")
  loadings <- subset(model_result$results, op == "=~")
  
  # Identify unique latent factors (from lhs) and observed variables (from rhs)
  latent_factors <- unique(loadings$lhs)
  observed_vars <- unique(loadings$rhs)
  
  # Begin constructing the DOT string
  dot_string <- "digraph pathDiagram {\n"
  dot_string <- paste0(
    dot_string,
    "  graph [layout = dot, rankdir = TB, overlap = false, ",
    "splines = true, ranksep = 2.0, nodesep = 0.5];\n",
    "  node [fontname = Helvetica, fontsize = 14];\n",
    "  edge [fontname = Helvetica, fontsize = 14];\n"
  )
  
  # Define latent factor nodes as ellipses
  for (lat in latent_factors) {
    dot_string <- paste0(
      dot_string,
      "  ", lat, " [shape = ellipse, style = filled, fillcolor = lightblue, ",
      "label = '", lat, "'];\n"
    )
  }
  
  # Define observed variable nodes as boxes
  for (obs in observed_vars) {
    dot_string <- paste0(
      dot_string,
      "  ", obs, " [shape = box, label = '", obs, "'];\n"
    )
  }
  
  # Add edges (latent -> observed) with standardized loadings as labels
  for (i in seq_len(nrow(loadings))) {
    edge_label <- sprintf("%.2f", loadings$STD_Genotype[i])
    edge_label <- gsub("^0\\.", ".", edge_label)
    # also take out leading 0 if its negative
    edge_label <- gsub("^-0\\.", "-", edge_label)
    dot_string <- paste0(
      dot_string,
      "  ", loadings$lhs[i], " -> ", loadings$rhs[i],
      " [label = '", edge_label, "'];\n"
    )
  }
  
  # ------------------------------------------------
  # ADD FACTOR CORRELATIONS (double-headed arrows)
  # ------------------------------------------------
  # Subset rows for factor correlations (op == "~~") among latent factors
  corrs <- subset(model_result$results, 
                  op == "~~" & lhs %in% latent_factors & rhs %in% latent_factors & lhs != rhs)
  
  # For each correlation, add a double-headed arrow (dir=both) between factors
  for (i in seq_len(nrow(corrs))) {
    # Use the standardized correlation if available, e.g. STD_Genotype or STD_All
    edge_label <- sprintf("%.2f", corrs$STD_Genotype[i])
    # take out the leading 0
    edge_label <- gsub("^0\\.", ".", edge_label)
    # also take out leading 0 if its negative
    edge_label <- gsub("^-0\\.", "-", edge_label)
    dot_string <- paste0(
      dot_string,
      "  ", corrs$lhs[i], " -> ", corrs$rhs[i],
      " [label = '", edge_label, "', dir=both, arrowtail=normal, arrowhead=normal, style=dashed];\n"
    )
  }
  
  # Force latent factors to appear at the top (rank=min)
  dot_string <- paste0(dot_string, "  { rank = min; ")
  for (lat in latent_factors) {
    dot_string <- paste0(dot_string, lat, "; ")
  }
  dot_string <- paste0(dot_string, "}\n")
  
  # Force observed variables to appear at the bottom (rank=max)
  dot_string <- paste0(dot_string, "  { rank = max; ")
  for (obs in observed_vars) {
    dot_string <- paste0(dot_string, obs, "; ")
  }
  dot_string <- paste0(dot_string, "}\n")
  
  dot_string <- paste0(dot_string, "}\n")
  
  # Render the diagram using DiagrammeR
  graph <- grViz(dot_string)
  
  # Convert to SVG string
  svg <- export_svg(graph)
  
  # Define filenames
  svg_file <- paste0("gsem/gsem_output/", model_name, "_pathDiagram.svg")
  png_file <- paste0("gsem/gsem_output/", model_name, "_pathDiagram.png")
  
  # Save as SVG
  cat(svg, file = svg_file)
  
  # Convert SVG to PNG
  rsvg::rsvg_png(svg_file, png_file)
  
  message("Saved ", model_name, " diagram as ", svg_file, " and ", png_file)
}

# combine the model result summaries and put them in a table
library(kableExtra)
library(dplyr)
library(tidyr)
library(knitr)


# Create a list of the model fit results
model_fits <- list(
  Two_Factor_Symptoms = cfa_results$two_fac_symptoms$modelfit,
  Three_Factor_Symptoms = cfa_results$three_fac_symptoms$modelfit,
  Three_Factor_Symptoms_Pared = cfa_results$three_fac_symptoms_pared$modelfit
  
)

# Convert each one to a row in a data frame
model_fit_df <- do.call(rbind, lapply(names(model_fits), function(name) {
  fit <- model_fits[[name]]
  data.frame(
    Model = name,
    Chisq = fit$chisq,
    df = fit$df,
    p_chisq = fit$p_chisq,
    AIC = fit$AIC,
    CFI = fit$CFI,
    SRMR = fit$SRMR
  )
}))

# round the values to 3 decimal places
model_fit_df[, c("Chisq", "df", "p_chisq", "AIC", "CFI", "SRMR")] <- round(
  model_fit_df[, c("Chisq", "df", "p_chisq", "AIC", "CFI", "SRMR")], 3
)

# View the result
print(model_fit_df)
# Save the model fit results to a CSV file
write.csv(model_fit_df, "gsem/gsem_output/model_fit_results_symptoms.csv", row.names = FALSE)

