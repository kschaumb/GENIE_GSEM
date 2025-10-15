
source("scripts/0.Packages.R")
load("gsem/gsem_settings/folderpaths.RData")
load("gsem/gsem_settings/sumstats_metadata.RData")
load("gsem/gsem_settings/variant_paths.RData")

library(reshape2)
library(ggplot2)
library(dplyr)
library(tibble)

# ---------- Load LDSC ----------
# Filter for European-ancestry traits
symptoms <- sumstats %>%
  filter(Symptom_RGs == 1)
symptoms <- symptoms$code

# Dynamically construct the expected mvLD filename
# (Assumes consistent naming scheme: "mvLD_<traits joined by _>.Rds")
ldsc_path <- file.path(folderpaths$gsem.mvLD,
  paste0("mvLD_", paste(symptoms, collapse = "_"), ".Rds")
)

# Confirm the path and read
if (!file.exists(ldsc_path)) {
  stop("❌ mvLD file not found: ", ldsc_path)
}

ldsc <- readRDS(ldsc_path)
ldsc <- ldsc$hm3

# Ensure S has proper dimnames and is symmetric
if (is.null(rownames(ldsc$S))) rownames(ldsc$S) <- colnames(ldsc$S)
ldsc$S <- (ldsc$S + t(ldsc$S)) / 2

trait_names <- colnames(ldsc$S)

# ---------- Genetic correlations ----------
gen_cor <- cov2cor(ldsc$S)
gen_cor[gen_cor > 1]  <- 1
gen_cor[gen_cor < -1] <- -1
gen_cor <- round(gen_cor, 2)
rownames(gen_cor) <- trait_names
colnames(gen_cor) <- trait_names

# ---------- Heritabilities ----------
# Prefer LDSC-provided per-trait h2 if available; else fall back to diag(S)
h2_vec <- if ("h2" %in% names(ldsc)) {
  # make sure it's in same order as S columns
  hv <- as.numeric(ldsc$h2)
  names(hv) <- trait_names
  hv
} else {
  message("NOTE: ldsc$h2 not found; using diag(ldsc$S) as an approximation.")
  setNames(diag(ldsc$S), trait_names)
}
h2_vec <- pmin(pmax(h2_vec, 0), 1)   # clamp for display
h2_df  <- data.frame(Trait = trait_names, h2 = round(h2_vec, 2))

# ---------- Cluster order ----------
hc <- hclust(as.dist(1 - abs(gen_cor)))
trait_order <- rownames(gen_cor)[hc$order]

# ---------- Long form; keep UPPER triangle (we'll flip y so triangle sits at bottom) ----------
cor_df <- reshape2::melt(gen_cor, varnames = c("Trait1","Trait2"), value.name = "rG")
cor_df$Trait1 <- factor(cor_df$Trait1, levels = trait_order)
cor_df$Trait2 <- factor(cor_df$Trait2, levels = trait_order)
cor_df <- subset(cor_df, as.numeric(Trait1) <= as.numeric(Trait2))

# ---------- Clean labels ----------
clean_names <- trait_order |>
  gsub("_Eur", "", x = _) |>
  gsub("_noAN", "", x = _) |>
  gsub("_StageMidAdol", " Stage", x = _) |>
  gsub("_", " ", x = _) %>% 
  recode("ANR" = "AN [Restriction]", "ANX FS" = "Anxiety Factor Score", "OCSymptoms" = "OC Symptoms", "Depressive Sx" = "Depressive Symptoms",
                 "WtGain Case" = "Weight Gain Case", "WtLoss Case" = "Weight Loss Case",
                 "BEBROAD" = "BE Broad", "BENARROW" = "BE Narrow", "Neurotic" = "Neuroticism", "FFM" = "Fat Free Mass", "BMI Childhood" = "Childhood BMI")
names(clean_names) <- trait_order

cor_df$Trait1 <- factor(cor_df$Trait1, levels = trait_order, labels = clean_names)
cor_df$Trait2 <- factor(cor_df$Trait2, levels = trait_order, labels = clean_names)
h2_df$Trait   <- factor(h2_df$Trait,   levels = trait_order, labels = clean_names)

# ---------- Plot (preserve your end-of-plot formatting) ----------
ggplot() +
  geom_tile(data = subset(cor_df, Trait1 != Trait2),
            aes(x = Trait1, y = Trait2, fill = rG), color = "white") +
  geom_tile(data = subset(cor_df, Trait1 == Trait2),
            aes(x = Trait1, y = Trait2), fill = "white", color = "white") +
  geom_text(
    data = subset(cor_df, Trait1 != Trait2) |> 
      mutate(label_fmt = sub("(?<=-|^)0\\.", ".", sprintf("%.2f", rG), perl = TRUE)),
    aes(x = Trait1, y = Trait2, label = label_fmt),
    size = 5
  ) +
  geom_text(
    data = h2_df,
    aes(x = Trait, y = Trait, label = sprintf("%.2f", h2)),
    size = 5, fontface = "bold", color = "#1A4F66"
  ) +
  scale_fill_gradient2(low = "#1A4F66", mid = "white", high = "#EF6C45",
                       midpoint = 0, limits = c(-1, 1),
                       name = expression(paste("Genetic ", italic("r")))) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(cor_df$Trait2)), expand = c(0, 0)) +
  coord_equal(clip = "off") +
  theme_minimal() +
  labs(title = "(b) Genetic Correlation Matrix - Key Traits and Symptoms",
       x = "Trait", y = "Trait") +
  embarktools::embark_theme_a +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# keep no legend

ggsave("gsem/gsem_output/Genetic_Correlation_Matrix_Symptoms.png",
       width = 10, height = 10, units = "in")

# ---- Sanity print (optional): quick table of diagonal labels being used ----
print(h2_df[order(as.numeric(h2_df$Trait)), ], row.names = FALSE)
