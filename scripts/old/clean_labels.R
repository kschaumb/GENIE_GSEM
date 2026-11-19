# scripts/labels_clean.R
# ------------------------------------------------------------
# Clean, human-readable labels for traits/symptoms/disorders
# - Works on any character vector of codes (e.g., colnames(S))
# - Keeps your previous gsub rules (remove _Eur, _noAN, etc.)
# - Normalizes common variants before mapping
# - Warns on unmapped codes so you can extend the map
# ------------------------------------------------------------

# 1) Optional pre-normalization of code variants before mapping
.normalize_codes <- function(x) {
  x <- gsub("\\s+", "_", x)                    # spaces -> underscore
  x <- gsub("\\bSx\\b", "Symptoms", x)         # Sx -> Symptoms
  x <- gsub("Depressive_Sx", "Depressive_Symptoms", x, ignore.case = TRUE)
  x <- gsub("Dep_Sx",        "Depressive_Symptoms", x, ignore.case = TRUE)
  x <- gsub("ANX_FS",        "ANX_FS", x, ignore.case = TRUE)  # ensure consistent
  x <- gsub("OCSymptoms",    "OC_Symptoms", x)
  x
}

# 2) Master label map (extend as needed)
.label_map <- c(
  # Core ED phenotypes
  "AN"                 = "Anorexia Nervosa",
  "ANR"                = "AN [Restriction]",
  "BEBROAD"            = "Binge Eating",
  "BENARROW"           = "Binge Eating (Narrow)",
  "BE"                 = "Binge Eating",
  "BE_noAN"            = "Binge Eating (not AN ascertained)",
  "BE_no_ascertained"  = "Binge Eating (not AN ascertained)",
  
  # Anxiety / OC / mood
  "ANX"                = "Anxiety Disorders",
  "ANX_FS"             = "Anxiety Factor Score",
  "OC_Symptoms"        = "OC Symptoms",
  "OCD"                = "Obsessive–Compulsive Disorder",
  "Depressive_Symptoms"= "Depressive Symptoms",
  "Neurotic"           = "Neuroticism",
  
  # MDD weight-change cases
  "MDD_WtGain_Case"    = "MDD Weight Gain",
  "MDD_WtLoss_Case"    = "MDD Weight Loss",
  
  # PTSD
  "PTSD"               = "PTSD",
  "PTSD_PCS"           = "PTSD Symptoms",
  
  # Anthropometry & activity
  "BMI"                = "BMI",
  "BMI_Childhood"      = "Childhood BMI",
  "BFP"                = "Body Fat Percentage",
  "FFM"                = "Fat-Free Mass",
  "PA"                 = "Physical Activity",
  
  # Puberty / timing
  "Age_Menarche"       = "Age at Menarche",
  "Tanner_StageMidAdol"= "Tanner Stage (Mid-Adolescence)"
)

# 3) Humanizing pass for leftovers (title-case, underscore->space, specific tweaks)
.prettify_fallback <- function(x) {
  y <- x
  y <- gsub("_", " ", y)
  y <- gsub("\\bEur\\b", "", y, ignore.case = TRUE)
  y <- gsub("\\bnoAN\\b", "(not AN ascertained)", y)
  y <- gsub("StageMidAdol", "Stage (Mid-Adolescence)", y)
  # Title case (simple)
  y <- tools::toTitleCase(trimws(y))
  # Keep common acronyms nice:
  y <- gsub("\\bOcd\\b", "OCD", y)
  y <- gsub("\\bPtsd\\b", "PTSD", y)
  y <- gsub("\\bBmi\\b", "BMI", y)
  y <- gsub("\\bBfp\\b", "BFP", y)
  y <- gsub("\\bFfm\\b", "FFM", y)
  y
}

# 4) Public function: clean_recode()
#    Input: character vector of codes (e.g., colnames(S))
#    Output: same-length character vector of pretty labels
clean_recode <- function(trait_codes) {
  stopifnot(is.character(trait_codes))
  raw <- trait_codes
  
  # strip standard suffixes first (as you do in your plots)
  base <- raw |>
    gsub("_Eur\\b", "", x = _, ignore.case = TRUE) |>
    gsub("_noAN\\b", "", x = _, ignore.case = TRUE)
  
  # normalize variants
  norm <- .normalize_codes(base)
  
  # map known labels
  mapped <- .label_map[norm]
  out <- mapped
  
  # unmapped -> fallback prettifier
  needs_fallback <- is.na(out)
  if (any(needs_fallback)) {
    out[needs_fallback] <- .prettify_fallback(norm[needs_fallback])
  }
  
  # optional: warn on truly unknown (i.e., didn’t hit explicit map and still looks rough)
  unknown <- setdiff(unique(norm[is.na(mapped)]), names(.label_map))
  if (length(unknown)) {
    msg <- paste0("Labels not in map (using fallback prettifier): ",
                  paste(sort(unknown), collapse = ", "))
    message(msg)
  }
  
  names(out) <- raw
  out
}

# 5) Helper to apply clean labels to a factor given a desired order
#    (useful before plotting or building diagrams)
apply_clean_labels <- function(levels_in_order) {
  lbls <- clean_recode(levels_in_order)
  stats::setNames(lbls, levels_in_order)
}
