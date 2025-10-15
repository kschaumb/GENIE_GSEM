# ============================================================
# Path diagrams for GenomicSEM CFA fits (DiagrammeR)
# ============================================================
# Prefers: gsem/gsem_output/disorder_cfa/cfa_DWLS_results.RDS (consolidated)
# Fallback: per-fit CSVs: *_loadings.csv, *_factor_corr.csv
# Writes : gsem/gsem_output/disorder_cfa/diagrams/<model>_pathDiagram.(svg|png)
# Usage  : source("scripts/cfa_make_diagrams.R")
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
  library(lavaan)  # for parameterEstimates if using fits$lavaan_output
})

# -----------------------
# Config / What to render
# -----------------------
BASE_DIR        <- "~/Library/CloudStorage/Box-Box/GENIE_GSEM"
OUT_DIR         <- file.path(BASE_DIR, "gsem/gsem_output/disorder_cfa")
DIAG_DIR        <- file.path(OUT_DIR, "diagrams")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

# Choose which models to render (names must match consolidated$fits / per-fit CSV stems)
# e.g., with floor 0.01 -> "floor001"
MODEL_SELECTION <- c("orig_DWLS_floor001", "constrained_DWLS_floor001")

# Edge rendering knobs
MIN_LOADING_TO_DRAW  <- 0.05   # hide arrows with |loading| below this
MIN_LABEL_TO_PRINT   <- 0.00   # set 0.01 to hide tiny labels
DIGITS               <- 2      # decimals in edge labels

# -------------------------------------------------
# Helpers
# -------------------------------------------------
format_coef <- function(x, digits = DIGITS, min_label = MIN_LABEL_TO_PRINT) {
  if (length(x) == 0 || is.na(x) || is.nan(x)) return("")
  if (abs(x) < min_label) return("")
  lbl <- sprintf(paste0("%.", digits, "f"), x)
  lbl <- sub("^0\\.", ".",  lbl)
  lbl <- sub("^-0\\.", "-.", lbl)
  lbl
}

pick_std_col <- function(df, prefer = c("STD_Genotype","std.all","Std.all","STD.ALL","std","estimate","est","Estimate","EST")){
  cand <- intersect(prefer, names(df))
  if (length(cand)) cand[1] else NA_character_
}

# Always ensure we keep exactly one correlation per latent pair (lhs<rhs)
canonicalize_corrs <- function(R){
  if (!nrow(R)) return(R[, c("lhs","rhs","std_corr")])
  R$l_can <- pmin(R$lhs, R$rhs)
  R$r_can <- pmax(R$lhs, R$rhs)
  R <- R[!duplicated(R[, c("l_can","r_can")]), ]
  R <- R[, c("l_can","r_can","std_corr")]
  names(R) <- c("lhs","rhs","std_corr")
  R
}

# Prefer consolidated RDS -> per_fit_tables -> fits$lavaan_output; fallback: CSVs
extract_tables_from_fit <- function(model_name,
                                    consolidated_path = file.path(OUT_DIR, "cfa_DWLS_results.RDS")) {
  # ---------- Try consolidated RDS ----------
  if (file.exists(consolidated_path)) {
    cons <- try(readRDS(consolidated_path), silent = TRUE)
    if (!inherits(cons, "try-error")) {
      # 1) Direct per-fit tables (saved by the pipeline)
      if (!is.null(cons$per_fit_tables) && !is.null(cons$per_fit_tables[[model_name]])) {
        pt <- cons$per_fit_tables[[model_name]]
        L <- as.data.frame(pt$loadings)
        R <- as.data.frame(pt$fac_corr)
        
        # Standardized loading column
        if (!"std_loading" %in% names(L)) {
          sc <- pick_std_col(L)
          if (is.na(sc)) stop("No standardized loading column found in per_fit_tables for ", model_name)
          L$std_loading <- as.numeric(L[[sc]])
        } else {
          L$std_loading <- as.numeric(L$std_loading)
        }
        
        # Standardized corr column
        if (nrow(R)) {
          if (!"std_corr" %in% names(R)) {
            rc <- pick_std_col(R)
            R$std_corr <- if (!is.na(rc)) as.numeric(R[[rc]]) else NA_real_
          } else {
            R$std_corr <- as.numeric(R$std_corr)
          }
          R <- canonicalize_corrs(R)
        } else {
          R <- data.frame(lhs=character(), rhs=character(), std_corr=numeric())
        }
        
        # Guarantee Puberty ~~ Internal appears
        R <- ensure_PI_edge(R, model_name, cons$fits[[model_name]])
        return(list(loadings = L[, c("lhs","rhs","std_loading")], corrs = R))
      }
      
      # 2) Derive from lavaan_output if available
      if (!is.null(cons$fits) && !is.null(cons$fits[[model_name]])) {
        f <- cons$fits[[model_name]]
        if (!is.null(f$lavaan_output)) {
          pe <- lavaan::parameterEstimates(f$lavaan_output, standardized = TRUE)
          if (!"STD_Genotype" %in% names(pe) && "std.all" %in% names(pe)) pe$STD_Genotype <- pe$std.all
          
          L <- subset(pe, op == "=~", select = c("lhs","rhs","STD_Genotype","est","estimate","std.all"))
          lc <- pick_std_col(L)
          L$std_loading <- as.numeric(L[[lc]])
          L <- L[, c("lhs","rhs","std_loading")]
          
          latents <- unique(subset(pe, op == "=~")$lhs)
          R <- subset(pe, op == "~~" & lhs %in% latents & rhs %in% latents & lhs != rhs,
                      select = c("lhs","rhs","STD_Genotype","est","estimate","std.all"))
          rc <- pick_std_col(R)
          if (nrow(R)) {
            R$std_corr <- as.numeric(R[[rc]])
            R <- canonicalize_corrs(R)
          } else {
            R <- data.frame(lhs=character(), rhs=character(), std_corr=numeric())
          }
          
          R <- ensure_PI_edge(R, model_name, f)
          return(list(loadings = L, corrs = R))
        }
      }
    }
  }
  
  # ---------- Fallback to CSVs ----------
  L_csv <- file.path(OUT_DIR, paste0(model_name, "_loadings.csv"))
  R_csv <- file.path(OUT_DIR, paste0(model_name, "_factor_corr.csv"))
  if (!file.exists(L_csv)) stop("Could not find consolidated RDS or ", L_csv, " for model ", model_name)
  
  L <- fread(L_csv)
  if (!"std_loading" %in% names(L)) {
    lc <- pick_std_col(L)
    if (is.na(lc)) stop("No standardized loading column in ", L_csv)
    L$std_loading <- as.numeric(L[[lc]])
  } else L$std_loading <- as.numeric(L$std_loading)
  
  if (file.exists(R_csv)) {
    R <- fread(R_csv)
    if (!"std_corr" %in% names(R)) {
      rc <- pick_std_col(R)
      R$std_corr <- if (!is.na(rc)) as.numeric(R[[rc]]) else NA_real_
    } else R$std_corr <- as.numeric(R$std_corr)
    R <- canonicalize_corrs(R)
  } else {
    R <- data.frame(lhs=character(), rhs=character(), std_corr=numeric())
  }
  
  # With CSVs we can't inspect constraints; just add the edge with NA (will print blank)
  list(loadings = L[, c("lhs","rhs","std_loading")], corrs = R)
}


build_dot <- function(L, R, model_name, min_loading = MIN_LOADING_TO_DRAW) {
  Lf <- L[!is.na(L$std_loading) & abs(L$std_loading) >= min_loading, , drop=FALSE]
  latents  <- unique(L$lhs)
  observed <- unique(L$rhs)
  
  dot <- "digraph pathDiagram {\n"
  dot <- paste0(dot,
                "  graph [layout = dot, rankdir = TB, overlap = false, splines = true, ranksep = 1.4, nodesep = 0.6];\n",
                "  node  [fontname = Helvetica, fontsize = 12];\n",
                "  edge  [fontname = Helvetica, fontsize = 11];\n"
  )
  
  # Latent nodes
  for (lat in latents) {
    dot <- paste0(dot, "  \"", lat, "\" [shape = ellipse, style = filled, fillcolor = lightblue, label = \"", lat, "\"];\n")
  }
  # Observed nodes
  for (obs in observed) {
    dot <- paste0(dot, "  \"", obs, "\" [shape = box, label = \"", obs, "\"];\n")
  }
  
  # Latent -> observed (loadings)
  for (i in seq_len(nrow(Lf))) {
    lbl <- format_coef(Lf$std_loading[i])
    dot <- paste0(dot, "  \"", Lf$lhs[i], "\" -> \"", Lf$rhs[i], "\" [label = \"", lbl, "\"];\n")
  }
  
  # Factor correlations (dashed, double-headed)
  if (!is.null(R) && nrow(R)) {
    for (i in seq_len(nrow(R))) {
      lbl <- format_coef(R$std_corr[i])
      # If it's exactly zero (or numerically tiny), show "0"
      if (lbl == "" && !is.na(R$std_corr[i]) && abs(R$std_corr[i]) < 1e-12) lbl <- "0"
      # If constrained model, annotate Puberty–Internal as fixed when zero-ish
      is_constrained <- grepl("constrained", model_name, ignore.case = TRUE)
      pi_pair <- (R$lhs[i] %in% c("Puberty","Internal")) && (R$rhs[i] %in% c("Puberty","Internal"))
      if (is_constrained && pi_pair && !is.na(R$std_corr[i]) && abs(R$std_corr[i]) < 1e-12) lbl <- "0 (fixed)"
      
      dot <- paste0(
        dot,
        "  \"", R$lhs[i], "\" -> \"", R$rhs[i],
        "\" [dir=both, arrowtail=normal, arrowhead=normal, style=dashed, label = \"", lbl, "\"];\n"
      )
    }
  }
  
  # Rank constraints
  dot <- paste0(dot, "  { rank = min; ", paste(sprintf("\"%s\"", latents), collapse="; "), " }\n")
  dot <- paste0(dot, "  { rank = max; ", paste(sprintf("\"%s\"", observed), collapse="; "), " }\n")
  dot <- paste0(dot, "}\n")
  dot
}

render_and_save <- function(model_name, L, R) {
  dot <- build_dot(L, R, model_name)
  g <- DiagrammeR::grViz(dot)
  svg_txt <- DiagrammeRsvg::export_svg(g)
  
  svg_file <- file.path(DIAG_DIR, paste0(model_name, "_pathDiagram.svg"))
  png_file <- file.path(DIAG_DIR, paste0(model_name, "_pathDiagram.png"))
  cat(svg_txt, file = svg_file)
  rsvg::rsvg_png(svg_file, png_file)
  message("Saved: ", basename(svg_file), " and ", basename(png_file))
}

discover_models <- function() {
  # Try consolidated RDS first
  cons_path <- file.path(OUT_DIR, "cfa_DWLS_results.RDS")
  rds_models <- character(0)
  if (file.exists(cons_path)) {
    cons <- try(readRDS(cons_path), silent = TRUE)
    if (!inherits(cons, "try-error")) {
      if (!is.null(cons$fits)) rds_models <- names(cons$fits)
      if (length(rds_models) == 0 && !is.null(cons$per_fit_tables)) rds_models <- names(cons$per_fit_tables)
    }
  }
  # Also include any CSV-derived models
  csv_models <- list.files(OUT_DIR, pattern = "_loadings\\.csv$", full.names = FALSE)
  csv_models <- sub("_loadings\\.csv$", "", csv_models)
  sort(unique(c(rds_models, csv_models)))
}

# -----------------------
# Main
# -----------------------
all_models <- if (length(MODEL_SELECTION)) MODEL_SELECTION else discover_models()
if (!length(all_models)) stop("No models found to render. Check outputs in: ", OUT_DIR)

for (mn in all_models) {
  tab <- extract_tables_from_fit(mn, consolidated_path = file.path(OUT_DIR, "cfa_DWLS_results.RDS"))
  render_and_save(mn, tab$loadings, tab$corrs)
}

message("All diagrams saved to: ", DIAG_DIR)
