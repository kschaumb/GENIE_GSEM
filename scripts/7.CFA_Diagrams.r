#!/usr/bin/env Rscript
# ============================================================
# cfa_make_diagrams_master.R
# Path diagrams for GenomicSEM CFA fits (DiagrammeR)
# - SYMPTOMS  -> ~/Library/CloudStorage/Box-Box/GENIE_GSEM/gsem/gsem_output/symptoms_cfa/diagrams/
# - DISORDERS -> ~/Library/CloudStorage/Box-Box/GENIE_GSEM/gsem/gsem_output/disorder_cfa/diagrams/
# Prefers consolidated RDS (cfa_DWLS_results.RDS); falls back to CSVs
# CLI:
#   Rscript scripts/cfa_make_diagrams_master.R                         # both
#   Rscript scripts/cfa_make_diagrams_master.R --which symptoms        # symptoms only
#   Rscript scripts/cfa_make_diagrams_master.R --which disorders       # disorders only
#   Rscript scripts/cfa_make_diagrams_master.R --models orig_DWLS_floor001,constrained_DWLS_floor001
#     (optional, applies to whichever pipelines you run)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
  library(lavaan)
})

# -----------------------
# Config (shared defaults)
# -----------------------
BASE_DIR <- "~/Library/CloudStorage/Box-Box/GENIE_GSEM"
DEFAULT_MODEL_SELECTION <- c("constrained_DWLS_floor001")

MIN_LOADING_TO_DRAW <- 0.05
MIN_CORR_TO_DRAW    <- 0.05   # << new
MIN_LABEL_TO_PRINT  <- 0.00
DIGITS              <- 2

stamp <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

source("scripts/clean_labels.R")

# If you have a diagrammer graph object 'g', replace node labels:
# g <- set_node_attrs(g, "label", value = pretty[match(get_node_attrs(g, "name"), names(pretty))])

# -----------------------
# CLI
# -----------------------
args <- commandArgs(trailingOnly = TRUE)
which_to_run <- c("symptoms", "disorders")
models_cli <- NULL
if (length(args)) {
  i <- 1
  while (i <= length(args)) {
    if (identical(args[[i]], "--which") && i < length(args)) {
      which_to_run <- strsplit(args[[i + 1]], ",", fixed = TRUE)[[1]]
      which_to_run <- intersect(which_to_run, c("symptoms","disorders"))
    } else if (identical(args[[i]], "--models") && i < length(args)) {
      models_cli <- strsplit(args[[i + 1]], ",", fixed = TRUE)[[1]]
    }
    i <- i + 1
  }
}
if (!length(which_to_run)) which_to_run <- c("symptoms","disorders")

# -----------------------
# Helpers (shared)
# -----------------------
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

canonicalize_corrs <- function(R){
  if (!nrow(R)) return(R[, c("lhs","rhs","std_corr")])
  R$l_can <- pmin(R$lhs, R$rhs)
  R$r_can <- pmax(R$lhs, R$rhs)
  R <- R[!duplicated(R[, c("l_can","r_can")]), ]
  R <- R[, c("l_can","r_can","std_corr")]
  names(R) <- c("lhs","rhs","std_corr")
  R
}



build_dot <- function(L, R, model_name, min_loading = MIN_LOADING_TO_DRAW, min_corr = MIN_CORR_TO_DRAW) {
  Lf <- L[!is.na(L$std_loading) & abs(L$std_loading) >= min_loading, , drop = FALSE]
  latents  <- unique(L$lhs)
  observed <- unique(L$rhs)
  
  # NEW: filter latent-latent correlations
  Rf <- R
  if (!is.null(Rf) && nrow(Rf)) {
    Rf <- Rf[!is.na(Rf$std_corr) & abs(Rf$std_corr) >= min_corr, , drop = FALSE]
  }
  
  dot <- "digraph pathDiagram {\n"
  dot <- paste0(dot,
                "  graph [layout = dot, rankdir = TB, overlap = false, splines = true, ranksep = 1.4, nodesep = 0.6];\n",
                "  node  [fontname = Helvetica, fontsize = 12];\n",
                "  edge  [fontname = Helvetica, fontsize = 11];\n"
  )
  
  # nodes ...
  for (lat in latents)   dot <- paste0(dot, "  \"", lat, "\" [shape = ellipse, style = filled, fillcolor = lightblue, label = \"", lat, "\"];\n")
  for (obs in observed)  dot <- paste0(dot, "  \"", obs, "\" [shape = box, label = \"", obs, "\"];\n")
  
  # loadings
  for (i in seq_len(nrow(Lf))) {
    lbl <- format_coef(Lf$std_loading[i])
    dot <- paste0(dot, "  \"", Lf$lhs[i], "\" -> \"", Lf$rhs[i], "\" [label = \"", lbl, "\"];\n")
  }
  
  # correlations (use Rf)
  if (!is.null(Rf) && nrow(Rf)) {
    for (i in seq_len(nrow(Rf))) {
      lbl <- format_coef(Rf$std_corr[i])
      dot <- paste0(
        dot, "  \"", Rf$lhs[i], "\" -> \"", Rf$rhs[i],
        "\" [dir=both, arrowtail=normal, arrowhead=normal, style=dashed, label = \"", lbl, "\"];\n"
      )
    }
  }
  
  dot <- paste0(dot, "  { rank = min; ", paste(sprintf("\"%s\"", latents), collapse = "; "), " }\n")
  dot <- paste0(dot, "  { rank = max; ", paste(sprintf("\"%s\"", observed), collapse = "; "), " }\n")
  dot <- paste0(dot, "}\n")
  dot
}

render_and_save <- function(model_name, L, R, diag_dir) {
  # Use your cleaner on *all* names that appear
  all_symbols <- unique(c(L$lhs, L$rhs, R$lhs, R$rhs))
  pretty <- clean_recode(all_symbols)
  
  # Relabel in-place
  L$lhs <- pretty[L$lhs]
  L$rhs <- pretty[L$rhs]
  if (nrow(R)) {
    R$lhs <- pretty[R$lhs]
    R$rhs <- pretty[R$rhs]
  }
  
  # Now build and render as usual
  dot <- build_dot(L, R, model_name)
  g <- DiagrammeR::grViz(dot)
  svg_txt <- DiagrammeRsvg::export_svg(g)
  
  svg_file <- file.path(diag_dir, paste0(model_name, "_pathDiagram.svg"))
  png_file <- file.path(diag_dir, paste0(model_name, "_pathDiagram.png"))
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  cat(svg_txt, file = svg_file)
  rsvg::rsvg_png(svg_file, png_file)
  stamp("Saved:", basename(svg_file), "and", basename(png_file))
}

# Prefer consolidated RDS -> per_fit_tables -> fits$lavaan_output; fallback: CSVs
extract_tables_from_fit <- function(model_name, out_dir) {
  consolidated_path <- file.path(out_dir, "cfa_DWLS_results.RDS")
  
  # ---------- Try consolidated RDS ----------
  if (file.exists(consolidated_path)) {
    cons <- try(readRDS(consolidated_path), silent = TRUE)
    if (!inherits(cons, "try-error")) {
      # 1) Direct per-fit tables (saved by pipeline)
      if (!is.null(cons$per_fit_tables) && !is.null(cons$per_fit_tables[[model_name]])) {
        pt <- cons$per_fit_tables[[model_name]]
        L <- as.data.frame(pt$loadings)
        R <- as.data.frame(pt$fac_corr)
        
        # Standardized loading column
        if (!"std_loading" %in% names(L)) {
          sc <- pick_std_col(L)
          if (is.na(sc)) stop("No standardized loading column found in per_fit_tables for ", model_name)
          L$std_loading <- as.numeric(L[[sc]])
        } else L$std_loading <- as.numeric(L$std_loading)
        
        # Standardized corr column
        if (nrow(R)) {
          if (!"std_corr" %in% names(R)) {
            rc <- pick_std_col(R)
            R$std_corr <- if (!is.na(rc)) as.numeric(R[[rc]]) else NA_real_
          } else R$std_corr <- as.numeric(R$std_corr)
          R <- canonicalize_corrs(R)
        } else {
          R <- data.frame(lhs = character(), rhs = character(), std_corr = numeric())
        }
        
        # Ensure Puberty ~~ Internal is present for legibility (esp. constrained)
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
            R <- data.frame(lhs = character(), rhs = character(), std_corr = numeric())
          }
          
          R <- ensure_edge(R, "Puberty", "Internal", 0)
          return(list(loadings = L, corrs = R))
        }
      }
    }
  }
  
  # ---------- Fallback to CSVs ----------
  L_csv <- file.path(out_dir, paste0(model_name, "_loadings.csv"))
  R_csv <- file.path(out_dir, paste0(model_name, "_factor_corr.csv"))
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
    R <- data.frame(lhs = character(), rhs = character(), std_corr = numeric())
  }
  
  R <- ensure_edge(R, "Puberty", "Internal", 0)
  list(loadings = L[, c("lhs","rhs","std_loading")], corrs = R)
}

# -----------------------
# Pipeline runners
# -----------------------
run_pipeline <- function(kind, models_override = NULL) {
  stopifnot(kind %in% c("symptoms","disorders"))
  out_dir <- file.path(BASE_DIR, "gsem/gsem_output", if (kind == "symptoms") "symptoms_cfa" else "disorder_cfa")
  diag_dir <- file.path(out_dir, "diagrams")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  
  # If user provided --models, use those; else default list. If empty, auto-discover.
  MODEL_SELECTION <- if (length(models_override)) models_override else DEFAULT_MODEL_SELECTION
  
  discover_models <- function() {
    cons_path <- file.path(out_dir, "cfa_DWLS_results.RDS")
    rds_models <- character(0)
    if (file.exists(cons_path)) {
      cons <- try(readRDS(cons_path), silent = TRUE)
      if (!inherits(cons, "try-error")) {
        if (!is.null(cons$fits)) rds_models <- names(cons$fits)
        if (length(rds_models) == 0 && !is.null(cons$per_fit_tables)) rds_models <- names(cons$per_fit_tables)
      }
    }
    csv_models <- list.files(out_dir, pattern = "_loadings\\.csv$", full.names = FALSE)
    csv_models <- sub("_loadings\\.csv$", "", csv_models)
    sort(unique(c(rds_models, csv_models)))
  }
  
  all_models <- if (length(MODEL_SELECTION)) MODEL_SELECTION else discover_models()
  if (!length(all_models)) stop("No models found to render. Check outputs in: ", out_dir)
  
  stamp("Rendering", kind, "diagrams for:", paste(all_models, collapse = ", "))
  for (mn in all_models) {
    tab <- extract_tables_from_fit(mn, out_dir)
    render_and_save(mn, tab$loadings, tab$corrs, diag_dir)
  }
  stamp("Done:", kind, "->", diag_dir)
}

# -----------------------
# Execute selection
# -----------------------
if ("symptoms" %in% which_to_run)  run_pipeline("symptoms", models_cli)
if ("disorders" %in% which_to_run) run_pipeline("disorders", models_cli)
stamp("All requested diagram runs finished.")
