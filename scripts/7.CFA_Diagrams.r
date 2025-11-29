#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
})

BASE_DIR <- "~/Library/CloudStorage/Box-Box/GENIE_GSEM"

DIGITS <- 2

`%||%` <- function(x, y) if (!is.null(x)) x else y
stamp <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")

# CLEAN LABEL FORMATTER (key improvement)
format_coef <- function(x, digits = DIGITS) {
  if (is.na(x)) return("")
  r <- round(x, digits)
  if (abs(r) < (10^(-digits))) return("")  # hide values that round to zero
  fmt <- sprintf(paste0("%.", digits, "f"), r)
  fmt <- sub("^0\\.", ".",  fmt)
  fmt <- sub("^-0\\.", "-.", fmt)
  fmt
}

# -------------------------------------------------------------------
# BUILD DIAGRAM (polished layout & styling)
# -------------------------------------------------------------------
build_dot <- function(L, R, E, model_name) {
  
  lat <- unique(L$lhs)
  obs <- unique(L$rhs)
  
  dot <- "digraph G {\n"
  dot <- paste0(dot,
                "  graph [rankdir=TB, splines=true, overlap=false, nodesep=0.5, ranksep=0.9];\n",
                "  node  [fontname=Helvetica, fontsize=12, penwidth=1.2];\n",
                "  edge  [fontname=Helvetica, fontsize=11, penwidth=1.2, color=\"#333333\"];\n"
  )
  
  for (f in lat)
    dot <- paste0(dot, "  \"", f, "\" [shape=ellipse, style=filled, fillcolor=\"#B7D4E8\", penwidth=1.6];\n")
  
  for (v in obs)
    dot <- paste0(dot, "  \"", v, "\" [shape=box, style=rounded, penwidth=1.3, fillcolor=\"white\"];\n")
  
  for (v in E$var)
    dot <- paste0(dot, "  \"e_", v, "\" [shape=circle, width=0.28, height=0.28, fixedsize=true, style=filled, fillcolor=\"#DDDDDD\", color=\"#888888\", penwidth=0.7, label=\"\"];\n")
  
  for (i in seq_len(nrow(L))) {
    lbl <- format_coef(L$std_loading[i])
    dot <- paste0(dot, "  \"", L$lhs[i], "\" -> \"", L$rhs[i], "\" [label=\"", lbl, "\", penwidth=1.6];\n")
  }
  
  for (i in seq_len(nrow(R))) {
    lbl <- format_coef(R$std_corr[i])
    dot <- paste0(dot, "  \"", R$lhs[i], "\" -> \"", R$rhs[i], "\" [dir=both, style=dashed, color=\"#666666\", arrowsize=0.8, label=\"", lbl, "\", constraint=false, minlen=2];\n")
  }
  
  for (i in seq_len(nrow(E))) {
    lbl <- format_coef(E$resid_var_std[i])
    v   <- E$var[i]
    dot <- paste0(dot, "  \"e_", v, "\" -> \"", v, "\" [label=\"", lbl, "\", color=\"#777777\", fontcolor=\"#555555\", arrowsize=0.6, penwidth=0.9, minlen=0.4, constraint=true];\n")
  }
  
  dot <- paste0(dot,
                "  { rank = same; ", paste(sprintf("\"%s\"", lat), collapse="; "), " }\n",
                "  { rank = same; ", paste(sprintf("\"%s\"", obs), collapse="; "), " }\n"
  )
  for (v in E$var)
    dot <- paste0(dot, "  { rank = same; \"e_", v, "\"; \"", v, "\" }\n")
  
  dot <- paste0(dot, "}\n")
  dot
}

# ---------------------------------------------------
# TABLE EXTRACTION (simple)
# ---------------------------------------------------
extract_model <- function(cons, model_name) {
  
  pt <- cons$per_fit_tables[[model_name]]
  
  L <- pt$loadings |>
    mutate(std_loading = STD_Genotype %||% std_loading) |>
    select(lhs, rhs, std_loading)
  
  R <- pt$fac_corr |>
    mutate(std_corr = STD_Genotype %||% std_corr) |>
    select(lhs, rhs, std_corr)
  
  TH <- pt$theta
  E <- data.frame(
    var = TH$lhs,
    resid_var_std = as.numeric(TH$std_resid %||% TH$STD_Genotype %||% TH$STD_All),
    stringsAsFactors = FALSE
  )
  
  list(L=L, R=R, E=E)
}

# ---------------------------------------------------
# RENDER
# ---------------------------------------------------
render_model <- function(kind, model_name) {
  
  out_dir <- file.path(BASE_DIR, "gsem/gsem_output", paste0(kind, "_cfa"))
  cons <- readRDS(file.path(out_dir, "cfa_DWLS_results.RDS"))
  
  d <- extract_model(cons, model_name)
  dot <- build_dot(d$L, d$R, d$E, model_name)
  g <- DiagrammeR::grViz(dot)
  svg <- DiagrammeRsvg::export_svg(g)
  
  diag_dir <- file.path(out_dir, "diagrams")
  dir.create(diag_dir, recursive=TRUE, showWarnings=FALSE)
  
  svg_file <- file.path(diag_dir, paste0(model_name, ".svg"))
  png_file <- file.path(diag_dir, paste0(model_name, ".png"))
  
  cat(svg, file = svg_file)
  rsvg::rsvg_png(svg_file, png_file)
  
  stamp("Saved", basename(svg_file), "and", basename(png_file))
}

# ---------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
which <- if ("--which" %in% args) args[which(args=="--which")+1] else "both"

if (which %in% c("both","symptoms"))  render_model("symptoms",  "constrained_DWLS_floor001")
if (which %in% c("both","disorders")) render_model("disorder", "constrained_DWLS_floor001")

stamp("Done.")
