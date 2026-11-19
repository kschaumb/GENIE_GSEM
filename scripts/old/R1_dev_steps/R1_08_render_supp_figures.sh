#!/usr/bin/env bash
# R&R step 08: render Supplementary Figures (S1-S4) with corrected captions.
# Output: gsem/gsem_output/Supplemental_Figures.R1.pdf  (original Supplemental_Figures.pdf untouched)
# Requires Quarto (https://quarto.org) and a LaTeX install (e.g., `quarto install tinytex`).
# Usage:   bash R1_08_render_supp_figures.sh
set -euo pipefail
cd "$HOME/Library/CloudStorage/Box-Box/EMBARK/GENIE_GSEM/gsem/gsem_output"
quarto render Supplemental_Figures.R1.qmd --to pdf 2>&1 | tee log_files/R1_08_render_supp_figures.log
ls -la Supplemental_Figures.R1.pdf
