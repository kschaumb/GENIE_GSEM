# ---- Small helpers ----
read_if <- function(fname) readr::read_csv(file.path(out_dir, fname), show_col_types = FALSE)

load_efa_csv <- function(file) {
  df <- read_if(file)
  if (is.na(names(df)[1]) || names(df)[1] %in% c("", "...1")) names(df)[1] <- "Trait"
  df %>% column_to_rownames("Trait")
}

round_df <- function(df, digits = 2) {
  num <- sapply(df, is.numeric); df[num] <- lapply(df[num], round, digits); df
}

order_by_top_factor <- function(L) {
  A <- as.matrix(L)
  ord <- order(apply(abs(A), 1, which.max), -apply(abs(A), 1, max), rownames(L))
  L[ord, , drop = FALSE]
}

fmt_sci2 <- function(x) if (is.numeric(x)) sprintf("%.2e", x) else x

add_sheet_simple <- function(wb, sheet, df) {
  addWorksheet(wb, sheet); writeData(wb, sheet, df, withFilter = TRUE)
  freezePane(wb, sheet, firstRow = TRUE); setColWidths(wb, sheet, 1:ncol(df), "auto")
}

write_block <- function(wb, sheet, title, df, start_row) {
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:(ncol(df) + 1))
  tbl <- data.frame(Trait = rownames(df), df, row.names = NULL, check.names = FALSE)
  writeData(wb, sheet, tbl, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  num_cols <- which(sapply(tbl, is.numeric))
  if (length(num_cols)) conditionalFormatting(
    wb, sheet, cols = num_cols,
    rows = (start_row + 2):(start_row + 1 + nrow(tbl)),
    style = c("#4575b4", "#ffffff", "#d73027"), type = "colorScale"
  )
  start_row + 1 + nrow(tbl) + 2
}

write_block_center0 <- function(wb, sheet, title, df, start_row) {
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:(ncol(df) + 1))
  tbl <- data.frame(Trait = rownames(df), df, row.names = NULL, check.names = FALSE)
  writeData(wb, sheet, tbl, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  
  num_cols <- which(sapply(tbl, is.numeric))
  if (length(num_cols)) {
    vals <- unlist(tbl[num_cols], use.names = FALSE)
    M <- max(abs(vals), na.rm = TRUE)
    if (is.finite(M) && M > 0) {
      conditionalFormatting(
        wb, sheet, cols = num_cols,
        rows = (start_row + 2):(start_row + 1 + nrow(tbl)),
        type  = "colorScale",
        style = c("#4575b4", "#ffffff", "#d73027"),
        rule  = c(-M, 0, M)
      )
    }
  }
  setColWidths(wb, sheet, 1:ncol(tbl), "auto")
  start_row + 1 + nrow(tbl) + 2
}

write_block_plain <- function(wb, sheet, title, df, start_row) {
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:ncol(df))
  writeData(wb, sheet, df, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  setColWidths(wb, sheet, 1:ncol(df), "auto")
  start_row + 1 + nrow(df) + 2
}

write_long_center0 <- function(wb, sheet, title, df, start_row) {
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:5)
  writeData(wb, sheet, df, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  if (nrow(df)) {
    M <- max(abs(df$Loading), na.rm = TRUE)
    if (is.finite(M) && M > 0) conditionalFormatting(
      wb, sheet, cols = 3, rows = (start_row + 2):(start_row + 1 + nrow(df)),
      type = "colorScale", style = c("#4575b4", "#ffffff", "#d73027"), rule = c(-M, 0, M)
    )
  }
  setColWidths(wb, sheet, 1:5, "auto")
  start_row + 1 + nrow(df) + 2
}

write_corr_center0 <- function(wb, sheet, title, Phi, start_row) {
  if (is.null(Phi)) {
    writeData(wb, sheet, paste0(title, " (not available)"), startCol = 1, startRow = start_row)
    return(start_row + 2)
  }
  Phi <- round(Phi, 3)
  if (is.null(colnames(Phi))) colnames(Phi) <- paste0("F", seq_len(ncol(Phi)))
  if (is.null(rownames(Phi))) rownames(Phi) <- colnames(Phi)
  df <- tibble::rownames_to_column(as.data.frame(Phi, check.names = FALSE), "Factor")
  addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = start_row, cols = 1)
  writeData(wb, sheet, title, startCol = 1, startRow = start_row)
  mergeCells(wb, sheet, rows = start_row, cols = 1:ncol(df))
  writeData(wb, sheet, df, startCol = 1, startRow = start_row + 1, withFilter = TRUE)
  conditionalFormatting(
    wb, sheet, cols = 2:ncol(df),
    rows = (start_row + 2):(start_row + 1 + nrow(df)),
    type = "colorScale", style = c("#4575b4", "#ffffff", "#d73027"), rule = c(-1, 0, 1)
  )
  setColWidths(wb, sheet, 1:ncol(df), "auto")
  start_row + 1 + nrow(df) + 2
}

# Minimal picker / null-coalescer
pick_col <- function(nm, choices) { hit <- intersect(choices, nm); if (length(hit)) hit[1] else NA_character_ }
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Choose best fit from $fits ----
.pick_fit <- function(nms) {
  for (pat in c("^constrained_DWLS_floor001$",
                "^constrained_DWLS_base$",
                "^constrained",
                "^orig_DWLS_floor001$",
                "^orig_DWLS_base$",
                "^orig")) {
    hit <- grep(pat, nms, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  nms[1]
}

# ---- Final constrained loader (reads from $fits; returns long 'results' if available) ----
final_constrained_loadings <- function(rds_path) {
  obj  <- readRDS(rds_path)
  fits <- if (is.list(obj$fits)) obj$fits else obj
  if (!length(fits)) return(data.frame())
  
  fit_name <- .pick_fit(names(fits))
  fit <- fits[[fit_name]]
  
  if (is.list(fit) && is.data.frame(fit$results) && all(c("lhs","op","rhs") %in% names(fit$results))) {
    return(fit$results)  # long parameter table; normalize later
  }
  
  L <- fit$loadings %||% fit$L %||% fit$std_loadings %||% fit$load
  if (!is.null(L)) return(as.data.frame(L, check.names = FALSE))
  
  L <- obj$loadings %||% obj$L %||% obj$std_loadings %||% obj$load
  as.data.frame(L, check.names = FALSE)
}

# ---- Φ extractor (reads from $fits; builds from latent '~~' in results) ----
extract_phi_matrix <- function(rds_path) {
  obj  <- readRDS(rds_path)
  fits <- if (is.list(obj$fits)) obj$fits else obj
  if (!length(fits)) return(NULL)
  
  fit_name <- .pick_fit(names(fits))
  fit <- fits[[fit_name]]
  
  if (is.matrix(fit$Phi)) return(fit$Phi)
  if (!is.null(fit$lavaan_output) && inherits(fit$lavaan_output, "lavaan")) {
    ph <- try(lavaan::inspect(fit$lavaan_output, "cor.lv"), silent = TRUE)
    if (!inherits(ph, "try-error") && is.matrix(ph)) return(ph)
  }
  
  if (!is.data.frame(fit$results) || !all(c("lhs","op","rhs") %in% names(fit$results))) return(NULL)
  res <- fit$results
  
  latents <- unique(res$lhs[res$op == "=~"])
  latents <- latents[!is.na(latents) & nzchar(latents)]
  if (!length(latents)) return(NULL)
  
  std_col <- pick_col(names(res), c("STD_Genotype","Std.all","STD_All","std.all"))
  if (is.na(std_col)) return(NULL)
  
  phi_long <- subset(res, op == "~~" & lhs %in% latents & rhs %in% latents,
                     select = c("lhs","rhs", std_col))
  names(phi_long)[3] <- "val"
  if (!nrow(phi_long)) return(NULL)
  
  Phi <- matrix(0, nrow = length(latents), ncol = length(latents),
                dimnames = list(latents, latents))
  diag(Phi) <- 1
  for (i in seq_len(nrow(phi_long))) {
    a <- phi_long$lhs[i]; b <- phi_long$rhs[i]
    v <- suppressWarnings(as.numeric(phi_long$val[i]))
    if (is.finite(v)) { Phi[a,b] <- v; Phi[b,a] <- v }
  }
  Phi
}

# ---- Normalize to 5 cols (robust SE/p; keeps strings like "< 5e-300") ----
# ---- Normalize to 5 cols (ONLY loadings; robust SE/p; keeps strings like "< 5e-300") ----
normalize_cfa_output <- function(df) {
  nm <- names(df)
  
  fmt_numeric_or_keep <- function(x, digits = 2) {
    if (is.null(x)) return(rep(NA_character_, nrow(df)))
    x_num <- suppressWarnings(as.numeric(x))
    if (all(is.na(x_num)) && is.character(x)) {
      x_num <- suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", x)))
    }
    out <- as.character(x)
    ok <- is.finite(x_num)
    out[ok] <- sprintf(paste0("%.", digits, "e"), x_num[ok])
    out
  }
  
  if (all(c("lhs","rhs") %in% nm)) {
    # LONG FORM: keep ONLY loadings rows
    df2 <- if ("op" %in% nm) df[df$op == "=~", , drop = FALSE] else df
    
    est <- pick_col(names(df2), c("std_loading","STD_Genotype","STD_All","Std.all","std.all","est","Estimate"))
    se  <- pick_col(names(df2), c("STD_Genotype_SE","SE","se","Unstand_SE"))
    pv  <- pick_col(names(df2), c("p_value","p","pvalue"))
    
    out <- df2[, c("rhs","lhs", est, se, pv)]
    names(out) <- c("Trait","Factor","Loading","SE","p")
    
    out$Loading <- suppressWarnings(as.numeric(out$Loading))
    out$SE <- fmt_numeric_or_keep(out$SE)
    out$p  <- fmt_numeric_or_keep(out$p)
    
  } else {
    # WIDE FORM
    if (!"Trait" %in% nm && !is.null(rownames(df))) df <- tibble::rownames_to_column(df, "Trait")
    num_cols <- names(df)[sapply(df, is.numeric)]
    num_cols <- setdiff(num_cols, c("SE","se","p","p_value","pvalue"))
    out <- tidyr::pivot_longer(df, cols = dplyr::all_of(num_cols),
                               names_to = "Factor", values_to = "Loading") %>%
      mutate(SE = NA_character_, p = NA_character_) %>%
      select(Trait, Factor, Loading, SE, p)
    out$Loading <- as.numeric(out$Loading)
  }
  
  # Order nicely within factor
  dplyr::arrange(out, Factor, dplyr::desc(abs(Loading)), Trait)
}


# ---- Residuals with Φ aligned to loading columns ----
compute_residuals <- function(loads_wide, phi = NULL) {
  if (!"Trait" %in% names(loads_wide)) stop("compute_residuals() needs a Trait column")
  
  factors <- setdiff(names(loads_wide), "Trait")
  A <- as.matrix(loads_wide[, factors, drop = FALSE])
  A[is.na(A)] <- 0
  
  # align/build Phi to match 'factors'
  if (is.null(phi)) {
    Phi <- diag(length(factors))
    dimnames(Phi) <- list(factors, factors)
  } else {
    Phi0 <- as.matrix(phi)
    Phi <- matrix(0, nrow = length(factors), ncol = length(factors),
                  dimnames = list(factors, factors))
    diag(Phi) <- 1
    rn <- rownames(Phi0); cn <- colnames(Phi0)
    if (!is.null(rn) && !is.null(cn)) {
      common <- intersect(factors, intersect(rn, cn))
      if (length(common)) {
        Phi[common, common] <- Phi0[common, common, drop = FALSE]
      }
    } else if (nrow(Phi0) == ncol(Phi0) && ncol(Phi0) == length(factors)) {
      Phi <- Phi0
      dimnames(Phi) <- list(factors, factors)
    }
  }
  
  h2 <- rowSums((A %*% Phi) * A)
  tibble::tibble(Trait = loads_wide$Trait, Residual = round(pmax(0, 1 - h2), 3))
}