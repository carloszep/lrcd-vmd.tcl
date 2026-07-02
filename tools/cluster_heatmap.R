#!/usr/bin/env Rscript
#
# cluster_heatmap.R
#
# Generate a clustered heatmap (hierarchical clustering + dendrogram) from a
# square matrix stored in a CSV or TSV file. R equivalent of cluster_heatmap_cli.py.
#
# Usage:
#   Rscript cluster_heatmap.R <input_file> [colorbar_label] [-o OUTPUT] [--sep SEP] [--cmap CMAP] [--method METHOD]
#
# Arguments:
#   input_file       Path to a .csv or .tsv file containing a square matrix.
#                     First row = column labels, first column = row labels.
#                     Delimiter is auto-detected from the file extension
#                     (.csv -> comma, .tsv/.txt -> tab) unless overridden
#                     with --sep.
#   colorbar_label   Label for the colorbar (default: "Value").
#
# Options:
#   -o, --output     Output PNG path (default: heatmap_cluster_<input_basename>.png,
#                     saved next to the input file).
#   --sep            Force a specific delimiter (e.g. ',' or '\t').
#   --cmap           Force a specific colormap, overriding automatic selection.
#                     One of: "viridis", "viridis_r" (or any palette name
#                     accepted by build_palette()).
#   --method         Linkage method for hierarchical clustering (default: average).
#
# Automatic colormap and clustering behaviour (same as the Python version):
#   If all values lie within [0, 1], the matrix is treated as a similarity-like
#   metric (higher = more similar) -> "viridis" (high values = bright/yellow).
#   For clustering, values are converted to distances via (1 - similarity) so
#   that higher similarity correctly produces shorter dendrogram branches.
#   If values fall outside [0, 1] (e.g. distances, RMSD, etc.), the matrix is
#   treated as a dissimilarity/distance metric (higher = more different) ->
#   "viridis_r" (low values = bright/yellow, i.e. "hot" = similar).
#   The display matrix (heatmap cell values) is never modified by this
#   conversion; only the internal clustering matrix is affected.
#
# Required packages: optparse, pheatmap, viridisLite
# Install with:
#   install.packages(c("optparse", "pheatmap", "viridisLite"))

suppressWarnings(suppressMessages({
  library(optparse)
  library(pheatmap)
  library(viridisLite)
}))

## ---------------------------------------------------------------------------
## Argument parsing
## ---------------------------------------------------------------------------

option_list <- list(
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output PNG path (default: next to input file)."),
  make_option(c("--sep"), type = "character", default = NULL,
              help = "Force a delimiter, e.g. ',' or '\\t'."),
  make_option(c("--cmap"), type = "character", default = NULL,
              help = "Force a colormap ('viridis' or 'viridis_r'), overriding automatic selection."),
  make_option(c("--method"), type = "character", default = "average",
              help = "Linkage method for hierarchical clustering (default: average).")
)

parser <- OptionParser(
  usage = "Rscript %prog <input_file> [colorbar_label] [options]",
  option_list = option_list
)

parsed <- parse_args(parser, positional_arguments = c(1, 2))
opts <- parsed$options
positional <- parsed$args

input_file <- positional[1]
colorbar_label <- if (length(positional) >= 2) positional[2] else "Value"

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------

detect_separator <- function(path, forced_sep = NULL) {
  if (!is.null(forced_sep)) return(forced_sep)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(",")
  if (ext %in% c("tsv", "txt")) return("\t")
  # Fallback: sniff the first line
  first_line <- readLines(path, n = 1)
  n_comma <- lengths(regmatches(first_line, gregexpr(",", first_line)))
  n_tab <- lengths(regmatches(first_line, gregexpr("\t", first_line)))
  if (n_comma >= n_tab) "," else "\t"
}

load_matrix <- function(path, forced_sep = NULL) {
  sep <- detect_separator(path, forced_sep)

  df <- read.delim(path, sep = sep, header = TRUE, row.names = 1,
                    check.names = FALSE, stringsAsFactors = FALSE)

  # Strip stray whitespace from labels (rows and columns)
  rownames(df) <- trimws(rownames(df))
  colnames(df) <- trimws(colnames(df))

  mat <- as.matrix(df)
  mat <- apply(mat, c(1, 2), as.numeric)
  rownames(mat) <- rownames(df)
  colnames(mat) <- colnames(df)

  if (nrow(mat) != ncol(mat)) {
    stop(sprintf("Input matrix is not square: %d rows vs %d columns.",
                  nrow(mat), ncol(mat)))
  }

  if (!isTRUE(all.equal(mat, t(mat), tolerance = 1e-3))) {
    warning("Matrix is not perfectly symmetric (tolerance 1e-3).")
  }

  mat
}

build_palette <- function(cmap_name, n = 256) {
  if (cmap_name == "viridis") {
    return(viridisLite::viridis(n))
  } else if (cmap_name == "viridis_r") {
    return(rev(viridisLite::viridis(n)))
  }
  # Allow falling back to any palette name viridisLite supports
  if (exists(cmap_name, where = asNamespace("viridisLite"), mode = "function")) {
    pal_fun <- get(cmap_name, asNamespace("viridisLite"))
    return(pal_fun(n))
  }
  stop(sprintf("Unsupported colormap: %s", cmap_name))
}

choose_colormap <- function(mat, forced_cmap = NULL) {
  if (!is.null(forced_cmap)) return(forced_cmap)
  finite_vals <- mat[is.finite(mat)]
  if (length(finite_vals) > 0 && min(finite_vals) >= 0 && max(finite_vals) <= 1) {
    # Similarity-like metric: higher = more similar -> bright at high values
    return("viridis")
  }
  # Distance-like metric: higher = more different -> bright at low values
  return("viridis_r")
}

build_title <- function(metric_label, input_path) {
  base <- tools::file_path_sans_ext(basename(input_path))
  sprintf("Hierarchical clustering heatmap (%s) - %s", metric_label, base)
}

## ---------------------------------------------------------------------------
## Main
## ---------------------------------------------------------------------------

mat_orig <- load_matrix(input_file, opts$sep)  # used for display (keeps real diagonal values)

# Detect metric type BEFORE modifying the matrix
finite_vals <- mat_orig[is.finite(mat_orig)]
is_similarity <- length(finite_vals) > 0 && min(finite_vals) >= 0 && max(finite_vals) <= 1

# Build matrix for clustering only (never used for display):
#   - similarity metrics: convert to distance (1 - sim) so that higher
#     similarity correctly maps to shorter linkage distance.
#   - distance metrics: use as-is; zero out diagonal for as.dist() safety.
mat_clust <- mat_orig
if (is_similarity) {
  mat_clust <- 1.0 - mat_clust
}
diag(mat_clust) <- 0

dist_obj <- as.dist(mat_clust)
hc <- hclust(dist_obj, method = opts$method)

cmap_name <- choose_colormap(mat_orig, opts$cmap)
palette <- build_palette(cmap_name)
title <- build_title(colorbar_label, input_file)

n <- nrow(mat_orig)
is_fractional <- max(mat_orig[is.finite(mat_orig)]) <= 1.0
number_format <- if (is_fractional) "%.2f" else "%.0f"

# Scale figure size and font sizes with n, mirroring the Python version
fig_width <- max(8, n * 0.55)
fig_height <- max(7, n * 0.5)
annot_size <- max(5, min(9, 90 %/% n))

if (is.null(opts$output)) {
  base <- tools::file_path_sans_ext(basename(input_file))
  input_dir <- dirname(normalizePath(input_file, mustWork = FALSE))
  output_path <- file.path(input_dir, sprintf("heatmap_cluster_%s.png", base))
} else {
  output_path <- opts$output
}

out_dir <- dirname(output_path)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pheatmap(
  mat_orig,
  color = palette,
  cluster_rows = hc,
  cluster_cols = hc,
  display_numbers = matrix(sprintf(number_format, mat_orig), nrow = n),
  number_color = "white",
  fontsize_number = annot_size,
  fontsize = 10,
  main = title,
  legend = TRUE,
  border_color = "grey80",
  cellwidth = NA,
  cellheight = NA,
  width = fig_width,
  height = fig_height,
  filename = output_path,
  silent = TRUE
)

cat(sprintf("Saved: %s\n", output_path))
cat(sprintf("Colormap used: %s\n", cmap_name))
cat(sprintf("Linkage method: %s\n", opts$method))

