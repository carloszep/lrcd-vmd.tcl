#!/usr/bin/env python3
"""
cluster_heatmap.py

Generate a clustered heatmap (hierarchical clustering + dendrogram) from a
square matrix stored in a CSV or TSV file.

Usage:
    python3 cluster_heatmap.py <input_file> [colorbar_label] [-o OUTPUT] [--cmap CMAP]

Arguments:
    input_file       Path to a .csv or .tsv file containing a square matrix.
                      First row = column labels, first column = row labels.
                      Delimiter is auto-detected from the file extension
                      (.csv -> comma, .tsv/.txt -> tab) unless overridden
                      with --sep.
    colorbar_label    Label for the colorbar (default: "Value").

Options:
    -o, --output      Output PNG path (default: derived from input filename).
    --sep             Force a specific delimiter (e.g. ',' or '\\t').
    --cmap            Force a specific matplotlib/seaborn colormap, overriding
                       automatic selection.
    --method          Linkage method for hierarchical clustering (default: average).

Automatic colormap and clustering behaviour:
    If all values lie within [0, 1], the matrix is treated as a similarity-like
    metric (higher = more similar). The colormap "viridis" is used (high values
    = bright/yellow). For clustering, values are converted to distances via
    (1 - similarity) so that higher similarity correctly produces shorter
    dendrogram branches (i.e. closer grouping).
    If values fall outside [0, 1] (e.g. distances, RMSD, etc.), the matrix is
    treated as a dissimilarity/distance metric (higher = more different) and
    the reversed colormap "viridis_r" is used (low values = bright/yellow).
    The display matrix (heatmap cell values) is never modified by this
    conversion; only the internal clustering matrix is affected.

Required libraries: pandas numpy seaborn matplotlib scipy
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


def detect_separator(path, forced_sep=None):
    if forced_sep is not None:
        return forced_sep
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        return ","
    if ext in (".tsv", ".txt"):
        return "\t"
    # Fallback: sniff the first line
    with open(path, "r") as f:
        first_line = f.readline()
    return "," if first_line.count(",") >= first_line.count("\t") else "\t"


def load_matrix(path, forced_sep=None):
    sep = detect_separator(path, forced_sep)
    df = pd.read_csv(path, sep=sep, index_col=0)

    # Strip stray whitespace from labels (rows and columns)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()

    if df.shape[0] != df.shape[1]:
        raise ValueError(
            f"Input matrix is not square: {df.shape[0]} rows vs {df.shape[1]} columns."
        )

    mat = df.values.astype(float)
    if not np.allclose(mat, mat.T, atol=1e-3, equal_nan=True):
        print("Warning: matrix is not perfectly symmetric (tolerance 1e-3).", file=sys.stderr)

    return df, mat


def choose_colormap(mat, forced_cmap=None):
    if forced_cmap is not None:
        return forced_cmap
    finite_vals = mat[np.isfinite(mat)]
    if finite_vals.size and finite_vals.min() >= 0.0 and finite_vals.max() <= 1.0:
        # Similarity-like metric: higher = more similar -> dark/hot at high values
        return "viridis"
    # Distance-like metric: higher = more different -> dark/hot at low values
    return "viridis_r"


def build_title(metric_label, input_path):
    base = os.path.splitext(os.path.basename(input_path))[0]
    return f"Hierarchical clustering heatmap ({metric_label}) - {base}"


def main():
    parser = argparse.ArgumentParser(
        description="Generate a clustered heatmap with dendrogram from a square matrix (CSV/TSV)."
    )
    parser.add_argument("input_file", help="Path to the input .csv or .tsv file.")
    parser.add_argument(
        "colorbar_label",
        nargs="?",
        default="Value",
        help='Label for the colorbar (default: "Value").',
    )
    parser.add_argument("-o", "--output", default=None, help="Output PNG path.")
    parser.add_argument("--sep", default=None, help="Force a delimiter, e.g. ',' or '\\t'.")
    parser.add_argument(
        "--cmap", default=None, help="Force a colormap, overriding automatic selection."
    )
    parser.add_argument(
        "--method",
        default="average",
        help="Linkage method for hierarchical clustering (default: average).",
    )
    args = parser.parse_args()

    df, mat = load_matrix(args.input_file, args.sep)

    # Detect metric type BEFORE modifying the matrix
    finite_vals = mat[np.isfinite(mat)]
    is_similarity = finite_vals.size > 0 and finite_vals.min() >= 0.0 and finite_vals.max() <= 1.0

    # Build matrix for clustering only (never used for display):
    #   - similarity metrics: convert to distance (1 - sim) so that higher
    #     similarity correctly maps to shorter linkage distance.
    #   - distance metrics: use as-is; zero out diagonal for squareform safety.
    if is_similarity:
        mat_clust = 1.0 - mat
    else:
        mat_clust = mat.copy()
    np.fill_diagonal(mat_clust, 0.0)

    condensed = squareform(mat_clust, checks=False)
    Z = linkage(condensed, method=args.method)

    cmap = choose_colormap(mat, args.cmap)
    title = build_title(args.colorbar_label, args.input_file)

    n = df.shape[0]
    figsize = (max(8, n * 0.55), max(7, n * 0.5))
    annot_size = max(5, min(9, 90 // n))

    g = sns.clustermap(
        df,
        row_linkage=Z,
        col_linkage=Z,
        cmap=cmap,
        figsize=figsize,
        annot=True,
        fmt=".2f" if mat[np.isfinite(mat)].max() <= 1.0 else ".0f",
        annot_kws={"size": annot_size},
        linewidths=0.4,
        cbar_kws={"label": args.colorbar_label},
        dendrogram_ratio=(0.18, 0.18),
        cbar_pos=(0.02, 0.83, 0.03, 0.15),
    )

    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.fig.suptitle(title, y=1.02, fontsize=13)

    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    if args.output:
        output_path = args.output
    else:
        base = os.path.splitext(os.path.basename(args.input_file))[0]
        input_dir = os.path.dirname(os.path.abspath(args.input_file))
        output_path = os.path.join(input_dir, f"heatmap_cluster_{base}.png")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    g.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {output_path}")
    print(f"Colormap used: {cmap}")
    print(f"Linkage method: {args.method}")


if __name__ == "__main__":
    main()

