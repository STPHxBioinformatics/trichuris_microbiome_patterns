#!/usr/bin/env python3
"""
Cluster patristic (cophenetic) distance matrices into discrete strain clusters
and summarize per-species clustering quality metrics - merged

Inputs:
  - A directory containing CSV/TSV square distance matrices (rows=cols=sample IDs).

Outputs (per input file <name>.*):
  - <name>_clusters.csv            : columns [sample_id, cluster_id]
  - (optional) <name>_leaforder.txt: dendrogram leaf order (one sample per line), if --write_leaf_order
  - clusters_summary.csv           : one-row-per-species summary with K, silhouette, intra/inter distances, etc.

Clustering:
  - Hierarchical clustering (SciPy) on precomputed distances.
  - Choose clusters via:
      * --k N                    : force K clusters
      * --height H               : cut tree at distance H
      * --auto-k                 : choose K by silhouette over [--k_min, --k_max]

Usage examples:
  # Fixed K
  python cluster_from_distance_matrices.py \
      --in_dir derived/cophenetic --glob "*_cophenetic.csv" --sep "," \
      --k 4 --out_dir derived/clusters

  # Height cutoff at patristic distance 0.02
  python cluster_from_distance_matrices.py \
      --in_dir derived/cophenetic --glob "*.csv" --height 0.02 \
      --out_dir derived/clusters

  # Auto-select K by silhouette (tries 2..12)
  python cluster_from_distance_matrices.py \
      --in_dir derived/cophenetic --glob "*.csv" --auto-k --k_min 2 --k_max 12 \
      --out_dir derived/clusters

Requirements:
  - pandas, numpy, scipy, scikit-learn
    (conda install -c conda-forge pandas numpy scipy scikit-learn)
"""

import argparse
import glob
import os
import sys
from typing import Optional, Tuple, List, Dict

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, cut_tree
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score


def _read_distance_matrix(path: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    if df.shape[0] != df.shape[1]:
        raise ValueError(f"Matrix is not square: {df.shape}")
    # Align columns to rows if sets equal but order differs
    if set(df.index) != set(df.columns):
        raise ValueError("Row/column labels differ.")
    df = df.loc[df.index, df.index]

    # Symmetrize + zero diagonal
    A = df.values.astype(float)
    A = (A + A.T) / 2.0
    np.fill_diagonal(A, 0.0)
    if np.isnan(A).any():
        raise ValueError("Matrix contains NaNs. Clean upstream or impute carefully.")

    df.loc[:, :] = A
    return df


def _linkage_from_distance(D: pd.DataFrame, method: str = "average") -> Tuple[np.ndarray, List[str]]:
    Y = squareform(D.values, checks=False)  # condensed vector
    Z = linkage(Y, method=method)
    return Z, list(D.index)


def _choose_k_by_silhouette(D: pd.DataFrame, Z: np.ndarray, k_min: int, k_max: int) -> Tuple[int, float]:
    """Pick K maximizing silhouette score using the precomputed distance matrix."""
    best_k, best_s = None, -1.0
    n = D.shape[0]
    k_max_eff = max(k_min, min(k_max, n - 1))
    for k in range(k_min, k_max_eff + 1):
        labels = cut_tree(Z, n_clusters=k).reshape(-1)
        try:
            s = silhouette_score(D.values, labels, metric="precomputed")
        except Exception:
            s = -1.0
        if s > best_s:
            best_k, best_s = k, s
    if best_k is None:
        best_k, best_s = 2, -1.0
    return best_k, best_s


def _cluster_labels(D: pd.DataFrame,
                    method: str = "average",
                    k: Optional[int] = None,
                    height: Optional[float] = None,
                    auto_k: bool = False,
                    k_min: int = 2,
                    k_max: int = 10) -> Tuple[pd.Series, np.ndarray, List[str], List[str], float]:
    """Return cluster labels (1..K), linkage, labels, leaf order, and silhouette score for final labels."""
    Z, labels = _linkage_from_distance(D, method=method)
    dend = dendrogram(Z, no_plot=True, labels=labels)
    leaf_order = dend["ivl"]

    silhouette = np.nan
    if auto_k and (k is None) and (height is None):
        k, silhouette = _choose_k_by_silhouette(D, Z, k_min, k_max)

    if k is not None and height is not None:
        raise ValueError("Specify only one of --k or --height (or use --auto-k).")

    if k is not None:
        labs = cut_tree(Z, n_clusters=k).reshape(-1)
        labs = pd.Series(labs, index=labels).astype(int) + 1
        if np.isnan(silhouette):
            try:
                silhouette = silhouette_score(D.values, labs.values, metric="precomputed")
            except Exception:
                silhouette = np.nan
        return labs, Z, labels, leaf_order, silhouette

    if height is not None:
        labs = fcluster(Z, t=height, criterion="distance")  # 1..K
        labs = pd.Series(labs, index=labels).astype(int)
        try:
            silhouette = silhouette_score(D.values, labs.values, metric="precomputed")
        except Exception:
            silhouette = np.nan
        return labs, Z, labels, leaf_order, silhouette

    # Default fallback: k=2
    labs = cut_tree(Z, n_clusters=2).reshape(-1)
    labs = pd.Series(labs, index=labels).astype(int) + 1
    try:
        silhouette = silhouette_score(D.values, labs.values, metric="precomputed")
    except Exception:
        silhouette = np.nan
    return labs, Z, labels, leaf_order, silhouette


def _pairwise_mask(tri: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Upper-triangle mask without diagonal; returns mask and indices."""
    n = tri.shape[0]
    iu1 = np.triu_indices(n, k=1)
    return iu1, iu1  # kept for clarity


def _cluster_metrics(D: pd.DataFrame, labels: pd.Series) -> Dict[str, float]:
    """Compute intra-/inter-cluster means and separation margin."""
    X = D.values
    n = X.shape[0]
    labs = labels.values
    # Upper triangle indices
    iu = np.triu_indices(n, 1)

    # For each pair, intra if same label else inter
    same = (labs[iu[0]] == labs[iu[1]])
    intra_vals = X[iu][same]
    inter_vals = X[iu][~same]

    mean_intra = float(np.nanmean(intra_vals)) if intra_vals.size else np.nan
    mean_inter = float(np.nanmean(inter_vals)) if inter_vals.size else np.nan
    ratio = mean_intra / mean_inter if (np.isfinite(mean_intra) and np.isfinite(mean_inter) and mean_inter > 0) else np.nan

    # Separation margin: min(inter) - max(intra)
    min_inter = float(np.nanmin(inter_vals)) if inter_vals.size else np.nan
    max_intra = float(np.nanmax(intra_vals)) if intra_vals.size else np.nan
    margin = min_inter - max_intra if (np.isfinite(min_inter) and np.isfinite(max_intra)) else np.nan

    return {
        "mean_intra": mean_intra,
        "mean_inter": mean_inter,
        "intra_over_inter": ratio,
        "min_inter_minus_max_intra": margin
    }


def main():
    ap = argparse.ArgumentParser(description="Hierarchical clustering of patristic distance matrices with per-species metrics.")
    ap.add_argument("--in_dir", required=True, help="Directory with distance matrices (CSV/TSV).")
    ap.add_argument("--out_dir", required=True, help="Directory to write cluster assignments and summary.")
    ap.add_argument("--glob", default="*.csv", help="Glob pattern (e.g., '*_cophenetic.csv').")
    ap.add_argument("--sep", default=",", choices=[",", "tab"], help="Field separator for input matrices.")
    ap.add_argument("--linkage", default="average", choices=["single", "complete", "average", "weighted", "ward"],
                    help="Linkage method (avoid 'ward' for precomputed patristic distances).")
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--k", type=int, help="Number of clusters.")
    group.add_argument("--height", type=float, help="Distance cutoff for clusters.")
    ap.add_argument("--auto-k", action="store_true", help="Select K by silhouette over [--k_min, --k_max].")
    ap.add_argument("--k_min", type=int, default=2, help="Minimum K for auto selection.")
    ap.add_argument("--k_max", type=int, default=10, help="Maximum K for auto selection.")
    ap.add_argument("--write_leaf_order", action="store_true", help="Also write dendrogram leaf order per species.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    args = ap.parse_args()

    in_dir = os.path.abspath(args.in_dir)
    out_dir = os.path.abspath(args.out_dir)
    sep = "\t" if args.sep == "tab" else ","
    pattern = os.path.join(in_dir, args.glob)

    os.makedirs(out_dir, exist_ok=True)
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[ERROR] No files match: {pattern}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Found {len(files)} matrix file(s). Writing to: {out_dir}")

    summary_rows = []
    failures = []

    for f in files:
        base = os.path.basename(f)
        stem = os.path.splitext(base)[0]
        out_clusters = os.path.join(out_dir, f"{stem}_clusters.csv")
        out_order = os.path.join(out_dir, f"{stem}_leaforder.txt")

        if os.path.exists(out_clusters) and not args.overwrite:
            print(f"[SKIP] {base} -> {out_clusters} (exists; use --overwrite)")
            continue

        try:
            D = _read_distance_matrix(f, sep=sep)
            labels, Z, lablist, leaf_order, sil = _cluster_labels(
                D,
                method=args.linkage,
                k=args.k,
                height=args.height,
                auto_k=args.auto_k,
                k_min=args.k_min,
                k_max=args.k_max
            )

            # Write clusters
            out_df = labels.reset_index()
            out_df.columns = ["sample_id", "cluster_id"]
            out_df.to_csv(out_clusters, index=False)

            # Optional leaf order
            if args.write_leaf_order:
                with open(out_order, "w") as fh:
                    for s in leaf_order:
                        fh.write(f"{s}\n")

            # Metrics
            mets = _cluster_metrics(D, labels)
            n = D.shape[0]
            K = int(labels.nunique())

            mode = "fixed_k" if args.k else ("height" if args.height is not None else ("auto_k" if args.auto_k else "default_k2"))
            height_used = args.height if args.height is not None else np.nan

            summary_rows.append({
                "matrix_file": base,
                "species_stem": stem,
                "n_samples": n,
                "K": K,
                "linkage": args.linkage,
                "selection_mode": mode,
                "height_cutoff": height_used,
                "silhouette": sil,
                **mets
            })

            print(f"[OK]   {base} -> {out_clusters}  (K={K}, n={n}, silhouette={sil if np.isfinite(sil) else 'NA'})")

        except Exception as e:
            print(f"[FAIL] {base}: {e}", file=sys.stderr)
            failures.append((base, str(e)))

    # Write summary
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(out_dir, "clusters_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"[SUMMARY] Wrote {len(summary_rows)} rows to {summary_path}")

    if failures:
        print("\n[SUMMARY] Some files failed:", file=sys.stderr)
        for fname, err in failures:
            print(f"  - {fname}: {err}", file=sys.stderr)
        sys.exit(2)

    print("[DONE] All cluster files and summary written.")


if __name__ == "__main__":
    main()
