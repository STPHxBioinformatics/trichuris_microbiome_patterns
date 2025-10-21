#!/usr/bin/env python3
"""
Compute country/region homogeneity of strain clusters per species.

Inputs:
  - A directory of per-species cluster files (from your previous step), each with:
        sample_id,cluster_id
    e.g., Streptococcus_gallolyticus_cophenetic_clusters.csv
  - A metadata CSV/TSV file with EXACT headers:
        sampleID,country,region

Outputs:
  - geo_structure_summary.csv  (one row per species with metrics & class)

Metrics:
  - H_country: size-weighted mean of per-cluster majority country proportion (0..1; higher = purer)
  - H_region : same at region level
  - ARI_country / ARI_region (optional): adjusted Rand index of (clusters vs labels)
  - structure_type: 'country-specific' (H_country>=thr) /
                    'region-specific'  (H_country<thr & H_region>=thr) /
                    'no-geographic-structure' (else)

Usage:
  python3 compute_geo_homogeneity.py \
    --clusters_dir clusters \
    --glob "*_clusters.csv" \
    --metadata Country_metadata.txt \
    --meta_sep tab \
    --out_dir geo_structure \
    --hom_threshold 0.8 \
    --compute_ari

Requires: pandas, numpy, scikit-learn
"""

import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd

try:
    from sklearn.metrics import adjusted_rand_score
    HAVE_ARI = True
except Exception:
    HAVE_ARI = False


def read_metadata(path: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=("\t" if sep == "tab" else sep))
    needed = {"sampleID", "country", "region"}
    if not needed.issubset(df.columns):
        raise ValueError(f"Metadata must have columns: {needed}. Found: {set(df.columns)}")
    # normalize types/whitespace
    for c in ["sampleID", "country", "region"]:
        df[c] = df[c].astype(str).str.strip()
    return df[["sampleID", "country", "region"]].drop_duplicates()


def majority_prop(series: pd.Series) -> float:
    """Proportion of the majority label in a vector (0..1)."""
    if series.empty:
        return np.nan
    p = series.value_counts(normalize=True, dropna=False)
    return float(p.iloc[0])


def homogeneity_for_species(clust_df: pd.DataFrame, meta: pd.DataFrame) -> dict:
    """
    clust_df: columns ['sample_id','cluster_id','species_stem']
    meta: columns ['sampleID','country','region']
    """
    m = clust_df.merge(meta, left_on="sample_id", right_on="sampleID", how="inner")

    if m.empty:
        raise ValueError("After merging with metadata, no samples remain for this species.")

    # Ensure clean labels
    m["country"] = m["country"].astype(str).str.strip()
    m["region"]  = m["region"].astype(str).str.strip()

    # Per-cluster homogeneity at both levels
    per_cluster = (
        m.groupby("cluster_id")
         .agg(
             n=("sample_id", "size"),
             H_country=("country", majority_prop),
             H_region=("region", majority_prop)
         )
         .reset_index()
    )

    # Weighted averages (weights = cluster sizes)
    H_country = float(np.average(per_cluster["H_country"], weights=per_cluster["n"]))
    H_region  = float(np.average(per_cluster["H_region"],  weights=per_cluster["n"]))

    # ARI vs labels (optional)
    if HAVE_ARI:
        # Align arrays
        labs = m.sort_values(["sample_id", "cluster_id"])
        ari_country = adjusted_rand_score(labs["cluster_id"].astype(str), labs["country"].astype(str))
        ari_region  = adjusted_rand_score(labs["cluster_id"].astype(str), labs["region"].astype(str))
    else:
        ari_country = np.nan
        ari_region  = np.nan

    n_samples = int(m["sample_id"].nunique())
    K = int(m["cluster_id"].nunique())

    return {
        "n_samples": n_samples,
        "K": K,
        "H_country": H_country,
        "H_region": H_region,
        "ARI_country": ari_country,
        "ARI_region": ari_region
    }


def infer_species_stem(filename: str) -> str:
    # drop a trailing suffix like "_clusters"
    stem = os.path.splitext(os.path.basename(filename))[0]
    if stem.endswith("_clusters"):
        stem = stem[:-9]  # remove "_clusters"
    return stem


def main():
    ap = argparse.ArgumentParser(description="Compute country/region homogeneity of strain clusters per species.")
    ap.add_argument("--clusters_dir", required=True, help="Directory containing per-species *_clusters.csv files.")
    ap.add_argument("--glob", default="*_clusters.csv", help="Glob for cluster files (default '*_clusters.csv').")
    ap.add_argument("--metadata", required=True, help="Metadata file with columns: sampleID,country,region.")
    ap.add_argument("--meta_sep", default="tab", choices=[",", "tab", ";"], help="Separator for metadata (default: tab).")
    ap.add_argument("--out_dir", required=True, help="Directory to write summary CSV.")
    ap.add_argument("--hom_threshold", type=float, default=0.8, help="Threshold to classify homogeneity (default 0.8).")
    ap.add_argument("--compute_ari", action="store_true", help="If set, require ARI calculation (needs scikit-learn).")
    ap.add_argument("--min_samples", type=int, default=4, help="Skip species with < min samples (default 4).")
    args = ap.parse_args()

    if args.compute_ari and not HAVE_ARI:
        print("[ERROR] scikit-learn not available for ARI. Install scikit-learn or omit --compute_ari.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)
    meta = read_metadata(args.metadata, sep=args.meta_sep)

    pattern = os.path.join(os.path.abspath(args.clusters_dir), args.glob)
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[ERROR] No cluster files matched: {pattern}", file=sys.stderr)
        sys.exit(1)

    rows = []
    failures = []
    for f in files:
        try:
            df = pd.read_csv(f)
            # Expect columns sample_id, cluster_id; tolerate case variants
            cols = {c.lower(): c for c in df.columns}
            if "sample_id" not in cols or "cluster_id" not in cols:
                raise ValueError(f"File {os.path.basename(f)} must have columns sample_id,cluster_id")
            # normalize
            df = df.rename(columns={cols["sample_id"]: "sample_id", cols["cluster_id"]: "cluster_id"})
            df["sample_id"] = df["sample_id"].astype(str).str.strip()
            df["cluster_id"] = df["cluster_id"].astype(str).str.strip()

            species = infer_species_stem(f)
            df["species_stem"] = species

            # Filter by min samples
            if df["sample_id"].nunique() < args.min_samples:
                print(f"[SKIP] {os.path.basename(f)}: n<{args.min_samples}")
                continue

            metrics = homogeneity_for_species(df, meta)
            # classify
            if metrics["H_country"] >= args.hom_threshold:
                structure = "country-specific"
            elif metrics["H_region"] >= args.hom_threshold:
                structure = "region-specific"
            else:
                structure = "no-geographic-structure"

            rows.append({
                "species_stem": species,
                **metrics,
                "hom_threshold": args.hom_threshold,
                "structure_type": structure,
                "clusters_file": os.path.basename(f)
            })
            print(f"[OK] {os.path.basename(f)}: K={metrics['K']} H_country={metrics['H_country']:.2f} "
                  f"H_region={metrics['H_region']:.2f} -> {structure}")

        except Exception as e:
            print(f"[FAIL] {os.path.basename(f)}: {e}", file=sys.stderr)
            failures.append((f, str(e)))

    if rows:
        out_path = os.path.join(args.out_dir, "geo_structure_summary.csv")
        pd.DataFrame(rows).sort_values("species_stem").to_csv(out_path, index=False)
        print(f"[SUMMARY] Wrote {len(rows)} species to {out_path}")
    else:
        print("[SUMMARY] No species processed.", file=sys.stderr)

    if failures:
        print("\n[WARN] Some files failed:", file=sys.stderr)
        for f, err in failures:
            print(f"  - {os.path.basename(f)}: {err}", file=sys.stderr)


if __name__ == "__main__":
    main()
