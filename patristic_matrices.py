#!/usr/bin/env python3
"""
Compute patristic (cophenetic) distance matrices from Newick trees using Biopython.

Works on Python 3.13+ (no ete3 dependency).

Usage:
  python3 patristic_matrices.py \
      --in_dir /path/to/trees \
      --out_dir /path/to/output \
      --glob "*.tre" \
      --sep "," \
      --suffix "_cophenetic.csv"
"""

import argparse
import glob
import os
import sys
from typing import Dict, List

import pandas as pd
from Bio import Phylo

def _tip_dict(tree) -> Dict[str, "Clade"]:
    tips = tree.get_terminals()
    names = [t.name for t in tips]
    if any(n is None or n == "" for n in names):
        raise RuntimeError("Empty tip labels found; all tips must have sample IDs.")
    # check duplicates
    seen = set()
    dups = []
    for n in names:
        if n in seen:
            dups.append(n)
        seen.add(n)
    if dups:
        raise RuntimeError("Duplicate tip labels found: " + ", ".join(sorted(set(dups))))
    return {t.name: t for t in tips}

def compute_patristic_matrix(tree_path: str) -> pd.DataFrame:
    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        raise RuntimeError(f"Failed to parse '{tree_path}': {e}")

    name_to_tip = _tip_dict(tree)
    tips: List[str] = sorted(name_to_tip.keys())

    # Biopython's tree.distance(a, b) sums branch lengths along the path; requires branch lengths.
    # If a branch has length None, Biopython treats it as 0.0.
    n = len(tips)
    D = pd.DataFrame(0.0, index=tips, columns=tips)

    # Pre-resolve clades for speed
    tips_clades = {name: name_to_tip[name] for name in tips}

    for i in range(n):
        ai = tips[i]
        for j in range(i + 1, n):
            bj = tips[j]
            dist = tree.distance(tips_clades[ai], tips_clades[bj])
            D.iat[i, j] = dist
            D.iat[j, i] = dist

    return D

def main():
    ap = argparse.ArgumentParser(description="Compute patristic distance matrices from Newick trees (Biopython).")
    ap.add_argument("--in_dir", required=True, help="Directory with Newick trees.")
    ap.add_argument("--out_dir", required=True, help="Directory to write matrices.")
    ap.add_argument("--glob", default="*.nwk", help="Glob for tree files (e.g., '*.nwk' or '*.treefile').")
    ap.add_argument("--sep", default=",", choices=[",", "tab"], help="CSV (',') or TSV ('tab').")
    ap.add_argument("--suffix", default="_cophenetic.csv", help="Output filename suffix.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    args = ap.parse_args()

    in_dir = os.path.abspath(args.in_dir)
    out_dir = os.path.abspath(args.out_dir)
    pattern = os.path.join(in_dir, args.glob)
    sep = "\t" if args.sep == "tab" else ","

    os.makedirs(out_dir, exist_ok=True)
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[ERROR] No files match: {pattern}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Found {len(files)} tree(s). Writing to: {out_dir}")

    failures = []
    for f in files:
        base = os.path.basename(f)
        species = os.path.splitext(base)[0]
        out_path = os.path.join(out_dir, f"{species}{args.suffix}")
        if os.path.exists(out_path) and not args.overwrite:
            print(f"[SKIP] {base} -> {out_path} (exists; use --overwrite)")
            continue
        try:
            D = compute_patristic_matrix(f)
            D.to_csv(out_path, sep=sep)
            print(f"[OK]   {base} -> {out_path} (n={D.shape[0]} tips)")
        except Exception as e:
            print(f"[FAIL] {base}: {e}", file=sys.stderr)
            failures.append((base, str(e)))

    if failures:
        print("\n[SUMMARY] Some files failed:", file=sys.stderr)
        for base, err in failures:
            print(f"  - {base}: {err}", file=sys.stderr)
        sys.exit(2)

    print("[DONE] All matrices written.")

if __name__ == "__main__":
    main()
