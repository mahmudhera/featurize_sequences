#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import pandas as pd


NONFEATURE_COLS = {
    "index", "name", "name.1", "X", "V1",
    "name_meta1", "name_meta2", "name_meta3",
}


def drop_nonfeatures(df: pd.DataFrame) -> pd.DataFrame:
    """Drop obvious metadata columns if present; keep sequence_name."""
    cols = [c for c in df.columns if c not in NONFEATURE_COLS]
    # also drop name_metaN pattern
    cols = [c for c in cols if not c.startswith("name_meta")]
    df = df.loc[:, cols]
    return df


def ensure_sequence_name(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure the ID column is named sequence_name."""
    if "sequence_name" in df.columns:
        return df
    if "name" in df.columns:
        # If there are multiple 'name' columns, pandas may auto-rename to name, name.1, ...
        # Use the first 'name' as sequence_name.
        df = df.rename(columns={"name": "sequence_name"})
        return df
    # Fallback: use first column as ID
    df = df.rename(columns={df.columns[0]: "sequence_name"})
    return df


def to_numeric_block(df: pd.DataFrame, cols: list[str], strict: bool) -> pd.DataFrame:
    """Convert selected columns to numeric; either fail or fill NaN with 0."""
    out = df.loc[:, cols].copy()
    for c in cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    if strict:
        # if any NaN introduced, fail with diagnostic
        bad = out.columns[out.isna().any()].tolist()
        if bad:
            raise ValueError(
                f"Non-numeric values detected after coercion in {len(bad)} columns; "
                f"example columns: {bad[:20]}"
            )
    # default behavior: set NaN to 0
    out = out.fillna(0.0)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Compute feature deltas: alt summary minus ref summary for one cell type."
    )
    ap.add_argument("--cell-dir", required=True,
                    help="Path to <CELL>/features directory containing ref/ and alt/")
    ap.add_argument("--out", default=None,
                    help="Output CSV path (default: <cell-dir>/delta_alt_minus_ref.csv)")
    ap.add_argument("--strict", action="store_true",
                    help="Fail if any feature column contains non-numeric values after coercion.")
    ap.add_argument("--sort", action="store_true",
                    help="Sort by sequence_name before comparing (recommended).")
    args = ap.parse_args()

    cell_dir = Path(args.cell_dir).resolve()
    if not cell_dir.exists():
        print(f"ERROR: --cell-dir does not exist: {cell_dir}", file=sys.stderr)
        return 2

    ref_path = cell_dir / "ref" / "summary_ref.csv"
    alt_path = cell_dir / "alt" / "summary_alt.csv"

    if not ref_path.exists():
        print(f"ERROR: missing: {ref_path}", file=sys.stderr)
        return 2
    if not alt_path.exists():
        print(f"ERROR: missing: {alt_path}", file=sys.stderr)
        return 2

    out_path = Path(args.out).resolve() if args.out else (cell_dir / "delta_alt_minus_ref.csv")

    # Read
    ref = pd.read_csv(ref_path)
    alt = pd.read_csv(alt_path)

    # Normalize ID + remove obvious junk columns if any
    ref = ensure_sequence_name(ref)
    alt = ensure_sequence_name(alt)
    ref = drop_nonfeatures(ref)
    alt = drop_nonfeatures(alt)

    if args.sort:
        ref = ref.sort_values("sequence_name").reset_index(drop=True)
        alt = alt.sort_values("sequence_name").reset_index(drop=True)

    if "sequence_name" not in ref.columns or "sequence_name" not in alt.columns:
        raise ValueError("sequence_name column missing after normalization.")

    if not ref["sequence_name"].equals(alt["sequence_name"]):
        # Provide a helpful mismatch report
        ref_set = set(ref["sequence_name"])
        alt_set = set(alt["sequence_name"])
        only_ref = list(ref_set - alt_set)[:10]
        only_alt = list(alt_set - ref_set)[:10]
        raise ValueError(
            "sequence_name mismatch between ref and alt summaries.\n"
            f"Examples only in ref: {only_ref}\n"
            f"Examples only in alt: {only_alt}\n"
            "Tip: ensure both summaries were generated from the same locs.fasta ordering."
        )

    common = sorted(set(ref.columns).intersection(set(alt.columns)))
    common = [c for c in common if c != "sequence_name"]

    if not common:
        raise ValueError("No common feature columns found between ref and alt summaries.")

    ref_num = to_numeric_block(ref, common, strict=args.strict)
    alt_num = to_numeric_block(alt, common, strict=args.strict)

    delta = alt_num.values - ref_num.values
    out = pd.DataFrame(delta, columns=common)
    out.insert(0, "sequence_name", ref["sequence_name"].values)

    out.to_csv(out_path, index=False)
    print(f"Wrote: {out_path}")
    print(f"n_rows={out.shape[0]}  n_features={out.shape[1]-1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

## USAGE:
## python /home/jl2791/scripts/featurizer/feature_diff.py \
#  --cell-dir #/projectsp/f_ak1833_1/data/Empirical_MPRA/processed/empirical_mpra/emvar_exports_encode/SKNSH/features \
#  --sort
