#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import pandas as pd


DROP_FIXED = {
    "sequence_name", "name", "name.1", "index", "X", "V1",
}
DROP_PREFIXES = ("name_meta",)


def drop_nonfeatures(df: pd.DataFrame) -> pd.DataFrame:
    keep = []
    for c in df.columns:
        if c in DROP_FIXED:
            continue
        if any(c.startswith(p) for p in DROP_PREFIXES):
            continue
        keep.append(c)
    return df.loc[:, keep]


def to_numeric(df: pd.DataFrame, strict: bool) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    if strict:
        bad = out.columns[out.isna().any()].tolist()
        if bad:
            raise ValueError(
                f"Non-numeric values detected after coercion in {len(bad)} columns; "
                f"example columns: {bad[:20]}"
            )
    return out.fillna(0.0)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Compute feature deltas (alt - ref) WITHOUT any key: subtract row-by-row by position."
    )
    ap.add_argument("--cell-dir", required=True,
                    help="Path to <CELL>/features directory containing ref/ and alt/")
    ap.add_argument("--out", default=None,
                    help="Output CSV path (default: <cell-dir>/delta_alt_minus_ref.csv)")
    ap.add_argument("--strict-rows", action="store_true",
                    help="Fail if ref and alt have different number of rows. Default: truncate to min rows.")
    ap.add_argument("--strict", action="store_true",
                    help="Fail if any feature column contains non-numeric values after coercion.")
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

    ref = pd.read_csv(ref_path)
    alt = pd.read_csv(alt_path)

    # Drop ID/metadata columns (including sequence_name)
    ref = drop_nonfeatures(ref)
    alt = drop_nonfeatures(alt)

    # Use only common columns in both
    common = sorted(set(ref.columns).intersection(set(alt.columns)))
    if not common:
        raise ValueError("No common feature columns found between ref and alt summaries (after dropping IDs).")

    ref = ref.loc[:, common]
    alt = alt.loc[:, common]

    # Row count handling
    n_ref, n_alt = ref.shape[0], alt.shape[0]
    if n_ref != n_alt:
        msg = f"Row count differs: ref={n_ref}, alt={n_alt}."
        if args.strict_rows:
            raise ValueError(msg + " Use consistent inputs or drop --strict-rows.")
        n = min(n_ref, n_alt)
        print(msg + f" Truncating both to first n={n} rows.", file=sys.stderr)
        ref = ref.iloc[:n, :].reset_index(drop=True)
        alt = alt.iloc[:n, :].reset_index(drop=True)

    # Numeric conversion + subtraction
    ref_num = to_numeric(ref, strict=args.strict)
    alt_num = to_numeric(alt, strict=args.strict)

    delta = alt_num.values - ref_num.values
    out = pd.DataFrame(delta, columns=common)
    out.insert(0, "row_index", range(out.shape[0]))

    out.to_csv(out_path, index=False)
    print(f"Wrote: {out_path}")
    print(f"n_rows={out.shape[0]}  n_features={out.shape[1]-1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
