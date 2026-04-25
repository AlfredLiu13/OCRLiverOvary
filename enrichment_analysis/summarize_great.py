#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

from enrichment_analysis.config import GREAT_OUTPUT_DIR, SUMMARY_OUTPUT_DIR


def find_go_file(folder: Path) -> Path | None:
    for f in folder.glob("*.csv"):
        if "GO_Biological_Process" in f.name:
            return f
    return None


def load_and_filter_go(file_path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(file_path)

    # find p-value column automatically
    pval_col = None
    for col in df.columns:
        if "binom" in col.lower():
            pval_col = col
            break

    if pval_col is None:
        print(f"[WARN] No binomial p-value column found in {file_path}")
        return pd.DataFrame()

    df[pval_col] = pd.to_numeric(df[pval_col], errors="coerce")

    df = df[df[pval_col] < 0.05].copy()

    if df.empty:
        return df

    df["dataset"] = label

    # keep only useful columns
    keep_cols = ["dataset", "ID", "Description", pval_col]
    keep_cols = [c for c in keep_cols if c in df.columns]

    return df[keep_cols]


def main():
    SUMMARY_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    for folder in sorted(GREAT_OUTPUT_DIR.iterdir()):
        if not folder.is_dir():
            continue

        label = folder.name
        go_file = find_go_file(folder)

        if go_file is None:
            print(f"[WARN] No GO BP file found for {label}")
            continue

        print(f"[LOAD] {label}")

        df = load_and_filter_go(go_file, label)

        if not df.empty:
            all_results.append(df)

    if not all_results:
        print("[WARN] No enrichment results found.")
        return

    merged = pd.concat(all_results, ignore_index=True)

    output_file = SUMMARY_OUTPUT_DIR / "great_summary.tsv"
    merged.to_csv(output_file, sep="\t", index=False)

    print(f"\n[OK] Summary written to: {output_file}")


if __name__ == "__main__":
    main()