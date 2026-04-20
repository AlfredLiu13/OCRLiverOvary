#!/usr/bin/env python3

"""
Analysis of rGREAT enrichment analysis results.

Generates:
    - summary.csv – Overview statistics for each GREAT file
    - *_filtered.tsv – All significant terms for each file
    - *_top15.tsv – Top 15 terms by fold enrichment
    - *_top_terms.png – Bar plot of top enriched terms
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

#Set visualization style for better-looking plots
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

#Results directory -> Contains all output TSV files from rGREAT
results_dir = Path("results/enrichment_results")
output_dir = Path("results/enrichment_results/results_summaries")
output_dir.mkdir(parents=True, exist_ok=True)

#Filtering thresholds
min_pval = 0.05
min_fold_enrichment = 2.0
min_observed_genes = 5
top_n = 15

print("GREAT Enrichment Analysis Results Summary")

#Load all TSV files (recursive search through all subfolders)
tsv_files = list(results_dir.glob("**/*.tsv"))
if not tsv_files:
    print(f"ERROR: No TSV files found in {results_dir}")
    exit(1)

results = {}
for tsv_file in tsv_files:
    try:
        df = pd.read_csv(tsv_file, sep="\t", low_memory=False)
        results[tsv_file.name] = df
        print(f"Loaded {tsv_file.name} ({len(df)} rows)")
    except Exception as e:
        print(f"Failed to load {tsv_file.name}: {e}")


def normalize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize column names for easier access.
    rGREAT uses various naming conventions; this normalizes them.
    """
    df_copy = df.copy()
    
    #Create a mapping of potential column names
    name_mapping = {
        'Binomial p-value': 'pval',
        'p-value': 'pval',
        'Binom_Raw_PValue': 'pval',
        'Hyper_Raw_PValue': 'pval',
        'FDR': 'fdr',
        'Adjusted p-value': 'fdr',
        'Binom_Adjp_BH': 'fdr',
        'Hyper_Adjp_BH': 'fdr',
        'Fold enrichment': 'fold_enrichment',
        'Fold Change': 'fold_enrichment',
        'Binom_Fold_Enrichment': 'fold_enrichment',
        'Hyper_Fold_Enrichment': 'fold_enrichment',
        'Observed genes': 'observed_genes',
        'Observed Genes': 'observed_genes',
        '# of Observed Genes': 'observed_genes',
        'Binom_Observed_Region_Hits': 'observed_genes',
        'Hyper_Observed_Gene_Hits': 'observed_genes',
        'Term': 'term',
        'Term name': 'term',
        'name': 'term',
        'Ontology ID': 'term_id',
        'Term ID': 'term_id',
        'ID': 'term_id',
    }
    
    #Apply mapping where applicable
    for old_name, new_name in name_mapping.items():
        if old_name in df_copy.columns and new_name not in df_copy.columns:
            df_copy.rename(columns={old_name: new_name}, inplace=True)
    
    return df_copy


def filter_results(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter enrichment table by significance and effect size.
    """
    df = normalize_column_names(df)
    
    #Determine which p-value column to use
    pval_col = 'pval' if 'pval' in df.columns else 'fdr'
    if pval_col not in df.columns:
        pval_cols = [c for c in df.columns if 'p-value' in c.lower() or 'pval' in c.lower()]
        pval_col = pval_cols[0] if pval_cols else df.select_dtypes(include=[np.number]).columns[0]
    
    fold_col = 'fold_enrichment' if 'fold_enrichment' in df.columns else None
    observed_col = 'observed_genes' if 'observed_genes' in df.columns else None
    
    initial_rows = len(df)
    
    #Filter by p-value
    df = df[pd.to_numeric(df[pval_col], errors='coerce') < min_pval]
    
    #Filter by fold enrichment (if available)
    if fold_col and fold_col in df.columns:
        df = df[pd.to_numeric(df[fold_col], errors='coerce') >= min_fold_enrichment]
    
    #Filter by minimum observed genes (if available)
    if observed_col and observed_col in df.columns:
        df = df[pd.to_numeric(df[observed_col], errors='coerce') >= min_observed_genes]
    
    filtered_rows = len(df)
    print(f"  Filtered {initial_rows} -> {filtered_rows} rows (p<{min_pval}, fold>{min_fold_enrichment}, genes>{min_observed_genes})")
    
    return df


def get_top_terms(df: pd.DataFrame, n: int = 15, sort_by: str = "fold_enrichment") -> pd.DataFrame:
    """Get top N enriched terms."""
    df = normalize_column_names(df)
    
    if sort_by == "fold_enrichment":
        col = 'fold_enrichment'
        ascending = False
    else:
        col = 'pval' if 'pval' in df.columns else 'fdr'
        ascending = True
    
    if col not in df.columns:
        col = df.select_dtypes(include=[np.number]).columns[0]
    
    df[col] = pd.to_numeric(df[col], errors='coerce')
    return df.nlargest(n, col) if not ascending else df.nsmallest(n, col)


def plot_top_terms(df: pd.DataFrame, filename: str, n: int = 15):
    """Generate a bar plot of top enriched terms."""
    df = normalize_column_names(df)
    top_terms = get_top_terms(df, n=n, sort_by="fold_enrichment")
    
    if len(top_terms) == 0:
        print(f"  Warning: No top terms to plot for {filename}")
        return
    
    #Prepare data
    term_col = 'term' if 'term' in top_terms.columns else top_terms.columns[0]
    fold_col = 'fold_enrichment' if 'fold_enrichment' in top_terms.columns else None
    
    if not fold_col:
        print(f"  Warning: No fold enrichment column in {filename}. Available: {list(top_terms.columns)}")
        return
    
    top_terms[fold_col] = pd.to_numeric(top_terms[fold_col], errors='coerce')
    plot_data = top_terms[[term_col, fold_col]].dropna().sort_values(fold_col)
    
    #Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(range(len(plot_data)), plot_data[fold_col], color='steelblue', edgecolor='black', alpha=0.7)
    ax.set_yticks(range(len(plot_data)))
    ax.set_yticklabels([str(t)[:50] for t in plot_data[term_col]], fontsize=9)
    ax.set_xlabel("Fold Enrichment", fontsize=11, fontweight='bold')
    ax.set_title(f"Top {n} Enriched Terms: {filename}", fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    output_path = output_dir / filename.replace(".tsv", "_top_terms.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved top_terms plot to {output_path}")


#GENERATE SUMMARY STATISTICS
print("\nSUMMARY STATISTICS")

summaries = []
for filename, df in results.items():
    df_filtered = filter_results(df)
    
    if len(df_filtered) == 0:
        summaries.append({
            "File": filename,
            "Total Terms": len(df),
            "Significant Terms": 0,
            "Mean Fold Enrichment": np.nan,
            "Median Fold Enrichment": np.nan,
            "Max Fold Enrichment": np.nan,
        })
    else:
        df_filtered = normalize_column_names(df_filtered)
        fold_col = 'fold_enrichment' if 'fold_enrichment' in df_filtered.columns else None
        
        fold_values = pd.Series([])
        if fold_col:
            fold_values = pd.to_numeric(df_filtered[fold_col], errors='coerce').dropna()
        
        summaries.append({
            "File": filename,
            "Total Terms": len(df),
            "Significant Terms": len(df_filtered),
            "Mean Fold Enrichment": fold_values.mean() if fold_values.size > 0 else np.nan,
            "Median Fold Enrichment": fold_values.median() if fold_values.size > 0 else np.nan,
            "Max Fold Enrichment": fold_values.max() if fold_values.size > 0 else np.nan,
        })

summary_df = pd.DataFrame(summaries)
summary_df.to_csv(output_dir / "summary.csv", index=False)
print(summary_df.to_string(index=False))


#PROCESS EACH FILE
print("\nPROCESSING RESULTS")

for filename, df in results.items():
    print(f"{filename}:")

    #Filter results
    df_filtered = filter_results(df)
    
    if len(df_filtered) == 0:
        print(f"No significant terms after filtering")
        continue
    
    #Get top terms
    top_terms = get_top_terms(df_filtered, n=top_n, sort_by="fold_enrichment")
    
    #Save filtered and top results
    output_basename = filename.replace(".tsv", "")
    filtered_output = output_dir / f"{output_basename}_filtered.tsv"
    top_output = output_dir / f"{output_basename}_top{top_n}.tsv"
    
    df_filtered.to_csv(filtered_output, sep="\t", index=False)
    top_terms.to_csv(top_output, sep="\t", index=False)
    
    print(f"Saved filtered results to {filtered_output}")
    print(f"Saved top {top_n} terms to {top_output}")
    
    #Generate plots
    plot_top_terms(df, filename, n=top_n)
    print(f"Saved plots: *_top_terms.png")


#SUMMARY
print("\nANALYSIS COMPLETE")
print("\nGenerated outputs:")
print("  1. summary.csv")
print("  2. *_filtered.tsv (all significant terms)")
print("  3. *_top15.tsv (top enriched terms)")
print("  4. *_top_terms.png (bar plots)")