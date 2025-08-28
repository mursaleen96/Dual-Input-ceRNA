# modules/qc_normalization.py

import os
import sys
import yaml
import pandas as pd
import numpy as np

# =========================
# Helper Functions
# =========================
def detect_id_type(sample_ids):
    """Detect if IDs are Ensembl, Entrez, or gene symbols."""
    ids = [str(x) for x in sample_ids]
    if any(x.startswith("ENSG") for x in ids):
        return "ensembl"
    # Heuristic for Entrez: mostly numeric IDs of length 5-10
    numeric_like = sum(x.replace(".", "").isdigit() and 5 <= len(x) <= 10 for x in ids)
    if numeric_like >= max(1, len(ids) // 2):
        return "entrez"
    return "symbol"

def clean_ensembl_ids(gene_ids_list):
    """Remove version suffix from Ensembl IDs (e.g., ENSG00000000005.6 -> ENSG00000000005)."""
    cleaned_list = [str(gid).split(".")[0] if str(gid).startswith("ENSG") and "." in str(gid) else str(gid) for gid in gene_ids_list]
    if len(cleaned_list) != len(gene_ids_list):
        raise ValueError("Length mismatch after cleaning Ensembl IDs")
    return cleaned_list

def load_gene_annotation(path="databases/gene_annotation.csv"):
    """Load gene annotation CSV and validate required columns."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Gene annotation file not found: {path}")
    df = pd.read_csv(path)
    required = {"ensembl_gene_id", "gene_symbol", "gene_biotype"}
    if not required.issubset(df.columns):
        raise ValueError(f"gene_annotation.csv must contain at least: {required}")
    # Clean Ensembl IDs in annotation (strip versions)
    df["ensembl_gene_id"] = [x.split(".")[0] if "." in str(x) else str(x) for x in df["ensembl_gene_id"]]
    return df

# =========================
# Expression Processing
# =========================
def cpm_normalization(df):
    """CPM normalization followed by log2 transformation."""
    counts = df.values
    library_sizes = counts.sum(axis=0)
    library_sizes[library_sizes == 0] = 1  # Avoid division by zero
    cpm = (counts / library_sizes) * 1e6
    return pd.DataFrame(np.log2(cpm + 1), index=df.index, columns=df.columns)

def filter_low_expression(df, min_count=1, min_samples=2):
    """Filter features with > min_count in >= min_samples."""
    keep = (df > min_count).sum(axis=1) >= min_samples
    return df.loc[keep]

def load_and_process_counts(path, label, annot=None):
    print(f"Loading {label} from {path}")
    df = pd.read_csv(path, index_col=0)
    print(f"Loaded: {df.shape[0]} features × {df.shape[1]} samples")

    if label == "RNA-seq":
        # Detect ID type
        head_ids = df.index.tolist()[:10]
        id_type = detect_id_type(head_ids)
        print(f"Detected ID type for RNA-seq: {id_type}")

        if id_type == "entrez":
            raise ValueError("Entrez gene IDs are not supported yet in this pipeline. Please convert to Ensembl or gene symbols.")

        if id_type == "ensembl" and annot is not None:
            # Clean Ensembl IDs in RNA-seq index
            df.index = clean_ensembl_ids(df.index.tolist())
            # Map to gene symbols
            common = set(df.index) & set(annot["ensembl_gene_id"])
            if not common:
                print("WARNING: No common IDs between RNA-seq and annotation! Proceeding without mapping.")
            else:
                df = df.loc[list(common)]
                symbol_map = dict(zip(annot["ensembl_gene_id"], annot["gene_symbol"]))
                df.index = [symbol_map.get(x, x) for x in df.index]
                df = df.groupby(df.index).sum()
                print(f"After Ensembl → gene symbol mapping: {df.shape[0]} genes")

        elif id_type == "symbol":
            print("Detected gene symbols; no mapping needed.")

    elif label == "miRNA-seq":
        # No changes to miRNA names, as per user request
        print("miRNA names already standardized; no changes applied.")

    return df

# =========================
# Sample Alignment
# =========================
def align_samples(rna_df, mirna_df):
    """Align to common samples between RNA-seq and miRNA-seq."""
    common = sorted(set(rna_df.columns) & set(mirna_df.columns))
    print(f"Common samples: {len(common)}")
    if len(common) == 0:
        raise ValueError("No common samples! Check column names in input CSVs.")
    return rna_df[common], mirna_df[common]

# =========================
# Main
# =========================
def main():
    print("=" * 50)
    print("STARTING QC AND NORMALIZATION")
    print("=" * 50)

    cfg_path = "config/config.yaml"
    cfg = yaml.safe_load(open(cfg_path)) if os.path.exists(cfg_path) else {}

    annot = load_gene_annotation(cfg.get("gene_anno_path", "databases/gene_annotation.csv"))

    rna_counts = load_and_process_counts("data/rna_counts.csv", "RNA-seq", annot)
    mirna_counts = load_and_process_counts("data/mirna_counts.csv", "miRNA-seq")

    rna_aligned, mirna_aligned = align_samples(rna_counts, mirna_counts)

    combined = pd.concat([rna_aligned, mirna_aligned], axis=0)
    print(f"Combined: {combined.shape}")

    # Filter low-expression (configurable)
    min_count = cfg.get("min_count", 1)
    min_samples = cfg.get("min_samples", 2)
    filtered = filter_low_expression(combined, min_count, min_samples)
    n_filtered = combined.shape[0] - filtered.shape[0]
    print(f"Filtered out {n_filtered} features (> {min_count} in >= {min_samples} samples)")

    norm = cpm_normalization(filtered)
    os.makedirs("results", exist_ok=True)
    norm.to_csv("results/norm_counts.csv")

    metadata = pd.DataFrame(index=filtered.columns)  # Use columns (samples) as index
    metadata.to_csv("results/sample_metadata.csv")

    print("=" * 50)
    print("COMPLETED")
    print(f"Final dataset: {norm.shape}")
    print("=" * 50)

if __name__ == "__main__":
    main()