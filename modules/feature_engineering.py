# modules/feature_engineering.py

import pandas as pd
import numpy as np
import pickle
from scipy.stats import pearsonr
from itertools import product
from collections import defaultdict
from sklearn.linear_model import LinearRegression
import yaml
import os  # FIXED: Added missing import for os.path.exists

def compute_partial_correlation(x, y, z):
    try:
        x = x.values.reshape(-1,1) if hasattr(x, 'values') else x.reshape(-1,1)
        y = y.values.reshape(-1,1) if hasattr(y, 'values') else y.reshape(-1,1)
        z = z.values.reshape(-1,1) if hasattr(z, 'values') else z.reshape(-1,1)
        
        model_x = LinearRegression().fit(z, x)
        model_y = LinearRegression().fit(z, y)
        
        x_res = x - model_x.predict(z)
        y_res = y - model_y.predict(z)
        
        r, p = pearsonr(x_res.flatten(), y_res.flatten())
        return r, p
    except Exception as e:
        print(f"    WARNING: Partial correlation failed: {e}")
        return np.nan, np.nan

def sponge_effect(lnc_expr, mirna_expr, mrna_expr):
    try:
        r_lnc_mrna, _ = pearsonr(lnc_expr, mrna_expr)
        
        def residuals(x, z):
            model = LinearRegression().fit(z.reshape(-1,1), x.reshape(-1,1))
            return x - model.predict(z.reshape(-1,1)).flatten()
        
        x_res = residuals(lnc_expr, mirna_expr)
        y_res = residuals(mrna_expr, mirna_expr)
        
        r_partial, _ = pearsonr(x_res, y_res)
        
        sensitivity = r_lnc_mrna - r_partial
        return sensitivity
    except Exception as e:
        print(f"    WARNING: SPONGE computation failed: {e}")
        return np.nan

def load_interaction_db(path):
    try:
        df = pd.read_csv(path, sep='\t', low_memory=False)
        if df.empty:
            print(f"Warning: Empty database file {path}")
            return {}
        
        mirna_to_targets = defaultdict(set)
        for _, row in df.iterrows():
            mirna_to_targets[row['miRNA']].add(row.iloc[1])
        return dict(mirna_to_targets)
    except Exception as e:
        print(f"Error loading {path}: {e}")
        return {}

def main():
    # Load config for thresholds
    cfg_path = "config/config.yaml"
    if os.path.exists(cfg_path):
        with open(cfg_path, "r") as f:
            cfg = yaml.safe_load(f)
    else:
        cfg = {}
        print("WARNING: config.yaml not found; using defaults.")

    corr_cutoff = cfg.get("corr_cutoff", 0.5)
    pval_cutoff = cfg.get("pval_cutoff", 0.01)
    sponge_cutoff = cfg.get("sponge_cutoff", 0.1)

    print("Loading normalized expression data...")
    norm_counts = pd.read_csv("results/norm_counts.csv", index_col=0)

    print("Loading miRNA - mRNA and miRNA - lncRNA interaction data...")
    mirna_mrna = load_interaction_db("databases/miRTarBase.txt")
    mirna_lncrna = load_interaction_db("databases/LncBase.txt")

    # Extract miRNA, lncRNA, mRNA indices
    genes = norm_counts.index.tolist()
    
    mirnas = set(list(mirna_mrna.keys()) + list(mirna_lncrna.keys()))
    mirnas = [m for m in mirnas if m in norm_counts.index]

    mRNAs = set()
    for targets in mirna_mrna.values():
        mRNAs.update(targets)
    mRNAs = [g for g in mRNAs if g in norm_counts.index]

    lncRNAs = set()
    for targets in mirna_lncrna.values():
        lncRNAs.update(targets)
    lncRNAs = [g for g in lncRNAs if g in norm_counts.index]

    print(f"Number of miRNAs: {len(mirnas)}")
    print(f"Number of mRNAs: {len(mRNAs)}")
    print(f"Number of lncRNAs: {len(lncRNAs)}")

    # Build triplets
    triplets = []
    for miRNA in mirnas:
        targets_mrna = mirna_mrna.get(miRNA, set())
        targets_lncrna = mirna_lncrna.get(miRNA, set())
        
        targets_mrna = set(t for t in targets_mrna if t in mRNAs)
        targets_lncrna = set(t for t in targets_lncrna if t in lncRNAs)
        
        for lnc in targets_lncrna:
            for mrna in targets_mrna:
                triplets.append((lnc, miRNA, mrna))

    print(f"Total candidate triplets: {len(triplets)}")

    if len(triplets) == 0:
        print("WARNING: No triplets found! Check database files and gene name consistency.")
        empty_features = pd.DataFrame(columns=[
            "lncRNA", "miRNA", "mRNA", "pearson_lnc_mrna", "pval_lnc_mrna",
            "pearson_lnc_mirna", "pval_lnc_mirna", "pearson_mrna_mirna", "pval_mrna_mirna",
            "partial_corr_lnc_mrna_mirna", "partial_corr_pval", "sponge_score",
            "mre_counts", "seed_match_energy", "cytoplasmic_localization"
        ])
        empty_features.to_pickle("results/features.pkl")
        print("Created empty features file.")
        return

    # Process triplets with filtering
    feature_rows = []
    skipped = 0
    for (lnc, miRNA, mrna) in triplets:
        try:
            lnc_expr = norm_counts.loc[lnc]
            miRNA_expr = norm_counts.loc[miRNA]
            mrna_expr = norm_counts.loc[mrna]

            common_samples = sorted(set(lnc_expr.index) & set(miRNA_expr.index) & set(mrna_expr.index))
            if not common_samples:
                skipped += 1
                continue

            lnc_vals = lnc_expr[common_samples].values
            miRNA_vals = miRNA_expr[common_samples].values
            mrna_vals = mrna_expr[common_samples].values

            if not (len(lnc_vals) == len(miRNA_vals) == len(mrna_vals)):
                skipped += 1
                continue

            # Pairwise correlations with filtering
            r_lncmrna, p_lncmrna = pearsonr(lnc_vals, mrna_vals)
            if abs(r_lncmrna) < corr_cutoff or p_lncmrna > pval_cutoff:
                skipped += 1
                continue

            r_lncmirna, p_lncmirna = pearsonr(lnc_vals, miRNA_vals)
            if abs(r_lncmirna) < corr_cutoff or p_lncmirna > pval_cutoff:
                skipped += 1
                continue

            r_mrnamirna, p_mrnamirna = pearsonr(mrna_vals, miRNA_vals)
            if abs(r_mrnamirna) < corr_cutoff or p_mrnamirna > pval_cutoff:
                skipped += 1
                continue

            # Partial correlation
            r_partial, p_partial = compute_partial_correlation(lnc_vals, mrna_vals, miRNA_vals)
            if p_partial > pval_cutoff:
                skipped += 1
                continue

            # SPONGE score with filter
            sponge_score = sponge_effect(lnc_vals, miRNA_vals, mrna_vals)
            if sponge_score < sponge_cutoff:
                skipped += 1
                continue

            # Placeholders for sequence features
            mre_counts = np.nan
            seed_match_energy = np.nan
            cytoplasmic_localization = np.nan

            feature_rows.append({
                "lncRNA": lnc,
                "miRNA": miRNA,
                "mRNA": mrna,
                "pearson_lnc_mrna": r_lncmrna,
                "pval_lnc_mrna": p_lncmrna,
                "pearson_lnc_mirna": r_lncmirna,
                "pval_lnc_mirna": p_lncmirna,
                "pearson_mrna_mirna": r_mrnamirna,
                "pval_mrna_mirna": p_mrnamirna,
                "partial_corr_lnc_mrna_mirna": r_partial,
                "partial_corr_pval": p_partial,
                "sponge_score": sponge_score,
                "mre_counts": mre_counts,
                "seed_match_energy": seed_match_energy,
                "cytoplasmic_localization": cytoplasmic_localization
            })
            
        except KeyError as e:
            print(f"Gene not found in expression data: {e}")
            skipped += 1
            continue
        except Exception as e:
            print(f"Error processing triplet {lnc}-{miRNA}-{mrna}: {e}")
            skipped += 1
            continue

    print(f"Processed {len(triplets)} triplets; skipped {skipped} due to filters or errors")

    features_df = pd.DataFrame(feature_rows)
    features_df.to_pickle("results/features.pkl")
    print(f"Feature engineering completed and saved to results/features.pkl")
    print(f"Generated features for {len(features_df)} triplets")

if __name__ == "__main__":
    main()
