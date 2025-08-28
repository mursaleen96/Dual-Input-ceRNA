# modules/statistical_validation.py

import pandas as pd
import statsmodels.api as sm
import numpy as np
import os
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

def calculate_adjusted_pvalues(pvalues, method='fdr_bh'):
    """
    Calculate adjusted p-values using multiple testing correction
    """
    try:
        # Remove NaN values for correction
        valid_mask = ~np.isnan(pvalues)
        valid_pvalues = pvalues[valid_mask]
        
        if len(valid_pvalues) == 0:
            return pvalues, pvalues  # Return original if all NaN
        
        # Apply correction
        rejected, pvals_corrected, alpha_sidak, alpha_bonf = multipletests(
            valid_pvalues, method=method, alpha=0.05
        )
        
        # Create full array with corrected values
        adjusted_pvalues = np.full_like(pvalues, np.nan)
        adjusted_pvalues[valid_mask] = pvals_corrected
        
        return adjusted_pvalues, rejected[valid_mask] if len(rejected) > 0 else np.array([])
        
    except Exception as e:
        print(f"Error in p-value adjustment: {e}")
        return pvalues, np.array([])

def rank_triplets(validated_df):
    """
    Rank validated triplets based on statistical significance and effect size
    """
    if validated_df.empty:
        return validated_df
    
    print("Ranking validated triplets...")
    
    # Calculate adjusted p-values
    pvalues = validated_df['mediation_pvalue'].values
    adjusted_pvals, rejected = calculate_adjusted_pvalues(pvalues, method='fdr_bh')
    
    validated_df['adjusted_pvalue'] = adjusted_pvals
    validated_df['significant_after_correction'] = False
    
    if len(rejected) > 0:
        valid_indices = ~np.isnan(pvalues)
        validated_df.loc[valid_indices, 'significant_after_correction'] = rejected
    
    # Calculate ranking score (lower is better)
    # Combines p-value significance with effect size
    validated_df['ranking_score'] = (
        -np.log10(validated_df['mediation_pvalue'] + 1e-100) * 0.6 +  # P-value contribution
        -np.log10(validated_df['adjusted_pvalue'] + 1e-100) * 0.3 +   # Adjusted p-value contribution
        np.abs(validated_df['indirect_effect']) * 10 * 0.1            # Effect size contribution
    )
    
    # Sort by ranking score (highest first)
    ranked_df = validated_df.sort_values([
        'ranking_score',           # Primary: combined ranking score
        'mediation_pvalue',        # Secondary: raw p-value
        'adjusted_pvalue'          # Tertiary: adjusted p-value
    ], ascending=[False, True, True]).reset_index(drop=True)
    
    # Add final rank
    ranked_df['final_rank'] = range(1, len(ranked_df) + 1)
    
    # Reorder columns for better readability
    column_order = [
        'final_rank', 'lncRNA', 'miRNA', 'mRNA', 'score',
        'mediation_pvalue', 'adjusted_pvalue', 'significant_after_correction',
        'indirect_effect', 'ranking_score', 'sensitivity', 'mediation_type'
    ]
    
    # Add any additional columns that exist
    additional_cols = [col for col in ranked_df.columns if col not in column_order]
    column_order.extend(additional_cols)
    
    # Filter to existing columns
    existing_cols = [col for col in column_order if col in ranked_df.columns]
    ranked_df = ranked_df[existing_cols]
    
    return ranked_df

def main():
    predictions_path = "results/predicted_triplets.csv"
    norm_counts_path = "results/norm_counts.csv"
    validated_path = "results/validated_triplets.csv"
    ranked_path = "results/validated_triplets_ranked.csv"

    # Load predictions and normalized counts
    predictions = pd.read_csv(predictions_path)
    norm_counts = pd.read_csv(norm_counts_path, index_col=0)

    if predictions.empty:
        print("No predicted triplets. Saving empty validated files.")
        empty_df = pd.DataFrame(columns=predictions.columns.tolist() + 
                               ['mediation_pvalue', 'adjusted_pvalue', 'indirect_effect', 
                                'sensitivity', 'mediation_type', 'ranking_score', 'final_rank'])
        empty_df.to_csv(validated_path, index=False)
        empty_df.to_csv(ranked_path, index=False)
        return

    validated = []
    processed = 0
    failed = 0
    
    print(f"Performing statistical validation on {len(predictions)} predicted triplets...")
    
    for _, row in predictions.iterrows():
        lnc = row['lncRNA']
        mir = row['miRNA']
        mrna = row['mRNA']
        
        processed += 1
        if processed % 100 == 0:
            print(f"Processed {processed}/{len(predictions)} triplets")

        try:
            if lnc in norm_counts.index and mir in norm_counts.index and mrna in norm_counts.index:
                lnc_expr = norm_counts.loc[lnc].values
                mir_expr = norm_counts.loc[mir].values
                mrna_expr = norm_counts.loc[mrna].values

                # Check for sufficient variance
                if any(np.var(expr) < 1e-10 for expr in [lnc_expr, mir_expr, mrna_expr]):
                    failed += 1
                    continue

                # Mediation analysis (using OLS)
                X = sm.add_constant(np.column_stack((lnc_expr, mir_expr)))
                
                # Path a: lncRNA -> miRNA
                med_model = sm.OLS(mir_expr, sm.add_constant(lnc_expr)).fit()
                
                # Path b: miRNA -> mRNA (controlled for lncRNA)
                out_model = sm.OLS(mrna_expr, X).fit()

                # Extract coefficients and standard errors
                a = med_model.params[1]
                sa = med_model.bse[1]
                b = out_model.params[2]
                sb = out_model.bse[2]

                # Mediation effect and p-value (Sobel test)
                mediation_effect = a * b
                se_med = np.sqrt((a**2 * sb**2) + (b**2 * sa**2))
                
                if se_med > 0:
                    z_med = mediation_effect / se_med
                    p_med = 2 * (1 - norm.cdf(abs(z_med)))
                else:
                    p_med = 1.0

                # Determine mediation type
                if p_med < 0.05:
                    # Check if direct effect is still significant
                    c_prime_pval = out_model.pvalues[1]  # lncRNA coefficient p-value
                    if c_prime_pval < 0.05:
                        mediation_type = 'partial_mediation'
                    else:
                        mediation_type = 'complete_mediation'
                else:
                    mediation_type = 'no_mediation'

                # Sensitivity measure
                sensitivity = abs(mediation_effect)

                validated_row = row.to_dict()
                validated_row['mediation_pvalue'] = p_med
                validated_row['indirect_effect'] = mediation_effect
                validated_row['sensitivity'] = sensitivity
                validated_row['mediation_type'] = mediation_type
                
                # Additional path information
                validated_row['path_a_coeff'] = a
                validated_row['path_a_pvalue'] = med_model.pvalues[1]
                validated_row['path_b_coeff'] = b
                validated_row['path_b_pvalue'] = out_model.pvalues[2]
                
                validated.append(validated_row)
                
            else:
                failed += 1
                
        except Exception as e:
            failed += 1
            if failed <= 5:  # Only print first 5 errors
                print(f"Error validating {lnc}-{mir}-{mrna}: {e}")

    # Create validated DataFrame
    validated_df = pd.DataFrame(validated)
    
    # Filter for significant mediation effects
    if not validated_df.empty:
        significant_mask = validated_df['mediation_pvalue'] < 0.05
        validated_df = validated_df[significant_mask]
        
        print(f"Found {len(validated_df)} triplets with significant mediation effects")
        
        # Rank the triplets
        if not validated_df.empty:
            ranked_df = rank_triplets(validated_df)
            
            # Save both versions
            validated_df.to_csv(validated_path, index=False)
            ranked_df.to_csv(ranked_path, index=False)
            
            # Print ranking summary
            print(f"\nRanking Summary:")
            print(f"  Total significant triplets: {len(ranked_df)}")
            if 'significant_after_correction' in ranked_df.columns:
                sig_after_correction = ranked_df['significant_after_correction'].sum()
                print(f"  Significant after FDR correction: {sig_after_correction}")
            
            if len(ranked_df) > 0:
                print(f"  Best p-value: {ranked_df['mediation_pvalue'].min():.2e}")
                print(f"  Median p-value: {ranked_df['mediation_pvalue'].median():.2e}")
                print(f"  Best adjusted p-value: {ranked_df['adjusted_pvalue'].min():.2e}")
                
                print(f"\nTop 5 ranked triplets:")
                top_5 = ranked_df.head(5)
                for _, row in top_5.iterrows():
                    print(f"  Rank {int(row['final_rank'])}: {row['lncRNA']}-{row['miRNA']}-{row['mRNA']} "
                          f"(p={row['mediation_pvalue']:.2e}, adj_p={row['adjusted_pvalue']:.2e})")
        else:
            # Empty dataframes
            pd.DataFrame().to_csv(validated_path, index=False)
            pd.DataFrame().to_csv(ranked_path, index=False)
    else:
        # Empty dataframes
        pd.DataFrame().to_csv(validated_path, index=False)
        pd.DataFrame().to_csv(ranked_path, index=False)
    
    print(f"\nStatistical validation completed:")
    print(f"  Processed: {processed} triplets")
    print(f"  Failed: {failed} triplets") 
    print(f"  Validated and ranked: {len(ranked_df) if not validated_df.empty else 0} triplets")
    print(f"  Results saved to:")
    print(f"    {validated_path}")
    print(f"    {ranked_path}")

if __name__ == "__main__":
    main()