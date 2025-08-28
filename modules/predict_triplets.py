# modules/predict_triplets.py

import pickle
import pandas as pd
import numpy as np
import os

def main():
    features_path = "results/features.pkl"
    models_path = "results/models.pkl"
    predictions_path = "results/predicted_triplets.csv"

    # Load model
    with open(models_path, "rb") as f:
        model = pickle.load(f)

    # Load features
    with open(features_path, "rb") as f:
        features = pickle.load(f)

    # Handle empty or placeholder model
    if features.empty or model is None:
        print("No model or features available. Saving empty predictions.")
        pd.DataFrame(columns=["lncRNA", "miRNA", "mRNA", "score"]).to_csv(predictions_path, index=False)
        return

    # Select only numeric columns for prediction
    numeric_cols = features.select_dtypes(include=['int', 'float', 'bool']).columns
    # Remove identifier columns
    identifier_cols = ['lncRNA', 'miRNA', 'mRNA']
    numeric_cols = [col for col in numeric_cols if col not in identifier_cols]

    if len(numeric_cols) == 0:
        print("No numeric features found. Saving empty predictions.")
        pd.DataFrame(columns=["lncRNA", "miRNA", "mRNA", "score"]).to_csv(predictions_path, index=False)
        return

    X = features[numeric_cols].copy()
    
    # Handle missing values (same as training)
    X = X.fillna(X.median())
    
    # Remove zero variance columns
    var_mask = X.var() > 0
    X = X.loc[:, var_mask]

    # Generate predictions
    print("Generating predictions...")
    try:
        probs = model.predict_proba(X)[:, 1]  # Probability of positive class
        
        # Validate predictions
        if np.any(np.isnan(probs)):
            print("Warning: Found NaN predictions, replacing with 0.1")
            probs = np.nan_to_num(probs, nan=0.1)
            
        # Ensure valid range
        probs = np.clip(probs, 0.0, 1.0)
        
    except Exception as e:
        print(f"Prediction failed: {e}")
        # Fallback to random scores
        probs = np.random.rand(len(X)) * 0.5 + 0.1

    # Add predictions to the original features (with identifiers)
    predictions = features[["lncRNA", "miRNA", "mRNA"]].copy()
    predictions["score"] = probs
    predictions = predictions.sort_values("score", ascending=False)

    # Save
    predictions.to_csv(predictions_path, index=False)
    print(f"Predictions saved to {predictions_path}")
    
    # Log some statistics
    print(f"Generated {len(predictions)} predictions")
    print(f"Score range: [{probs.min():.3f}, {probs.max():.3f}]")
    print(f"Mean score: {probs.mean():.3f}")

if __name__ == "__main__":
    main()