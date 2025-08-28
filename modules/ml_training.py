# modules/ml_training.py

import pickle
import pandas as pd
import numpy as np
import os
from xgboost import XGBClassifier

def main():
    features_path = "results/features.pkl"
    models_path = "results/models.pkl"

    # Load features
    with open(features_path, "rb") as f:
        features = pickle.load(f)

    print(f"Loaded {features.shape[0]} feature rows")

    # Handle empty features
    if features.empty:
        print("No features available. Saving placeholder model.")
        with open(models_path, "wb") as f:
            pickle.dump(None, f)
        return

    # Prepare labels (example; adjust threshold/column as needed)
    if 'sponge_score' not in features.columns:
        print("Warning: sponge_score column not found. Using random labels for demonstration.")
        y = np.random.randint(0, 2, size=len(features))
    else:
        # Use top 30% of SPONGE scores as positive class
        threshold = features['sponge_score'].quantile(0.7)
        y = (features['sponge_score'] > threshold).astype(int)

    print(f"Label distribution: {np.bincount(y)}")

    # Select only numeric columns for features
    numeric_cols = features.select_dtypes(include=['int', 'float', 'bool']).columns
    # Remove identifier columns
    identifier_cols = ['lncRNA', 'miRNA', 'mRNA']
    numeric_cols = [col for col in numeric_cols if col not in identifier_cols]

    if len(numeric_cols) == 0:
        print("No numeric features found. Saving placeholder model.")
        with open(models_path, "wb") as f:
            pickle.dump(None, f)
        return

    X = features[numeric_cols].copy()
    
    # Handle missing values
    X = X.fillna(X.median())
    
    # Remove columns with zero variance
    var_mask = X.var() > 0
    X = X.loc[:, var_mask]
    
    if X.shape[1] == 0:
        print("No features with variance. Saving placeholder model.")
        with open(models_path, "wb") as f:
            pickle.dump(None, f)
        return

    print(f"Training with {X.shape[1]} features")

    # Train model with basic parameters
    model = XGBClassifier(
        n_estimators=100,
        max_depth=6,
        learning_rate=0.1,
        random_state=42,
        eval_metric='logloss'
    )
    
    try:
        model.fit(X, y)
        print("Model training completed")
        
        # Basic feature importance
        feature_importance = pd.DataFrame({
            'feature': X.columns,
            'importance': model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        print("Top 5 most important features:")
        print(feature_importance.head())
        
    except Exception as e:
        print(f"Model training failed: {e}")
        model = None

    # Save model
    with open(models_path, "wb") as f:
        pickle.dump(model, f)

    print(f"Model saved to {models_path}")

if __name__ == "__main__":
    main()