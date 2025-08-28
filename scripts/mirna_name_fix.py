# scripts/mirna_name_fix.py

import pandas as pd

def standardize_miRNA_counts(in_path, out_path):
    """
    Standardize miRNA names in count matrix
    
    Converts formats like:
    - MIR3935 -> hsa-miR-3935
    - mir-21 -> hsa-miR-21
    - MIR21A -> hsa-miR-21a
    """
    df = pd.read_csv(in_path, index_col=0)
    print(f"Loaded {df.shape[0]} genes across {df.shape[1]} samples")

    # Standardize miRNA rows
    def standardize(name):
        name_upper = str(name).upper()
        
        # Check if already standardized
        if name.startswith('hsa-miR-') or name.startswith('hsa-mir-'):
            return name
        
        # Handle various miRNA naming formats
        if name_upper.startswith("MIR"):
            # Remove MIR prefix and add standard prefix
            stem = name[3:].lower()  # Remove 'MIR' and lowercase
            return f"hsa-miR-{stem}"
        elif name_upper.startswith("MIRN"):
            # Handle MIRN prefix
            stem = name[4:].lower()
            return f"hsa-miR-{stem}"
        elif name.lower().startswith("mir-"):
            # Handle mir- prefix
            stem = name[4:]  # Remove 'mir-'
            return f"hsa-miR-{stem}"
        elif name.lower().startswith("mir"):
            # Handle mir prefix (without dash)
            stem = name[3:].lower()
            return f"hsa-miR-{stem}"
        else:
            # Leave non-miRNA rows unchanged
            return name

    # Apply standardization
    original_index = df.index.tolist()
    standardized_index = [standardize(idx) for idx in original_index]
    
    # Count changes
    changes = sum(1 for orig, std in zip(original_index, standardized_index) if orig != std)
    print(f"Standardized {changes} miRNA names")

    # Update index
    df.index = standardized_index

    # Drop duplicates if any after mapping
    initial_count = len(df)
    df = df[~df.index.duplicated(keep='first')]
    final_count = len(df)
    
    if initial_count != final_count:
        print(f"Removed {initial_count - final_count} duplicate entries after standardization")

    # Save result
    df.to_csv(out_path)
    print(f"Standardized data saved to {out_path}")
    print(f"Final dataset: {df.shape[0]} genes across {df.shape[1]} samples")

def main():
    """
    Main function - modify these paths for your data
    """
    # MODIFY THESE PATHS FOR YOUR DATA
    input_file = "your_input_file.csv"  # Change this
    output_file = "your_standardized_file.csv"  # Change this
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found!")
        print("Please modify the input_file and output_file paths in this script.")
        return
    
    try:
        standardize_miRNA_counts(input_file, output_file)
        print("miRNA name standardization completed successfully!")
    except Exception as e:
        print(f"Error during standardization: {e}")

if __name__ == "__main__":
    import os
    main()