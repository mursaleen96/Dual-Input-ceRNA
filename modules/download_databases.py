# modules/download_databases.py

import os
import pandas as pd
import shutil

def process_mirtarbase(src, dest):
    """Process miRTarBase txt file ('miRNA', 'Target'), rename and export."""
    print(f"Processing miRTarBase from {src} ...")
    
    if not os.path.exists(src):
        print(f"Warning: {src} not found. Please download miRTarBase data.")
        # Create empty file as placeholder
        pd.DataFrame(columns=['miRNA', 'mRNA']).to_csv(dest, sep='\t', index=False)
        return
    
    try:
        df = pd.read_csv(src, sep='\t', low_memory=False)
        df = df[['miRNA', 'Target']].drop_duplicates()
        df = df.rename(columns={'Target': 'mRNA'})
        df.to_csv(dest, sep='\t', index=False)
        print(f"Processed miRTarBase to {dest}")
    except Exception as e:
        print(f"Error processing miRTarBase: {e}")
        # Create empty file as placeholder
        pd.DataFrame(columns=['miRNA', 'mRNA']).to_csv(dest, sep='\t', index=False)

def process_starbase(src, dest):
    """Process starBase txt file ('miRNA', 'lncRNA'), export as is."""
    print(f"Processing starBase miRNA–lncRNA from {src} ...")
    
    if not os.path.exists(src):
        print(f"Warning: {src} not found. Please download starBase data.")
        # Create empty file as placeholder
        pd.DataFrame(columns=['miRNA', 'lncRNA']).to_csv(dest, sep='\t', index=False)
        return
    
    try:
        df = pd.read_csv(src, sep='\t', low_memory=False)
        df = df[['miRNA', 'lncRNA']].drop_duplicates()
        df.to_csv(dest, sep='\t', index=False)
        print(f"Processed starBase to {dest}")
    except Exception as e:
        print(f"Error processing starBase: {e}")
        # Create empty file as placeholder
        pd.DataFrame(columns=['miRNA', 'lncRNA']).to_csv(dest, sep='\t', index=False)

def copy_annotation(src, dest):
    """Copy your pre-prepared gene_annotation.csv into the pipeline databases folder."""
    print(f"Copying gene annotation from {src} to {dest} ...")
    
    if os.path.exists(src):
        shutil.copy(src, dest)
        print(f"Copied annotation to {dest}")
    else:
        print(f"Warning: {src} not found. Creating placeholder annotation file.")
        # Create minimal annotation file
        pd.DataFrame({
            'gene_id': ['ENSG00000000001', 'ENSG00000000002'],
            'gene_name': ['GENE1', 'GENE2'],
            'gene_biotype': ['protein_coding', 'lncRNA']
        }).to_csv(dest, index=False)

def main():
    # Ensure the output folder exists
    os.makedirs("databases", exist_ok=True)

    # Input files (must exist before running)
    mirtarbase_src = "databases/miRTarBase_MTI.txt"
    starbase_src = "databases/starBase_miRNA_lncRNA.txt"
    annotation_src = "gene_annotation.csv"

    # Output files
    mirtarbase_dest = "databases/miRTarBase.txt"
    starbase_dest = "databases/LncBase.txt"
    annotation_dest = "databases/gene_annotation.csv"

    # Process each file
    process_mirtarbase(mirtarbase_src, mirtarbase_dest)
    process_starbase(starbase_src, starbase_dest)
    copy_annotation(annotation_src, annotation_dest)

    print("Finished processing all interaction databases.")

if __name__ == "__main__":
    main()