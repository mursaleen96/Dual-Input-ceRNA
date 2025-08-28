# convert_gtf.py
import pandas as pd
import re

def gtf_to_csv(gtf_file, output_csv):
    genes = []
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'gene':
                attributes = parts[8]
                gene_id = re.search(r'gene_id "([^"]+)"', attributes)
                gene_name = re.search(r'gene_name "([^"]+)"', attributes)
                gene_biotype = re.search(r'gene_biotype "([^"]+)"', attributes)
                
                if gene_id:
                    genes.append({
                        'gene_id': gene_id.group(1),
                        'gene_name': gene_name.group(1) if gene_name else '',
                        'gene_biotype': gene_biotype.group(1) if gene_biotype else ''
                    })
    
    df = pd.DataFrame(genes)
    df.to_csv(output_csv, index=False)
    print(f"Converted {len(genes)} genes to {output_csv}")

if __name__ == "__main__":
    gtf_to_csv("Homo_sapiens.GRCh38.114.gtf", "gene_annotation.csv")
