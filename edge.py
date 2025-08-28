import pandas as pd

# === CONFIGURATION ===
input_file = "validated_triplets.csv"     # Your original triplets file
output_file = "edges_for_cytoscape.csv"   # Output file for Cytoscape

# === LOAD DATA ===
df = pd.read_csv(input_file)

# === CREATE EDGES ===
edges = []

for _, row in df.iterrows():
    # Edge 1: lncRNA -> miRNA
    edges.append({
        "Source": row["lncRNA"],
        "Target": row["miRNA"],
        "score": row["score"],
        "mediation_pvalue": row["mediation_pvalue"],
        "indirect_effect": row["indirect_effect"],
        "sensitivity": row["sensitivity"],
        "mediation_type": row["mediation_type"]
    })

    # Edge 2: miRNA -> mRNA
    edges.append({
        "Source": row["miRNA"],
        "Target": row["mRNA"],
        "score": row["score"],
        "mediation_pvalue": row["mediation_pvalue"],
        "indirect_effect": row["indirect_effect"],
        "sensitivity": row["sensitivity"],
        "mediation_type": row["mediation_type"]
    })

# === SAVE EDGES ===
edges_df = pd.DataFrame(edges)
edges_df.to_csv(output_file, index=False)

print(f"Edge list saved to: {output_file}")