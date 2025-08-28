# modules/network_analysis.py

import pandas as pd
import networkx as nx
import os

def main():
    # Use the ranked triplets file
    validated_path = "results/validated_triplets_ranked.csv"
    network_path = "results/cerna_network.graphml"
    cytoscape_path = "results/cerna_network.sif"
    centrality_path = "results/centrality_scores.csv"
    nodes_path = "results/cerna_network_nodes.csv"
    edges_path = "results/cerna_network_edges.csv"

    # Load validated triplets
    try:
        validated = pd.read_csv(validated_path)
        print(f"Loaded {len(validated)} ranked validated triplets")
    except:
        validated = pd.DataFrame()
        print("No validated triplets found")
    
    if validated.empty:
        print("No validated triplets found. Creating empty network files.")
        # Create empty network
        G = nx.Graph()
        nx.write_graphml(G, network_path)
        
        # Create empty files
        pd.DataFrame(columns=['gene', 'degree_centrality']).to_csv(centrality_path, index=False)
        pd.DataFrame(columns=['gene', 'attributes']).to_csv(nodes_path, index=False)
        pd.DataFrame(columns=['source', 'target', 'attributes']).to_csv(edges_path, index=False)
        
        with open(cytoscape_path, 'w') as f:
            f.write("# Empty network\n")
        
        print("Empty network files created.")
        return

    # Build network using top-ranked triplets
    G = nx.Graph()
    edge_data = {}  # Store edge information
    
    for _, row in validated.iterrows():
        lnc = row['lncRNA']
        mrna = row['mRNA']
        mir = row['miRNA']
        score = row['score']
        rank = row.get('final_rank', 0)
        p_val = row['mediation_pvalue']
        adj_p_val = row.get('adjusted_pvalue', p_val)
        
        # Create edge key
        edge_key = tuple(sorted([lnc, mrna]))
        
        # If edge exists, keep the better-ranked triplet
        if edge_key in edge_data:
            if rank < edge_data[edge_key]['rank']:  # Lower rank is better
                edge_data[edge_key] = {
                    'miRNA': mir,
                    'score': score,
                    'rank': rank,
                    'mediation_pvalue': p_val,
                    'adjusted_pvalue': adj_p_val
                }
        else:
            edge_data[edge_key] = {
                'miRNA': mir,
                'score': score,
                'rank': rank,
                'mediation_pvalue': p_val,
                'adjusted_pvalue': adj_p_val
            }

    # Add edges to network
    for (node1, node2), data in edge_data.items():
        G.add_edge(node1, node2, **data)

    # Add gene names as node attributes
    all_genes = set(validated['lncRNA']) | set(validated['mRNA'])
    for gene in all_genes:
        if gene in G.nodes:
            G.nodes[gene]['name'] = gene
            # Determine gene type
            if gene in validated['lncRNA'].values:
                G.nodes[gene]['type'] = 'lncRNA'
            else:
                G.nodes[gene]['type'] = 'mRNA'
            
            # Add node statistics
            gene_triplets = validated[
                (validated['lncRNA'] == gene) | (validated['mRNA'] == gene)
            ]
            G.nodes[gene]['triplet_count'] = len(gene_triplets)
            G.nodes[gene]['best_rank'] = gene_triplets['final_rank'].min() if len(gene_triplets) > 0 else float('inf')
            G.nodes[gene]['avg_pvalue'] = gene_triplets['mediation_pvalue'].mean() if len(gene_triplets) > 0 else 1.0

    print(f"Network constructed: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Save as GraphML
    nx.write_graphml(G, network_path)
    print(f"Network saved to {network_path}")

    # Save as SIF for Cytoscape
    with open(cytoscape_path, "w") as f:
        for u, v, data in G.edges(data=True):
            mir = data.get('miRNA', 'unknown')
            rank = data.get('rank', 'unknown')
            f.write(f"{u}\tceRNA_interaction\t{v}\n")
    print(f"Cytoscape SIF file saved to {cytoscape_path}")

    # Compute and save centrality scores
    if G.number_of_nodes() > 0:
        centrality_data = []
        
        # Basic centrality measures
        degree_cent = nx.degree_centrality(G)
        
        # Additional measures for smaller networks
        if G.number_of_nodes() < 1000:
            try:
                betweenness_cent = nx.betweenness_centrality(G)
            except:
                betweenness_cent = {node: 0 for node in G.nodes()}
                
            try:
                if nx.is_connected(G):
                    closeness_cent = nx.closeness_centrality(G)
                else:
                    closeness_cent = {node: 0 for node in G.nodes()}
            except:
                closeness_cent = {node: 0 for node in G.nodes()}
        else:
            betweenness_cent = {node: 0 for node in G.nodes()}
            closeness_cent = {node: 0 for node in G.nodes()}
        
        # Combine centrality measures
        for node in G.nodes():
            node_data = G.nodes[node]
            centrality_data.append({
                'gene': node,
                'gene_type': node_data.get('type', 'unknown'),
                'degree_centrality': degree_cent.get(node, 0),
                'betweenness_centrality': betweenness_cent.get(node, 0),
                'closeness_centrality': closeness_cent.get(node, 0),
                'triplet_count': node_data.get('triplet_count', 0),
                'best_rank': node_data.get('best_rank', float('inf')),
                'avg_pvalue': node_data.get('avg_pvalue', 1.0)
            })
        
        centrality_df = pd.DataFrame(centrality_data)
        centrality_df = centrality_df.sort_values('degree_centrality', ascending=False)
        centrality_df.to_csv(centrality_path, index=False)
        print(f"Centrality scores saved to {centrality_path}")
    else:
        pd.DataFrame(columns=['gene', 'degree_centrality']).to_csv(centrality_path, index=False)

    # Export nodes to CSV
    try:
        nodes_data = []
        for node, data in G.nodes(data=True):
            node_record = {'gene': node}
            node_record.update(data)
            nodes_data.append(node_record)
        
        nodes_df = pd.DataFrame(nodes_data)
        nodes_df.to_csv(nodes_path, index=False)
        print(f"Nodes exported to {nodes_path}")
    except Exception as e:
        print(f"Error exporting nodes: {e}")
        pd.DataFrame(columns=['gene']).to_csv(nodes_path, index=False)

    # Export edges to CSV
    try:
        edges_data = []
        for u, v, data in G.edges(data=True):
            edge_record = {
                'source': u,
                'target': v,
                'source_type': G.nodes[u].get('type', 'unknown'),
                'target_type': G.nodes[v].get('type', 'unknown')
            }
            edge_record.update(data)
            edges_data.append(edge_record)
        
        edges_df = pd.DataFrame(edges_data)
        edges_df = edges_df.sort_values('rank') if 'rank' in edges_df.columns else edges_df
        edges_df.to_csv(edges_path, index=False)
        print(f"Edges exported to {edges_path}")
    except Exception as e:
        print(f"Error exporting edges: {e}")
        pd.DataFrame(columns=['source', 'target']).to_csv(edges_path, index=False)

    # Print network summary
    if G.number_of_nodes() > 0:
        print(f"\nNetwork Summary:")
        print(f"  Nodes: {G.number_of_nodes()}")
        print(f"  Edges: {G.number_of_edges()}")
        print(f"  Density: {nx.density(G):.4f}")
        print(f"  Connected components: {nx.number_connected_components(G)}")
        
        # Top hub genes
        if len(centrality_df) > 0:
            print(f"\nTop 5 hub genes (by degree centrality):")
            for _, row in centrality_df.head(5).iterrows():
                print(f"    {row['gene']} ({row['gene_type']}): {row['degree_centrality']:.3f}")

if __name__ == "__main__":
    main()