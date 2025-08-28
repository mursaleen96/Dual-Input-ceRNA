# modules/generate_report.py

import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from plotly.offline import plot
import os
import numpy as np

def main():
    # Use ranked triplets file
    validated_path = "results/validated_triplets_ranked.csv"
    network_path = "results/cerna_network.graphml"
    report_path = "results/cerna_analysis_report.html"

    # Load validated triplets
    try:
        validated = pd.read_csv(validated_path)
        print(f"Loaded {len(validated)} ranked validated triplets")
    except:
        validated = pd.DataFrame()
        print("No validated triplets found")

    # Load network
    try:
        G = nx.read_graphml(network_path)
        print(f"Loaded network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    except:
        G = nx.Graph()
        print("No network found, using empty graph")

    # Generate network visualization
    if G.number_of_nodes() > 0:
        network_html = create_network_visualization(G)
    else:
        network_html = "<p>No network data available for visualization.</p>"

    # Generate summary statistics
    summary_stats = generate_summary_stats(validated, G)

    # Create ranking summary table
    ranking_table_html = create_ranking_table(validated)

    # Create HTML report
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>ceRNA Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
            .header {{ background-color: #f4f4f4; padding: 20px; border-radius: 5px; margin-bottom: 30px; }}
            .summary {{ background-color: #e9f5ff; padding: 20px; border-radius: 5px; margin-bottom: 30px; }}
            .section {{ margin-bottom: 40px; }}
            .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }}
            .stat-box {{ background: white; border: 1px solid #ddd; padding: 15px; border-radius: 5px; text-align: center; }}
            .stat-number {{ font-size: 2em; font-weight: bold; color: #2c3e50; }}
            .stat-label {{ color: #7f8c8d; font-size: 0.9em; }}
            table {{ width: 100%; border-collapse: collapse; margin-top: 20px; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            .network-container {{ margin: 20px 0; min-height: 500px; border: 1px solid #ddd; border-radius: 5px; }}
            .highlight {{ background-color: #fff3cd; }}
            .scientific {{ font-family: 'Courier New', monospace; }}
        </style>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
        <div class="header">
            <h1>ceRNA Network Analysis Report</h1>
            <p><strong>Generated:</strong> {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Analysis:</strong> Competing Endogenous RNA Network Discovery</p>
            <p><strong>Input:</strong> Dual RNA-seq and miRNA-seq datasets</p>
        </div>

        <div class="summary">
            <h2>Executive Summary</h2>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="stat-number">{summary_stats['total_triplets']}</div>
                    <div class="stat-label">Validated & Ranked Triplets</div>
                </div>
                <div class="stat-box">
                    <div class="stat-number">{summary_stats['network_nodes']}</div>
                    <div class="stat-label">Network Nodes</div>
                </div>
                <div class="stat-box">
                    <div class="stat-number">{summary_stats['network_edges']}</div>
                    <div class="stat-label">Network Edges</div>
                </div>
                <div class="stat-box">
                    <div class="stat-number">{summary_stats['sig_after_correction']}</div>
                    <div class="stat-label">FDR-Significant</div>
                </div>
            </div>
            
            {f'<div class="highlight"><h3>Key Findings</h3><ul>' +
             f'<li><strong>Best P-value:</strong> {summary_stats["best_pvalue"]}</li>' +
             f'<li><strong>Best Adjusted P-value:</strong> {summary_stats["best_adj_pvalue"]}</li>' +
             f'<li><strong>Median Effect Size:</strong> {summary_stats["median_effect"]:.4f}</li>' +
             f'<li><strong>Data Sources:</strong> Combined RNA-seq and miRNA-seq</li></ul></div>' 
             if summary_stats['total_triplets'] > 0 else ''}
        </div>

        <div class="section">
            <h2>Top-Ranked ceRNA Triplets</h2>
            {ranking_table_html}
        </div>

        <div class="section">
            <h2>ceRNA Interaction Network</h2>
            <div class="network-container">
                {network_html}
            </div>
            <p><em>Network shows validated ceRNA interactions ranked by statistical significance. 
            Node colors represent gene types (blue: lncRNA, red: mRNA).</em></p>
        </div>

        <div class="section">
            <h2>Analysis Methods</h2>
            <h3>Enhanced Pipeline Features</h3>
            <ol>
                <li><strong>Dual Input Processing:</strong> Separate RNA-seq and miRNA-seq datasets for enhanced accuracy</li>
                <li><strong>Gene ID Conversion:</strong> Automatic conversion of Ensembl/Entrez IDs to gene names</li>
                <li><strong>Statistical Ranking:</strong> FDR-corrected p-values and comprehensive ranking system</li>
                <li><strong>Mediation Analysis:</strong> Sobel test for indirect effect significance</li>
                <li><strong>Network Construction:</strong> Top-ranked interactions form the final network</li>
            </ol>
            
            <h3>Statistical Measures</h3>
            <ul>
                <li><strong>Mediation P-value:</strong> Sobel test significance for indirect effect</li>
                <li><strong>Adjusted P-value:</strong> FDR-corrected for multiple testing</li>
                <li><strong>Ranking Score:</strong> Combined metric of significance and effect size</li>
                <li><strong>Indirect Effect:</strong> Magnitude of ceRNA-mediated regulation</li>
            </ul>
            
            <h3>Quality Control</h3>
            <ul>
                <li>Sample alignment between RNA-seq and miRNA-seq datasets</li>
                <li>Low-expression filtering and CPM normalization</li>
                <li>Multiple testing correction (Benjamini-Hochberg FDR)</li>
                <li>Comprehensive validation with mediation analysis</li>
            </ul>
        </div>

        <footer style="margin-top: 50px; padding-top: 20px; border-top: 1px solid #ddd; color: #666; font-size: 0.9em;">
            <p>Generated by Enhanced ceRNA Pipeline v2.0 | Dual RNA-seq + miRNA-seq Analysis | 
            For questions or issues, please contact the development team.</p>
        </footer>
    </body>
    </html>
    """

    # Write report
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Enhanced HTML report generated: {report_path}")

def create_ranking_table(validated):
    """Create HTML table for top-ranked triplets"""
    if validated.empty:
        return "<p>No validated triplets found.</p>"
    
    # Show top 20 or all if fewer
    display_count = min(20, len(validated))
    top_triplets = validated.head(display_count)
    
    # Format for display
    display_df = top_triplets[['final_rank', 'lncRNA', 'miRNA', 'mRNA', 
                               'mediation_pvalue', 'adjusted_pvalue', 
                               'indirect_effect', 'significant_after_correction']].copy()
    
    # Format scientific notation
    display_df['mediation_pvalue'] = display_df['mediation_pvalue'].apply(lambda x: f"{x:.2e}")
    display_df['adjusted_pvalue'] = display_df['adjusted_pvalue'].apply(lambda x: f"{x:.2e}")
    display_df['indirect_effect'] = display_df['indirect_effect'].apply(lambda x: f"{x:.4f}")
    
    # Rename columns
    display_df.columns = ['Rank', 'lncRNA', 'miRNA', 'mRNA', 'P-value', 'Adj. P-value', 
                         'Indirect Effect', 'FDR Significant']
    
    table_html = display_df.to_html(classes='table', escape=False, index=False)
    
    if len(validated) > display_count:
        table_html += f'<p><em>Showing top {display_count} of {len(validated)} total validated triplets. ' \
                     f'See results/validated_triplets_ranked.csv for complete results.</em></p>'
    
    return table_html

def create_network_visualization(G):
    """Create interactive network visualization with ranking information"""
    try:
        # Use spring layout for node positioning
        pos = nx.spring_layout(G, k=1/np.sqrt(G.number_of_nodes()), iterations=50)
        
        # Create edge traces with ranking-based styling
        edge_x = []
        edge_y = []
        edge_colors = []
        edge_widths = []
        
        for edge in G.edges(data=True):
            u, v = edge[0], edge[1]
            data = edge[2]
            
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            
            # Color and width based on rank (lower rank = better)
            rank = data.get('rank', float('inf'))
            if rank <= 10:
                edge_colors.extend(['red', 'red', None])
                edge_widths.extend([3, 3, None])
            elif rank <= 50:
                edge_colors.extend(['orange', 'orange', None])
                edge_widths.extend([2, 2, None])
            else:
                edge_colors.extend(['gray', 'gray', None])
                edge_widths.extend([1, 1, None])

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=2, color='rgba(125, 125, 125, 0.5)'),
            hoverinfo='none',
            mode='lines'
        )

        # Create node traces
        node_x = []
        node_y = []
        node_text = []
        node_info = []
        node_colors = []
        node_sizes = []

        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(node)
            
            # Node info for hover
            node_data = G.nodes[node]
            degree = G.degree(node)
            node_type = node_data.get('type', 'unknown')
            best_rank = node_data.get('best_rank', 'N/A')
            triplet_count = node_data.get('triplet_count', 0)
            
            info = f"<b>{node}</b><br>Type: {node_type}<br>Degree: {degree}<br>Best Rank: {best_rank}<br>Triplets: {triplet_count}"
            node_info.append(info)
            
            # Color by type
            if node_type == 'lncRNA':
                node_colors.append('lightblue')
            elif node_type == 'mRNA':
                node_colors.append('lightcoral')
            else:
                node_colors.append('lightgray')
            
            # Size based on degree and best rank
            base_size = 15
            size_bonus = min(degree * 3, 20) + (10 if best_rank != 'N/A' and best_rank <= 10 else 0)
            node_sizes.append(base_size + size_bonus)

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            hovertemplate='%{hovertext}<extra></extra>',
            hovertext=node_info,
            text=node_text,
            textposition='top center',
            textfont=dict(size=10),
            marker=dict(
                size=node_sizes,
                color=node_colors,
                line=dict(width=2, color='black')
            )
        )

        # Create figure
        fig = go.Figure(data=[edge_trace, node_trace])
        fig.update_layout(
            title=f"ceRNA Interaction Network ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)<br><sub>Colored by gene type, sized by importance</sub>",
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=60),
            annotations=[
                dict(
                    text="Red edges: Top 10 ranked • Orange edges: Top 50 ranked • Blue nodes: lncRNAs • Red nodes: mRNAs",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002,
                    xanchor='left', yanchor='bottom',
                    font=dict(color='gray', size=11)
                )
            ],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            height=600
        )

        return plot(fig, output_type='div', include_plotlyjs=True)
        
    except Exception as e:
        print(f"Error creating network visualization: {e}")
        return "<p>Error generating network visualization.</p>"

def generate_summary_stats(validated, G):
    """Generate enhanced summary statistics"""
    stats = {
        'total_triplets': len(validated) if not validated.empty else 0,
        'network_nodes': G.number_of_nodes(),
        'network_edges': G.number_of_edges(),
        'sig_after_correction': 0,
        'best_pvalue': 'N/A',
        'best_adj_pvalue': 'N/A',
        'median_effect': 0.0
    }
    
    if not validated.empty:
        if 'significant_after_correction' in validated.columns:
            stats['sig_after_correction'] = validated['significant_after_correction'].sum()
        
        if 'mediation_pvalue' in validated.columns:
            stats['best_pvalue'] = f"{validated['mediation_pvalue'].min():.2e}"
        
        if 'adjusted_pvalue' in validated.columns:
            stats['best_adj_pvalue'] = f"{validated['adjusted_pvalue'].min():.2e}"
        
        if 'indirect_effect' in validated.columns:
            stats['median_effect'] = validated['indirect_effect'].abs().median()
    
    return stats

if __name__ == "__main__":
    main()