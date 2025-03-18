#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import scanpy as sc
from collections import defaultdict
from sklearn.neighbors import KNeighborsClassifier
import warnings
warnings.filterwarnings('ignore')

# Set working directory
os.chdir('/Users/jacob/Desktop/astrocytes-ontology-research')

print("Loading interaction sites data...")
# Load interaction sites identified in the previous analysis
interaction_sites = pd.read_csv('interaction_sites.csv')

# Summarize interaction pairs
pair_counts = interaction_sites['pair'].value_counts()
print("\nInteraction pairs summary:")
print(pair_counts)

print("\nLoading single-cell data to get cell type information...")
# Load the single-cell data with cell type annotations
adata = sc.read_h5ad('raw/CTR081_Fron.h5ad')

# Get cell type information
print("\nCell types in the dataset:")
cell_types = adata.obs['celltype'].unique()
print(sorted(cell_types))

# Function to extract ligand and receptor from pair string
def split_lr_pair(pair):
    return tuple(pair.split('-'))

# Extract unique ligands and receptors from the interaction pairs
unique_pairs = interaction_sites['pair'].unique()
unique_ligands = set()
unique_receptors = set()

for pair in unique_pairs:
    ligand, receptor = split_lr_pair(pair)
    unique_ligands.add(ligand)
    unique_receptors.add(receptor)

print(f"\nFound {len(unique_ligands)} unique ligands and {len(unique_receptors)} unique receptors in interactions")

# Create a mapping between genes and cell types based on expression
print("\nAnalyzing gene expression patterns across cell types...")

# Get average expression of each gene in each cell type
cell_type_expr = {}

for cell_type in cell_types:
    # Subset cells of this type
    cells_mask = adata.obs['celltype'] == cell_type
    
    # Calculate average expression for ligands and receptors
    gene_expr = {}
    
    # For ligands
    for gene in unique_ligands:
        if gene in adata.var_names:
            expr_values = adata[cells_mask, gene].X
            if hasattr(expr_values, "toarray"):
                expr_values = expr_values.toarray()
            gene_expr[gene] = float(np.mean(expr_values))
        else:
            gene_expr[gene] = 0.0
            
    # For receptors
    for gene in unique_receptors:
        if gene in adata.var_names:
            expr_values = adata[cells_mask, gene].X
            if hasattr(expr_values, "toarray"):
                expr_values = expr_values.toarray()
            gene_expr[gene] = float(np.mean(expr_values))
        else:
            gene_expr[gene] = 0.0
    
    cell_type_expr[cell_type] = gene_expr

# Create DataFrame for ligand expression by cell type
ligand_expr_df = pd.DataFrame.from_dict({ct: {g: cell_type_expr[ct][g] for g in unique_ligands} 
                                         for ct in cell_types}).T
print("\nLigand expression by cell type (first 5 ligands):")
print(ligand_expr_df.iloc[:, :5])

# Create DataFrame for receptor expression by cell type
receptor_expr_df = pd.DataFrame.from_dict({ct: {g: cell_type_expr[ct][g] for g in unique_receptors} 
                                          for ct in cell_types}).T
print("\nReceptor expression by cell type (first 5 receptors):")
print(receptor_expr_df.iloc[:, :5])

# Save these expression matrices
ligand_expr_df.to_csv('ligand_expression_by_cell_type.csv')
receptor_expr_df.to_csv('receptor_expression_by_cell_type.csv')

# For each ligand-receptor pair, identify the most likely producing and receiving cell types
print("\nPredicting cell types involved in interactions...")

interaction_cell_types = []

for pair in unique_pairs:
    ligand, receptor = split_lr_pair(pair)
    
    # Get expression of this ligand across cell types
    ligand_expr = {ct: cell_type_expr[ct][ligand] for ct in cell_types}
    receptor_expr = {ct: cell_type_expr[ct][receptor] for ct in cell_types}
    
    # Get cell type with highest expression
    ligand_cell_type = max(ligand_expr.items(), key=lambda x: x[1])[0]
    receptor_cell_type = max(receptor_expr.items(), key=lambda x: x[1])[0]
    
    # Add to results
    interaction_cell_types.append({
        'pair': pair,
        'ligand': ligand,
        'receptor': receptor,
        'ligand_cell_type': ligand_cell_type,
        'receptor_cell_type': receptor_cell_type,
        'ligand_expression': ligand_expr[ligand_cell_type],
        'receptor_expression': receptor_expr[receptor_cell_type]
    })

# Convert to DataFrame and save
interaction_cell_types_df = pd.DataFrame(interaction_cell_types)
interaction_cell_types_df.to_csv('interaction_cell_types.csv', index=False)

print("\nPredicted cell types for each ligand-receptor pair:")
print(interaction_cell_types_df[['pair', 'ligand_cell_type', 'receptor_cell_type']])

# Find the most common interaction patterns
print("\nMost common cell-cell interaction patterns:")
interaction_patterns = interaction_cell_types_df.groupby(['ligand_cell_type', 'receptor_cell_type']).size().reset_index(name='count')
interaction_patterns = interaction_patterns.sort_values('count', ascending=False)
print(interaction_patterns.head(10))

# Create a heatmap of interaction patterns
plt.figure(figsize=(12, 10))
interaction_matrix = pd.crosstab(
    interaction_cell_types_df['ligand_cell_type'],
    interaction_cell_types_df['receptor_cell_type']
)
sns.heatmap(interaction_matrix, annot=True, cmap='viridis', fmt='d')
plt.title('Predicted Cell-Cell Interaction Patterns')
plt.xlabel('Receptor-Expressing Cell Type')
plt.ylabel('Ligand-Producing Cell Type')
plt.tight_layout()
plt.savefig('cell_interaction_patterns.png', dpi=300)
print("Cell interaction patterns heatmap saved to cell_interaction_patterns.png")

# Create a circos-like plot for cell-cell interactions
try:
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    
    def draw_chord(ax, start, end, radius=1.0, width=0.1, alpha=0.5, color='blue'):
        # Draw a chord between two points on a circle
        # Adapted from matplotlib documentation
        rad = np.pi/180.
        start_angle = start * rad
        end_angle = end * rad
        
        # Create the bezier curve
        n_points = 100
        theta1 = np.linspace(start_angle, end_angle, n_points)
        theta2 = np.linspace(end_angle, start_angle, n_points)
        
        # Fix for chord drawing
        if abs(end_angle - start_angle) > np.pi:
            if end_angle > start_angle:
                theta1 = np.linspace(start_angle, end_angle - 2*np.pi, n_points)
                theta2 = np.linspace(end_angle - 2*np.pi, start_angle, n_points)
            else:
                theta1 = np.linspace(start_angle, end_angle + 2*np.pi, n_points)
                theta2 = np.linspace(end_angle + 2*np.pi, start_angle, n_points)
        
        # Path for the chord
        r1 = radius + width
        r2 = radius
        
        verts = np.vstack([
            np.column_stack([r1 * np.cos(theta1), r1 * np.sin(theta1)]),
            np.column_stack([r2 * np.cos(theta2), r2 * np.sin(theta2)]),
            np.column_stack([r1 * np.cos(theta1[0]), r1 * np.sin(theta1[0])])
        ])
        
        codes = [Path.MOVETO] + [Path.LINETO] * (len(verts) - 2) + [Path.CLOSEPOLY]
        path = Path(verts, codes)
        patch = PathPatch(path, facecolor=color, alpha=alpha, edgecolor='none')
        ax.add_patch(patch)
    
    # Get unique cell types
    unique_cell_types = np.unique(np.concatenate([
        interaction_cell_types_df['ligand_cell_type'].unique(),
        interaction_cell_types_df['receptor_cell_type'].unique()
    ]))
    
    # Create a mapping between cell types and angles
    n_cell_types = len(unique_cell_types)
    angle_map = {ct: i * 360 / n_cell_types for i, ct in enumerate(unique_cell_types)}
    
    # Create a color map for cell types
    color_map = dict(zip(unique_cell_types, plt.cm.tab20(np.linspace(0, 1, n_cell_types))))
    
    # Prepare the figure
    plt.figure(figsize=(14, 14))
    ax = plt.subplot(111, polar=True)
    
    # Draw the cell type segments
    for i, ct in enumerate(unique_cell_types):
        start_angle = i * 2 * np.pi / n_cell_types
        end_angle = (i + 1) * 2 * np.pi / n_cell_types
        
        # Draw segment
        ax.bar(start_angle, 1.0, width=(end_angle - start_angle), bottom=0.0,
               color=color_map[ct], edgecolor='white', alpha=0.8)
        
        # Add label
        mid_angle = (start_angle + end_angle) / 2
        x = 1.3 * np.cos(mid_angle)
        y = 1.3 * np.sin(mid_angle)
        
        # Rotate text label based on position
        if 0 <= mid_angle < np.pi:
            plt.text(mid_angle, 1.3, ct, ha='center', va='bottom', rotation=mid_angle * 180 / np.pi - 90)
        else:
            plt.text(mid_angle, 1.3, ct, ha='center', va='top', rotation=mid_angle * 180 / np.pi + 90)
    
    # Draw connections for each interaction
    for _, row in interaction_patterns.iterrows():
        source = row['ligand_cell_type']
        target = row['receptor_cell_type']
        weight = row['count']
        
        if source == target:
            continue  # Skip self-loops for simplicity
        
        source_angle = angle_map[source]
        target_angle = angle_map[target]
        
        alpha = min(0.8, 0.2 + 0.6 * weight / interaction_patterns['count'].max())
        
        # Get source cell type color
        color = color_map[source]
        
        draw_chord(ax, source_angle, target_angle, radius=0.7, width=0.05, 
                  alpha=alpha, color=color)
    
    # Remove axis elements
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['polar'].set_visible(False)
    
    plt.title('Cell-Cell Interaction Network', size=20, pad=20)
    plt.tight_layout()
    plt.savefig('cell_interaction_network.png', dpi=300, bbox_inches='tight')
    print("Cell interaction network plot saved to cell_interaction_network.png")
except Exception as e:
    print(f"Error creating circos plot: {e}")

# Identify specific regions with high interaction density
print("\nIdentifying interaction hotspots in the tissue...")

# Use DBSCAN clustering to identify high-density regions
try:
    from sklearn.cluster import DBSCAN
    
    # Extract coordinates from interaction sites
    coords = interaction_sites[['x', 'y']].values
    
    # Run DBSCAN clustering
    dbscan = DBSCAN(eps=200, min_samples=5)  # Parameters can be adjusted
    interaction_sites['cluster'] = dbscan.fit_predict(coords)
    
    # Get statistics for clusters
    cluster_stats = interaction_sites.groupby('cluster').agg({
        'pair': 'count',
        'x': 'mean',
        'y': 'mean',
        'distance': 'mean'
    }).rename(columns={'pair': 'interaction_count'})
    
    # Filter out noise points (cluster = -1)
    valid_clusters = cluster_stats[cluster_stats.index >= 0]
    
    print(f"Identified {len(valid_clusters)} interaction hotspots in the tissue")
    print("Hotspot statistics:")
    print(valid_clusters)
    
    # Save cluster information
    valid_clusters.to_csv('interaction_hotspots.csv')
    
    # Visualize the hotspots
    plt.figure(figsize=(12, 10))
    
    # Plot all interaction sites, colored by cluster
    scatter = plt.scatter(interaction_sites['x'], interaction_sites['y'], 
                         c=interaction_sites['cluster'], cmap='viridis', 
                         s=20, alpha=0.7)
    
    # Highlight cluster centers
    for idx, row in valid_clusters.iterrows():
        plt.scatter(row['x'], row['y'], s=300, 
                   facecolors='none', edgecolors='red', 
                   linewidths=2, alpha=0.8)
        plt.text(row['x'], row['y'], f"Cluster {idx}\n({int(row['interaction_count'])} interactions)", 
                ha='center', va='center', fontsize=8)
    
    plt.colorbar(scatter, label='Cluster ID (-1 is noise)')
    plt.title('Interaction Hotspots in Tissue')
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.tight_layout()
    plt.savefig('interaction_hotspots.png', dpi=300)
    print("Interaction hotspots visualization saved to interaction_hotspots.png")
    
    # For each hotspot, analyze the cell types involved
    if len(valid_clusters) > 0:
        hotspot_analysis = []
        
        for cluster_id in valid_clusters.index:
            # Get interactions in this cluster
            cluster_interactions = interaction_sites[interaction_sites['cluster'] == cluster_id]
            
            # Get unique interaction pairs in this cluster
            cluster_pairs = cluster_interactions['pair'].unique()
            
            # For each pair, get the cell types involved
            cluster_cell_types = set()
            for pair in cluster_pairs:
                pair_data = interaction_cell_types_df[interaction_cell_types_df['pair'] == pair]
                if not pair_data.empty:
                    cluster_cell_types.add(pair_data.iloc[0]['ligand_cell_type'])
                    cluster_cell_types.add(pair_data.iloc[0]['receptor_cell_type'])
            
            # Add to analysis
            hotspot_analysis.append({
                'cluster_id': cluster_id,
                'interaction_count': len(cluster_interactions),
                'unique_pairs': len(cluster_pairs),
                'pairs': ', '.join(cluster_pairs),
                'cell_types': ', '.join(cluster_cell_types),
                'x': valid_clusters.loc[cluster_id, 'x'],
                'y': valid_clusters.loc[cluster_id, 'y']
            })
        
        # Convert to DataFrame and save
        hotspot_analysis_df = pd.DataFrame(hotspot_analysis)
        hotspot_analysis_df.to_csv('hotspot_cell_types.csv', index=False)
        
        print("\nHotspot cell type analysis:")
        print(hotspot_analysis_df[['cluster_id', 'interaction_count', 'cell_types']])
except Exception as e:
    print(f"Error in hotspot analysis: {e}")

print("\nIntegrated spatial analysis completed successfully!") 