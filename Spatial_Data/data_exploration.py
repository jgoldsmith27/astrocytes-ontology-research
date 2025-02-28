import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.stats import pearsonr

# === Step 1: Load the H5AD File ===
file_path = "/Users/jacob/Desktop/Karolinska Research/astrocytes-ontology-research/single_cell_data_cleaned.h5ad"
adata = sc.read_h5ad(file_path)

# === Step 2: Basic Inspection ===
print("Dataset Summary:")
print(adata)

# Inspect cell metadata
if adata.obs.shape[1] > 0:
    print("\nCell Metadata Columns:")
    print(adata.obs.head())
else:
    print("\nNo cell metadata available.")

# Inspect gene metadata
if adata.var.shape[1] > 0:
    print("\nGene Metadata Columns:")
    print(adata.var.head())
else:
    print("\nNo gene metadata available.")

# Check for spatial data
if 'spatial' in adata.obsm.keys():
    print("\nSpatial coordinates are present in the dataset.")
else:
    print("\nNo spatial coordinates found.")

# Convert expression data to dense format
dense_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
genes = adata.var_names
cell_types = adata.obs['cell_type'].unique() if 'cell_type' in adata.obs.columns else None

print(f"\nNumber of genes: {len(genes)}")
if cell_types is not None:
    print(f"Unique Cell Types: {cell_types}")

# === Step 3: Gene Expression Distribution ===
plt.figure(figsize=(10, 5))
sns.histplot(dense_matrix.flatten(), bins=50, kde=True)
plt.title("Gene Expression Distribution")
plt.xlabel("Expression Level")
plt.ylabel("Frequency")
plt.show()

# === Step 4: Spatial Gene Expression (if available) ===
if 'spatial' in adata.obsm.keys():
    spatial_coords = adata.obsm['spatial']
    x, y = spatial_coords[:, 0], spatial_coords[:, 1]
    
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, c=np.mean(dense_matrix, axis=1), cmap='viridis', alpha=0.7)
    plt.title("Spatial Distribution of Gene Expression")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.colorbar(label="Mean Expression")
    plt.show()
else:
    print("Skipping spatial visualization: No spatial coordinates found.")

# === Step 5: Gene-Gene Coexpression Analysis ===
expr_df = pd.DataFrame(dense_matrix, columns=genes)

# Compute Pearson correlation matrix
correlation_matrix = expr_df.corr()

# Display heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(correlation_matrix, cmap="coolwarm", center=0)
plt.title("Gene-Gene Coexpression Heatmap")
plt.show()

# === Step 6: Construct Gene Coexpression Network ===
G = nx.Graph()

# Add genes as nodes
for gene in genes:
    G.add_node(gene, type="gene")

# Add edges based on coexpression threshold
threshold = 0.7  # Adjust based on significance
for gene1 in genes:
    for gene2 in genes:
        if gene1 != gene2 and correlation_matrix.loc[gene1, gene2] > threshold:
            G.add_edge(gene1, gene2, weight=correlation_matrix.loc[gene1, gene2])

# Visualize the coexpression network
plt.figure(figsize=(10, 6))
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=False, node_size=20, edge_color="gray")
plt.title("Gene-Gene Coexpression Network")
plt.show()

# === Step 7: Extract Data for Ontology Planning ===
ontology_data = {
    "num_genes": len(genes),
    "num_cells": adata.n_obs,
    "has_spatial": 'spatial' in adata.obsm.keys(),
    "cell_types": list(cell_types) if cell_types is not None else None,
    "highly_coexpressed_genes": correlation_matrix.stack()[
        correlation_matrix.stack() > 0.7].index.tolist(),
}