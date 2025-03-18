#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Set working directory
os.chdir('/Users/jacob/Desktop/astrocytes-ontology-research')

# Load the h5ad file
print("Loading h5ad file...")
adata = sc.read_h5ad('raw/CTR081_Fron.h5ad')

# Print basic information
print("\n=== Basic Information ===")
print(f"Dataset shape: {adata.shape}")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")

# Check available annotations
print("\n=== Available Annotations ===")
print("Observation (cell) annotations:")
for key in adata.obs.columns:
    print(f"  - {key}")

print("\nVariable (gene) annotations:")
for key in adata.var.columns:
    print(f"  - {key}")

# Check if cell type annotations are available
if 'cell_type' in adata.obs.columns:
    print("\n=== Cell Type Statistics ===")
    cell_type_counts = adata.obs['cell_type'].value_counts()
    print(cell_type_counts)
    
    # Export cell type counts to CSV
    cell_type_counts.to_csv('cell_type_counts.csv')
    print("Cell type counts saved to cell_type_counts.csv")
elif 'leiden' in adata.obs.columns:
    print("\n=== Leiden Cluster Statistics ===")
    leiden_counts = adata.obs['leiden'].value_counts()
    print(leiden_counts)
    
    # Export leiden counts to CSV
    leiden_counts.to_csv('leiden_counts.csv')
    print("Leiden cluster counts saved to leiden_counts.csv")

# Look for ligand and receptor genes in the dataset
print("\n=== Loading CellChat data ===")
try:
    # Try loading CellChat ligand-receptor pairs using rpy2
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    
    ro.r('''
    load("/Users/jacob/Desktop/astrocytes-ontology-research/cell-chat-data/CellChatDB.human.rda")
    ligands <- unique(CellChatDB.human$interaction$ligand)
    receptors <- unique(CellChatDB.human$interaction$receptor)
    write.csv(ligands, "cellchat_ligands.csv", row.names=FALSE)
    write.csv(receptors, "cellchat_receptors.csv", row.names=FALSE)
    ''')
    
    ligands = pd.read_csv('cellchat_ligands.csv')['x'].tolist()
    receptors = pd.read_csv('cellchat_receptors.csv')['x'].tolist()
    
    print(f"Loaded {len(ligands)} ligands and {len(receptors)} receptors from CellChatDB")
    
    # Check which ligands and receptors are present in the dataset
    genes_in_data = adata.var_names.tolist()
    
    ligands_in_data = [gene for gene in ligands if gene in genes_in_data]
    receptors_in_data = [gene for gene in receptors if gene in genes_in_data]
    
    print(f"Found {len(ligands_in_data)} ligands and {len(receptors_in_data)} receptors in the dataset")
    
    # Export the lists
    pd.DataFrame(ligands_in_data, columns=['gene']).to_csv('ligands_in_data.csv', index=False)
    pd.DataFrame(receptors_in_data, columns=['gene']).to_csv('receptors_in_data.csv', index=False)
    
    print("Ligands and receptors in the dataset saved to ligands_in_data.csv and receptors_in_data.csv")
    
    # Print some examples
    print("\nExample ligands found in the dataset:")
    print(ligands_in_data[:10])
    
    print("\nExample receptors found in the dataset:")
    print(receptors_in_data[:10])

except Exception as e:
    print(f"Error loading CellChat data with rpy2: {e}")
    print("Continuing without ligand-receptor analysis...")

# Save a overview plot of cell types
if 'cell_type' in adata.obs.columns:
    try:
        plt.figure(figsize=(10, 8))
        sc.pl.umap(adata, color='cell_type', legend_loc='on data', legend_fontsize=8, save='_cell_types.pdf')
        print("\nCell type UMAP plot saved to figures/umap_cell_types.pdf")
    except Exception as e:
        print(f"Error creating UMAP plot: {e}")

print("\nExploration completed successfully!") 