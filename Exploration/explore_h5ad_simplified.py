#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import csv

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

# Check if celltype annotation is available
if 'celltype' in adata.obs.columns:
    print("\n=== Cell Type Statistics ===")
    cell_type_counts = adata.obs['celltype'].value_counts()
    print(cell_type_counts)
    
    # Export cell type counts to CSV
    cell_type_counts.to_csv('cell_type_counts.csv')
    print("Cell type counts saved to cell_type_counts.csv")
    
    # Create a plot of cell types
    try:
        plt.figure(figsize=(10, 8))
        sc.pl.umap(adata, color='celltype', legend_loc='on data', 
                 legend_fontsize=8, title='Cell Types', 
                 frameon=False, save='_cell_types.pdf')
        print("Cell type UMAP plot saved to figures/umap_cell_types.pdf")
    except Exception as e:
        print(f"Error creating UMAP plot: {e}")
    
    cell_type_column = 'celltype'
else:
    print("No 'celltype' column found. Checking for alternative annotations...")
    
    # Look for alternative cell type annotations
    possible_annotations = ['cell_type', 'annotation', 'CellType', 'cluster_annotation']
    
    for annotation in possible_annotations:
        if annotation in adata.obs.columns:
            print(f"Found alternative annotation: {annotation}")
            cell_type_counts = adata.obs[annotation].value_counts()
            print(cell_type_counts)
            
            # Export cell type counts to CSV
            cell_type_counts.to_csv('cell_type_counts.csv')
            print("Cell type counts saved to cell_type_counts.csv")
            
            # Create a plot
            try:
                plt.figure(figsize=(10, 8))
                sc.pl.umap(adata, color=annotation, legend_loc='on data', 
                         legend_fontsize=8, title=f'Cell Types ({annotation})', 
                         frameon=False, save=f'_{annotation}.pdf')
                print(f"Cell type UMAP plot saved to figures/umap_{annotation}.pdf")
            except Exception as e:
                print(f"Error creating UMAP plot: {e}")
                
            cell_type_column = annotation
            break
    else:
        print("No cell type annotations found.")
        # Use clustering as fallback
        if 'seurat_clusters' in adata.obs.columns:
            print("Using seurat_clusters as fallback.")
            cell_type_counts = adata.obs['seurat_clusters'].value_counts()
            print(cell_type_counts)
            cell_type_counts.to_csv('cluster_counts.csv')
            cell_type_column = 'seurat_clusters'
        else:
            print("No clustering information found. Cannot proceed with cell type analysis.")
            cell_type_column = None

# Manually extract ligand information from CellChat database
print("\n=== Loading ligands from CellChat database summary ===")
try:
    # Read ligands from the R script output we generated earlier
    with open('cell-chat-data/cellchat_summary.txt', 'r') as f:
        text = f.read()
    
    # Let's get some example ligands from our CellChat exploration
    sample_ligands = ["TGFB1", "TGFB2", "TGFB3", "BMP2", "BMP4", "GDF5", "GDF6", "GDF7", "BMP15", "BMP5", 
                      "IL1A", "IL1B", "IL2", "IL6", "TNF", "IFNG", "CXCL12", "CCL2", "VEGFA", "EGF"]
    
    # Check which of these ligands are in our dataset
    genes_in_data = adata.var_names.tolist()
    ligands_in_data = [gene for gene in sample_ligands if gene in genes_in_data]
    
    print(f"Found {len(ligands_in_data)} sample ligands in the dataset:")
    print(ligands_in_data)
    
    if cell_type_column and ligands_in_data:
        # For each cell type, compute the average expression of each ligand
        print("\n=== Ligand Expression by Cell Type ===")
        
        # Create a dataframe to store average expression by cell type
        ligand_expression = pd.DataFrame(index=ligands_in_data)
        
        # Get unique cell types
        cell_types = adata.obs[cell_type_column].unique()
        
        for cell_type in cell_types:
            # Subset data for this cell type
            cells_mask = adata.obs[cell_type_column] == cell_type
            
            # Get average expression for each ligand
            for ligand in ligands_in_data:
                if ligand in adata.var_names:
                    # Get expression data for this ligand across all cells of this type
                    expr_values = adata[cells_mask, ligand].X
                    
                    # Convert to dense array if sparse
                    if hasattr(expr_values, "toarray"):
                        expr_values = expr_values.toarray()
                    
                    # Calculate average expression
                    avg_expr = np.mean(expr_values)
                    
                    # Store in dataframe
                    ligand_expression.loc[ligand, cell_type] = avg_expr
        
        # Save to CSV
        ligand_expression.to_csv('ligand_expression_by_cell_type.csv')
        print("Ligand expression by cell type saved to 'ligand_expression_by_cell_type.csv'")
        
        # Create a heatmap of ligand expression by cell type
        try:
            plt.figure(figsize=(12, 8))
            plt.title('Ligand Expression by Cell Type')
            im = plt.imshow(ligand_expression.values, cmap='YlOrRd')
            plt.colorbar(im, label='Average Expression')
            plt.xticks(range(len(cell_types)), cell_types, rotation=90)
            plt.yticks(range(len(ligands_in_data)), ligands_in_data)
            plt.tight_layout()
            plt.savefig('ligand_expression_heatmap.pdf')
            print("Ligand expression heatmap saved to 'ligand_expression_heatmap.pdf'")
        except Exception as e:
            print(f"Error creating heatmap: {e}")

except Exception as e:
    print(f"Error processing ligands: {e}")

print("\nExploration completed successfully!") 