#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
import os
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap
from collections import defaultdict

# Set working directory
os.chdir('/Users/jacob/Desktop/astrocytes-ontology-research')

# Define ligand-receptor pairs of interest based on the CellChat database
# This is a subset of well-known interactions
LR_PAIRS = [
    ("TGFB1", "TGFBR1"),  # TGF-beta signaling
    ("TGFB2", "TGFBR1"),
    ("IL1B", "IL1R1"),    # Inflammatory signaling
    ("IL6", "IL6R"),
    ("TNF", "TNFRSF1A"),
    ("BDNF", "NTRK2"),    # Neurotrophic signaling
    ("NGF", "NTRK1"),
    ("VEGFA", "FLT1"),    # Vascular signaling
    ("VEGFA", "KDR"),
    ("CXCL12", "CXCR4"),  # Chemokine signaling
    ("CCL2", "CCR2"),
    ("EGF", "EGFR"),      # Growth factor signaling
    ("PDGFB", "PDGFRB"),  # Pericyte signaling
    ("APOE", "LRP1"),     # Astrocyte-microglia signaling
    ("CX3CL1", "CX3CR1")  # Neuronal-microglia signaling
]

print("Loading spatial data...")
# Load the spatial data
spatial_data = pd.read_csv('raw/1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv')

print("Data loaded! Basic statistics:")
print(f"Total rows: {len(spatial_data)}")
print(f"Columns: {spatial_data.columns.tolist()}")

# Get count of unique genes
unique_genes = spatial_data['geneID'].unique()
print(f"Number of unique genes: {len(unique_genes)}")

# Check if our ligands and receptors of interest are in the data
lr_genes = set([gene for pair in LR_PAIRS for gene in pair])
found_genes = set(unique_genes).intersection(lr_genes)
print(f"Found {len(found_genes)} out of {len(lr_genes)} ligands/receptors of interest:")
print(sorted(found_genes))

# Initialize a dictionary to store proximity data for each L-R pair
proximity_counts = defaultdict(int)
interaction_points = []

# Define interaction distance threshold (in spatial units)
INTERACTION_THRESHOLD = 50  # Can be adjusted based on tissue scale

print("\nAnalyzing potential interaction sites...")
# For each ligand-receptor pair, find potential interaction sites
for ligand, receptor in LR_PAIRS:
    if ligand not in unique_genes or receptor not in unique_genes:
        print(f"Skipping {ligand}-{receptor} pair: One or both genes not found in data")
        continue
        
    print(f"Analyzing {ligand}-{receptor} interactions...")
    
    # Get coordinates for ligand expressions
    ligand_points = spatial_data[spatial_data['geneID'] == ligand][['x', 'y']].values
    
    # Get coordinates for receptor expressions
    receptor_points = spatial_data[spatial_data['geneID'] == receptor][['x', 'y']].values
    
    if len(ligand_points) == 0 or len(receptor_points) == 0:
        print(f"  No spatial data for {ligand}-{receptor} pair")
        continue
        
    print(f"  Found {len(ligand_points)} ligand and {len(receptor_points)} receptor points")
    
    # Use KD-tree for efficient nearest neighbor search
    tree = cKDTree(receptor_points)
    
    # For each ligand point, find the nearest receptor point
    for ligand_point in ligand_points:
        dist, idx = tree.query(ligand_point, k=1)
        
        # If the nearest receptor is within the threshold distance
        if dist <= INTERACTION_THRESHOLD:
            pair_key = f"{ligand}-{receptor}"
            proximity_counts[pair_key] += 1
            
            # Store the midpoint between ligand and receptor
            midpoint = (ligand_point + receptor_points[idx]) / 2
            interaction_points.append({
                'pair': pair_key,
                'x': midpoint[0],
                'y': midpoint[1],
                'distance': dist
            })

# Convert interaction_points to DataFrame
interaction_df = pd.DataFrame(interaction_points)

# Save interaction data
if not interaction_df.empty:
    interaction_df.to_csv('interaction_sites.csv', index=False)
    print(f"\nFound {len(interaction_df)} potential interaction sites")
    print("Interaction sites saved to interaction_sites.csv")
    
    # Summary of interactions by pair
    print("\nInteraction summary by ligand-receptor pair:")
    for pair, count in sorted(proximity_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {pair}: {count} potential interaction sites")
else:
    print("\nNo potential interaction sites found with the current threshold")

# Create visualization of the spatial data and potential interaction sites
print("\nCreating visualization...")
try:
    # Plot overall spatial data distribution (subsample if too large)
    plt.figure(figsize=(12, 10))
    
    if len(spatial_data) > 10000:
        sample_size = 10000
        sample_data = spatial_data.sample(sample_size)
    else:
        sample_data = spatial_data
    
    # Create a background plot of all gene expressions
    plt.scatter(sample_data['x'], sample_data['y'], s=1, alpha=0.1, color='lightgray')
    
    # Add interaction sites
    if not interaction_df.empty:
        # Create a colormap for different L-R pairs
        pair_types = interaction_df['pair'].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(pair_types)))
        pair_color_map = {pair: colors[i] for i, pair in enumerate(pair_types)}
        
        # Plot each interaction site
        for pair in pair_types:
            pair_data = interaction_df[interaction_df['pair'] == pair]
            plt.scatter(pair_data['x'], pair_data['y'], 
                       s=20, alpha=0.7, label=pair, 
                       color=pair_color_map[pair])
    
    plt.title('Spatial Distribution of Gene Expression and Potential Interaction Sites')
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('spatial_interactions.png', dpi=300)
    print("Spatial visualization saved to spatial_interactions.png")
    
    # Generate heatmap of interaction density
    if not interaction_df.empty and len(interaction_df) > 10:
        plt.figure(figsize=(10, 8))
        
        # Create 2D histogram (heatmap) of interaction sites
        heatmap, xedges, yedges = np.histogram2d(
            interaction_df['x'], interaction_df['y'], 
            bins=50, 
            range=[[min(spatial_data['x']), max(spatial_data['x'])], 
                   [min(spatial_data['y']), max(spatial_data['y'])]]
        )
        
        # Smooth the heatmap
        heatmap = np.ma.masked_where(heatmap == 0, heatmap)
        
        plt.imshow(heatmap.T, cmap='viridis', origin='lower', 
                  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], 
                  aspect='auto', interpolation='bilinear')
        
        plt.colorbar(label='Interaction Density')
        plt.title('Heatmap of Potential Cell-Cell Interaction Density')
        plt.xlabel('X coordinate')
        plt.ylabel('Y coordinate')
        plt.tight_layout()
        plt.savefig('interaction_heatmap.png', dpi=300)
        print("Interaction density heatmap saved to interaction_heatmap.png")
        
except Exception as e:
    print(f"Error creating visualization: {e}")

print("\nSpatial interaction analysis completed!") 