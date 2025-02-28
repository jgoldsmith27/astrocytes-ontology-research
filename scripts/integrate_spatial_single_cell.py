#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, XSD, OWL
import scanpy as sc
import anndata
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import argparse
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
from scipy.stats import pearsonr, spearmanr
import networkx as nx

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
SC = Namespace("http://example.org/single-cell#")

# Define astrocyte marker genes for different types
# Based on literature and single-cell data
ASTROCYTE_MARKERS = {
    'Protoplasmic': ['GFAP', 'S100B', 'AQP4', 'GJA1', 'SLC1A2', 'SLC1A3'],
    'Fibrous': ['GFAP', 'VIM', 'ALDH1L1', 'CD44', 'CRYAB', 'HOPX'],
    'Reactive': ['GFAP', 'VIM', 'SERPINA3', 'C3', 'CXCL10', 'STAT3', 'LCN2']
}

class SpatialSingleCellIntegrator:
    """
    A class to integrate spatial transcriptomics data with single-cell RNA sequencing data
    for comprehensive analysis of astrocyte heterogeneity.
    """
    
    def __init__(self, spatial_ttl, single_cell_h5ad, output_dir="../output/integrated_analysis"):
        """
        Initialize the integrator with paths to spatial and single-cell data files.
        
        Parameters:
        -----------
        spatial_ttl : str
            Path to the TURTLE file containing the spatial ontology
        single_cell_h5ad : str
            Path to the H5AD file containing the single-cell data
        output_dir : str
            Path to the directory where analysis results will be saved
        """
        self.spatial_ttl = spatial_ttl
        self.single_cell_h5ad = single_cell_h5ad
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize data structures
        self.spatial_graph = None
        self.single_cell_data = None
        self.gene_coordinates = {}
        self.astrocyte_spatial_data = None
        self.astrocyte_sc_data = None
        self.integrated_data = None
        
        # Load data
        self.load_data()
    
    def load_data(self):
        """
        Load spatial and single-cell data.
        """
        print("Loading data...")
        
        # Load spatial ontology
        print(f"Loading spatial ontology from {self.spatial_ttl}")
        self.spatial_graph = Graph()
        self.spatial_graph.parse(self.spatial_ttl, format="turtle")
        print(f"Spatial graph loaded with {len(self.spatial_graph)} triples")
        
        # Extract gene coordinates from spatial data
        self.gene_coordinates = self.get_gene_coordinates()
        
        # Load single-cell data
        print(f"Loading single-cell data from {self.single_cell_h5ad}")
        self.single_cell_data = sc.read_h5ad(self.single_cell_h5ad)
        print(f"Single-cell data loaded with {self.single_cell_data.shape[0]} cells and {self.single_cell_data.shape[1]} genes")
        
        # Process data
        self.process_spatial_data()
        self.process_single_cell_data()
    
    def get_gene_coordinates(self):
        """
        Extract coordinates for all genes from the spatial ontology.
        
        Returns:
        --------
        dict
            Dictionary mapping gene IDs to lists of (x, y, expression_level) tuples
        """
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?geneID ?x ?y ?expressionLevel
        WHERE {
            ?gene rdf:type :Gene ;
                  :expressedAt ?point ;
                  :hasExpressionLevel ?expressionLevel .
            
            ?point :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        """
        
        results = self.spatial_graph.query(query)
        
        gene_coords = {}
        for row in results:
            gene_id = str(row.geneID)
            x = float(row.x)
            y = float(row.y)
            expr = float(row.expressionLevel)
            
            if gene_id not in gene_coords:
                gene_coords[gene_id] = []
            
            gene_coords[gene_id].append((x, y, expr))
        
        return gene_coords
    
    def process_spatial_data(self):
        """
        Process spatial data to identify astrocyte types.
        """
        print("Processing spatial data...")
        
        # Create a spatial grid to aggregate gene expression
        grid_size = 50  # Size of grid cells
        grid_data = {}
        
        # Aggregate gene expression in grid cells
        for gene, coords_list in self.gene_coordinates.items():
            for x, y, expr in coords_list:
                # Discretize coordinates to grid cells
                grid_x = int(x // grid_size)
                grid_y = int(y // grid_size)
                grid_key = (grid_x, grid_y)
                
                if grid_key not in grid_data:
                    grid_data[grid_key] = {
                        'grid_x': grid_x,
                        'grid_y': grid_y,
                        'center_x': grid_x * grid_size + grid_size / 2,
                        'center_y': grid_y * grid_size + grid_size / 2
                    }
                
                # Store expression value in grid
                grid_data[grid_key][gene] = max(grid_data[grid_key].get(gene, 0), expr)
        
        # Convert to DataFrame
        grid_df = pd.DataFrame(list(grid_data.values()))
        
        # Calculate scores for each astrocyte type
        for cell_type, markers in ASTROCYTE_MARKERS.items():
            # Initialize score column
            grid_df[f'{cell_type}_score'] = 0
            
            # Sum expression of marker genes (if present in data)
            for marker in markers:
                if marker in grid_df.columns:
                    grid_df[f'{cell_type}_score'] += grid_df[marker].fillna(0)
            
            # Normalize by number of markers found
            markers_found = sum(1 for marker in markers if marker in grid_df.columns)
            if markers_found > 0:
                grid_df[f'{cell_type}_score'] /= markers_found
        
        # Determine dominant astrocyte type for each grid cell
        type_columns = [f'{cell_type}_score' for cell_type in ASTROCYTE_MARKERS.keys()]
        grid_df['dominant_type'] = grid_df[type_columns].idxmax(axis=1)
        grid_df['dominant_type'] = grid_df['dominant_type'].str.replace('_score', '')
        
        # Store results
        self.astrocyte_spatial_data = grid_df
        
        print(f"Identified {len(grid_df)} grid cells with astrocyte expression")
        for cell_type in ASTROCYTE_MARKERS.keys():
            count = sum(grid_df['dominant_type'] == cell_type)
            print(f"  - {cell_type}: {count} cells")
    
    def process_single_cell_data(self):
        """
        Process single-cell data to identify astrocyte types.
        """
        print("Processing single-cell data...")
        
        # Filter for astrocytes if cell type annotations are available
        if 'cell_type' in self.single_cell_data.obs.columns:
            astrocyte_mask = self.single_cell_data.obs['cell_type'].str.contains('astrocyte', case=False)
            self.astrocyte_sc_data = self.single_cell_data[astrocyte_mask].copy()
        else:
            # If no cell type annotations, use all cells
            self.astrocyte_sc_data = self.single_cell_data.copy()
        
        # Calculate scores for each astrocyte type
        for cell_type, markers in ASTROCYTE_MARKERS.items():
            # Initialize score column
            self.astrocyte_sc_data.obs[f'{cell_type}_score'] = 0
            
            # Sum expression of marker genes (if present in data)
            for marker in markers:
                if marker in self.astrocyte_sc_data.var_names:
                    self.astrocyte_sc_data.obs[f'{cell_type}_score'] += self.astrocyte_sc_data[:, marker].X.toarray().flatten()
            
            # Normalize by number of markers found
            markers_found = sum(1 for marker in markers if marker in self.astrocyte_sc_data.var_names)
            if markers_found > 0:
                self.astrocyte_sc_data.obs[f'{cell_type}_score'] /= markers_found
        
        # Determine dominant astrocyte type for each cell
        type_columns = [f'{cell_type}_score' for cell_type in ASTROCYTE_MARKERS.keys()]
        self.astrocyte_sc_data.obs['dominant_type'] = self.astrocyte_sc_data.obs[type_columns].idxmax(axis=1)
        self.astrocyte_sc_data.obs['dominant_type'] = self.astrocyte_sc_data.obs['dominant_type'].str.replace('_score', '')
        
        print(f"Processed {self.astrocyte_sc_data.shape[0]} astrocytes from single-cell data")
        for cell_type in ASTROCYTE_MARKERS.keys():
            count = sum(self.astrocyte_sc_data.obs['dominant_type'] == cell_type)
            print(f"  - {cell_type}: {count} cells")
    
    def integrate_data(self):
        """
        Integrate spatial and single-cell data.
        """
        print("Integrating spatial and single-cell data...")
        
        # Create a common gene set
        spatial_genes = set(self.astrocyte_spatial_data.columns) - set(['grid_x', 'grid_y', 'center_x', 'center_y', 'dominant_type'] + 
                                                                      [f'{cell_type}_score' for cell_type in ASTROCYTE_MARKERS.keys()])
        sc_genes = set(self.astrocyte_sc_data.var_names)
        common_genes = list(spatial_genes.intersection(sc_genes))
        
        print(f"Found {len(common_genes)} genes in common between spatial and single-cell data")
        
        # Create expression matrices for common genes
        spatial_expr = self.astrocyte_spatial_data[common_genes].values
        sc_expr = self.astrocyte_sc_data[:, common_genes].X.toarray()
        
        # Calculate correlation between spatial and single-cell expression patterns
        gene_correlations = []
        for i, gene in enumerate(common_genes):
            # Calculate correlation between spatial and average single-cell expression
            spatial_gene_expr = spatial_expr[:, i]
            sc_gene_expr = sc_expr[:, i]
            
            # Calculate average expression per cell type in single-cell data
            sc_type_expr = {}
            for cell_type in ASTROCYTE_MARKERS.keys():
                mask = self.astrocyte_sc_data.obs['dominant_type'] == cell_type
                if sum(mask) > 0:
                    sc_type_expr[cell_type] = np.mean(sc_gene_expr[mask])
                else:
                    sc_type_expr[cell_type] = 0
            
            # Calculate average expression per cell type in spatial data
            spatial_type_expr = {}
            for cell_type in ASTROCYTE_MARKERS.keys():
                mask = self.astrocyte_spatial_data['dominant_type'] == cell_type
                if sum(mask) > 0:
                    spatial_type_expr[cell_type] = np.mean(spatial_gene_expr[mask])
                else:
                    spatial_type_expr[cell_type] = 0
            
            # Calculate correlation between cell type expressions
            cell_types = list(ASTROCYTE_MARKERS.keys())
            sc_values = [sc_type_expr[ct] for ct in cell_types]
            spatial_values = [spatial_type_expr[ct] for ct in cell_types]
            
            if len(set(sc_values)) > 1 and len(set(spatial_values)) > 1:
                corr, p_value = pearsonr(sc_values, spatial_values)
            else:
                corr, p_value = 0, 1
            
            gene_correlations.append({
                'gene': gene,
                'correlation': corr,
                'p_value': p_value,
                'sc_protoplasmic': sc_type_expr.get('Protoplasmic', 0),
                'sc_fibrous': sc_type_expr.get('Fibrous', 0),
                'sc_reactive': sc_type_expr.get('Reactive', 0),
                'spatial_protoplasmic': spatial_type_expr.get('Protoplasmic', 0),
                'spatial_fibrous': spatial_type_expr.get('Fibrous', 0),
                'spatial_reactive': spatial_type_expr.get('Reactive', 0)
            })
        
        # Convert to DataFrame
        self.integrated_data = pd.DataFrame(gene_correlations)
        
        # Sort by correlation
        self.integrated_data = self.integrated_data.sort_values('correlation', ascending=False)
        
        print(f"Integrated data created with {len(self.integrated_data)} genes")
        
        # Save integrated data
        output_file = os.path.join(self.output_dir, "integrated_gene_data.csv")
        self.integrated_data.to_csv(output_file, index=False)
        print(f"Integrated data saved to {output_file}")
    
    def visualize_integrated_data(self):
        """
        Create visualizations of the integrated data.
        """
        print("Creating visualizations...")
        
        # Create PDF for all visualizations
        pdf_file = os.path.join(self.output_dir, "integrated_visualizations.pdf")
        pdf = PdfPages(pdf_file)
        
        # 1. Visualize astrocyte type distribution in spatial data
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Define colors for each astrocyte type
        colors = {
            'Protoplasmic': 'blue',
            'Fibrous': 'red',
            'Reactive': 'green'
        }
        
        # Plot each astrocyte type
        for cell_type, color in colors.items():
            # Filter data for this cell type
            type_data = self.astrocyte_spatial_data[self.astrocyte_spatial_data['dominant_type'] == cell_type]
            
            if not type_data.empty:
                # Plot points
                ax.scatter(
                    type_data['center_x'], 
                    type_data['center_y'],
                    c=color,
                    s=type_data[f'{cell_type}_score'] * 20,  # Size based on score
                    alpha=0.7,
                    label=f'{cell_type} Astrocytes'
                )
        
        # Add legend and labels
        ax.legend(fontsize=12)
        ax.set_xlabel('X Coordinate', fontsize=14)
        ax.set_ylabel('Y Coordinate', fontsize=14)
        ax.set_title('Spatial Distribution of Astrocyte Types', fontsize=16)
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Save to PDF
        pdf.savefig(fig)
        plt.close(fig)
        
        # 2. Visualize astrocyte type distribution in single-cell data (UMAP)
        sc.pp.neighbors(self.astrocyte_sc_data)
        sc.tl.umap(self.astrocyte_sc_data)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot each astrocyte type
        for cell_type, color in colors.items():
            # Filter data for this cell type
            mask = self.astrocyte_sc_data.obs['dominant_type'] == cell_type
            
            if sum(mask) > 0:
                # Plot points
                ax.scatter(
                    self.astrocyte_sc_data.obsm['X_umap'][mask, 0],
                    self.astrocyte_sc_data.obsm['X_umap'][mask, 1],
                    c=color,
                    s=50,
                    alpha=0.7,
                    label=f'{cell_type} Astrocytes'
                )
        
        # Add legend and labels
        ax.legend(fontsize=12)
        ax.set_xlabel('UMAP 1', fontsize=14)
        ax.set_ylabel('UMAP 2', fontsize=14)
        ax.set_title('UMAP of Astrocyte Types in Single-Cell Data', fontsize=16)
        
        # Save to PDF
        pdf.savefig(fig)
        plt.close(fig)
        
        # 3. Visualize top correlated genes
        top_genes = self.integrated_data.head(20)
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create bar plot of correlations
        bars = ax.bar(top_genes['gene'], top_genes['correlation'], color='skyblue')
        
        # Add labels and title
        ax.set_xlabel('Gene', fontsize=14)
        ax.set_ylabel('Correlation between Spatial and Single-Cell', fontsize=14)
        ax.set_title('Top 20 Genes with Highest Correlation', fontsize=16)
        
        # Rotate x-axis labels
        plt.xticks(rotation=90)
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7, axis='y')
        
        # Tight layout
        plt.tight_layout()
        
        # Save to PDF
        pdf.savefig(fig)
        plt.close(fig)
        
        # 4. Visualize expression patterns of top marker genes
        for cell_type, markers in ASTROCYTE_MARKERS.items():
            # Find markers that are in both datasets
            common_markers = [m for m in markers if m in self.integrated_data['gene'].values]
            
            if common_markers:
                # Create figure with subplots for each marker
                fig, axes = plt.subplots(len(common_markers), 2, figsize=(15, 5 * len(common_markers)))
                
                if len(common_markers) == 1:
                    axes = np.array([axes])
                
                fig.suptitle(f'Expression Patterns of {cell_type} Astrocyte Markers', fontsize=18)
                
                for i, marker in enumerate(common_markers):
                    # Spatial visualization
                    if marker in self.gene_coordinates:
                        coords = self.gene_coordinates[marker]
                        x_vals = [x for x, y, expr in coords]
                        y_vals = [y for x, y, expr in coords]
                        expr_vals = [expr for x, y, expr in coords]
                        
                        scatter = axes[i, 0].scatter(
                            x_vals, 
                            y_vals, 
                            c=expr_vals, 
                            cmap='viridis', 
                            s=30, 
                            alpha=0.8
                        )
                        
                        # Add colorbar
                        cbar = plt.colorbar(scatter, ax=axes[i, 0])
                        cbar.set_label(f'{marker} Expression Level', fontsize=10)
                        
                        axes[i, 0].set_title(f'Spatial Expression of {marker}', fontsize=12)
                        axes[i, 0].set_xlabel('X Coordinate', fontsize=10)
                        axes[i, 0].set_ylabel('Y Coordinate', fontsize=10)
                    
                    # Single-cell visualization
                    if marker in self.astrocyte_sc_data.var_names:
                        sc.pl.umap(self.astrocyte_sc_data, color=marker, ax=axes[i, 1], show=False, title=f'Single-Cell Expression of {marker}')
                
                # Tight layout
                plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust for suptitle
                
                # Save to PDF
                pdf.savefig(fig)
                plt.close(fig)
        
        # 5. Visualize correlation heatmap between spatial and single-cell data
        # Prepare data for heatmap
        cell_types = list(ASTROCYTE_MARKERS.keys())
        sc_cols = [f'sc_{ct.lower()}' for ct in cell_types]
        spatial_cols = [f'spatial_{ct.lower()}' for ct in cell_types]
        
        # Get top 50 correlated genes
        top50_genes = self.integrated_data.head(50)
        
        # Create correlation matrix
        corr_matrix = np.zeros((len(cell_types), len(cell_types)))
        
        for i, sc_type in enumerate(cell_types):
            for j, spatial_type in enumerate(cell_types):
                sc_col = f'sc_{sc_type.lower()}'
                spatial_col = f'spatial_{spatial_type.lower()}'
                
                if sc_col in top50_genes.columns and spatial_col in top50_genes.columns:
                    corr, _ = pearsonr(top50_genes[sc_col], top50_genes[spatial_col])
                    corr_matrix[i, j] = corr
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))
        
        sns.heatmap(
            corr_matrix,
            annot=True,
            fmt=".2f",
            cmap="coolwarm",
            xticklabels=cell_types,
            yticklabels=cell_types,
            ax=ax
        )
        
        # Add labels and title
        ax.set_xlabel('Spatial Data Astrocyte Types', fontsize=14)
        ax.set_ylabel('Single-Cell Data Astrocyte Types', fontsize=14)
        ax.set_title('Correlation Between Astrocyte Types in Spatial and Single-Cell Data', fontsize=16)
        
        # Save to PDF
        pdf.savefig(fig)
        plt.close(fig)
        
        # Close PDF
        pdf.close()
        
        print(f"Visualizations saved to {pdf_file}")
    
    def run_analysis(self):
        """
        Run the full analysis pipeline.
        """
        self.integrate_data()
        self.visualize_integrated_data()
        
        print("Analysis complete!")

def main():
    """
    Main function to run the integrator.
    """
    parser = argparse.ArgumentParser(description="Integrate Spatial and Single-Cell Data for Astrocyte Analysis")
    parser.add_argument("--spatial", type=str, required=True, help="Path to the TURTLE file containing the spatial ontology")
    parser.add_argument("--single-cell", type=str, required=True, help="Path to the H5AD file containing the single-cell data")
    parser.add_argument("--output", type=str, default="../output/integrated_analysis", help="Path to the output directory")
    
    args = parser.parse_args()
    
    integrator = SpatialSingleCellIntegrator(
        spatial_ttl=args.spatial,
        single_cell_h5ad=args.single_cell,
        output_dir=args.output
    )
    
    integrator.run_analysis()

if __name__ == "__main__":
    main() 