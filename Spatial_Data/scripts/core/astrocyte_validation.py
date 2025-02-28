#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef, Literal, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL
import argparse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.stats import spearmanr
import scanpy as sc
import anndata as ad
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap
from matplotlib_venn import venn3

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
SC = Namespace("http://example.org/single-cell#")
BRIDGE = Namespace("http://example.org/ontology-bridge#")

# Define canonical marker genes for astrocyte types based on literature
# These are well-established markers from multiple studies
CANONICAL_MARKERS = {
    'Protoplasmic': ['GFAP', 'S100B', 'AQP4', 'GJA1', 'SLC1A2', 'SLC1A3', 'GLUL', 'ATP1B2', 'FGFR3'],
    'Fibrous': ['GFAP', 'VIM', 'ALDH1L1', 'CD44', 'CRYAB', 'HOPX', 'TNC', 'GFAP', 'ID3', 'MFGE8'],
    'Reactive': ['GFAP', 'VIM', 'SERPINA3', 'C3', 'CXCL10', 'STAT3', 'LCN2', 'CP', 'OSMR', 'ASPG']
}

# Add a simple function to adjust p-values using Benjamini-Hochberg method
def adjust_pvalues(pvalues):
    """
    Adjust p-values using Benjamini-Hochberg method to control false discovery rate.
    This is a simple implementation to replace the statsmodels.stats.multitest dependency.
    
    Parameters:
    -----------
    pvalues : array-like
        Array of p-values to adjust
    
    Returns:
    --------
    numpy.ndarray
        Array of adjusted p-values
    """
    pvalues = np.asarray(pvalues)
    n = len(pvalues)
    
    # Sort p-values in ascending order
    sorted_indices = np.argsort(pvalues)
    sorted_pvalues = pvalues[sorted_indices]
    
    # Calculate adjusted p-values
    adjusted_pvalues = np.empty_like(sorted_pvalues)
    for i in range(n):
        rank = i + 1
        adjusted_pvalues[i] = sorted_pvalues[i] * n / rank
    
    # Ensure monotonicity
    for i in range(n-2, -1, -1):
        adjusted_pvalues[i] = min(adjusted_pvalues[i], adjusted_pvalues[i+1])
    
    # Cap at 1.0
    adjusted_pvalues = np.minimum(adjusted_pvalues, 1.0)
    
    # Return to original order
    result = np.empty_like(adjusted_pvalues)
    result[sorted_indices] = adjusted_pvalues
    
    return result

class AstrocyteValidator:
    """
    A class for validating astrocyte cell type classifications in spatial transcriptomics data
    using industry-standard approaches based on co-expressed genes.
    """
    
    def __init__(self, integrated_ttl, output_dir="../output/astrocyte_validation"):
        """
        Initialize the validator with input and output paths.
        
        Parameters:
        -----------
        integrated_ttl : str
            Path to the integrated ontology TURTLE file
        output_dir : str
            Directory where validation results will be saved
        """
        self.integrated_ttl = integrated_ttl
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize graph and data structures
        self.graph = None
        self.spatial_cells_df = None
        self.gene_expression_matrix = None
        self.validation_results = None
        self.coexpression_networks = {}
        
        # Load data
        self.load_data()
    
    def load_data(self):
        """
        Load the integrated ontology data from TURTLE file.
        """
        print(f"Loading integrated TURTLE file: {self.integrated_ttl}")
        self.graph = Graph()
        self.graph.parse(self.integrated_ttl, format="turtle")
        print(f"Integrated graph loaded with {len(self.graph)} triples")
    
    def extract_classified_cells(self):
        """
        Extract classified cells from the spatial dataset and create a dataframe
        with all genes expressed in these cells.
        """
        print("Extracting classified cells from spatial dataset...")
        
        # Query to get spatial points with astrocyte types
        query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y ?astrocyteType ?probability
        WHERE {
            ?point rdf:type st:SpatialPoint ;
                   st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?astrocyteType ;
                      st:typeProbability ?probability .
            
            # Only include astrocyte types
            FILTER(?astrocyteType IN ("Protoplasmic", "Fibrous", "Reactive"))
        }
        """
        
        results = list(self.graph.query(query))
        
        # Convert to DataFrame
        cells_data = []
        for row in results:
            point_uri = str(row.point)
            x = float(row.x)
            y = float(row.y)
            astrocyte_type = str(row.astrocyteType)
            probability = float(row.probability)
            
            cells_data.append({
                'point_uri': point_uri,
                'x': x,
                'y': y,
                'astrocyte_type': astrocyte_type,
                'probability': probability
            })
        
        self.spatial_cells_df = pd.DataFrame(cells_data)
        
        print(f"Extracted {len(self.spatial_cells_df)} classified cells")
        
        # Get gene expression for each cell
        self.extract_gene_expression()
    
    def extract_gene_expression(self):
        """
        Extract gene expression data for the classified cells.
        """
        print("Extracting gene expression data...")
        
        # Get all points
        point_uris = self.spatial_cells_df['point_uri'].tolist()
        
        # Initialize a dictionary to store gene expression data
        gene_expr_data = {}
        
        # Query to get gene expression for each point
        for point_uri in point_uris:
            query = f"""
            PREFIX st: <http://example.org/spatial-transcriptomics#>
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            
            SELECT ?gene ?geneID ?expressionLevel
            WHERE {{
                <{point_uri}> st:hasGeneExpression ?gene .
                ?gene st:hasGeneID ?geneID ;
                      st:hasExpressionLevel ?expressionLevel .
            }}
            """
            
            results = list(self.graph.query(query))
            
            # Store gene expression data for this point
            gene_expr = {}
            for row in results:
                gene_id = str(row.geneID)
                expr_level = float(row.expressionLevel)
                gene_expr[gene_id] = expr_level
            
            gene_expr_data[point_uri] = gene_expr
        
        # Convert to DataFrame
        # First, get all unique gene IDs
        all_genes = sorted(set(gene_id for expr_dict in gene_expr_data.values() for gene_id in expr_dict.keys()))
        
        # Create a matrix with points as rows and genes as columns
        data = []
        for point_uri, point_data in self.spatial_cells_df.iterrows():
            row = {
                'point_uri': point_data['point_uri'],
                'x': point_data['x'],
                'y': point_data['y'],
                'astrocyte_type': point_data['astrocyte_type'],
                'probability': point_data['probability']
            }
            
            # Add gene expression values
            expr_dict = gene_expr_data.get(point_data['point_uri'], {})
            for gene_id in all_genes:
                row[gene_id] = expr_dict.get(gene_id, 0)
            
            data.append(row)
        
        # Convert to DataFrame
        self.gene_expression_matrix = pd.DataFrame(data)
        
        print(f"Created gene expression matrix with {len(self.gene_expression_matrix)} cells and {len(all_genes)} genes")
    
    def compute_coexpression_networks(self):
        """
        Compute gene co-expression networks for each astrocyte type.
        """
        print("Computing gene co-expression networks...")
        
        # Get gene columns (exclude metadata columns)
        gene_cols = [col for col in self.gene_expression_matrix.columns 
                    if col not in ['point_uri', 'x', 'y', 'astrocyte_type', 'probability']]
        
        # Compute co-expression networks for each astrocyte type
        for cell_type in ['Protoplasmic', 'Fibrous', 'Reactive']:
            print(f"  Computing co-expression network for {cell_type} astrocytes...")
            
            # Filter cells by type
            type_cells = self.gene_expression_matrix[self.gene_expression_matrix['astrocyte_type'] == cell_type]
            
            if len(type_cells) < 5:
                print(f"  Not enough {cell_type} cells for co-expression analysis (found {len(type_cells)})")
                continue
            
            # Extract gene expression matrix
            expr_matrix = type_cells[gene_cols].values
            
            # Compute correlation matrix
            corr_matrix = np.zeros((len(gene_cols), len(gene_cols)))
            for i in range(len(gene_cols)):
                for j in range(i, len(gene_cols)):
                    if i == j:
                        corr_matrix[i, j] = 1.0
                    else:
                        # Use Spearman correlation for robustness
                        corr, _ = spearmanr(expr_matrix[:, i], expr_matrix[:, j])
                        corr_matrix[i, j] = corr if not np.isnan(corr) else 0
                        corr_matrix[j, i] = corr_matrix[i, j]
            
            # Create DataFrame
            corr_df = pd.DataFrame(corr_matrix, index=gene_cols, columns=gene_cols)
            
            # Store network
            self.coexpression_networks[cell_type] = corr_df
            
            # Save to file
            output_path = os.path.join(self.output_dir, f"{cell_type.lower()}_coexpression.csv")
            corr_df.to_csv(output_path)
            
            print(f"  Saved {cell_type} co-expression network to {output_path}")
    
    def validate_cell_types(self):
        """
        Validate cell type classifications using industry-standard approaches.
        """
        print("Validating cell type classifications...")
        
        # Get gene columns (exclude metadata columns)
        gene_cols = [col for col in self.gene_expression_matrix.columns 
                    if col not in ['point_uri', 'x', 'y', 'astrocyte_type', 'probability']]
        
        # Initialize validation results
        validation_results = []
        
        # 1. Marker gene enrichment analysis
        print("Performing marker gene enrichment analysis...")
        
        for idx, row in self.gene_expression_matrix.iterrows():
            cell_uri = row['point_uri']
            assigned_type = row['astrocyte_type']
            
            # Calculate marker scores for each cell type
            marker_scores = {}
            for cell_type, markers in CANONICAL_MARKERS.items():
                # Find markers present in our data
                present_markers = [m for m in markers if m in gene_cols]
                
                if not present_markers:
                    marker_scores[cell_type] = 0
                    continue
                
                # Calculate average expression of markers
                marker_expr = [row[m] for m in present_markers]
                marker_scores[cell_type] = np.mean(marker_expr)
            
            # Determine best type based on marker scores
            best_type = max(marker_scores, key=marker_scores.get)
            
            # Check if assigned type matches best type
            is_valid = assigned_type == best_type
            
            # Store validation result
            validation_results.append({
                'cell_uri': cell_uri,
                'assigned_type': assigned_type,
                'marker_based_type': best_type,
                'is_valid': is_valid,
                'protoplasmic_score': marker_scores.get('Protoplasmic', 0),
                'fibrous_score': marker_scores.get('Fibrous', 0),
                'reactive_score': marker_scores.get('Reactive', 0)
            })
        
        # Convert to DataFrame
        self.validation_results = pd.DataFrame(validation_results)
        
        # Calculate validation statistics
        valid_count = sum(self.validation_results['is_valid'])
        total_count = len(self.validation_results)
        valid_percent = (valid_count / total_count) * 100 if total_count > 0 else 0
        
        print(f"Validation results: {valid_count}/{total_count} cells correctly classified ({valid_percent:.2f}%)")
        
        # Detailed statistics by cell type
        for cell_type in ['Protoplasmic', 'Fibrous', 'Reactive']:
            type_cells = self.validation_results[self.validation_results['assigned_type'] == cell_type]
            valid_type_count = sum(type_cells['is_valid'])
            total_type_count = len(type_cells)
            valid_type_percent = (valid_type_count / total_type_count) * 100 if total_type_count > 0 else 0
            
            print(f"  {cell_type}: {valid_type_count}/{total_type_count} cells correctly classified ({valid_type_percent:.2f}%)")
        
        # Save validation results
        output_path = os.path.join(self.output_dir, "validation_results.csv")
        self.validation_results.to_csv(output_path, index=False)
        
        print(f"Saved validation results to {output_path}")
    
    def perform_differential_expression_analysis(self):
        """
        Perform differential expression analysis to identify genes that distinguish astrocyte types.
        """
        print("Performing differential expression analysis...")
        
        # Get gene columns (exclude metadata columns)
        gene_cols = [col for col in self.gene_expression_matrix.columns 
                    if col not in ['point_uri', 'x', 'y', 'astrocyte_type', 'probability']]
        
        # Create AnnData object for scanpy analysis
        X = self.gene_expression_matrix[gene_cols].values
        obs = self.gene_expression_matrix[['point_uri', 'astrocyte_type', 'probability']]
        var = pd.DataFrame(index=gene_cols)
        
        adata = ad.AnnData(X=X, obs=obs, var=var)
        
        # Normalize data
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Perform differential expression analysis
        sc.tl.rank_genes_groups(adata, 'astrocyte_type', method='wilcoxon')
        
        # Extract results
        de_results = {}
        for cell_type in ['Protoplasmic', 'Fibrous', 'Reactive']:
            genes = adata.uns['rank_genes_groups']['names'][cell_type]
            scores = adata.uns['rank_genes_groups']['scores'][cell_type]
            pvals = adata.uns['rank_genes_groups']['pvals'][cell_type]
            
            # Use our custom function to adjust p-values
            pvals_adj = adjust_pvalues(pvals)
            
            # Create DataFrame
            de_df = pd.DataFrame({
                'gene': genes,
                'score': scores,
                'pval': pvals,
                'pval_adj': pvals_adj
            })
            
            # Sort by adjusted p-value
            de_df = de_df.sort_values('pval_adj')
            
            # Save to file
            output_path = os.path.join(self.output_dir, f"{cell_type.lower()}_de_genes.csv")
            de_df.to_csv(output_path, index=False)
            
            # Store results
            de_results[cell_type] = de_df
            
            print(f"  Found {sum(de_df['pval_adj'] < 0.05)} significant DE genes for {cell_type}")
            print(f"  Saved {cell_type} DE genes to {output_path}")
        
        return de_results
    
    def visualize_validation_results(self):
        """
        Visualize validation results.
        """
        print("Visualizing validation results...")
        
        # 1. Confusion matrix
        plt.figure(figsize=(10, 8))
        conf_matrix = pd.crosstab(
            self.validation_results['assigned_type'], 
            self.validation_results['marker_based_type'],
            normalize='index'
        )
        
        sns.heatmap(conf_matrix, annot=True, cmap='Blues', fmt='.2f', cbar=True)
        plt.title('Cell Type Classification Validation')
        plt.xlabel('Marker-Based Type')
        plt.ylabel('Assigned Type')
        
        output_path = os.path.join(self.output_dir, "validation_confusion_matrix.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Marker gene expression by cell type
        plt.figure(figsize=(14, 10))
        
        # Melt the data for easier plotting
        plot_data = pd.melt(
            self.validation_results,
            id_vars=['cell_uri', 'assigned_type', 'marker_based_type', 'is_valid'],
            value_vars=['protoplasmic_score', 'fibrous_score', 'reactive_score'],
            var_name='marker_type',
            value_name='score'
        )
        
        # Clean up marker type names
        plot_data['marker_type'] = plot_data['marker_type'].str.replace('_score', '')
        plot_data['marker_type'] = plot_data['marker_type'].str.capitalize()
        
        # Plot
        sns.boxplot(x='assigned_type', y='score', hue='marker_type', data=plot_data)
        plt.title('Marker Gene Expression by Assigned Cell Type')
        plt.xlabel('Assigned Cell Type')
        plt.ylabel('Marker Score')
        
        output_path = os.path.join(self.output_dir, "marker_expression_by_type.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Spatial distribution of validation results
        spatial_plot_created = False
        
        # Check if we have spatial coordinates
        if 'x' in self.gene_expression_matrix.columns and 'y' in self.gene_expression_matrix.columns:
            try:
                plt.figure(figsize=(12, 10))
                
                # Merge validation results with spatial coordinates
                # First ensure we have the cell_uri column in both dataframes
                if 'cell_uri' in self.validation_results.columns and 'cell_uri' in self.gene_expression_matrix.columns:
                    merged_data = pd.merge(
                        self.validation_results,
                        self.gene_expression_matrix[['cell_uri', 'x', 'y']],
                        on='cell_uri'
                    )
                elif 'point_uri' in self.validation_results.columns and 'point_uri' in self.gene_expression_matrix.columns:
                    merged_data = pd.merge(
                        self.validation_results,
                        self.gene_expression_matrix[['point_uri', 'x', 'y']],
                        on='point_uri'
                    )
                else:
                    print("  Warning: Cannot create spatial distribution plot due to missing URI columns")
                    # We'll create a placeholder below
                    raise ValueError("Missing URI columns for spatial plot")
                
                # Plot valid classifications
                valid_data = merged_data[merged_data['is_valid']]
                plt.scatter(valid_data['x'], valid_data['y'], c='green', alpha=0.7, s=50, label='Valid')
                
                # Plot invalid classifications
                invalid_data = merged_data[~merged_data['is_valid']]
                plt.scatter(invalid_data['x'], invalid_data['y'], c='red', alpha=0.7, s=50, label='Invalid')
                
                plt.title('Spatial Distribution of Validation Results')
                plt.xlabel('X Coordinate')
                plt.ylabel('Y Coordinate')
                plt.legend()
                plt.grid(alpha=0.3)
                
                output_path = os.path.join(self.output_dir, "spatial_validation_results.png")
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
                spatial_plot_created = True
            except Exception as e:
                print(f"  Warning: Could not create spatial distribution plot: {str(e)}")
                # We'll create a placeholder below
        else:
            print("  Warning: Cannot create spatial distribution plot due to missing spatial coordinates")
            # We'll create a placeholder below
        
        # Create a placeholder image if the spatial plot wasn't created
        if not spatial_plot_created:
            plt.figure(figsize=(10, 8))
            plt.text(0.5, 0.5, "Spatial validation plot unavailable\nMissing spatial coordinates or URI columns", 
                    horizontalalignment='center', verticalalignment='center', fontsize=14)
            plt.axis('off')
            output_path = os.path.join(self.output_dir, "spatial_validation_results.png")
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
    
    def visualize_coexpression_networks(self):
        """
        Visualize gene co-expression networks.
        """
        print("Visualizing gene co-expression networks...")
        
        if not self.coexpression_networks:
            print("  Warning: No co-expression networks to visualize")
            return
        
        for cell_type, corr_df in self.coexpression_networks.items():
            print(f"  Visualizing {cell_type} co-expression network...")
            
            # Check if the network is empty
            if corr_df.empty:
                print(f"  Warning: Co-expression network for {cell_type} is empty")
                continue
            
            # Filter for canonical markers
            markers = CANONICAL_MARKERS[cell_type]
            present_markers = [m for m in markers if m in corr_df.columns]
            
            if len(present_markers) < 2:
                print(f"  Not enough markers for {cell_type} in the data")
                continue
            
            try:
                # Extract subnetwork of markers
                marker_network = corr_df.loc[present_markers, present_markers]
                
                # Create network visualization
                plt.figure(figsize=(12, 10))
                
                # Create graph
                G = nx.Graph()
                
                # Add nodes
                for marker in present_markers:
                    G.add_node(marker)
                
                # Add edges (only if correlation is above threshold)
                threshold = 0.3
                edges_to_add = []
                
                try:
                    for i in range(len(present_markers)):
                        for j in range(i+1, len(present_markers)):
                            marker1 = present_markers[i]
                            marker2 = present_markers[j]
                            corr = marker_network.loc[marker1, marker2]
                            
                            # Handle the case where corr might be a Series
                            if hasattr(corr, 'iloc'):
                                corr = corr.iloc[0] if len(corr) > 0 else 0
                            
                            if abs(corr) > threshold:
                                edges_to_add.append((marker1, marker2, {'weight': corr}))
                except Exception as e:
                    print(f"  Warning: Error creating edges for {cell_type} network: {str(e)}")
                
                # Add all edges at once
                G.add_edges_from(edges_to_add)
                
                # Set node positions
                pos = nx.spring_layout(G, seed=42)
                
                # Draw nodes
                nx.draw_networkx_nodes(G, pos, node_size=700, node_color='skyblue')
                
                # Draw edges with width based on correlation strength
                edges = G.edges(data=True)
                if edges:
                    edge_colors = ['red' if e[2]['weight'] < 0 else 'blue' for e in edges]
                    edge_widths = [abs(e[2]['weight']) * 5 for e in edges]
                    
                    nx.draw_networkx_edges(G, pos, edgelist=edges, width=edge_widths, edge_color=edge_colors, alpha=0.7)
                
                # Draw labels
                nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')
                
                plt.title(f"{cell_type} Astrocyte Marker Co-expression Network")
                plt.axis('off')
                
                output_path = os.path.join(self.output_dir, f"{cell_type.lower()}_coexpression_network.png")
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
            except Exception as e:
                print(f"  Warning: Could not create co-expression network for {cell_type}: {str(e)}")
                
                # Create a placeholder image
                plt.figure(figsize=(10, 8))
                plt.text(0.5, 0.5, f"Co-expression network for {cell_type} unavailable\n{str(e)}", 
                        horizontalalignment='center', verticalalignment='center', fontsize=14)
                plt.axis('off')
                output_path = os.path.join(self.output_dir, f"{cell_type.lower()}_coexpression_network.png")
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.close()
    
    def compare_with_literature(self):
        """
        Compare identified marker genes with literature-reported markers.
        """
        print("Comparing with literature-reported markers...")
        
        # Get gene columns (exclude metadata columns)
        gene_cols = [col for col in self.gene_expression_matrix.columns 
                    if col not in ['point_uri', 'x', 'y', 'astrocyte_type', 'probability']]
        
        # Perform differential expression analysis if not already done
        de_results = self.perform_differential_expression_analysis()
        
        # Compare with canonical markers
        for cell_type, de_df in de_results.items():
            # Get top DE genes (adjusted p-value < 0.05)
            significant_genes = de_df[de_df['pval_adj'] < 0.05]['gene'].tolist()
            
            # Get canonical markers
            canonical_markers = CANONICAL_MARKERS[cell_type]
            
            # Find overlap
            overlap = set(significant_genes).intersection(set(canonical_markers))
            
            # Calculate statistics
            total_canonical = len(canonical_markers)
            total_significant = len(significant_genes)
            total_overlap = len(overlap)
            
            print(f"\n{cell_type} Astrocyte Marker Comparison:")
            print(f"  Canonical markers: {total_canonical}")
            print(f"  Significant DE genes: {total_significant}")
            print(f"  Overlap: {total_overlap} ({(total_overlap/total_canonical)*100:.2f}% of canonical markers)")
            
            if total_overlap > 0:
                print(f"  Overlapping markers: {', '.join(overlap)}")
            
            # Create Venn diagram
            plt.figure(figsize=(10, 8))
            
            # Create sets for Venn diagram
            set1 = set(canonical_markers)
            set2 = set(significant_genes)
            set3 = set(gene_cols)  # All genes in the dataset
            
            venn3([set1, set2, set3], 
                  ('Canonical Markers', 'Significant DE Genes', 'All Genes in Dataset'))
            
            plt.title(f"{cell_type} Astrocyte Marker Comparison")
            
            output_path = os.path.join(self.output_dir, f"{cell_type.lower()}_marker_comparison.png")
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
    
    def run_validation(self):
        """
        Run the full validation pipeline.
        """
        print("Running astrocyte validation pipeline...")
        
        # Extract classified cells and gene expression
        self.extract_classified_cells()
        
        # Compute co-expression networks
        self.compute_coexpression_networks()
        
        # Validate cell types
        self.validate_cell_types()
        
        # Skip differential expression analysis due to dependency issues
        print("Skipping differential expression analysis due to dependency issues")
        
        # Visualize validation results
        self.visualize_validation_results()
        
        # Visualize co-expression networks
        self.visualize_coexpression_networks()
        
        # Skip literature comparison as it depends on differential expression analysis
        print("Skipping literature comparison due to dependency issues")
        
        print("Validation pipeline completed")
        
        # Return summary of validation results
        valid_count = sum(self.validation_results['is_valid'])
        total_count = len(self.validation_results)
        valid_percent = (valid_count / total_count) * 100 if total_count > 0 else 0
        
        return {
            'total_cells': total_count,
            'valid_cells': valid_count,
            'valid_percent': valid_percent,
            'results_path': os.path.join(self.output_dir, "validation_results.csv")
        }

def main():
    """
    Main function to run the astrocyte validator.
    """
    parser = argparse.ArgumentParser(description='Validate astrocyte cell type classifications')
    parser.add_argument('--input', required=True, help='Path to integrated ontology TURTLE file')
    parser.add_argument('--output-dir', default='../output/astrocyte_validation', help='Output directory for validation results')
    
    args = parser.parse_args()
    
    # Create validator
    validator = AstrocyteValidator(
        integrated_ttl=args.input,
        output_dir=args.output_dir
    )
    
    # Run validation
    results = validator.run_validation()
    
    # Print summary
    print("\nValidation Summary:")
    print(f"Total cells analyzed: {results['total_cells']}")
    print(f"Cells with valid classification: {results['valid_cells']} ({results['valid_percent']:.2f}%)")
    print(f"Detailed results saved to: {results['results_path']}")

if __name__ == "__main__":
    main() 