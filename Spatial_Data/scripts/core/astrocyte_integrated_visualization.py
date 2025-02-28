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
import networkx as nx
from collections import defaultdict
import argparse
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import scanpy as sc
import anndata
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import base64
import io
from PIL import Image

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
SC = Namespace("http://example.org/single-cell#")
CELL = Namespace("http://example.org/cell-ontology#")

# Define astrocyte marker genes for different types
ASTROCYTE_MARKERS = {
    'Protoplasmic': ['GFAP', 'S100B', 'AQP4', 'GJA1', 'SLC1A2', 'SLC1A3'],
    'Fibrous': ['GFAP', 'VIM', 'ALDH1L1', 'CD44', 'CRYAB', 'HOPX'],
    'Reactive': ['GFAP', 'VIM', 'SERPINA3', 'C3', 'CXCL10', 'STAT3', 'LCN2']
}

class AstrocyteIntegratedVisualization:
    """
    A class for integrating and visualizing spatial and single-cell data
    for astrocyte analysis.
    """
    
    def __init__(self, spatial_ttl, single_cell_ttl, output_dir="../output/integrated_visualization"):
        """
        Initialize the visualization tool with input and output file paths.
        
        Parameters:
        -----------
        spatial_ttl : str
            Path to the TURTLE file containing the spatial ontology
        single_cell_ttl : str
            Path to the TURTLE file containing the single-cell ontology
        output_dir : str
            Path to the directory where visualization results will be saved
        """
        self.spatial_ttl = spatial_ttl
        self.single_cell_ttl = single_cell_ttl
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize data structures
        self.spatial_graph = None
        self.single_cell_graph = None
        self.spatial_data = None
        self.single_cell_data = None
        self.integrated_data = None
        self.astrocyte_types = None
        self.gene_expression_patterns = None
        
        # Load data
        self.load_data()
    
    def load_data(self):
        """
        Load the spatial and single-cell ontology data from TURTLE files.
        """
        print(f"Loading spatial TURTLE file: {self.spatial_ttl}")
        self.spatial_graph = Graph()
        self.spatial_graph.parse(self.spatial_ttl, format="turtle")
        print(f"Spatial graph loaded with {len(self.spatial_graph)} triples")
        
        print(f"Loading single-cell TURTLE file: {self.single_cell_ttl}")
        self.single_cell_graph = Graph()
        self.single_cell_graph.parse(self.single_cell_ttl, format="turtle")
        print(f"Single-cell graph loaded with {len(self.single_cell_graph)} triples")
        
        # Extract data from ontologies
        self.extract_spatial_data()
        self.extract_single_cell_data()
        
    def extract_spatial_data(self):
        """
        Extract spatial data from the spatial ontology.
        """
        print("Extracting spatial data...")
        
        # Query for astrocyte types and their spatial locations
        query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX cell: <http://example.org/cell-ontology#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y ?astrocyteType ?probability
        WHERE {
            ?point rdf:type st:SpatialPoint ;
                   st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?astrocyteType ;
                      st:typeProbability ?probability .
        }
        """
        
        results = self.spatial_graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'point_id': str(row.point).split('#')[-1],
                'x': float(row.x),
                'y': float(row.y),
                'astrocyte_type': str(row.astrocyteType),
                'probability': float(row.probability)
            })
        
        self.spatial_data = pd.DataFrame(data)
        print(f"Extracted {len(self.spatial_data)} spatial points with astrocyte type information")
    
    def extract_single_cell_data(self):
        """
        Extract single-cell data from the single-cell ontology.
        """
        print("Extracting single-cell data...")
        
        # Query for cells, their types, and gene expression
        query = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?cell ?cellType ?gene ?expressionLevel
        WHERE {
            ?cell rdf:type sc:Cell ;
                  sc:hasCellType ?cellType .
            
            ?expr sc:belongsToCell ?cell ;
                  sc:forGene ?gene ;
                  sc:hasExpressionLevel ?expressionLevel .
        }
        """
        
        results = self.single_cell_graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'cell_id': str(row.cell).split('#')[-1],
                'cell_type': str(row.cellType),
                'gene': str(row.gene).split('#')[-1],
                'expression_level': float(row.expressionLevel)
            })
        
        self.single_cell_data = pd.DataFrame(data)
        print(f"Extracted {len(self.single_cell_data)} gene expression records from single-cell data")
    
    def integrate_data(self):
        """
        Integrate spatial and single-cell data for comprehensive analysis.
        """
        print("Integrating spatial and single-cell data...")
        
        if self.spatial_data is None or self.single_cell_data is None:
            print("Error: Data not loaded properly")
            return
        
        # Get average gene expression for each astrocyte type from single-cell data
        cell_type_gene_expr = self.single_cell_data.pivot_table(
            index='cell_type', 
            columns='gene', 
            values='expression_level', 
            aggfunc='mean'
        ).fillna(0)
        
        # Map astrocyte types between datasets
        type_mapping = {
            'Protoplasmic': 'Protoplasmic',
            'Fibrous': 'Fibrous',
            'Reactive': 'Reactive'
        }
        
        # Create integrated dataset
        integrated_data = []
        
        for _, point in self.spatial_data.iterrows():
            point_data = {
                'point_id': point['point_id'],
                'x': point['x'],
                'y': point['y'],
                'astrocyte_type': point['astrocyte_type'],
                'probability': point['probability']
            }
            
            # Add gene expression data from single-cell for this astrocyte type
            if point['astrocyte_type'] in type_mapping:
                sc_type = type_mapping[point['astrocyte_type']]
                if sc_type in cell_type_gene_expr.index:
                    for gene, expr in cell_type_gene_expr.loc[sc_type].items():
                        point_data[f'gene_{gene}'] = expr
            
            integrated_data.append(point_data)
        
        self.integrated_data = pd.DataFrame(integrated_data)
        print(f"Created integrated dataset with {len(self.integrated_data)} points")
    
    def analyze_gene_expression_patterns(self):
        """
        Analyze gene expression patterns for each astrocyte type.
        """
        print("Analyzing gene expression patterns...")
        
        if self.single_cell_data is None:
            print("Error: Single-cell data not loaded")
            return
        
        # Pivot data to get gene expression by cell type
        gene_expr_by_type = self.single_cell_data.pivot_table(
            index='cell_type', 
            columns='gene', 
            values='expression_level', 
            aggfunc='mean'
        ).fillna(0)
        
        # Get top expressed genes for each astrocyte type
        top_genes = {}
        for cell_type in gene_expr_by_type.index:
            # Sort genes by expression level
            sorted_genes = gene_expr_by_type.loc[cell_type].sort_values(ascending=False)
            # Get top 20 genes
            top_genes[cell_type] = sorted_genes.head(20)
        
        self.gene_expression_patterns = {
            'by_type': gene_expr_by_type,
            'top_genes': top_genes
        }
        
        print("Gene expression patterns analyzed")
    
    def visualize_spatial_distribution(self):
        """
        Visualize the spatial distribution of astrocyte types.
        """
        print("Visualizing spatial distribution of astrocyte types...")
        
        if self.spatial_data is None:
            print("Error: Spatial data not loaded")
            return
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Define colors for each astrocyte type
        colors = {
            'Protoplasmic': 'blue',
            'Fibrous': 'red',
            'Reactive': 'green'
        }
        
        # Plot each astrocyte type
        for cell_type, color in colors.items():
            # Filter data for this cell type
            type_data = self.spatial_data[self.spatial_data['astrocyte_type'] == cell_type]
            
            if not type_data.empty:
                # Plot points
                scatter = ax.scatter(
                    type_data['x'], 
                    type_data['y'],
                    c=color,
                    s=type_data['probability'] * 100,  # Size based on probability
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
        
        # Save the figure
        fig_file = os.path.join(self.output_dir, "astrocyte_spatial_distribution.png")
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        print(f"Visualization saved to {fig_file}")
        
        return fig_file
    
    def visualize_gene_expression_heatmap(self):
        """
        Create a heatmap of gene expression for different astrocyte types.
        """
        print("Creating gene expression heatmap...")
        
        if self.gene_expression_patterns is None:
            self.analyze_gene_expression_patterns()
        
        gene_expr = self.gene_expression_patterns['by_type']
        
        # Select marker genes for visualization
        all_markers = []
        for markers in ASTROCYTE_MARKERS.values():
            all_markers.extend(markers)
        all_markers = list(set(all_markers))  # Remove duplicates
        
        # Filter for marker genes that exist in our data
        existing_markers = [gene for gene in all_markers if gene in gene_expr.columns]
        
        if not existing_markers:
            print("No marker genes found in the data")
            return
        
        # Create heatmap
        plt.figure(figsize=(14, 8))
        sns.set(font_scale=1.2)
        
        # Plot heatmap
        heatmap = sns.heatmap(
            gene_expr[existing_markers],
            cmap="YlGnBu",
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            cbar_kws={"label": "Expression Level"}
        )
        
        plt.title("Gene Expression Levels by Astrocyte Type", fontsize=16)
        plt.ylabel("Astrocyte Type", fontsize=14)
        plt.xlabel("Gene", fontsize=14)
        plt.xticks(rotation=45, ha='right')
        
        # Save the figure
        fig_file = os.path.join(self.output_dir, "gene_expression_heatmap.png")
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Heatmap saved to {fig_file}")
        
        return fig_file
    
    def visualize_marker_expression_in_space(self, marker_gene):
        """
        Visualize the expression of a specific marker gene in spatial context.
        
        Parameters:
        -----------
        marker_gene : str
            The marker gene to visualize
        """
        print(f"Visualizing spatial expression of {marker_gene}...")
        
        if self.integrated_data is None:
            self.integrate_data()
        
        gene_col = f'gene_{marker_gene}'
        
        if gene_col not in self.integrated_data.columns:
            print(f"Gene {marker_gene} not found in integrated data")
            return
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Create scatter plot with expression levels as colors
        scatter = ax.scatter(
            self.integrated_data['x'],
            self.integrated_data['y'],
            c=self.integrated_data[gene_col],
            cmap='viridis',
            s=50,
            alpha=0.8
        )
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(f'{marker_gene} Expression Level', fontsize=12)
        
        # Add labels and title
        ax.set_xlabel('X Coordinate', fontsize=14)
        ax.set_ylabel('Y Coordinate', fontsize=14)
        ax.set_title(f'Spatial Distribution of {marker_gene} Expression', fontsize=16)
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Save the figure
        fig_file = os.path.join(self.output_dir, f"spatial_{marker_gene}_expression.png")
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        print(f"Visualization saved to {fig_file}")
        
        return fig_file
    
    def create_interactive_dashboard(self):
        """
        Create an interactive dashboard for exploring astrocyte data.
        """
        print("Creating interactive dashboard...")
        
        if self.integrated_data is None:
            self.integrate_data()
        
        if self.gene_expression_patterns is None:
            self.analyze_gene_expression_patterns()
        
        # Initialize Dash app
        app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
        
        # Define layout
        app.layout = dbc.Container([
            dbc.Row([
                dbc.Col([
                    html.H1("Astrocyte Integrated Analysis Dashboard", 
                            className="text-center my-4")
                ], width=12)
            ]),
            
            dbc.Row([
                dbc.Col([
                    html.H3("Visualization Controls", className="mb-3"),
                    html.Label("Select Visualization Type:"),
                    dcc.Dropdown(
                        id='viz-type',
                        options=[
                            {'label': 'Spatial Distribution of Astrocyte Types', 'value': 'spatial_types'},
                            {'label': 'Gene Expression Heatmap', 'value': 'gene_heatmap'},
                            {'label': 'Marker Gene Spatial Expression', 'value': 'marker_spatial'}
                        ],
                        value='spatial_types'
                    ),
                    
                    html.Div(id='marker-selection', className="mt-3", children=[
                        html.Label("Select Marker Gene:"),
                        dcc.Dropdown(
                            id='marker-gene',
                            options=[{'label': gene, 'value': gene} for gene in 
                                    sum(ASTROCYTE_MARKERS.values(), [])],
                            value='GFAP'
                        )
                    ], style={'display': 'none'})
                ], width=3),
                
                dbc.Col([
                    html.Div(id='visualization-output', className="border p-3")
                ], width=9)
            ]),
            
            dbc.Row([
                dbc.Col([
                    html.H3("Statistics", className="my-3"),
                    html.Div(id='statistics-output')
                ], width=12)
            ])
        ], fluid=True)
        
        # Define callbacks
        @app.callback(
            [Output('marker-selection', 'style'),
             Output('visualization-output', 'children'),
             Output('statistics-output', 'children')],
            [Input('viz-type', 'value'),
             Input('marker-gene', 'value')]
        )
        def update_visualization(viz_type, marker_gene):
            # Show/hide marker selection based on visualization type
            marker_style = {'display': 'block'} if viz_type == 'marker_spatial' else {'display': 'none'}
            
            # Generate visualization
            if viz_type == 'spatial_types':
                # Generate spatial distribution visualization
                fig_file = self.visualize_spatial_distribution()
                img_element = html.Img(src=app.get_asset_url(os.path.basename(fig_file)), 
                                      style={'width': '100%'})
                
                # Statistics
                type_counts = self.spatial_data['astrocyte_type'].value_counts()
                stats = html.Div([
                    html.H4("Astrocyte Type Distribution"),
                    dbc.Table.from_dataframe(
                        pd.DataFrame({
                            'Astrocyte Type': type_counts.index,
                            'Count': type_counts.values,
                            'Percentage': (type_counts.values / type_counts.sum() * 100).round(2)
                        }),
                        striped=True, bordered=True, hover=True
                    )
                ])
                
            elif viz_type == 'gene_heatmap':
                # Generate gene expression heatmap
                fig_file = self.visualize_gene_expression_heatmap()
                img_element = html.Img(src=app.get_asset_url(os.path.basename(fig_file)), 
                                      style={'width': '100%'})
                
                # Statistics
                top_genes = self.gene_expression_patterns['top_genes']
                stats_children = [html.H4("Top Expressed Genes by Astrocyte Type")]
                
                for cell_type, genes in top_genes.items():
                    stats_children.append(html.H5(f"{cell_type} Astrocytes"))
                    stats_children.append(
                        dbc.Table.from_dataframe(
                            pd.DataFrame({
                                'Gene': genes.index[:10],  # Top 10 genes
                                'Expression Level': genes.values[:10].round(2)
                            }),
                            striped=True, bordered=True, hover=True, size='sm'
                        )
                    )
                
                stats = html.Div(stats_children)
                
            elif viz_type == 'marker_spatial':
                # Generate marker gene spatial expression
                fig_file = self.visualize_marker_expression_in_space(marker_gene)
                img_element = html.Img(src=app.get_asset_url(os.path.basename(fig_file)), 
                                      style={'width': '100%'})
                
                # Statistics
                gene_col = f'gene_{marker_gene}'
                if gene_col in self.integrated_data.columns:
                    expr_by_type = self.integrated_data.groupby('astrocyte_type')[gene_col].mean()
                    stats = html.Div([
                        html.H4(f"{marker_gene} Expression Statistics"),
                        dbc.Table.from_dataframe(
                            pd.DataFrame({
                                'Astrocyte Type': expr_by_type.index,
                                'Mean Expression': expr_by_type.values.round(2)
                            }),
                            striped=True, bordered=True, hover=True
                        ),
                        html.P(f"Overall mean expression: {self.integrated_data[gene_col].mean().round(2)}")
                    ])
                else:
                    stats = html.Div([
                        html.H4(f"{marker_gene} Expression Statistics"),
                        html.P("Expression data not available for this gene")
                    ])
            
            return marker_style, img_element, stats
        
        # Run the app
        app.run_server(debug=True, port=8050)
        
        print("Interactive dashboard running at http://127.0.0.1:8050/")
    
    def run_sparql_queries(self):
        """
        Run SPARQL queries to analyze gene expression patterns within each astrocyte type.
        """
        print("Running SPARQL queries for gene expression analysis...")
        
        # Query 1: Find genes that are highly expressed in specific astrocyte types
        query1 = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?astrocyteType ?gene (AVG(?expressionLevel) as ?avgExpression)
        WHERE {
            # Get spatial points with astrocyte type
            ?point rdf:type st:SpatialPoint ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?astrocyteType ;
                      st:typeProbability ?probability .
            
            # Only consider high probability assignments
            FILTER(?probability > 0.7)
            
            # Get gene expression from single-cell data for this type
            ?cell sc:hasCellType ?astrocyteType .
            ?expr sc:belongsToCell ?cell ;
                  sc:forGene ?gene ;
                  sc:hasExpressionLevel ?expressionLevel .
        }
        GROUP BY ?astrocyteType ?gene
        HAVING (AVG(?expressionLevel) > 1.0)
        ORDER BY ?astrocyteType DESC(?avgExpression)
        LIMIT 100
        """
        
        # Execute query on the merged graph
        merged_graph = Graph()
        merged_graph += self.spatial_graph
        merged_graph += self.single_cell_graph
        
        results1 = merged_graph.query(query1)
        
        # Convert to DataFrame
        data1 = []
        for row in results1:
            data1.append({
                'astrocyte_type': str(row.astrocyteType),
                'gene': str(row.gene).split('#')[-1],
                'avg_expression': float(row.avgExpression)
            })
        
        df1 = pd.DataFrame(data1)
        
        # Save results
        results_file = os.path.join(self.output_dir, "high_expression_genes_by_type.csv")
        df1.to_csv(results_file, index=False)
        print(f"Query results saved to {results_file}")
        
        # Query 2: Find co-expressed genes in specific astrocyte types
        query2 = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?astrocyteType ?gene1 ?gene2 (COUNT(?cell) as ?coexpressionCount)
        WHERE {
            # Get cells of specific astrocyte type
            ?cell sc:hasCellType ?astrocyteType .
            
            # Get expression for gene1
            ?expr1 sc:belongsToCell ?cell ;
                   sc:forGene ?gene1 ;
                   sc:hasExpressionLevel ?expr1Level .
            
            # Get expression for gene2
            ?expr2 sc:belongsToCell ?cell ;
                   sc:forGene ?gene2 ;
                   sc:hasExpressionLevel ?expr2Level .
            
            # Only consider high expression
            FILTER(?expr1Level > 1.0 && ?expr2Level > 1.0)
            
            # Avoid duplicate pairs
            FILTER(STR(?gene1) < STR(?gene2))
        }
        GROUP BY ?astrocyteType ?gene1 ?gene2
        HAVING (COUNT(?cell) > 5)
        ORDER BY ?astrocyteType DESC(?coexpressionCount)
        LIMIT 100
        """
        
        results2 = merged_graph.query(query2)
        
        # Convert to DataFrame
        data2 = []
        for row in results2:
            data2.append({
                'astrocyte_type': str(row.astrocyteType),
                'gene1': str(row.gene1).split('#')[-1],
                'gene2': str(row.gene2).split('#')[-1],
                'coexpression_count': int(row.coexpressionCount)
            })
        
        df2 = pd.DataFrame(data2)
        
        # Save results
        results_file = os.path.join(self.output_dir, "coexpressed_genes_by_type.csv")
        df2.to_csv(results_file, index=False)
        print(f"Query results saved to {results_file}")
        
        return df1, df2
    
    def run_pipeline(self):
        """
        Run the full visualization pipeline.
        """
        print("Running integrated visualization pipeline...")
        
        # Load and integrate data
        self.integrate_data()
        
        # Analyze gene expression patterns
        self.analyze_gene_expression_patterns()
        
        # Create visualizations
        self.visualize_spatial_distribution()
        self.visualize_gene_expression_heatmap()
        
        # Visualize marker genes
        for cell_type, markers in ASTROCYTE_MARKERS.items():
            for marker in markers[:2]:  # Visualize first two markers of each type
                self.visualize_marker_expression_in_space(marker)
        
        # Run SPARQL queries
        self.run_sparql_queries()
        
        print("Visualization pipeline completed")

def main():
    """
    Main function to run the integrated visualization tool.
    """
    parser = argparse.ArgumentParser(description='Integrated visualization of astrocyte types')
    parser.add_argument('--spatial', required=True, help='Path to spatial ontology TURTLE file')
    parser.add_argument('--single-cell', required=True, help='Path to single-cell ontology TURTLE file')
    parser.add_argument('--output', default='../output/integrated_visualization', help='Output directory')
    parser.add_argument('--interactive', action='store_true', help='Launch interactive dashboard')
    
    args = parser.parse_args()
    
    # Create visualization tool
    viz_tool = AstrocyteIntegratedVisualization(
        spatial_ttl=args.spatial,
        single_cell_ttl=args.single_cell,
        output_dir=args.output
    )
    
    # Run pipeline
    viz_tool.run_pipeline()
    
    # Launch interactive dashboard if requested
    if args.interactive:
        viz_tool.create_interactive_dashboard()

if __name__ == "__main__":
    main() 