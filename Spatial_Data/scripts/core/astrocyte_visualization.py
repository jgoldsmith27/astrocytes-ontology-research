#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, XSD
import argparse
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
SC = Namespace("http://example.org/single-cell#")
BRIDGE = Namespace("http://example.org/ontology-bridge#")

class AstrocyteVisualizer:
    """
    A class for visualizing astrocyte data from the integrated ontology.
    """
    
    def __init__(self, integrated_ttl, output_dir="../output/visualizations"):
        """
        Initialize the visualizer with input and output paths.
        
        Parameters:
        -----------
        integrated_ttl : str
            Path to the integrated ontology TURTLE file
        output_dir : str
            Directory where visualization outputs will be saved
        """
        self.integrated_ttl = integrated_ttl
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize graph
        self.graph = None
        
        # Load data
        self.load_data()
        
        # Set up color maps for astrocyte types
        self.astrocyte_colors = {
            "Protoplasmic": "#4287f5",  # Blue
            "Fibrous": "#f54242",       # Red
            "Reactive": "#42f54e"       # Green
        }
        
        # Set up color maps for brain regions
        self.region_colors = {
            "PrefrontalCortex": "#e6194B",
            "MotorCortex": "#3cb44b",
            "Hippocampus": "#ffe119",
            "Cerebellum": "#4363d8",
            "Thalamus": "#f58231"
        }
    
    def load_data(self):
        """
        Load the integrated ontology data from TURTLE file.
        """
        print(f"Loading integrated TURTLE file: {self.integrated_ttl}")
        self.graph = Graph()
        self.graph.parse(self.integrated_ttl, format="turtle")
        print(f"Integrated graph loaded with {len(self.graph)} triples")
    
    def get_spatial_points_data(self):
        """
        Extract spatial points data from the ontology.
        
        Returns:
        --------
        pandas.DataFrame
            DataFrame containing spatial points data
        """
        # Query to get spatial points with astrocyte types and coordinates
        query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y ?astrocyteType ?probability ?region
        WHERE {
            ?point rdf:type st:SpatialPoint ;
                   st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?astrocyteType ;
                      st:typeProbability ?probability .
            
            # Only include astrocyte types
            FILTER(?astrocyteType IN ("Protoplasmic", "Fibrous", "Reactive"))
            
            # Get brain region if available
            OPTIONAL {
                ?point st:inBrainRegion ?regionNode .
                ?regionNode rdfs:label ?region .
            }
        }
        """
        
        results = list(self.graph.query(query))
        
        # Convert to DataFrame
        data = []
        for row in results:
            point_uri = str(row.point)
            x = float(row.x)
            y = float(row.y)
            astrocyte_type = str(row.astrocyteType)
            probability = float(row.probability)
            region = str(row.region) if row.region else "Unknown"
            
            data.append({
                'point_uri': point_uri,
                'x': x,
                'y': y,
                'astrocyte_type': astrocyte_type,
                'probability': probability,
                'region': region
            })
        
        return pd.DataFrame(data)
    
    def get_gene_expression_data(self, gene_id=None):
        """
        Extract gene expression data from the ontology.
        
        Parameters:
        -----------
        gene_id : str, optional
            Specific gene ID to extract data for. If None, extract data for all genes.
        
        Returns:
        --------
        pandas.DataFrame
            DataFrame containing gene expression data
        """
        # Query to get gene expression data
        query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y ?geneID ?expressionLevel ?astrocyteType
        WHERE {
            ?point rdf:type st:SpatialPoint ;
                   st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?astrocyteType .
            
            # Only include astrocyte types
            FILTER(?astrocyteType IN ("Protoplasmic", "Fibrous", "Reactive"))
            
            # Get gene expression
            ?gene st:expressedAt ?point ;
                  st:hasExpressionLevel ?expressionLevel ;
                  st:hasGeneID ?geneID .
        """
        
        if gene_id:
            query += f"""
            # Filter for specific gene
            FILTER(?geneID = "{gene_id}")
            """
        
        query += """
        }
        """
        
        results = list(self.graph.query(query))
        
        # Convert to DataFrame
        data = []
        for row in results:
            point_uri = str(row.point)
            x = float(row.x)
            y = float(row.y)
            gene_id = str(row.geneID)
            expression_level = float(row.expressionLevel)
            astrocyte_type = str(row.astrocyteType)
            
            data.append({
                'point_uri': point_uri,
                'x': x,
                'y': y,
                'gene_id': gene_id,
                'expression_level': expression_level,
                'astrocyte_type': astrocyte_type
            })
        
        return pd.DataFrame(data)
    
    def get_single_cell_data(self):
        """
        Extract single-cell data from the ontology.
        
        Returns:
        --------
        pandas.DataFrame
            DataFrame containing single-cell data
        """
        # Query to get single-cell data with spatial mappings
        query = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?cell ?cellType ?x ?y ?region
        WHERE {
            ?cell rdf:type sc:Cell ;
                  sc:hasCellType ?cellType ;
                  bridge:hasSpatialRepresentation ?point .
            
            # Only include astrocyte cell types
            FILTER(?cellType IN ("Protoplasmic", "Fibrous", "Reactive"))
            
            # Get spatial coordinates
            ?point st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y .
            
            # Get brain region if available
            OPTIONAL {
                ?point st:inBrainRegion ?regionNode .
                ?regionNode rdfs:label ?region .
            }
        }
        """
        
        results = list(self.graph.query(query))
        
        # Convert to DataFrame
        data = []
        for row in results:
            cell_uri = str(row.cell)
            cell_type = str(row.cellType)
            x = float(row.x)
            y = float(row.y)
            region = str(row.region) if row.region else "Unknown"
            
            data.append({
                'cell_uri': cell_uri,
                'cell_type': cell_type,
                'x': x,
                'y': y,
                'region': region
            })
        
        return pd.DataFrame(data)
    
    def visualize_astrocyte_distribution(self):
        """
        Visualize the spatial distribution of astrocyte types.
        """
        print("Visualizing astrocyte distribution...")
        
        # Get spatial points data
        df = self.get_spatial_points_data()
        
        if df.empty:
            print("No astrocyte data found for visualization")
            return
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Plot points colored by astrocyte type
        for astrocyte_type, color in self.astrocyte_colors.items():
            subset = df[df['astrocyte_type'] == astrocyte_type]
            if not subset.empty:
                plt.scatter(subset['x'], subset['y'], c=color, alpha=0.7, 
                           s=subset['probability'] * 100, label=astrocyte_type)
        
        # Add legend
        plt.legend(title="Astrocyte Type", fontsize=12)
        
        # Add labels and title
        plt.xlabel("X Coordinate", fontsize=14)
        plt.ylabel("Y Coordinate", fontsize=14)
        plt.title("Spatial Distribution of Astrocyte Types", fontsize=16)
        
        # Add grid
        plt.grid(alpha=0.3)
        
        # Save figure
        output_path = os.path.join(self.output_dir, "astrocyte_distribution.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Astrocyte distribution visualization saved to {output_path}")
    
    def visualize_brain_regions(self):
        """
        Visualize astrocytes by brain region.
        """
        print("Visualizing astrocytes by brain region...")
        
        # Get spatial points data
        df = self.get_spatial_points_data()
        
        if df.empty:
            print("No astrocyte data found for visualization")
            return
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Plot points colored by brain region
        for region in df['region'].unique():
            if region == "Unknown":
                continue
                
            subset = df[df['region'] == region]
            if not subset.empty:
                color = self.region_colors.get(region, "#999999")  # Default gray if region not in color map
                plt.scatter(subset['x'], subset['y'], c=color, alpha=0.7, 
                           s=subset['probability'] * 100, label=region)
        
        # Add legend
        plt.legend(title="Brain Region", fontsize=12)
        
        # Add labels and title
        plt.xlabel("X Coordinate", fontsize=14)
        plt.ylabel("Y Coordinate", fontsize=14)
        plt.title("Astrocytes by Brain Region", fontsize=16)
        
        # Add grid
        plt.grid(alpha=0.3)
        
        # Save figure
        output_path = os.path.join(self.output_dir, "astrocytes_by_region.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Brain region visualization saved to {output_path}")
    
    def visualize_gene_expression(self, gene_id="GENE1"):
        """
        Visualize gene expression levels across astrocyte types.
        
        Parameters:
        -----------
        gene_id : str
            Gene ID to visualize expression for
        """
        print(f"Visualizing expression of gene {gene_id}...")
        
        # Get gene expression data
        df = self.get_gene_expression_data(gene_id)
        
        if df.empty:
            print(f"No expression data found for gene {gene_id}")
            return
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Create a custom colormap for expression levels
        cmap = plt.cm.viridis
        
        # Plot points with color intensity based on expression level
        scatter = plt.scatter(df['x'], df['y'], c=df['expression_level'], 
                             cmap=cmap, alpha=0.8, s=80, 
                             edgecolor='black', linewidth=0.5)
        
        # Add colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('Expression Level', fontsize=12)
        
        # Add markers for different astrocyte types
        for astrocyte_type, marker in zip(["Protoplasmic", "Fibrous", "Reactive"], ["o", "s", "^"]):
            subset = df[df['astrocyte_type'] == astrocyte_type]
            if not subset.empty:
                plt.scatter(subset['x'], subset['y'], s=120, facecolors='none', 
                           edgecolors=self.astrocyte_colors[astrocyte_type], 
                           linewidth=2, marker=marker, label=astrocyte_type)
        
        # Add legend
        plt.legend(title="Astrocyte Type", fontsize=12)
        
        # Add labels and title
        plt.xlabel("X Coordinate", fontsize=14)
        plt.ylabel("Y Coordinate", fontsize=14)
        plt.title(f"Expression of {gene_id} Across Astrocyte Types", fontsize=16)
        
        # Add grid
        plt.grid(alpha=0.3)
        
        # Save figure
        output_path = os.path.join(self.output_dir, f"gene_expression_{gene_id}.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Gene expression visualization saved to {output_path}")
    
    def visualize_astrocyte_connectivity(self):
        """
        Visualize connectivity between astrocytes.
        """
        print("Visualizing astrocyte connectivity...")
        
        # Query to get connectivity data
        query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point1 ?x1 ?y1 ?type1 ?point2 ?x2 ?y2 ?type2 ?strength ?connType
        WHERE {
            # Get first point
            ?point1 rdf:type st:SpatialPoint ;
                    st:hasXCoordinate ?x1 ;
                    st:hasYCoordinate ?y1 ;
                    st:hasAstrocyteType ?typeNode1 ;
                    st:connectsTo ?point2 .
            
            ?typeNode1 st:astrocyteType ?type1 .
            
            # Get connection details
            ?point1 st:hasConnection ?connNode .
            ?connNode st:toPoint ?point2 ;
                      st:connectionStrength ?strength .
            
            OPTIONAL {
                ?connNode st:connectionType ?connType .
            }
            
            # Get second point
            ?point2 st:hasXCoordinate ?x2 ;
                    st:hasYCoordinate ?y2 ;
                    st:hasAstrocyteType ?typeNode2 .
            
            ?typeNode2 st:astrocyteType ?type2 .
            
            # Only include astrocyte types
            FILTER(?type1 IN ("Protoplasmic", "Fibrous", "Reactive"))
            FILTER(?type2 IN ("Protoplasmic", "Fibrous", "Reactive"))
        }
        LIMIT 100
        """
        
        results = list(self.graph.query(query))
        
        # Convert to DataFrame
        data = []
        for row in results:
            x1 = float(row.x1)
            y1 = float(row.y1)
            type1 = str(row.type1)
            x2 = float(row.x2)
            y2 = float(row.y2)
            type2 = str(row.type2)
            strength = float(row.strength) if row.strength else 0.5
            conn_type = str(row.connType) if row.connType else "unknown"
            
            data.append({
                'x1': x1, 'y1': y1, 'type1': type1,
                'x2': x2, 'y2': y2, 'type2': type2,
                'strength': strength, 'conn_type': conn_type
            })
        
        df = pd.DataFrame(data)
        
        if df.empty:
            print("No connectivity data found for visualization")
            return
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Plot points
        for astrocyte_type, color in self.astrocyte_colors.items():
            # Plot first points
            subset1 = df[df['type1'] == astrocyte_type]
            if not subset1.empty:
                plt.scatter(subset1['x1'], subset1['y1'], c=color, alpha=0.7, s=80, label=astrocyte_type)
            
            # Plot second points (without adding to legend)
            subset2 = df[df['type2'] == astrocyte_type]
            if not subset2.empty:
                plt.scatter(subset2['x2'], subset2['y2'], c=color, alpha=0.7, s=80)
        
        # Plot connections
        for _, row in df.iterrows():
            # Determine line style based on connection type
            if row['conn_type'] == 'homotypic':
                linestyle = '-'
                alpha = 0.7
            else:
                linestyle = '--'
                alpha = 0.5
            
            # Determine line width based on connection strength
            linewidth = row['strength'] * 2
            
            # Determine color based on connection type
            if row['conn_type'] == 'homotypic':
                color = self.astrocyte_colors.get(row['type1'], 'gray')
            else:
                color = 'gray'
            
            # Plot connection line
            plt.plot([row['x1'], row['x2']], [row['y1'], row['y2']], 
                    linestyle=linestyle, color=color, alpha=alpha, linewidth=linewidth)
        
        # Add legend for astrocyte types
        plt.legend(title="Astrocyte Type", fontsize=12)
        
        # Add legend for connection types
        conn_legend_elements = [
            Patch(facecolor='gray', edgecolor='black', label='Homotypic', alpha=0.7),
            Patch(facecolor='gray', edgecolor='black', label='Heterotypic', alpha=0.5, hatch='/')
        ]
        plt.legend(handles=conn_legend_elements, title="Connection Type", fontsize=12, loc='upper right')
        
        # Add labels and title
        plt.xlabel("X Coordinate", fontsize=14)
        plt.ylabel("Y Coordinate", fontsize=14)
        plt.title("Astrocyte Connectivity Network", fontsize=16)
        
        # Add grid
        plt.grid(alpha=0.3)
        
        # Save figure
        output_path = os.path.join(self.output_dir, "astrocyte_connectivity.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Astrocyte connectivity visualization saved to {output_path}")
    
    def visualize_cross_modality_integration(self):
        """
        Visualize the integration between spatial and single-cell data.
        """
        print("Visualizing cross-modality integration...")
        
        # Get spatial points data
        spatial_df = self.get_spatial_points_data()
        
        # Get single-cell data
        sc_df = self.get_single_cell_data()
        
        if spatial_df.empty or sc_df.empty:
            print("Insufficient data for cross-modality visualization")
            return
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Plot spatial points as circles
        for astrocyte_type, color in self.astrocyte_colors.items():
            subset = spatial_df[spatial_df['astrocyte_type'] == astrocyte_type]
            if not subset.empty:
                plt.scatter(subset['x'], subset['y'], c=color, alpha=0.3, 
                           s=subset['probability'] * 100, marker='o', 
                           label=f"Spatial {astrocyte_type}")
        
        # Plot single-cell points as triangles
        for astrocyte_type, color in self.astrocyte_colors.items():
            subset = sc_df[sc_df['cell_type'] == astrocyte_type]
            if not subset.empty:
                plt.scatter(subset['x'], subset['y'], c=color, alpha=0.7, 
                           s=80, marker='^', edgecolor='black', linewidth=0.5,
                           label=f"Single-cell {astrocyte_type}")
        
        # Add legend
        plt.legend(title="Data Source & Type", fontsize=10)
        
        # Add labels and title
        plt.xlabel("X Coordinate", fontsize=14)
        plt.ylabel("Y Coordinate", fontsize=14)
        plt.title("Integration of Spatial and Single-cell Astrocyte Data", fontsize=16)
        
        # Add grid
        plt.grid(alpha=0.3)
        
        # Save figure
        output_path = os.path.join(self.output_dir, "cross_modality_integration.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Cross-modality integration visualization saved to {output_path}")
    
    def run_visualizations(self):
        """
        Run all visualizations.
        """
        print("Running all visualizations...")
        
        # Visualize astrocyte distribution
        self.visualize_astrocyte_distribution()
        
        # Visualize brain regions
        self.visualize_brain_regions()
        
        # Visualize gene expression for each gene
        for gene_id in ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]:
            self.visualize_gene_expression(gene_id)
        
        # Visualize astrocyte connectivity
        self.visualize_astrocyte_connectivity()
        
        # Visualize cross-modality integration
        self.visualize_cross_modality_integration()
        
        print("All visualizations completed")

def main():
    """
    Main function to run the astrocyte visualizer.
    """
    parser = argparse.ArgumentParser(description='Visualize astrocyte data from integrated ontology')
    parser.add_argument('--input', required=True, help='Path to integrated ontology TURTLE file')
    parser.add_argument('--output-dir', default='output/visualizations', help='Output directory for visualizations')
    
    args = parser.parse_args()
    
    # Create visualizer
    visualizer = AstrocyteVisualizer(
        integrated_ttl=args.input,
        output_dir=args.output_dir
    )
    
    # Run visualizations
    visualizer.run_visualizations()

if __name__ == "__main__":
    main() 