#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef, Literal, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL
import networkx as nx
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime, timedelta

class TemporalSpatialAnalyzer:
    """
    A class to simulate and analyze temporal changes in spatial gene expression patterns.
    """
    
    def __init__(self, ttl_file, output_dir="../output/temporal_analysis"):
        """
        Initialize the analyzer with the path to the TURTLE file and output directory.
        
        Parameters:
        -----------
        ttl_file : str
            Path to the TURTLE file containing the enhanced spatial ontology
        output_dir : str
            Path to the directory where analysis results will be saved
        """
        self.ttl_file = ttl_file
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Load the ontology
        self.graph = Graph()
        print(f"Loading TURTLE file: {ttl_file}")
        self.graph.parse(ttl_file, format="turtle")
        print(f"Graph loaded with {len(self.graph)} triples")
        
        # Define namespaces
        self.ns = Namespace("http://example.org/spatial-transcriptomics#")
        
        # Create a new graph for the temporal extension
        self.temporal_graph = Graph()
        # Copy all triples from the original graph
        for s, p, o in self.graph:
            self.temporal_graph.add((s, p, o))
        
        # Number of time points to simulate
        self.num_timepoints = 5
    
    def get_spatial_points(self):
        """
        Get all spatial points from the ontology.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame containing spatial points with their coordinates
        """
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y
        WHERE {
            ?point rdf:type :SpatialPoint ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
        }
        """
        
        results = self.graph.query(query)
        
        data = []
        for row in results:
            data.append({
                'point': str(row.point),
                'x': float(str(row.x)),
                'y': float(str(row.y))
            })
        
        return pd.DataFrame(data)
    
    def get_genes(self):
        """
        Get all genes from the ontology.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame containing genes with their expression points
        """
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?point ?x ?y ?expressionLevel
        WHERE {
            ?gene rdf:type :Gene ;
                  :expressedAt ?point ;
                  :hasExpressionLevel ?expressionLevel .
                  
            ?point :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
        }
        """
        
        results = self.graph.query(query)
        
        data = []
        for row in results:
            gene_uri = str(row.gene)
            gene_symbol = gene_uri.split('_')[1] if '_' in gene_uri else gene_uri.split('#')[-1]
            
            data.append({
                'gene': gene_symbol,
                'gene_uri': gene_uri,
                'point': str(row.point),
                'x': float(str(row.x)),
                'y': float(str(row.y)),
                'expression_level': float(str(row.expressionLevel))
            })
        
        return pd.DataFrame(data)
    
    def simulate_temporal_expression(self):
        """
        Simulate temporal changes in gene expression levels.
        
        Returns:
        --------
        dict
            Dictionary mapping time points to DataFrames containing gene expression data
        """
        print("\nSimulating temporal gene expression changes...")
        
        # Get genes and their expression points
        genes_df = self.get_genes()
        
        if genes_df.empty:
            print("No gene expression data found")
            return {}
        
        # Group by gene
        gene_groups = genes_df.groupby('gene')
        
        # Create time points (e.g., days)
        time_points = [datetime.now() + timedelta(days=i) for i in range(self.num_timepoints)]
        time_point_strs = [tp.strftime("%Y-%m-%d") for tp in time_points]
        
        # Dictionary to store expression data for each time point
        temporal_data = {}
        
        # Seed for reproducibility
        np.random.seed(42)
        
        # For each time point, simulate expression changes
        for i, time_point in enumerate(time_points):
            time_point_str = time_point_strs[i]
            print(f"Simulating expression for time point: {time_point_str}")
            
            # Copy the original data
            time_data = genes_df.copy()
            
            # Add time point column
            time_data['time_point'] = time_point_str
            
            # Simulate expression changes for each gene
            for gene, group in gene_groups:
                # Different genes have different temporal patterns
                pattern_type = hash(gene) % 4
                
                if pattern_type == 0:
                    # Gradual increase
                    factor = 1.0 + (i * 0.2)
                elif pattern_type == 1:
                    # Gradual decrease
                    factor = 1.0 - (i * 0.15)
                elif pattern_type == 2:
                    # Oscillating
                    factor = 1.0 + 0.3 * np.sin(i * np.pi / 2)
                else:
                    # Random fluctuation
                    factor = 1.0 + 0.2 * np.random.randn()
                
                # Ensure factor is positive
                factor = max(0.1, factor)
                
                # Apply factor to expression levels for this gene
                mask = time_data['gene'] == gene
                time_data.loc[mask, 'expression_level'] = time_data.loc[mask, 'expression_level'] * factor
                
                # Add some random noise
                noise = 0.1 * np.random.randn(sum(mask))
                time_data.loc[mask, 'expression_level'] += noise
                
                # Ensure expression levels are positive
                time_data.loc[mask, 'expression_level'] = time_data.loc[mask, 'expression_level'].clip(lower=0.1)
            
            temporal_data[time_point_str] = time_data
            
            # Add temporal data to the RDF graph
            self.add_temporal_data_to_graph(time_data, time_point_str)
        
        return temporal_data
    
    def add_temporal_data_to_graph(self, time_data, time_point_str):
        """
        Add temporal gene expression data to the RDF graph.
        
        Parameters:
        -----------
        time_data : pd.DataFrame
            DataFrame containing gene expression data for a specific time point
        time_point_str : str
            String representation of the time point
        """
        # Create a time point node
        time_point_uri = URIRef(f"{self.ns}TimePoint_{time_point_str}")
        
        # Add time point to the graph
        self.temporal_graph.add((time_point_uri, RDF.type, self.ns.TimePoint))
        self.temporal_graph.add((time_point_uri, self.ns.hasTimeValue, Literal(time_point_str, datatype=XSD.date)))
        
        # Add expression data for each gene at this time point
        for _, row in time_data.iterrows():
            gene_uri = URIRef(row['gene_uri'])
            point_uri = URIRef(row['point'])
            
            # Create a temporal expression node
            temp_expr_node = BNode()
            
            # Add temporal expression data
            self.temporal_graph.add((temp_expr_node, RDF.type, self.ns.TemporalExpression))
            self.temporal_graph.add((temp_expr_node, self.ns.hasGene, gene_uri))
            self.temporal_graph.add((temp_expr_node, self.ns.hasSpatialPoint, point_uri))
            self.temporal_graph.add((temp_expr_node, self.ns.hasTimePoint, time_point_uri))
            self.temporal_graph.add((temp_expr_node, self.ns.hasExpressionLevel, 
                                    Literal(float(row['expression_level']), datatype=XSD.float)))
    
    def analyze_temporal_patterns(self, temporal_data):
        """
        Analyze temporal patterns in gene expression.
        
        Parameters:
        -----------
        temporal_data : dict
            Dictionary mapping time points to DataFrames containing gene expression data
        """
        print("\nAnalyzing temporal expression patterns...")
        
        if not temporal_data:
            print("No temporal data to analyze")
            return
        
        # Get all unique genes
        all_genes = set()
        for time_point, data in temporal_data.items():
            all_genes.update(data['gene'].unique())
        
        # Create a DataFrame to store temporal patterns for each gene
        pattern_data = []
        
        # Analyze each gene's temporal pattern
        for gene in all_genes:
            # Extract expression levels for this gene across time points
            gene_data = []
            for time_point, data in temporal_data.items():
                gene_expr = data[data['gene'] == gene]
                if not gene_expr.empty:
                    avg_expr = gene_expr['expression_level'].mean()
                    gene_data.append({
                        'gene': gene,
                        'time_point': time_point,
                        'avg_expression': avg_expr
                    })
            
            if gene_data:
                gene_df = pd.DataFrame(gene_data)
                
                # Calculate temporal pattern metrics
                expr_values = gene_df['avg_expression'].values
                
                # Trend: positive slope = increasing, negative = decreasing
                if len(expr_values) > 1:
                    slope = np.polyfit(range(len(expr_values)), expr_values, 1)[0]
                else:
                    slope = 0
                
                # Variability: standard deviation
                variability = np.std(expr_values)
                
                # Peak time point
                peak_idx = np.argmax(expr_values)
                peak_time = gene_df.iloc[peak_idx]['time_point']
                
                # Pattern classification
                if slope > 0.1:
                    pattern = "Increasing"
                elif slope < -0.1:
                    pattern = "Decreasing"
                elif variability > 0.5:
                    pattern = "Oscillating"
                else:
                    pattern = "Stable"
                
                pattern_data.append({
                    'gene': gene,
                    'pattern': pattern,
                    'slope': slope,
                    'variability': variability,
                    'peak_time': peak_time
                })
        
        # Convert to DataFrame
        patterns_df = pd.DataFrame(pattern_data)
        
        if not patterns_df.empty:
            print("\nTemporal pattern classification:")
            pattern_counts = patterns_df['pattern'].value_counts()
            print(pattern_counts)
            
            # Save pattern data
            output_file = os.path.join(self.output_dir, "temporal_patterns.csv")
            patterns_df.to_csv(output_file, index=False)
            print(f"Temporal pattern data saved to {output_file}")
            
            # Visualize pattern distribution
            plt.figure(figsize=(10, 6))
            sns.countplot(x='pattern', data=patterns_df)
            plt.title('Distribution of Temporal Expression Patterns')
            plt.xlabel('Pattern Type')
            plt.ylabel('Number of Genes')
            plt.tight_layout()
            
            output_file = os.path.join(self.output_dir, "pattern_distribution.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Pattern distribution visualization saved to {output_file}")
            plt.close()
            
            # Visualize example genes for each pattern
            self.visualize_example_patterns(patterns_df, temporal_data)
        
        return patterns_df
    
    def visualize_example_patterns(self, patterns_df, temporal_data, num_examples=3):
        """
        Visualize example genes for each temporal pattern.
        
        Parameters:
        -----------
        patterns_df : pd.DataFrame
            DataFrame containing temporal pattern classifications
        temporal_data : dict
            Dictionary mapping time points to DataFrames containing gene expression data
        num_examples : int
            Number of example genes to visualize for each pattern
        """
        # Get unique patterns
        patterns = patterns_df['pattern'].unique()
        
        # Time points in order
        time_points = sorted(temporal_data.keys())
        
        # Create figure with subplots for each pattern
        fig, axes = plt.subplots(len(patterns), 1, figsize=(12, 4 * len(patterns)))
        
        if len(patterns) == 1:
            axes = [axes]
        
        for i, pattern in enumerate(patterns):
            # Get example genes for this pattern
            pattern_genes = patterns_df[patterns_df['pattern'] == pattern]['gene'].head(num_examples).tolist()
            
            if not pattern_genes:
                continue
            
            # Plot each example gene
            for gene in pattern_genes:
                # Extract expression levels for this gene across time points
                expr_values = []
                for tp in time_points:
                    gene_expr = temporal_data[tp][temporal_data[tp]['gene'] == gene]
                    if not gene_expr.empty:
                        avg_expr = gene_expr['expression_level'].mean()
                        expr_values.append(avg_expr)
                    else:
                        expr_values.append(np.nan)
                
                # Plot expression over time
                axes[i].plot(time_points, expr_values, marker='o', label=gene)
            
            # Set title and labels
            axes[i].set_title(f'{pattern} Expression Pattern')
            axes[i].set_xlabel('Time Point')
            axes[i].set_ylabel('Average Expression Level')
            axes[i].legend()
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, "example_patterns.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Example patterns visualization saved to {output_file}")
        plt.close()
    
    def create_expression_animation(self, temporal_data, gene):
        """
        Create an animation of expression changes over time for a specific gene.
        
        Parameters:
        -----------
        temporal_data : dict
            Dictionary mapping time points to DataFrames containing gene expression data
        gene : str
            Gene symbol to visualize
        """
        print(f"\nCreating animation for gene: {gene}")
        
        # Time points in order
        time_points = sorted(temporal_data.keys())
        
        # Extract data for this gene
        gene_data = {}
        for tp in time_points:
            gene_expr = temporal_data[tp][temporal_data[tp]['gene'] == gene]
            if not gene_expr.empty:
                gene_data[tp] = gene_expr
        
        if not gene_data:
            print(f"No expression data found for gene: {gene}")
            return
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Function to update the plot for each frame
        def update(frame):
            ax.clear()
            
            tp = time_points[frame]
            data = gene_data[tp]
            
            # Create scatter plot of expression
            scatter = ax.scatter(data['x'], data['y'], 
                               c=data['expression_level'], cmap='viridis', 
                               s=100, alpha=0.7)
            
            # Add colorbar
            if frame == 0:
                self.cbar = plt.colorbar(scatter, ax=ax)
                self.cbar.set_label('Expression Level')
            
            # Set title and labels
            ax.set_title(f'Expression of {gene} at {tp}')
            ax.set_xlabel('X Coordinate')
            ax.set_ylabel('Y Coordinate')
            ax.grid(True, alpha=0.3)
            
            return scatter,
        
        # Create animation
        ani = FuncAnimation(fig, update, frames=len(time_points), blit=True)
        
        # Save animation
        output_file = os.path.join(self.output_dir, f"animation_{gene}.gif")
        ani.save(output_file, writer='pillow', fps=1)
        print(f"Animation saved to {output_file}")
        plt.close()
    
    def create_spatial_temporal_heatmap(self, temporal_data):
        """
        Create a heatmap showing spatial-temporal expression patterns.
        
        Parameters:
        -----------
        temporal_data : dict
            Dictionary mapping time points to DataFrames containing gene expression data
        """
        print("\nCreating spatial-temporal heatmap...")
        
        # Time points in order
        time_points = sorted(temporal_data.keys())
        
        if not time_points:
            print("No temporal data available")
            return
        
        # Get all unique genes
        all_genes = set()
        for tp in time_points:
            all_genes.update(temporal_data[tp]['gene'].unique())
        
        # Select a subset of genes for visualization (to avoid overcrowding)
        if len(all_genes) > 20:
            selected_genes = list(all_genes)[:20]
        else:
            selected_genes = list(all_genes)
        
        # Create a matrix of expression values: genes x time points
        expr_matrix = np.zeros((len(selected_genes), len(time_points)))
        
        for i, gene in enumerate(selected_genes):
            for j, tp in enumerate(time_points):
                gene_expr = temporal_data[tp][temporal_data[tp]['gene'] == gene]
                if not gene_expr.empty:
                    expr_matrix[i, j] = gene_expr['expression_level'].mean()
        
        # Create heatmap
        plt.figure(figsize=(12, 10))
        sns.heatmap(expr_matrix, cmap='viridis', 
                   xticklabels=time_points, yticklabels=selected_genes,
                   annot=False, fmt='.2f')
        
        plt.title('Spatial-Temporal Expression Heatmap')
        plt.xlabel('Time Point')
        plt.ylabel('Gene')
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, "spatial_temporal_heatmap.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved to {output_file}")
        plt.close()
    
    def save_temporal_ontology(self):
        """
        Save the temporal ontology to a TURTLE file.
        """
        # Add temporal ontology classes and properties
        self.temporal_graph.add((self.ns.TimePoint, RDF.type, OWL.Class))
        self.temporal_graph.add((self.ns.hasTimeValue, RDF.type, OWL.DatatypeProperty))
        self.temporal_graph.add((self.ns.TemporalExpression, RDF.type, OWL.Class))
        self.temporal_graph.add((self.ns.hasTimePoint, RDF.type, OWL.ObjectProperty))
        self.temporal_graph.add((self.ns.hasSpatialPoint, RDF.type, OWL.ObjectProperty))
        
        # Save to file
        output_file = os.path.join(self.output_dir, "temporal_spatial_ontology.ttl")
        self.temporal_graph.serialize(destination=output_file, format="turtle")
        print(f"\nTemporal spatial ontology saved to {output_file}")
    
    def run_analysis(self):
        """
        Run the temporal spatial analysis.
        """
        # Simulate temporal expression changes
        temporal_data = self.simulate_temporal_expression()
        
        # Analyze temporal patterns
        patterns_df = self.analyze_temporal_patterns(temporal_data)
        
        # Create animations for example genes from each pattern
        if patterns_df is not None and not patterns_df.empty:
            for pattern in patterns_df['pattern'].unique():
                example_gene = patterns_df[patterns_df['pattern'] == pattern]['gene'].iloc[0]
                self.create_expression_animation(temporal_data, example_gene)
        
        # Create spatial-temporal heatmap
        self.create_spatial_temporal_heatmap(temporal_data)
        
        # Save the temporal ontology
        self.save_temporal_ontology()
        
        print("\nTemporal spatial analysis completed!")

def main():
    parser = argparse.ArgumentParser(description='Temporal Spatial Analyzer')
    parser.add_argument('--ttl', type=str, default='../data/enhanced_spatial_ontology.ttl',
                        help='Path to the TURTLE file containing the enhanced spatial ontology')
    parser.add_argument('--output', type=str, default='../output/temporal_analysis',
                        help='Path to the output directory')
    
    args = parser.parse_args()
    
    analyzer = TemporalSpatialAnalyzer(args.ttl, args.output)
    analyzer.run_analysis()
    
    print("\nAnalysis completed successfully!")

if __name__ == "__main__":
    main() 