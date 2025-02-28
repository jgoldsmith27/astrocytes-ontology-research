#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF, RDFS

class SpatialGeneRelationshipAnalyzer:
    """
    A class to analyze spatial relationships between genes based on the enhanced ontology.
    """
    
    def __init__(self, ttl_file, output_dir="../output/gene_relationships"):
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
    
    def find_genes_with_spatial_relationship(self, relationship_property):
        """
        Find genes that are expressed at points with a specific spatial relationship.
        
        Parameters:
        -----------
        relationship_property : str
            The spatial relationship property to query for (e.g., 'directlyAbove', 'moderatelyCloseToPoint')
            
        Returns:
        --------
        pd.DataFrame
            DataFrame containing gene pairs with the specified spatial relationship
        """
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene1 ?gene2 ?x1 ?y1 ?x2 ?y2
        WHERE {{
            ?gene1 rdf:type :Gene ;
                   :expressedAt ?point1 .
                   
            ?gene2 rdf:type :Gene ;
                   :expressedAt ?point2 .
                   
            ?point1 :{relationship_property} ?point2 ;
                    :hasXCoordinate ?x1 ;
                    :hasYCoordinate ?y1 .
                    
            ?point2 :hasXCoordinate ?x2 ;
                    :hasYCoordinate ?y2 .
                    
            # Ensure we get distinct gene pairs
            FILTER(?gene1 != ?gene2)
        }}
        LIMIT 20
        """
        
        results = self.graph.query(query)
        
        data = []
        for row in results:
            # Extract gene symbols from URIs
            gene1_uri = str(row.gene1)
            gene2_uri = str(row.gene2)
            
            gene1_symbol = gene1_uri.split('_')[1] if '_' in gene1_uri else gene1_uri.split('#')[-1]
            gene2_symbol = gene2_uri.split('_')[1] if '_' in gene2_uri else gene2_uri.split('#')[-1]
            
            data.append({
                'gene1': gene1_symbol,
                'gene2': gene2_symbol,
                'x1': float(str(row.x1)),
                'y1': float(str(row.y1)),
                'x2': float(str(row.x2)),
                'y2': float(str(row.y2))
            })
        
        return pd.DataFrame(data)
    
    def find_genes_in_same_region(self, region_name):
        """
        Find genes expressed in the same spatial region.
        
        Parameters:
        -----------
        region_name : str
            The name of the region to query for (e.g., 'Region_1_1')
            
        Returns:
        --------
        pd.DataFrame
            DataFrame containing genes expressed in the specified region
        """
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?x ?y
        WHERE {{
            ?gene rdf:type :Gene ;
                  :expressedAt ?point .
                  
            ?point :locatedInRegion ?region ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
                   
            ?region :hasRegionName "{region_name}" .
        }}
        """
        
        results = self.graph.query(query)
        
        data = []
        for row in results:
            # Extract gene symbol from URI
            gene_uri = str(row.gene)
            gene_symbol = gene_uri.split('_')[1] if '_' in gene_uri else gene_uri.split('#')[-1]
            
            data.append({
                'gene': gene_symbol,
                'x': float(str(row.x)),
                'y': float(str(row.y))
            })
        
        return pd.DataFrame(data)
    
    def find_co_expressed_genes_by_distance(self, distance_property):
        """
        Find co-expressed genes based on distance relationships between their expression points.
        
        Parameters:
        -----------
        distance_property : str
            The distance property to query for (e.g., 'veryCloseToPoint', 'moderatelyCloseToPoint')
            
        Returns:
        --------
        pd.DataFrame
            DataFrame containing co-expressed gene pairs
        """
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene1 ?gene2 (COUNT(*) as ?count)
        WHERE {{
            ?gene1 rdf:type :Gene ;
                   :expressedAt ?point1 .
                   
            ?gene2 rdf:type :Gene ;
                   :expressedAt ?point2 .
                   
            ?point1 :{distance_property} ?point2 .
                    
            # Ensure we get distinct gene pairs
            FILTER(?gene1 != ?gene2)
        }}
        GROUP BY ?gene1 ?gene2
        ORDER BY DESC(?count)
        LIMIT 20
        """
        
        results = self.graph.query(query)
        
        data = []
        for row in results:
            # Extract gene symbols from URIs
            gene1_uri = str(row.gene1)
            gene2_uri = str(row.gene2)
            
            gene1_symbol = gene1_uri.split('_')[1] if '_' in gene1_uri else gene1_uri.split('#')[-1]
            gene2_symbol = gene2_uri.split('_')[1] if '_' in gene2_uri else gene2_uri.split('#')[-1]
            
            data.append({
                'gene1': gene1_symbol,
                'gene2': gene2_symbol,
                'count': int(row[2])  # Use index to access count
            })
        
        return pd.DataFrame(data)
    
    def visualize_gene_pairs_with_relationship(self, df, relationship_property):
        """
        Visualize gene pairs with a specific spatial relationship.
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame containing gene pairs with spatial relationships
        relationship_property : str
            The spatial relationship property being visualized
        """
        if df.empty:
            print(f"No gene pairs found with relationship: {relationship_property}")
            return
        
        plt.figure(figsize=(12, 10))
        
        # Plot points and arrows
        for _, row in df.iterrows():
            plt.scatter(row['x1'], row['y1'], color='blue', s=50)
            plt.scatter(row['x2'], row['y2'], color='red', s=50)
            
            # Draw arrow from point1 to point2
            plt.arrow(row['x1'], row['y1'], 
                     row['x2'] - row['x1'], row['y2'] - row['y1'],
                     head_width=20, head_length=20, fc='black', ec='black',
                     length_includes_head=True, alpha=0.6)
            
            # Add gene labels
            plt.text(row['x1'], row['y1'] - 30, row['gene1'], 
                    horizontalalignment='center', verticalalignment='center',
                    fontsize=8)
            plt.text(row['x2'], row['y2'] + 30, row['gene2'], 
                    horizontalalignment='center', verticalalignment='center',
                    fontsize=8)
        
        plt.title(f'Gene Pairs with {relationship_property} Relationship')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.grid(True, alpha=0.3)
        
        # Save figure
        output_file = os.path.join(self.output_dir, f"genes_{relationship_property}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
        plt.close()
    
    def visualize_genes_in_region(self, df, region_name):
        """
        Visualize genes expressed in a specific spatial region.
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame containing genes in the region
        region_name : str
            The name of the region being visualized
        """
        if df.empty:
            print(f"No genes found in region: {region_name}")
            return
        
        plt.figure(figsize=(12, 10))
        
        # Plot points
        plt.scatter(df['x'], df['y'], color='green', s=50, alpha=0.7)
        
        # Add gene labels
        for _, row in df.iterrows():
            plt.text(row['x'], row['y'] + 15, row['gene'], 
                    horizontalalignment='center', verticalalignment='center',
                    fontsize=8)
        
        plt.title(f'Genes Expressed in {region_name}')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.grid(True, alpha=0.3)
        
        # Save figure
        output_file = os.path.join(self.output_dir, f"genes_in_{region_name}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
        plt.close()
    
    def visualize_co_expressed_genes(self, df, distance_property):
        """
        Visualize co-expressed gene pairs based on distance relationships.
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame containing co-expressed gene pairs
        distance_property : str
            The distance property being visualized
        """
        if df.empty:
            print(f"No co-expressed gene pairs found with relationship: {distance_property}")
            return
        
        plt.figure(figsize=(12, 8))
        
        # Create a horizontal bar plot
        sns.barplot(x='count', y='gene1', data=df.head(10), palette='viridis')
        
        # Add gene2 labels to the bars
        for i, row in enumerate(df.head(10).itertuples()):
            plt.text(row.count + 0.5, i, f"with {row.gene2}", 
                    verticalalignment='center', fontsize=10)
        
        plt.title(f'Top Co-expressed Gene Pairs ({distance_property})')
        plt.xlabel('Number of Point Pairs')
        plt.ylabel('Gene')
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, f"co_expressed_genes_{distance_property}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
        plt.close()
    
    def analyze_spatial_gene_relationships(self):
        """
        Analyze spatial relationships between genes and generate visualizations.
        """
        print("\n=== Analyzing Directional Relationships Between Genes ===")
        for direction in ['directlyAbove', 'directlyBelow', 'directlyLeftOf', 'directlyRightOf', 'diagonallyRelatedTo']:
            print(f"\nFinding gene pairs with {direction} relationship:")
            df = self.find_genes_with_spatial_relationship(direction)
            if not df.empty:
                print(f"Found {len(df)} gene pairs with {direction} relationship")
                print(df.head())
                self.visualize_gene_pairs_with_relationship(df, direction)
            else:
                print(f"No gene pairs found with {direction} relationship")
        
        print("\n=== Analyzing Genes in Spatial Regions ===")
        for region in ['Region_1_1', 'Region_2_2', 'Region_3_3']:
            print(f"\nFinding genes in {region}:")
            df = self.find_genes_in_same_region(region)
            if not df.empty:
                print(f"Found {len(df)} genes in {region}")
                print(df.head())
                self.visualize_genes_in_region(df, region)
            else:
                print(f"No genes found in {region}")
        
        print("\n=== Analyzing Co-expressed Genes by Distance ===")
        for distance in ['moderatelyCloseToPoint', 'farFromPoint']:
            print(f"\nFinding co-expressed gene pairs with {distance} relationship:")
            df = self.find_co_expressed_genes_by_distance(distance)
            if not df.empty:
                print(f"Found {len(df)} co-expressed gene pairs with {distance} relationship")
                print(df.head())
                self.visualize_co_expressed_genes(df, distance)
            else:
                print(f"No co-expressed gene pairs found with {distance} relationship")
        
        print("\nSpatial gene relationship analysis completed!")

def main():
    parser = argparse.ArgumentParser(description='Analyze spatial relationships between genes')
    parser.add_argument('--ttl', type=str, default='../data/enhanced_spatial_ontology.ttl',
                        help='Path to the TURTLE file containing the enhanced spatial ontology')
    parser.add_argument('--output', type=str, default='../output/gene_relationships',
                        help='Path to the output directory')
    
    args = parser.parse_args()
    
    analyzer = SpatialGeneRelationshipAnalyzer(args.ttl, args.output)
    analyzer.analyze_spatial_gene_relationships()
    
    print("\nAnalysis completed successfully!")

if __name__ == "__main__":
    main() 