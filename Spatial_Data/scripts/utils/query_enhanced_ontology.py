#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF, RDFS

class EnhancedOntologyQuerier:
    """
    A class to query the enhanced spatial ontology for specific spatial relationships
    and extract valuable insights.
    """
    
    def __init__(self, ttl_file, output_dir="../output/query_results"):
        """
        Initialize the querier with the path to the TURTLE file and output directory.
        
        Parameters:
        -----------
        ttl_file : str
            Path to the TURTLE file containing the enhanced spatial ontology
        output_dir : str
            Path to the directory where query results will be saved
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
    
    def find_points_with_directional_relationships(self, direction_property, limit=10):
        """
        Find spatial points that have a specific directional relationship with other points.
        
        Parameters:
        -----------
        direction_property : str
            The directional property to query for (e.g., 'directlyAbove', 'directlyRightOf')
        limit : int
            Maximum number of results to return
            
        Returns:
        --------
        pd.DataFrame
            DataFrame containing point pairs with the specified directional relationship
        """
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point1 ?point2 ?x1 ?y1 ?x2 ?y2
        WHERE {{
            ?point1 :{direction_property} ?point2 ;
                    :hasXCoordinate ?x1 ;
                    :hasYCoordinate ?y1 .
                    
            ?point2 :hasXCoordinate ?x2 ;
                    :hasYCoordinate ?y2 .
        }}
        LIMIT {limit}
        """
        
        results = self.graph.query(query)
        
        data = []
        for row in results:
            data.append({
                'point1': str(row.point1),
                'point2': str(row.point2),
                'x1': float(str(row.x1)),
                'y1': float(str(row.y1)),
                'x2': float(str(row.x2)),
                'y2': float(str(row.y2))
            })
        
        return pd.DataFrame(data)
    
    def find_points_in_region(self, region_name):
        """
        Find spatial points in a specific spatial region.
        
        Parameters:
        -----------
        region_name : str
            The name of the region to query for (e.g., 'Region_1_1')
            
        Returns:
        --------
        pd.DataFrame
            DataFrame containing points in the specified region
        """
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y
        WHERE {{
            ?point rdf:type :SpatialPoint ;
                   :locatedInRegion ?region ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
                   
            ?region :hasRegionName "{region_name}" .
        }}
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
    
    def find_points_with_distance_relationship(self, distance_property):
        """
        Find spatial points that have a specific distance relationship with other points.
        
        Parameters:
        -----------
        distance_property : str
            The distance property to query for (e.g., 'veryCloseToPoint', 'moderatelyCloseToPoint')
            
        Returns:
        --------
        pd.DataFrame
            DataFrame containing point pairs with the specified distance relationship
        """
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point1 ?point2 ?x1 ?y1 ?x2 ?y2
        WHERE {{
            ?point1 :{distance_property} ?point2 ;
                    :hasXCoordinate ?x1 ;
                    :hasYCoordinate ?y1 .
                    
            ?point2 :hasXCoordinate ?x2 ;
                    :hasYCoordinate ?y2 .
        }}
        LIMIT 20
        """
        
        results = self.graph.query(query)
        
        data = []
        for row in results:
            data.append({
                'point1': str(row.point1),
                'point2': str(row.point2),
                'x1': float(str(row.x1)),
                'y1': float(str(row.y1)),
                'x2': float(str(row.x2)),
                'y2': float(str(row.y2))
            })
        
        return pd.DataFrame(data)
    
    def visualize_directional_point_pairs(self, df, direction_property):
        """
        Visualize point pairs with a specific directional relationship.
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame containing point pairs with directional relationships
        direction_property : str
            The directional property being visualized
        """
        if df.empty:
            print(f"No point pairs found with relationship: {direction_property}")
            return
        
        plt.figure(figsize=(10, 8))
        
        # Plot points and arrows
        for _, row in df.iterrows():
            plt.scatter(row['x1'], row['y1'], color='blue', s=50)
            plt.scatter(row['x2'], row['y2'], color='red', s=50)
            
            # Draw arrow from point1 to point2
            plt.arrow(row['x1'], row['y1'], 
                     row['x2'] - row['x1'], row['y2'] - row['y1'],
                     head_width=20, head_length=20, fc='black', ec='black',
                     length_includes_head=True, alpha=0.6)
        
        plt.title(f'Point Pairs with {direction_property} Relationship')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.grid(True, alpha=0.3)
        
        # Save figure
        output_file = os.path.join(self.output_dir, f"{direction_property}_point_pairs.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
        plt.close()
    
    def visualize_points_in_region(self, df, region_name):
        """
        Visualize points in a specific spatial region.
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame containing points in the region
        region_name : str
            The name of the region being visualized
        """
        if df.empty:
            print(f"No points found in region: {region_name}")
            return
        
        plt.figure(figsize=(10, 8))
        
        # Plot points
        plt.scatter(df['x'], df['y'], color='green', s=50, alpha=0.7)
        
        plt.title(f'Points in {region_name}')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.grid(True, alpha=0.3)
        
        # Save figure
        output_file = os.path.join(self.output_dir, f"points_in_{region_name}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
        plt.close()
    
    def run_example_queries(self):
        """
        Run example queries to demonstrate the capabilities of the enhanced ontology.
        """
        print("\n=== Querying for Directional Relationships ===")
        for direction in ['directlyAbove', 'directlyBelow', 'directlyLeftOf', 'directlyRightOf', 'diagonallyRelatedTo']:
            print(f"\nFinding point pairs with {direction} relationship:")
            df = self.find_points_with_directional_relationships(direction, limit=5)
            if not df.empty:
                print(f"Found {len(df)} point pairs with {direction} relationship")
                print(df.head())
                self.visualize_directional_point_pairs(df, direction)
            else:
                print(f"No point pairs found with {direction} relationship")
        
        print("\n=== Querying for Points in Spatial Regions ===")
        for region in ['Region_1_1', 'Region_2_2', 'Region_3_3']:
            print(f"\nFinding points in {region}:")
            df = self.find_points_in_region(region)
            if not df.empty:
                print(f"Found {len(df)} points in {region}")
                print(df.head())
                self.visualize_points_in_region(df, region)
            else:
                print(f"No points found in {region}")
        
        print("\n=== Querying for Distance Relationships ===")
        for distance in ['veryCloseToPoint', 'moderatelyCloseToPoint', 'farFromPoint']:
            print(f"\nFinding point pairs with {distance} relationship:")
            df = self.find_points_with_distance_relationship(distance)
            if not df.empty:
                print(f"Found {len(df)} point pairs with {distance} relationship")
                print(df.head())
            else:
                print(f"No point pairs found with {distance} relationship")

def main():
    parser = argparse.ArgumentParser(description='Query enhanced spatial ontology')
    parser.add_argument('--ttl', type=str, default='../data/enhanced_spatial_ontology.ttl',
                        help='Path to the TURTLE file containing the enhanced spatial ontology')
    parser.add_argument('--output', type=str, default='../output/query_results',
                        help='Path to the output directory')
    
    args = parser.parse_args()
    
    querier = EnhancedOntologyQuerier(args.ttl, args.output)
    querier.run_example_queries()
    
    print("\nQuery execution completed successfully!")

if __name__ == "__main__":
    main() 