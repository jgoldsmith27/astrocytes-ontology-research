import rdflib
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, XSD, OWL
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns
import networkx as nx
from collections import Counter, defaultdict
import argparse
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")

class SpatialPatternAnalyzer:
    """A class to analyze and visualize enhanced spatial relationships in the ontology"""
    
    def __init__(self, ttl_file):
        """Initialize with a TURTLE file"""
        self.ttl_file = ttl_file
        self.graph = self.load_turtle_file(ttl_file)
        self.output_dir = "../output/spatial_analysis"
        os.makedirs(self.output_dir, exist_ok=True)
    
    def load_turtle_file(self, ttl_file):
        """Load a TURTLE file into an RDF graph"""
        print(f"Loading TURTLE file: {ttl_file}")
        g = Graph()
        g.parse(ttl_file, format="turtle")
        print(f"Graph loaded with {len(g)} triples")
        return g
    
    def analyze_directional_relationships(self):
        """Analyze directional relationships between spatial points"""
        print("\nAnalyzing directional relationships...")
        
        # Check if directional properties exist in the ontology
        has_directional = False
        for s, p, o in self.graph.triples((None, RDF.type, OWL.ObjectProperty)):
            if s in [ST.directlyAbove, ST.directlyBelow, ST.directlyLeftOf, ST.directlyRightOf]:
                has_directional = True
                break
        
        if not has_directional:
            print("No directional relationships found in the ontology")
            return None
        
        # Query for directional relationships
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?relation (COUNT(*) as ?count)
        WHERE {
            ?point1 ?relation ?point2 .
            
            # Filter for directional relations
            FILTER(?relation IN (
                :directlyAbove, :directlyBelow, 
                :directlyLeftOf, :directlyRightOf,
                :diagonallyRelatedTo
            ))
        }
        GROUP BY ?relation
        ORDER BY DESC(?count)
        """
        
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            relation = str(row.relation).split('#')[-1]
            count = int(row[1])  # Access count using index instead of attribute
            data.append({
                'relation': relation,
                'count': count
            })
        
        df = pd.DataFrame(data)
        
        if not df.empty:
            print("Directional relationship counts:")
            print(df)
            
            # Visualize relationship counts
            plt.figure(figsize=(10, 6))
            sns.barplot(x='relation', y='count', data=df)
            plt.title('Directional Relationship Counts')
            plt.xlabel('Relationship Type')
            plt.ylabel('Count')
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            output_file = os.path.join(self.output_dir, "directional_relationships.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Directional relationship visualization saved to {output_file}")
            plt.close()
        
        return df
    
    def analyze_spatial_regions(self):
        """Analyze spatial regions and their gene distributions"""
        print("\nAnalyzing spatial regions...")
        
        # Check if spatial regions exist in the ontology
        has_regions = False
        for s, p, o in self.graph.triples((None, RDF.type, ST.SpatialRegion)):
            has_regions = True
            break
        
        if not has_regions:
            print("No spatial regions found in the ontology")
            return None
        
        # Query for regions and their boundaries
        query_regions = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?region ?name ?minX ?maxX ?minY ?maxY
        WHERE {
            ?region rdf:type :SpatialRegion ;
                    :hasRegionName ?name ;
                    :hasMinX ?minX ;
                    :hasMaxX ?maxX ;
                    :hasMinY ?minY ;
                    :hasMaxY ?maxY .
        }
        ORDER BY ?name
        """
        
        results_regions = self.graph.query(query_regions)
        
        # Convert to DataFrame
        regions = []
        for row in results_regions:
            regions.append({
                'uri': row.region,
                'name': str(row.name),
                'min_x': float(str(row.minX)),
                'max_x': float(str(row.maxX)),
                'min_y': float(str(row.minY)),
                'max_y': float(str(row.maxY))
            })
        
        regions_df = pd.DataFrame(regions)
        
        if regions_df.empty:
            print("No region data found")
            return None
        
        print(f"Found {len(regions_df)} spatial regions")
        
        # Query for points in regions
        query_points = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?region ?regionName ?x ?y
        WHERE {
            ?point rdf:type :SpatialPoint ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y ;
                   :locatedInRegion ?region .
            
            ?region :hasRegionName ?regionName .
        }
        """
        
        results_points = self.graph.query(query_points)
        
        # Convert to DataFrame
        points = []
        for row in results_points:
            points.append({
                'point': row.point,
                'region': row.region,
                'region_name': str(row.regionName),
                'x': float(str(row.x)),
                'y': float(str(row.y))
            })
        
        points_df = pd.DataFrame(points)
        
        if points_df.empty:
            print("No points associated with regions")
            return regions_df, None
        
        print(f"Found {len(points_df)} points associated with regions")
        
        # Count points per region
        region_counts = points_df['region_name'].value_counts().reset_index()
        region_counts.columns = ['region_name', 'point_count']
        
        # Query for genes in regions
        query_genes = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?geneID ?region ?regionName
        WHERE {
            ?gene rdf:type :Gene ;
                  :expressedAt ?point .
            
            ?point :locatedInRegion ?region .
            
            ?region :hasRegionName ?regionName .
            
            # Extract gene ID
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        """
        
        results_genes = self.graph.query(query_genes)
        
        # Group genes by region
        region_genes = defaultdict(set)
        for row in results_genes:
            region_name = str(row.regionName)
            gene_id = str(row.geneID)
            region_genes[region_name].add(gene_id)
        
        # Create a DataFrame with gene counts per region
        gene_counts = []
        for region, genes in region_genes.items():
            gene_counts.append({
                'region_name': region,
                'gene_count': len(genes),
                'genes': ', '.join(sorted(list(genes))[:5]) + ('...' if len(genes) > 5 else '')
            })
        
        gene_counts_df = pd.DataFrame(gene_counts)
        
        # Merge with region counts
        if not gene_counts_df.empty:
            region_data = pd.merge(region_counts, gene_counts_df, on='region_name', how='left')
            region_data['gene_count'] = region_data['gene_count'].fillna(0).astype(int)
            
            print("\nRegion statistics:")
            print(region_data)
            
            # Visualize regions and their gene distributions
            self.visualize_spatial_regions(regions_df, points_df, region_data)
        
        return regions_df, points_df
    
    def visualize_spatial_regions(self, regions_df, points_df, region_data):
        """Visualize spatial regions and their gene distributions"""
        if regions_df.empty:
            return
        
        # Create figure
        plt.figure(figsize=(14, 10))
        
        # Plot regions as polygons
        patches = []
        region_colors = []
        
        # Normalize gene counts for color mapping
        max_gene_count = region_data['gene_count'].max() if not region_data.empty else 1
        norm = plt.Normalize(0, max_gene_count)
        
        for _, region in regions_df.iterrows():
            # Create polygon for region
            polygon = Polygon([
                (region['min_x'], region['min_y']),
                (region['max_x'], region['min_y']),
                (region['max_x'], region['max_y']),
                (region['min_x'], region['max_y'])
            ], closed=True)
            
            patches.append(polygon)
            
            # Get gene count for color
            if not region_data.empty:
                gene_count = region_data.loc[region_data['region_name'] == region['name'], 'gene_count'].values
                gene_count = gene_count[0] if len(gene_count) > 0 else 0
            else:
                gene_count = 0
            
            region_colors.append(gene_count)
        
        # Create patch collection
        p = PatchCollection(patches, alpha=0.4, cmap='viridis')
        p.set_array(np.array(region_colors))
        plt.gca().add_collection(p)
        
        # Add colorbar
        cbar = plt.colorbar(p)
        cbar.set_label('Number of Genes')
        
        # Plot points if available
        if points_df is not None and not points_df.empty:
            plt.scatter(points_df['x'], points_df['y'], c='red', s=10, alpha=0.5)
        
        # Add region labels
        for _, region in regions_df.iterrows():
            center_x = (region['min_x'] + region['max_x']) / 2
            center_y = (region['min_y'] + region['max_y']) / 2
            plt.text(center_x, center_y, region['name'], 
                    horizontalalignment='center', verticalalignment='center',
                    fontsize=12, fontweight='bold')
        
        # Set axis limits
        plt.xlim(regions_df['min_x'].min() - 100, regions_df['max_x'].max() + 100)
        plt.ylim(regions_df['min_y'].min() - 100, regions_df['max_y'].max() + 100)
        
        # Set labels and title
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.title('Spatial Regions and Gene Distribution')
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, "spatial_regions.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Spatial regions visualization saved to {output_file}")
        plt.close()
    
    def analyze_distance_relationships(self):
        """Analyze distance-based relationships between points"""
        print("\nAnalyzing distance-based relationships...")
        
        # Check if distance properties exist in the ontology
        has_distance = False
        for s, p, o in self.graph.triples((None, RDF.type, OWL.ObjectProperty)):
            if s in [ST.veryCloseToPoint, ST.moderatelyCloseToPoint, ST.farFromPoint]:
                has_distance = True
                break
        
        if not has_distance:
            print("No distance-based relationships found in the ontology")
            return None
        
        # Query for distance relationships
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?relation (COUNT(*) as ?count)
        WHERE {
            ?point1 ?relation ?point2 .
            
            # Filter for distance relations
            FILTER(?relation IN (
                :veryCloseToPoint, :moderatelyCloseToPoint, :farFromPoint
            ))
        }
        GROUP BY ?relation
        ORDER BY DESC(?count)
        """
        
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            relation = str(row.relation).split('#')[-1]
            count = int(row[1])  # Access count using index instead of attribute
            data.append({
                'relation': relation,
                'count': count
            })
        
        df = pd.DataFrame(data)
        
        if not df.empty:
            print("Distance relationship counts:")
            print(df)
            
            # Visualize relationship counts
            plt.figure(figsize=(10, 6))
            sns.barplot(x='relation', y='count', data=df)
            plt.title('Distance Relationship Counts')
            plt.xlabel('Relationship Type')
            plt.ylabel('Count')
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            output_file = os.path.join(self.output_dir, "distance_relationships.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Distance relationship visualization saved to {output_file}")
            plt.close()
        
        return df
    
    def analyze_spatial_patterns(self):
        """Analyze spatial patterns of gene expression"""
        print("\nAnalyzing spatial expression patterns...")
        
        # Check if pattern classes exist in the ontology
        has_patterns = False
        for s, p, o in self.graph.triples((None, RDF.type, OWL.Class)):
            if s in [ST.SpatialPattern, ST.ClusteredPattern, ST.DispersedPattern, ST.GradientPattern]:
                has_patterns = True
                break
        
        if not has_patterns:
            print("No spatial pattern classes found in the ontology")
            return None
        
        # Query for genes and their spatial patterns
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?geneID ?pattern ?patternType
        WHERE {
            ?gene rdf:type :Gene ;
                  :exhibitsSpatialPattern ?pattern .
            
            ?pattern rdf:type ?patternType .
            
            # Extract gene ID
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
            
            # Filter for pattern types
            FILTER(?patternType IN (
                :ClusteredPattern, :DispersedPattern, :GradientPattern
            ))
        }
        """
        
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            gene_id = str(row.geneID)
            pattern_type = str(row.patternType).split('#')[-1]
            data.append({
                'gene': gene_id,
                'pattern_type': pattern_type
            })
        
        df = pd.DataFrame(data)
        
        if df.empty:
            print("No genes with spatial patterns found")
            return None
        
        print(f"Found {len(df)} genes with classified spatial patterns")
        
        # Count genes by pattern type
        pattern_counts = df['pattern_type'].value_counts().reset_index()
        pattern_counts.columns = ['pattern_type', 'gene_count']
        
        print("Pattern type counts:")
        print(pattern_counts)
        
        # Visualize pattern counts
        plt.figure(figsize=(10, 6))
        sns.barplot(x='pattern_type', y='gene_count', data=pattern_counts)
        plt.title('Gene Counts by Spatial Pattern Type')
        plt.xlabel('Pattern Type')
        plt.ylabel('Number of Genes')
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, "spatial_pattern_counts.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Spatial pattern counts visualization saved to {output_file}")
        plt.close()
        
        # Visualize example genes for each pattern type
        self.visualize_pattern_examples(df)
        
        return df
    
    def visualize_pattern_examples(self, pattern_df):
        """Visualize example genes for each pattern type"""
        if pattern_df.empty:
            return
        
        # Get example genes for each pattern type
        pattern_types = pattern_df['pattern_type'].unique()
        
        for pattern_type in pattern_types:
            # Get up to 3 example genes for this pattern
            example_genes = pattern_df[pattern_df['pattern_type'] == pattern_type]['gene'].head(3).tolist()
            
            if not example_genes:
                continue
            
            # Create figure with subplots for each example gene
            fig, axes = plt.subplots(1, len(example_genes), figsize=(15, 5))
            if len(example_genes) == 1:
                axes = [axes]
            
            for i, gene in enumerate(example_genes):
                # Query for gene expression points
                query = f"""
                PREFIX : <http://example.org/spatial-transcriptomics#>
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                
                SELECT ?x ?y ?expressionLevel
                WHERE {{
                    ?geneURI rdf:type :Gene ;
                            :expressedAt ?point ;
                            :hasExpressionLevel ?expressionLevel .
                    
                    ?point :hasXCoordinate ?x ;
                           :hasYCoordinate ?y .
                    
                    # Extract gene ID
                    BIND(REPLACE(STR(?geneURI), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
                    
                    # Filter for specific gene
                    FILTER(?geneID = "{gene}")
                }}
                """
                
                results = self.graph.query(query)
                
                # Extract coordinates and expression levels
                x_coords = []
                y_coords = []
                expr_levels = []
                
                for row in results:
                    x_coords.append(int(row.x))
                    y_coords.append(int(row.y))
                    expr_levels.append(int(row.expressionLevel))
                
                if not x_coords:
                    axes[i].text(0.5, 0.5, f"No data for {gene}", 
                                ha='center', va='center', transform=axes[i].transAxes)
                    continue
                
                # Plot expression points
                scatter = axes[i].scatter(x_coords, y_coords, c=expr_levels, 
                                         cmap='viridis', s=50, alpha=0.7)
                
                # Add colorbar
                cbar = plt.colorbar(scatter, ax=axes[i])
                cbar.set_label('Expression Level')
                
                # Set title and labels
                axes[i].set_title(f"Gene {gene}")
                axes[i].set_xlabel('X Coordinate')
                axes[i].set_ylabel('Y Coordinate')
            
            # Set overall title
            fig.suptitle(f'Example Genes with {pattern_type}', fontsize=16)
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            
            # Save figure
            output_file = os.path.join(self.output_dir, f"pattern_examples_{pattern_type}.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Example visualization for {pattern_type} saved to {output_file}")
            plt.close()
    
    def analyze_enhanced_ontology(self):
        """Analyze the enhanced spatial ontology"""
        print("\nAnalyzing enhanced spatial ontology...")
        
        # Analyze directional relationships
        directional_df = self.analyze_directional_relationships()
        
        # Analyze spatial regions
        regions_results = self.analyze_spatial_regions()
        
        # Analyze distance relationships
        distance_df = self.analyze_distance_relationships()
        
        # Analyze spatial patterns
        patterns_df = self.analyze_spatial_patterns()
        
        # Generate summary report
        self.generate_summary_report(directional_df, regions_results, distance_df, patterns_df)
    
    def generate_summary_report(self, directional_df, regions_results, distance_df, patterns_df):
        """Generate a summary report of the enhanced spatial ontology analysis"""
        print("\nGenerating summary report...")
        
        report = ["# Enhanced Spatial Ontology Analysis Report\n"]
        
        # Add directional relationships summary
        report.append("## Directional Relationships\n")
        if directional_df is not None and not directional_df.empty:
            report.append("The ontology contains the following directional relationships:\n")
            for _, row in directional_df.iterrows():
                report.append(f"- {row['relation']}: {row['count']} instances\n")
        else:
            report.append("No directional relationships found in the ontology.\n")
        
        # Add spatial regions summary
        report.append("\n## Spatial Regions\n")
        if regions_results is not None and regions_results[0] is not None:
            regions_df = regions_results[0]
            report.append(f"The ontology defines {len(regions_df)} spatial regions.\n")
            
            # If we have point data
            if regions_results[1] is not None:
                points_df = regions_results[1]
                region_counts = points_df['region_name'].value_counts()
                
                report.append("Region statistics:\n")
                for region in regions_df['name']:
                    count = region_counts.get(region, 0)
                    report.append(f"- {region}: {count} points\n")
        else:
            report.append("No spatial regions found in the ontology.\n")
        
        # Add distance relationships summary
        report.append("\n## Distance Relationships\n")
        if distance_df is not None and not distance_df.empty:
            report.append("The ontology contains the following distance-based relationships:\n")
            for _, row in distance_df.iterrows():
                report.append(f"- {row['relation']}: {row['count']} instances\n")
        else:
            report.append("No distance-based relationships found in the ontology.\n")
        
        # Add spatial patterns summary
        report.append("\n## Spatial Expression Patterns\n")
        if patterns_df is not None and not patterns_df.empty:
            pattern_counts = patterns_df['pattern_type'].value_counts()
            report.append(f"The ontology classifies {len(patterns_df)} genes into spatial expression patterns:\n")
            for pattern, count in pattern_counts.items():
                report.append(f"- {pattern}: {count} genes\n")
            
            # List some example genes for each pattern
            for pattern in pattern_counts.index:
                example_genes = patterns_df[patterns_df['pattern_type'] == pattern]['gene'].head(5).tolist()
                if example_genes:
                    report.append(f"\nExample genes with {pattern}:\n")
                    for gene in example_genes:
                        report.append(f"- {gene}\n")
        else:
            report.append("No spatial expression patterns found in the ontology.\n")
        
        # Write report to file
        report_file = os.path.join(self.output_dir, "spatial_ontology_report.md")
        with open(report_file, 'w') as f:
            f.writelines(report)
        
        print(f"Summary report saved to {report_file}")

def main():
    parser = argparse.ArgumentParser(description='Spatial Pattern Analyzer')
    parser.add_argument('--ttl', default="../data/enhanced_spatial_ontology.ttl",
                        help='Path to enhanced TURTLE file')
    args = parser.parse_args()
    
    analyzer = SpatialPatternAnalyzer(args.ttl)
    analyzer.analyze_enhanced_ontology()

if __name__ == "__main__":
    main() 