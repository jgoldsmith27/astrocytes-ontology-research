import rdflib
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, XSD, OWL
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display, HTML
import textwrap

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")

class SpatialQueryBuilder:
    """A class to build and execute SPARQL queries for spatial gene relationships"""
    
    def __init__(self, ttl_file):
        """Initialize with a TURTLE file"""
        self.ttl_file = ttl_file
        self.graph = self.load_turtle_file(ttl_file)
        self.output_dir = "../output/results"
    
    def load_turtle_file(self, ttl_file):
        """Load a TURTLE file into an RDF graph"""
        print(f"Loading TURTLE file: {ttl_file}")
        g = Graph()
        g.parse(ttl_file, format="turtle")
        print(f"Graph loaded with {len(g)} triples")
        return g
    
    def get_gene_list(self):
        """Get a list of all genes in the ontology"""
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT DISTINCT ?geneID
        WHERE {
            ?gene rdf:type :Gene .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        ORDER BY ?geneID
        """
        
        results = self.graph.query(query)
        return [str(row.geneID) for row in results]
    
    def find_genes_within_distance(self, target_gene, distance_threshold=200):
        """Find all genes within a specific distance of the target gene"""
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        
        SELECT ?targetGene ?nearbyGene ?targetX ?targetY ?nearbyX ?nearbyY ?distance
        WHERE {{
            # Target gene info
            ?targetGeneURI rdf:type :Gene ;
                           :expressedAt ?targetPoint .
            
            ?targetPoint :hasXCoordinate ?targetX ;
                         :hasYCoordinate ?targetY .
            
            # Nearby gene info
            ?nearbyGeneURI rdf:type :Gene ;
                           :expressedAt ?nearbyPoint .
            
            ?nearbyPoint :hasXCoordinate ?nearbyX ;
                         :hasYCoordinate ?nearbyY .
            
            # Extract gene IDs
            BIND(REPLACE(STR(?targetGeneURI), "^.*Gene_([^_]+)_.*$", "$1") AS ?targetGene)
            BIND(REPLACE(STR(?nearbyGeneURI), "^.*Gene_([^_]+)_.*$", "$1") AS ?nearbyGene)
            
            # Calculate distance
            BIND(SQRT(POW(?targetX - ?nearbyX, 2) + POW(?targetY - ?nearbyY, 2)) AS ?distance)
            
            # Filter for target gene and distance
            FILTER(?targetGene = "{target_gene}" && 
                   ?targetGeneURI != ?nearbyGeneURI &&
                   ?distance <= {distance_threshold})
        }}
        ORDER BY ?distance
        """
        
        print(f"Executing query to find genes within {distance_threshold} units of {target_gene}...")
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'target_gene': str(row.targetGene),
                'nearby_gene': str(row.nearbyGene),
                'target_x': int(row.targetX),
                'target_y': int(row.targetY),
                'nearby_x': int(row.nearbyX),
                'nearby_y': int(row.nearbyY),
                'distance': float(row.distance)
            })
        
        df = pd.DataFrame(data)
        
        if not df.empty:
            # Get unique gene pairs with minimum distance
            df = df.loc[df.groupby(['target_gene', 'nearby_gene'])['distance'].idxmin()]
        
        print(f"Found {len(df)} genes within {distance_threshold} units of {target_gene}")
        return df
    
    def find_co_expressed_genes(self, target_gene, min_count=1):
        """Find genes that are co-expressed with the target gene"""
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene2ID (COUNT(*) AS ?coExpressionCount)
        WHERE {{
            ?gene1 rdf:type :Gene ;
                   :expressedAt ?point1 .
            
            ?gene2 rdf:type :Gene ;
                   :expressedAt ?point2 .
            
            ?point1 :locatedNear ?point2 .
            
            # Extract gene IDs from URIs
            BIND(REPLACE(STR(?gene1), "^.*Gene_([^_]+)_.*$", "$1") AS ?gene1ID)
            BIND(REPLACE(STR(?gene2), "^.*Gene_([^_]+)_.*$", "$1") AS ?gene2ID)
            
            # Filter for target gene
            FILTER(?gene1ID = "{target_gene}" && ?gene1 != ?gene2)
        }}
        GROUP BY ?gene2ID
        HAVING (COUNT(*) >= {min_count})
        ORDER BY DESC(?coExpressionCount)
        """
        
        print(f"Executing query to find genes co-expressed with {target_gene}...")
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'gene': str(row.gene2ID),
                'co_expression_count': int(row.coExpressionCount)
            })
        
        df = pd.DataFrame(data)
        print(f"Found {len(df)} genes co-expressed with {target_gene}")
        return df
    
    def find_genes_in_same_cluster(self, target_gene):
        """Find genes that belong to the same cluster as the target gene"""
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?clusterID ?geneID
        WHERE {{
            # Get the cluster of the target gene
            ?targetGene rdf:type :Gene ;
                        :belongsToCluster ?cluster .
            
            # Get other genes in the same cluster
            ?gene rdf:type :Gene ;
                  :belongsToCluster ?cluster .
            
            # Extract IDs
            BIND(REPLACE(STR(?targetGene), "^.*Gene_([^_]+)_.*$", "$1") AS ?targetGeneID)
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
            BIND(REPLACE(STR(?cluster), "^.*Cluster_(.*)$", "$1") AS ?clusterID)
            
            # Filter for target gene
            FILTER(?targetGeneID = "{target_gene}" && ?targetGene != ?gene)
        }}
        ORDER BY ?geneID
        """
        
        print(f"Executing query to find genes in the same cluster as {target_gene}...")
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'cluster_id': str(row.clusterID),
                'gene': str(row.geneID)
            })
        
        df = pd.DataFrame(data)
        print(f"Found {len(df)} genes in the same cluster as {target_gene}")
        return df
    
    def find_genes_with_expression_pattern(self, min_x, max_x, min_y, max_y, min_expression=1):
        """Find genes expressed in a specific spatial region"""
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?geneID (COUNT(*) AS ?expressionCount) (AVG(?exprLevel) AS ?avgExpression)
        WHERE {{
            ?gene rdf:type :Gene ;
                  :expressedAt ?point ;
                  :hasExpressionLevel ?exprLevel .
            
            ?point :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
            
            # Filter for spatial region
            FILTER(?x >= {min_x} && ?x <= {max_x} && 
                   ?y >= {min_y} && ?y <= {max_y} &&
                   ?exprLevel >= {min_expression})
        }}
        GROUP BY ?geneID
        ORDER BY DESC(?expressionCount)
        """
        
        print(f"Executing query to find genes expressed in region ({min_x},{min_y}) to ({max_x},{max_y})...")
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'gene': str(row.geneID),
                'expression_count': int(row.expressionCount),
                'avg_expression': float(row.avgExpression)
            })
        
        df = pd.DataFrame(data)
        print(f"Found {len(df)} genes expressed in the specified region")
        return df
    
    def find_spatial_gene_relationships(self, gene1, gene2):
        """Find spatial relationships between two specific genes"""
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene1Point ?gene2Point ?x1 ?y1 ?x2 ?y2 
               (SQRT(POW(?x1 - ?x2, 2) + POW(?y1 - ?y2, 2)) AS ?distance)
        WHERE {{
            # Gene 1 info
            ?gene1URI rdf:type :Gene ;
                      :expressedAt ?gene1Point .
            
            ?gene1Point :hasXCoordinate ?x1 ;
                        :hasYCoordinate ?y1 .
            
            # Gene 2 info
            ?gene2URI rdf:type :Gene ;
                      :expressedAt ?gene2Point .
            
            ?gene2Point :hasXCoordinate ?x2 ;
                        :hasYCoordinate ?y2 .
            
            # Extract gene IDs
            BIND(REPLACE(STR(?gene1URI), "^.*Gene_([^_]+)_.*$", "$1") AS ?gene1ID)
            BIND(REPLACE(STR(?gene2URI), "^.*Gene_([^_]+)_.*$", "$1") AS ?gene2ID)
            
            # Filter for specific genes
            FILTER(?gene1ID = "{gene1}" && ?gene2ID = "{gene2}")
        }}
        ORDER BY ?distance
        LIMIT 10
        """
        
        print(f"Executing query to find spatial relationships between {gene1} and {gene2}...")
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'gene1_point': str(row.gene1Point),
                'gene2_point': str(row.gene2Point),
                'x1': int(row.x1),
                'y1': int(row.y1),
                'x2': int(row.x2),
                'y2': int(row.y2),
                'distance': float(row.distance)
            })
        
        df = pd.DataFrame(data)
        print(f"Found {len(df)} spatial relationships between {gene1} and {gene2}")
        return df
    
    def find_genes_by_expression_level(self, min_expression=5):
        """Find genes with high expression levels"""
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?geneID (MAX(?exprLevel) AS ?maxExpression) (COUNT(*) AS ?expressionCount)
        WHERE {{
            ?gene rdf:type :Gene ;
                  :hasExpressionLevel ?exprLevel .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
            
            # Filter for high expression
            FILTER(?exprLevel >= {min_expression})
        }}
        GROUP BY ?geneID
        ORDER BY DESC(?maxExpression) DESC(?expressionCount)
        """
        
        print(f"Executing query to find genes with expression level >= {min_expression}...")
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'gene': str(row.geneID),
                'max_expression': int(row.maxExpression),
                'expression_count': int(row.expressionCount)
            })
        
        df = pd.DataFrame(data)
        print(f"Found {len(df)} genes with expression level >= {min_expression}")
        return df
    
    def generate_custom_query(self, query_template, params=None):
        """Execute a custom SPARQL query with parameter substitution"""
        if params:
            query = query_template.format(**params)
        else:
            query = query_template
        
        print("Executing custom query...")
        results = self.graph.query(query)
        
        # Convert to DataFrame (assuming results have consistent columns)
        if results:
            data = []
            for row in results:
                row_dict = {}
                for var in results.vars:
                    var_name = str(var)
                    var_value = getattr(row, var_name)
                    
                    # Try to convert to appropriate type
                    try:
                        # Try as int
                        row_dict[var_name] = int(var_value)
                    except (ValueError, TypeError):
                        try:
                            # Try as float
                            row_dict[var_name] = float(var_value)
                        except (ValueError, TypeError):
                            # Keep as string
                            row_dict[var_name] = str(var_value)
                
                data.append(row_dict)
            
            df = pd.DataFrame(data)
            print(f"Query returned {len(df)} results")
            return df
        else:
            print("Query returned no results")
            return pd.DataFrame()
    
    def save_query_results(self, df, filename):
        """Save query results to CSV"""
        if df is not None and not df.empty:
            output_file = os.path.join(self.output_dir, filename)
            df.to_csv(output_file, index=False)
            print(f"Results saved to {output_file}")
    
    def visualize_query_results(self, df, plot_type='bar', x=None, y=None, title=None, filename=None):
        """Visualize query results"""
        if df is None or df.empty:
            print("No data to visualize")
            return
        
        plt.figure(figsize=(12, 8))
        
        if plot_type == 'bar':
            if x is None:
                x = df.columns[0]
            if y is None:
                y = df.columns[1] if len(df.columns) > 1 else df.columns[0]
            
            sns.barplot(x=x, y=y, data=df)
            
        elif plot_type == 'scatter':
            if x is None:
                x = df.columns[0]
            if y is None:
                y = df.columns[1] if len(df.columns) > 1 else df.columns[0]
            
            sns.scatterplot(x=x, y=y, data=df)
            
        elif plot_type == 'heatmap':
            # Assumes df is already in pivot format
            sns.heatmap(df, annot=True, cmap='YlGnBu')
        
        if title:
            plt.title(title)
        
        plt.tight_layout()
        
        if filename:
            output_file = os.path.join(self.output_dir, filename)
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Visualization saved to {output_file}")
        
        plt.close()
    
    def display_query_examples(self):
        """Display examples of SPARQL queries for spatial gene relationships"""
        examples = [
            {
                "title": "Find genes within a specific distance of a target gene",
                "query": """
                PREFIX : <http://example.org/spatial-transcriptomics#>
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                
                SELECT ?targetGene ?nearbyGene 
                       ?targetX ?targetY 
                       ?nearbyX ?nearbyY 
                       (SQRT(POW(?targetX - ?nearbyX, 2) + POW(?targetY - ?nearbyY, 2)) AS ?distance)
                WHERE {
                    # Target gene info
                    ?targetGeneURI rdf:type :Gene ;
                                   :expressedAt ?targetPoint .
                    
                    ?targetPoint :hasXCoordinate ?targetX ;
                                 :hasYCoordinate ?targetY .
                    
                    # Nearby gene info
                    ?nearbyGeneURI rdf:type :Gene ;
                                   :expressedAt ?nearbyPoint .
                    
                    ?nearbyPoint :hasXCoordinate ?nearbyX ;
                                 :hasYCoordinate ?nearbyY .
                    
                    # Extract gene IDs
                    BIND(REPLACE(STR(?targetGeneURI), "^.*Gene_([^_]+)_.*$", "$1") AS ?targetGene)
                    BIND(REPLACE(STR(?nearbyGeneURI), "^.*Gene_([^_]+)_.*$", "$1") AS ?nearbyGene)
                    
                    # Filter for target gene and distance
                    FILTER(?targetGene = "MBP" && 
                           ?targetGeneURI != ?nearbyGeneURI &&
                           SQRT(POW(?targetX - ?nearbyX, 2) + POW(?targetY - ?nearbyY, 2)) <= 200)
                }
                ORDER BY ?distance
                """
            },
            {
                "title": "Find genes co-expressed in nearby spatial points",
                "query": """
                PREFIX : <http://example.org/spatial-transcriptomics#>
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                
                SELECT ?gene1ID ?gene2ID (COUNT(*) AS ?coExpressionCount)
                WHERE {
                    ?gene1 rdf:type :Gene ;
                           :expressedAt ?point1 .
                    
                    ?gene2 rdf:type :Gene ;
                           :expressedAt ?point2 .
                    
                    ?point1 :locatedNear ?point2 .
                    
                    # Extract gene IDs from URIs
                    BIND(REPLACE(STR(?gene1), "^.*Gene_([^_]+)_.*$", "$1") AS ?gene1ID)
                    BIND(REPLACE(STR(?gene2), "^.*Gene_([^_]+)_.*$", "$1") AS ?gene2ID)
                    
                    # Ensure we don't count the same gene twice
                    FILTER(?gene1 != ?gene2)
                }
                GROUP BY ?gene1ID ?gene2ID
                HAVING (COUNT(*) >= 3)
                ORDER BY DESC(?coExpressionCount)
                """
            },
            {
                "title": "Find genes expressed in a specific spatial region",
                "query": """
                PREFIX : <http://example.org/spatial-transcriptomics#>
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                
                SELECT ?geneID (COUNT(*) AS ?expressionCount) (AVG(?exprLevel) AS ?avgExpression)
                WHERE {
                    ?gene rdf:type :Gene ;
                          :expressedAt ?point ;
                          :hasExpressionLevel ?exprLevel .
                    
                    ?point :hasXCoordinate ?x ;
                           :hasYCoordinate ?y .
                    
                    # Extract gene ID from URI
                    BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
                    
                    # Filter for spatial region
                    FILTER(?x >= 5000 && ?x <= 6000 && 
                           ?y >= 8000 && ?y <= 9000 &&
                           ?exprLevel >= 1)
                }
                GROUP BY ?geneID
                ORDER BY DESC(?expressionCount)
                """
            },
            {
                "title": "Find genes that belong to the same cluster",
                "query": """
                PREFIX : <http://example.org/spatial-transcriptomics#>
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                
                SELECT ?clusterID (COUNT(?gene) AS ?geneCount)
                WHERE {
                    ?cluster rdf:type :GeneCluster .
                    ?gene rdf:type :Gene ;
                          :belongsToCluster ?cluster .
                    
                    # Extract cluster ID
                    BIND(REPLACE(STR(?cluster), "^.*Cluster_(.*)$", "$1") AS ?clusterID)
                }
                GROUP BY ?clusterID
                ORDER BY DESC(?geneCount)
                """
            }
        ]
        
        print("\n=== SPARQL Query Examples for Spatial Gene Relationships ===\n")
        
        for i, example in enumerate(examples, 1):
            print(f"{i}. {example['title']}")
            print("-" * 80)
            print(textwrap.dedent(example['query']).strip())
            print("\n" + "=" * 80 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Spatial Query Builder')
    parser.add_argument('--ttl', default="../data/spatial_transcriptomics_advanced.ttl",
                        help='Path to TURTLE file')
    parser.add_argument('--gene', help='Target gene for queries')
    parser.add_argument('--distance', type=float, default=200, help='Distance threshold for spatial queries')
    parser.add_argument('--examples', action='store_true', help='Display example SPARQL queries')
    args = parser.parse_args()
    
    # Initialize query builder
    query_builder = SpatialQueryBuilder(args.ttl)
    
    # Display examples if requested
    if args.examples:
        query_builder.display_query_examples()
        return
    
    # If gene is specified, run some example queries
    if args.gene:
        print(f"\nRunning example queries for gene: {args.gene}\n")
        
        # Find nearby genes
        nearby_genes = query_builder.find_genes_within_distance(args.gene, args.distance)
        if not nearby_genes.empty:
            query_builder.save_query_results(nearby_genes, f"{args.gene}_nearby_genes.csv")
            query_builder.visualize_query_results(
                nearby_genes.head(15), 
                plot_type='bar', 
                x='nearby_gene', 
                y='distance',
                title=f'Genes within {args.distance} units of {args.gene}',
                filename=f"{args.gene}_nearby_genes.png"
            )
        
        # Find co-expressed genes
        co_expressed = query_builder.find_co_expressed_genes(args.gene)
        if not co_expressed.empty:
            query_builder.save_query_results(co_expressed, f"{args.gene}_co_expressed.csv")
            query_builder.visualize_query_results(
                co_expressed.head(15), 
                plot_type='bar', 
                x='gene', 
                y='co_expression_count',
                title=f'Genes co-expressed with {args.gene}',
                filename=f"{args.gene}_co_expressed.png"
            )
        
        # Find genes in same cluster
        cluster_genes = query_builder.find_genes_in_same_cluster(args.gene)
        if not cluster_genes.empty:
            query_builder.save_query_results(cluster_genes, f"{args.gene}_cluster_genes.csv")
            print(f"\nGenes in the same cluster as {args.gene}:")
            print(cluster_genes.head(10))
    else:
        # If no gene specified, just list available genes
        genes = query_builder.get_gene_list()
        print("\nAvailable genes in the dataset:")
        for i, gene in enumerate(genes[:20]):
            print(f"  {gene}", end=", " if (i+1) % 5 != 0 else "\n")
        
        if len(genes) > 20:
            print(f"... and {len(genes) - 20} more")
        
        print("\nUse --gene parameter to specify a target gene for analysis")
        print("Use --examples to see example SPARQL queries")

if __name__ == "__main__":
    main() 