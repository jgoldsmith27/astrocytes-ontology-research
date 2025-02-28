import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.spatial import KDTree, distance
import rdflib
from rdflib import Graph, Namespace, URIRef, Literal, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL
import seaborn as sns
from sklearn.cluster import DBSCAN
from collections import Counter
import networkx as nx
import argparse

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
GO = Namespace("http://purl.obolibrary.org/obo/GO_")
UBERON = Namespace("http://purl.obolibrary.org/obo/UBERON_")

class SpatialTranscriptomicsOntology:
    def __init__(self, csv_file=None, ttl_file=None):
        """
        Initialize the Spatial Transcriptomics Ontology
        
        Args:
            csv_file: Path to the CSV file with spatial transcriptomics data
            ttl_file: Path to an existing TURTLE file to load
        """
        self.csv_file = csv_file
        self.ttl_file = ttl_file
        self.df = None
        self.graph = Graph()
        
        # Bind namespaces
        self.graph.bind("", ST)
        self.graph.bind("rdf", RDF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind("owl", OWL)
        self.graph.bind("xsd", XSD)
        self.graph.bind("go", GO)
        self.graph.bind("uberon", UBERON)
        
        # Initialize ontology if loading from CSV
        if csv_file:
            self._init_ontology()
        
        # Load existing TURTLE file if provided
        if ttl_file:
            self.load_turtle(ttl_file)
    
    def _init_ontology(self):
        """Initialize the ontology with classes and properties"""
        # Add classes
        self.graph.add((ST.Gene, RDF.type, OWL.Class))
        self.graph.add((ST.SpatialPoint, RDF.type, OWL.Class))
        self.graph.add((ST.GeneCluster, RDF.type, OWL.Class))
        self.graph.add((ST.SpatialRegion, RDF.type, OWL.Class))
        
        # Add object properties
        self.graph.add((ST.expressedAt, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.expressedAt, RDFS.domain, ST.Gene))
        self.graph.add((ST.expressedAt, RDFS.range, ST.SpatialPoint))
        
        self.graph.add((ST.locatedNear, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.locatedNear, RDFS.domain, ST.SpatialPoint))
        self.graph.add((ST.locatedNear, RDFS.range, ST.SpatialPoint))
        
        self.graph.add((ST.belongsToCluster, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.belongsToCluster, RDFS.domain, ST.Gene))
        self.graph.add((ST.belongsToCluster, RDFS.range, ST.GeneCluster))
        
        self.graph.add((ST.locatedIn, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.locatedIn, RDFS.domain, ST.SpatialPoint))
        self.graph.add((ST.locatedIn, RDFS.range, ST.SpatialRegion))
        
        # Add datatype properties
        self.graph.add((ST.hasExpressionLevel, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasExpressionLevel, RDFS.domain, ST.Gene))
        self.graph.add((ST.hasExpressionLevel, RDFS.range, XSD.integer))
        
        self.graph.add((ST.hasExonCount, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasExonCount, RDFS.domain, ST.Gene))
        self.graph.add((ST.hasExonCount, RDFS.range, XSD.integer))
        
        self.graph.add((ST.hasIntronCount, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasIntronCount, RDFS.domain, ST.Gene))
        self.graph.add((ST.hasIntronCount, RDFS.range, XSD.integer))
        
        self.graph.add((ST.hasXCoordinate, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasXCoordinate, RDFS.domain, ST.SpatialPoint))
        self.graph.add((ST.hasXCoordinate, RDFS.range, XSD.integer))
        
        self.graph.add((ST.hasYCoordinate, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasYCoordinate, RDFS.domain, ST.SpatialPoint))
        self.graph.add((ST.hasYCoordinate, RDFS.range, XSD.integer))
        
        self.graph.add((ST.hasClusterSize, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasClusterSize, RDFS.domain, ST.GeneCluster))
        self.graph.add((ST.hasClusterSize, RDFS.range, XSD.integer))
    
    def load_data(self, sample_size=None):
        """
        Load data from CSV file
        
        Args:
            sample_size: Number of rows to sample (None for all)
        """
        print(f"Loading data from {self.csv_file}...")
        
        if sample_size:
            # Sample the dataset
            self.df = pd.read_csv(self.csv_file, skiprows=lambda x: x > 0 and np.random.random() > sample_size/9000000)
            if len(self.df) > sample_size:
                self.df = self.df.sample(sample_size, random_state=42)
        else:
            # Load the entire dataset
            self.df = pd.read_csv(self.csv_file)
        
        print(f"Loaded {len(self.df)} records")
        print(f"Columns: {self.df.columns.tolist()}")
        
        return self.df
    
    def load_turtle(self, ttl_file):
        """Load an existing TURTLE file"""
        print(f"Loading TURTLE file: {ttl_file}")
        self.graph.parse(ttl_file, format="turtle")
        print(f"Loaded {len(self.graph)} triples")
    
    def save_turtle(self, output_file):
        """Save the graph to a TURTLE file"""
        print(f"Saving {len(self.graph)} triples to {output_file}...")
        self.graph.serialize(destination=output_file, format="turtle")
        print(f"TURTLE file saved: {output_file}")
    
    def add_gene_instances(self):
        """Add gene instances to the graph"""
        if self.df is None:
            raise ValueError("No data loaded. Call load_data() first.")
        
        print("Adding gene instances to the graph...")
        
        # Dictionary to store points and their coordinates
        points = {}
        
        # Add gene and point instances
        for index, row in tqdm(self.df.iterrows(), total=len(self.df)):
            gene_id = row['geneID']
            x, y = int(row['x']), int(row['y'])
            expr_level = int(row['MIDCount'])
            exon_count = int(row['ExonCount'])
            intron_count = int(row['IntronCount'])
            
            # Create URIs
            gene_uri = URIRef(f"{ST}Gene_{gene_id}_{index}")
            point_uri = URIRef(f"{ST}Point_{x}_{y}")
            
            # Store point coordinates
            points[point_uri] = (x, y)
            
            # Add gene instance
            self.graph.add((gene_uri, RDF.type, ST.Gene))
            self.graph.add((gene_uri, ST.expressedAt, point_uri))
            self.graph.add((gene_uri, ST.hasExpressionLevel, Literal(expr_level, datatype=XSD.integer)))
            self.graph.add((gene_uri, ST.hasExonCount, Literal(exon_count, datatype=XSD.integer)))
            self.graph.add((gene_uri, ST.hasIntronCount, Literal(intron_count, datatype=XSD.integer)))
            
            # Add point instance if not already added
            if not (point_uri, RDF.type, ST.SpatialPoint) in self.graph:
                self.graph.add((point_uri, RDF.type, ST.SpatialPoint))
                self.graph.add((point_uri, ST.hasXCoordinate, Literal(x, datatype=XSD.integer)))
                self.graph.add((point_uri, ST.hasYCoordinate, Literal(y, datatype=XSD.integer)))
        
        return points
    
    def add_spatial_relationships(self, points, proximity_threshold=100, max_neighbors=3):
        """
        Add spatial relationships between points
        
        Args:
            points: Dictionary of point URIs and their coordinates
            proximity_threshold: Distance threshold for "locatedNear" relationship
            max_neighbors: Maximum number of neighbors to include for each point
        """
        print("Adding spatial relationships...")
        
        # Extract coordinates and point URIs
        coords = np.array(list(points.values()))
        point_uris = list(points.keys())
        
        # Build KDTree for efficient nearest neighbor search
        tree = KDTree(coords)
        
        # Query for neighbors within threshold
        for i, point_uri in enumerate(tqdm(point_uris)):
            # Find neighbors within threshold
            neighbors = tree.query_ball_point(coords[i], proximity_threshold)
            
            # Remove self from neighbors
            neighbors = [j for j in neighbors if j != i]
            
            # Limit number of neighbors
            neighbors = neighbors[:max_neighbors]
            
            # Add relationships to graph
            for neighbor_idx in neighbors:
                neighbor_uri = point_uris[neighbor_idx]
                self.graph.add((point_uri, ST.locatedNear, neighbor_uri))
    
    def identify_gene_clusters(self, eps=200, min_samples=3):
        """
        Identify clusters of genes using DBSCAN
        
        Args:
            eps: The maximum distance between two samples for them to be considered as in the same neighborhood
            min_samples: The number of samples in a neighborhood for a point to be considered as a core point
        """
        print("Identifying gene clusters...")
        
        # Get all spatial points
        points_query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y
        WHERE {
            ?point rdf:type :SpatialPoint ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
        }
        """
        
        results = self.graph.query(points_query)
        
        # Extract coordinates
        points = []
        point_uris = []
        
        for row in results:
            x, y = int(row.x), int(row.y)
            points.append([x, y])
            point_uris.append(row.point)
        
        # Convert to numpy array
        points = np.array(points)
        
        # Apply DBSCAN clustering
        clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(points)
        labels = clustering.labels_
        
        # Count clusters
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        print(f"Found {n_clusters} clusters")
        
        # Create cluster instances
        cluster_counts = Counter(labels)
        
        for cluster_id in set(labels):
            if cluster_id == -1:  # Skip noise points
                continue
            
            # Create cluster URI
            cluster_uri = URIRef(f"{ST}Cluster_{cluster_id}")
            
            # Add cluster instance
            self.graph.add((cluster_uri, RDF.type, ST.GeneCluster))
            self.graph.add((cluster_uri, ST.hasClusterSize, Literal(cluster_counts[cluster_id], datatype=XSD.integer)))
            
            # Add points to cluster
            for i, label in enumerate(labels):
                if label == cluster_id:
                    point_uri = point_uris[i]
                    
                    # Find genes at this point
                    gene_query = f"""
                    PREFIX : <http://example.org/spatial-transcriptomics#>
                    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                    
                    SELECT ?gene
                    WHERE {{
                        ?gene rdf:type :Gene ;
                              :expressedAt <{point_uri}> .
                    }}
                    """
                    
                    gene_results = self.graph.query(gene_query)
                    
                    # Add genes to cluster
                    for gene_row in gene_results:
                        self.graph.add((gene_row.gene, ST.belongsToCluster, cluster_uri))
        
        return labels, points, point_uris
    
    def visualize_clusters(self, labels, points, output_file=None):
        """
        Visualize gene clusters
        
        Args:
            labels: Cluster labels from DBSCAN
            points: Array of point coordinates
            output_file: Path to save the visualization
        """
        plt.figure(figsize=(12, 10))
        
        # Create a colormap for clusters
        unique_labels = set(labels)
        colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_labels)))
        
        # Plot points colored by cluster
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise
                col = [0, 0, 0, 1]
            
            class_member_mask = (labels == k)
            
            xy = points[class_member_mask]
            plt.scatter(xy[:, 0], xy[:, 1], s=30, c=[col], label=f"Cluster {k}" if k != -1 else "Noise")
        
        plt.title(f'Spatial Gene Clusters (DBSCAN)')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Cluster visualization saved to {output_file}")
        
        plt.show()
    
    def visualize_gene_network(self, output_file=None, max_points=100):
        """
        Visualize the gene expression network
        
        Args:
            output_file: Path to save the visualization
            max_points: Maximum number of points to visualize
        """
        # Query for spatial relationships
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y ?nearPoint ?nearX ?nearY
        WHERE {
            ?point rdf:type :SpatialPoint ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y ;
                   :locatedNear ?nearPoint .
            
            ?nearPoint :hasXCoordinate ?nearX ;
                       :hasYCoordinate ?nearY .
        }
        """
        
        results = self.graph.query(query)
        
        # Create a network graph
        G = nx.Graph()
        
        # Add nodes and edges
        for row in results:
            x1, y1 = int(row.x), int(row.y)
            x2, y2 = int(row.nearX), int(row.nearY)
            
            # Add nodes with position attributes
            G.add_node(f"({x1},{y1})", pos=(x1, y1))
            G.add_node(f"({x2},{y2})", pos=(x2, y2))
            
            # Add edge
            G.add_edge(f"({x1},{y1})", f"({x2},{y2})")
        
        # Limit the number of nodes if needed
        if len(G) > max_points:
            # Sample nodes
            nodes = list(G.nodes())
            sampled_nodes = np.random.choice(nodes, max_points, replace=False)
            G = G.subgraph(sampled_nodes)
        
        # Get node positions
        pos = nx.get_node_attributes(G, 'pos')
        
        # Create the plot
        plt.figure(figsize=(12, 10))
        
        # Draw the network
        nx.draw(G, pos, node_size=30, node_color='blue', alpha=0.7, 
                with_labels=False, width=0.5, edge_color='gray', edge_alpha=0.3)
        
        plt.title(f'Spatial Gene Expression Network (showing {len(G)} points)')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Network visualization saved to {output_file}")
        
        plt.show()
    
    def analyze_gene_distribution(self, output_file=None):
        """
        Analyze the distribution of genes
        
        Args:
            output_file: Path to save the visualization
        """
        # Query for genes and their expression levels
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?geneID (COUNT(?gene) as ?count) (SUM(?expressionLevel) as ?totalExpression)
        WHERE {
            ?gene rdf:type :Gene ;
                  :hasExpressionLevel ?expressionLevel .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        GROUP BY ?geneID
        ORDER BY DESC(?totalExpression)
        """
        
        results = self.graph.query(query)
        
        # Convert to DataFrame
        gene_data = []
        for row in results:
            gene_data.append({
                'geneID': str(row.geneID),
                'count': int(row.count),
                'totalExpression': int(row.totalExpression)
            })
        
        gene_df = pd.DataFrame(gene_data)
        
        # Plot the distribution
        plt.figure(figsize=(14, 8))
        
        # Top 20 genes by total expression
        top_genes = gene_df.sort_values('totalExpression', ascending=False).head(20)
        
        sns.barplot(x='geneID', y='totalExpression', data=top_genes)
        plt.title('Top 20 Genes by Total Expression Level')
        plt.xlabel('Gene ID')
        plt.ylabel('Total Expression Level')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Gene distribution visualization saved to {output_file}")
        
        plt.show()
        
        return gene_df

def main():
    parser = argparse.ArgumentParser(description='Spatial Transcriptomics Ontology Builder')
    parser.add_argument('--csv', type=str, help='Path to the CSV file')
    parser.add_argument('--ttl', type=str, help='Path to the TURTLE file (for loading or saving)')
    parser.add_argument('--sample', type=int, default=200, help='Number of rows to sample')
    parser.add_argument('--proximity', type=int, default=100, help='Distance threshold for spatial relationships')
    parser.add_argument('--neighbors', type=int, default=3, help='Maximum number of neighbors per point')
    parser.add_argument('--eps', type=int, default=200, help='DBSCAN eps parameter')
    parser.add_argument('--min_samples', type=int, default=3, help='DBSCAN min_samples parameter')
    parser.add_argument('--output_dir', type=str, help='Directory to save output files')
    
    args = parser.parse_args()
    
    # Set default CSV file if not provided
    if not args.csv:
        args.csv = "/Users/jacob/Desktop/Karolinska Research/astrocytes-ontology-research/Spatial_Data/1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv"
    
    # Set default output directory if not provided
    if not args.output_dir:
        args.output_dir = os.path.dirname(args.csv)
    
    # Set default TURTLE file if not provided
    if not args.ttl:
        args.ttl = os.path.join(args.output_dir, "spatial_transcriptomics_advanced.ttl")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize ontology
    ontology = SpatialTranscriptomicsOntology(csv_file=args.csv)
    
    # Load data
    ontology.load_data(sample_size=args.sample)
    
    # Add gene instances
    points = ontology.add_gene_instances()
    
    # Add spatial relationships
    ontology.add_spatial_relationships(points, proximity_threshold=args.proximity, max_neighbors=args.neighbors)
    
    # Identify gene clusters
    labels, points_array, point_uris = ontology.identify_gene_clusters(eps=args.eps, min_samples=args.min_samples)
    
    # Save ontology to TURTLE file
    ontology.save_turtle(args.ttl)
    
    # Visualize clusters
    cluster_viz_file = os.path.join(args.output_dir, "gene_clusters.png")
    ontology.visualize_clusters(labels, points_array, output_file=cluster_viz_file)
    
    # Visualize gene network
    network_viz_file = os.path.join(args.output_dir, "gene_network.png")
    ontology.visualize_gene_network(output_file=network_viz_file)
    
    # Analyze gene distribution
    dist_viz_file = os.path.join(args.output_dir, "gene_distribution.png")
    gene_df = ontology.analyze_gene_distribution(output_file=dist_viz_file)
    
    print("Analysis complete!")

if __name__ == "__main__":
    main() 