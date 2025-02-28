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

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")

def load_turtle_file(ttl_file):
    """Load a TURTLE file into an RDF graph"""
    print(f"Loading TURTLE file: {ttl_file}")
    g = Graph()
    g.parse(ttl_file, format="turtle")
    print(f"Graph loaded with {len(g)} triples")
    return g

def analyze_clusters(g):
    """Analyze gene clusters"""
    print("\nAnalyzing gene clusters...")
    
    # Query for clusters and their sizes
    query = """
    PREFIX : <http://example.org/spatial-transcriptomics#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?cluster ?clusterSize
    WHERE {
        ?cluster rdf:type :GeneCluster ;
                 :hasClusterSize ?clusterSize .
    }
    ORDER BY DESC(?clusterSize)
    """
    
    results = g.query(query)
    
    # Convert to DataFrame
    clusters = []
    for row in results:
        cluster_id = str(row.cluster).split('_')[-1]
        clusters.append({
            'cluster_id': cluster_id,
            'size': int(row.clusterSize)
        })
    
    clusters_df = pd.DataFrame(clusters)
    
    if not clusters_df.empty:
        print(f"Found {len(clusters_df)} clusters")
        print("Top 5 clusters by size:")
        print(clusters_df.head(5))
        
        # Plot cluster sizes
        plt.figure(figsize=(12, 6))
        sns.barplot(x='cluster_id', y='size', data=clusters_df)
        plt.title('Gene Cluster Sizes')
        plt.xlabel('Cluster ID')
        plt.ylabel('Number of Points')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        output_file = os.path.join(os.path.dirname(ttl_file), "cluster_sizes.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Cluster size visualization saved to {output_file}")
        plt.close()
    else:
        print("No clusters found in the data")
    
    return clusters_df

def analyze_gene_distribution_by_cluster(g):
    """Analyze the distribution of genes within clusters"""
    print("\nAnalyzing gene distribution by cluster...")
    
    # Query for genes and their clusters
    query = """
    PREFIX : <http://example.org/spatial-transcriptomics#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?geneID ?cluster
    WHERE {
        ?gene rdf:type :Gene ;
              :belongsToCluster ?cluster .
        
        # Extract gene ID from URI
        BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
    }
    """
    
    results = g.query(query)
    
    # Group genes by cluster
    cluster_genes = defaultdict(list)
    for row in results:
        cluster_id = str(row.cluster).split('_')[-1]
        gene_id = str(row.geneID)
        cluster_genes[cluster_id].append(gene_id)
    
    # Count genes per cluster
    gene_counts = {}
    for cluster_id, genes in cluster_genes.items():
        gene_counts[cluster_id] = Counter(genes)
    
    # Print summary
    print(f"Found genes in {len(cluster_genes)} clusters")
    
    # Find the most common genes in each cluster
    for cluster_id, gene_counter in gene_counts.items():
        if gene_counter:
            most_common = gene_counter.most_common(3)
            print(f"Cluster {cluster_id}: Most common genes - {', '.join([f'{gene} ({count})' for gene, count in most_common])}")
    
    # Create a heatmap of top genes across clusters
    all_genes = set()
    for counter in gene_counts.values():
        all_genes.update(counter.keys())
    
    # Get top 20 genes by total count
    gene_totals = Counter()
    for counter in gene_counts.values():
        gene_totals.update(counter)
    
    top_genes = [gene for gene, _ in gene_totals.most_common(20)]
    
    # Create matrix for heatmap
    heatmap_data = []
    for cluster_id in sorted(cluster_genes.keys()):
        row = []
        counter = gene_counts[cluster_id]
        for gene in top_genes:
            row.append(counter.get(gene, 0))
        heatmap_data.append(row)
    
    if heatmap_data and top_genes:
        # Plot heatmap
        plt.figure(figsize=(14, 10))
        sns.heatmap(heatmap_data, 
                   xticklabels=top_genes, 
                   yticklabels=sorted(cluster_genes.keys()),
                   cmap="YlGnBu",
                   annot=True,
                   fmt="d")
        plt.title('Gene Distribution Across Clusters')
        plt.xlabel('Gene')
        plt.ylabel('Cluster ID')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        output_file = os.path.join(os.path.dirname(ttl_file), "gene_cluster_heatmap.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Gene cluster heatmap saved to {output_file}")
        plt.close()
    
    return gene_counts

def analyze_spatial_network(g):
    """Analyze the spatial network of gene expression"""
    print("\nAnalyzing spatial network...")
    
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
    
    results = g.query(query)
    
    # Create network graph
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
    
    # Calculate network metrics
    print(f"Network has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    # Calculate connected components
    components = list(nx.connected_components(G))
    print(f"Network has {len(components)} connected components")
    
    # Calculate degree distribution
    degrees = [d for _, d in G.degree()]
    degree_counts = Counter(degrees)
    
    # Plot degree distribution
    plt.figure(figsize=(10, 6))
    plt.bar(degree_counts.keys(), degree_counts.values())
    plt.title('Degree Distribution of Spatial Network')
    plt.xlabel('Degree')
    plt.ylabel('Count')
    plt.tight_layout()
    
    output_file = os.path.join(os.path.dirname(ttl_file), "degree_distribution.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Degree distribution saved to {output_file}")
    plt.close()
    
    return G

def find_co_expressed_genes(g):
    """Find genes that are co-expressed at nearby spatial points"""
    print("\nFinding co-expressed genes...")
    
    # Query for genes expressed at points that are near each other
    query = """
    PREFIX : <http://example.org/spatial-transcriptomics#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?gene1ID ?gene2ID ?point1 ?point2
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
    """
    
    results = g.query(query)
    
    # Count co-expressions
    co_expressions = Counter()
    for row in results:
        gene_pair = tuple(sorted([str(row.gene1ID), str(row.gene2ID)]))
        co_expressions[gene_pair] += 1
    
    # Print top co-expressed gene pairs
    print("Top co-expressed gene pairs:")
    for (gene1, gene2), count in co_expressions.most_common(10):
        print(f"{gene1} and {gene2}: {count} co-expressions")
    
    # Create co-expression network
    if co_expressions:
        G = nx.Graph()
        
        # Add edges with weights based on co-expression count
        for (gene1, gene2), count in co_expressions.items():
            G.add_edge(gene1, gene2, weight=count)
        
        # Get top 30 genes by degree for visualization
        degrees = dict(G.degree())
        top_genes = sorted(degrees.keys(), key=lambda x: degrees[x], reverse=True)[:30]
        
        # Create subgraph with only top genes
        H = G.subgraph(top_genes)
        
        # Plot co-expression network
        plt.figure(figsize=(12, 12))
        
        # Use spring layout for positioning
        pos = nx.spring_layout(H, seed=42)
        
        # Get edge weights for line thickness
        edge_weights = [G[u][v]['weight'] for u, v in H.edges()]
        
        # Normalize weights for visualization
        max_weight = max(edge_weights)
        normalized_weights = [2 + 8 * (w / max_weight) for w in edge_weights]
        
        # Draw the network
        nx.draw_networkx_nodes(H, pos, node_size=500, node_color='lightblue', alpha=0.8)
        nx.draw_networkx_edges(H, pos, width=normalized_weights, alpha=0.5, edge_color='gray')
        nx.draw_networkx_labels(H, pos, font_size=10, font_family='sans-serif')
        
        plt.title('Gene Co-expression Network')
        plt.axis('off')
        plt.tight_layout()
        
        output_file = os.path.join(os.path.dirname(ttl_file), "co_expression_network.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Co-expression network saved to {output_file}")
        plt.close()
    
    return co_expressions

def find_expression_patterns(g):
    """Find patterns in gene expression across spatial regions"""
    print("\nAnalyzing spatial expression patterns...")
    
    # Query for genes, their expression levels, and coordinates
    query = """
    PREFIX : <http://example.org/spatial-transcriptomics#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?geneID ?expressionLevel ?x ?y
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
    
    results = g.query(query)
    
    # Convert to DataFrame
    expression_data = []
    for row in results:
        expression_data.append({
            'geneID': str(row.geneID),
            'expressionLevel': int(row.expressionLevel),
            'x': int(row.x),
            'y': int(row.y)
        })
    
    df = pd.DataFrame(expression_data)
    
    if not df.empty:
        # Get top 5 genes by total expression
        top_genes = df.groupby('geneID')['expressionLevel'].sum().nlargest(5).index.tolist()
        
        # Create a spatial heatmap for each top gene
        for gene in top_genes:
            gene_df = df[df['geneID'] == gene]
            
            plt.figure(figsize=(10, 8))
            
            # Create a 2D histogram
            heatmap, xedges, yedges = np.histogram2d(
                gene_df['x'], 
                gene_df['y'], 
                bins=20, 
                weights=gene_df['expressionLevel']
            )
            
            # Plot heatmap
            plt.imshow(heatmap.T, origin='lower', aspect='auto', 
                      extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                      cmap='viridis')
            
            plt.colorbar(label='Expression Level')
            plt.title(f'Spatial Expression Pattern of {gene}')
            plt.xlabel('X Coordinate')
            plt.ylabel('Y Coordinate')
            plt.tight_layout()
            
            output_file = os.path.join(os.path.dirname(ttl_file), f"expression_pattern_{gene}.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Expression pattern for {gene} saved to {output_file}")
            plt.close()
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Advanced Queries for Spatial Transcriptomics')
    parser.add_argument('--ttl', default="../data/spatial_transcriptomics_advanced.ttl",
                        help='Path to TURTLE file')
    args = parser.parse_args()
    
    # Load the TURTLE file
    g = load_turtle_file(args.ttl)
    
    # Analyze clusters
    clusters_df = analyze_clusters(g)
    
    # Analyze gene distribution by cluster
    gene_counts = analyze_gene_distribution_by_cluster(g)
    
    # Analyze spatial network
    G = analyze_spatial_network(g)
    
    # Find co-expressed genes
    co_expressions = find_co_expressed_genes(g)
    
    # Find expression patterns
    expression_df = find_expression_patterns(g)
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main() 