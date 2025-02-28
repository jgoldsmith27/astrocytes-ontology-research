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
from sklearn.cluster import DBSCAN
from scipy.spatial import distance
import community as community_louvain
import matplotlib.cm as cm

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
GO = Namespace("http://purl.obolibrary.org/obo/GO_")

class SpatialGenePathwayAnalyzer:
    def __init__(self, ttl_file):
        """Initialize the Spatial Gene Pathway Analyzer"""
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
    
    def identify_spatial_domains(self, eps=200, min_samples=3):
        """Identify spatial domains based on gene expression patterns"""
        print("\nIdentifying spatial domains...")
        
        # Query for all spatial points
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
        
        # Extract coordinates
        points = []
        point_uris = []
        for row in results:
            x = int(row.x)
            y = int(row.y)
            points.append([x, y])
            point_uris.append(str(row.point))
        
        # Convert to numpy array
        points = np.array(points)
        
        if len(points) == 0:
            print("No spatial points found in the data")
            return None
        
        # Apply DBSCAN clustering
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        clusters = dbscan.fit_predict(points)
        
        # Create a DataFrame with results
        domains_df = pd.DataFrame({
            'point_uri': point_uris,
            'x': points[:, 0],
            'y': points[:, 1],
            'domain': clusters
        })
        
        # Count points per domain
        domain_counts = domains_df['domain'].value_counts()
        print(f"Identified {len(domain_counts)} spatial domains")
        print("Domain sizes:")
        print(domain_counts)
        
        # Visualize domains
        plt.figure(figsize=(12, 10))
        
        # Plot points colored by domain
        scatter = plt.scatter(domains_df['x'], domains_df['y'], 
                             c=domains_df['domain'], cmap='viridis', 
                             s=30, alpha=0.7)
        
        # Add colorbar
        plt.colorbar(scatter, label='Domain ID')
        
        # Set labels and title
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.title('Spatial Domains Identified by Gene Expression Patterns')
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, "spatial_domains.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Spatial domains visualization saved to {output_file}")
        plt.close()
        
        return domains_df
    
    def analyze_gene_distribution_by_domain(self, domains_df):
        """Analyze the distribution of genes across spatial domains"""
        print("\nAnalyzing gene distribution across spatial domains...")
        
        if domains_df is None or domains_df.empty:
            print("No domain information available")
            return None
        
        # Create a dictionary mapping points to domains
        point_to_domain = {}
        for _, row in domains_df.iterrows():
            point_uri = row['point_uri']
            domain = row['domain']
            point_to_domain[point_uri] = domain
        
        # Query for genes and their expression points
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?geneID ?point
        WHERE {
            ?gene rdf:type :Gene ;
                  :expressedAt ?point .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        """
        
        results = self.graph.query(query)
        
        # Group genes by domain
        gene_domain_counts = defaultdict(lambda: defaultdict(int))
        for row in results:
            gene_id = str(row.geneID)
            point_uri = str(row.point)
            
            if point_uri in point_to_domain:
                domain = point_to_domain[point_uri]
                gene_domain_counts[gene_id][domain] += 1
        
        # Convert to DataFrame for analysis
        gene_domain_data = []
        for gene_id, domain_counts in gene_domain_counts.items():
            for domain, count in domain_counts.items():
                gene_domain_data.append({
                    'gene': gene_id,
                    'domain': domain,
                    'count': count
                })
        
        gene_domain_df = pd.DataFrame(gene_domain_data)
        
        if gene_domain_df.empty:
            print("No gene-domain associations found")
            return None
        
        # Find domain-specific genes (genes predominantly expressed in one domain)
        domain_specific_genes = []
        
        for gene in gene_domain_df['gene'].unique():
            gene_data = gene_domain_df[gene_domain_df['gene'] == gene]
            total_expr = gene_data['count'].sum()
            max_domain = gene_data.loc[gene_data['count'].idxmax()]
            
            # If >70% of expression is in one domain, consider it domain-specific
            if max_domain['count'] / total_expr > 0.7:
                domain_specific_genes.append({
                    'gene': gene,
                    'primary_domain': max_domain['domain'],
                    'specificity': max_domain['count'] / total_expr
                })
        
        domain_specific_df = pd.DataFrame(domain_specific_genes)
        
        if not domain_specific_df.empty:
            print(f"Found {len(domain_specific_df)} domain-specific genes")
            
            # Group by domain and count genes
            domain_gene_counts = domain_specific_df.groupby('primary_domain').size().reset_index(name='gene_count')
            
            # Plot domain-specific gene counts
            plt.figure(figsize=(12, 6))
            sns.barplot(x='primary_domain', y='gene_count', data=domain_gene_counts)
            plt.title('Domain-Specific Genes')
            plt.xlabel('Domain ID')
            plt.ylabel('Number of Domain-Specific Genes')
            plt.tight_layout()
            
            output_file = os.path.join(self.output_dir, "domain_specific_genes.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Domain-specific genes visualization saved to {output_file}")
            plt.close()
            
            # Create heatmap of top genes across domains
            # Get top 20 genes by total count
            gene_totals = gene_domain_df.groupby('gene')['count'].sum().nlargest(20)
            top_genes = gene_totals.index.tolist()
            
            # Pivot data for heatmap
            pivot_df = gene_domain_df[gene_domain_df['gene'].isin(top_genes)].pivot_table(
                index='gene', columns='domain', values='count', fill_value=0
            )
            
            # Normalize by row (gene) for better visualization
            normalized_pivot = pivot_df.div(pivot_df.sum(axis=1), axis=0)
            
            # Plot heatmap
            plt.figure(figsize=(14, 10))
            sns.heatmap(normalized_pivot, cmap='YlGnBu', annot=False)
            plt.title('Gene Expression Distribution Across Spatial Domains')
            plt.xlabel('Domain ID')
            plt.ylabel('Gene')
            plt.tight_layout()
            
            output_file = os.path.join(self.output_dir, "gene_domain_heatmap.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Gene-domain heatmap saved to {output_file}")
            plt.close()
        
        return gene_domain_df, domain_specific_df
    
    def build_gene_co_expression_network(self):
        """Build a gene co-expression network based on spatial proximity"""
        print("\nBuilding gene co-expression network...")
        
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
        
        results = self.graph.query(query)
        
        # Count co-expressions
        co_expressions = Counter()
        for row in results:
            gene_pair = tuple(sorted([str(row.gene1ID), str(row.gene2ID)]))
            co_expressions[gene_pair] += 1
        
        # Create network graph
        G = nx.Graph()
        
        # Add edges with weights based on co-expression count
        for (gene1, gene2), count in co_expressions.items():
            G.add_edge(gene1, gene2, weight=count)
        
        print(f"Co-expression network has {G.number_of_nodes()} genes and {G.number_of_edges()} connections")
        
        # Identify communities using Louvain method
        partition = community_louvain.best_partition(G)
        
        # Count genes per community
        community_counts = Counter(partition.values())
        print(f"Identified {len(community_counts)} gene communities")
        print("Top 5 communities by size:")
        for community, count in community_counts.most_common(5):
            print(f"Community {community}: {count} genes")
        
        # Visualize network with communities
        plt.figure(figsize=(14, 14))
        
        # Limit to top genes by degree for better visualization
        if G.number_of_nodes() > 100:
            degrees = dict(G.degree())
            top_genes = sorted(degrees.keys(), key=lambda x: degrees[x], reverse=True)[:100]
            G = G.subgraph(top_genes)
            # Update partition to only include top genes
            partition = {node: partition[node] for node in G.nodes()}
        
        # Position nodes using spring layout
        pos = nx.spring_layout(G, k=0.3, iterations=50, seed=42)
        
        # Color nodes by community
        cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos, 
                              node_size=80,
                              cmap=cmap, 
                              node_color=list(partition.values()),
                              alpha=0.8)
        
        # Draw edges with transparency based on weight
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        max_weight = max(edge_weights)
        edge_alphas = [0.1 + 0.6 * (w / max_weight) for w in edge_weights]
        
        for (u, v, w), alpha in zip(G.edges(data='weight'), edge_alphas):
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=1.0, alpha=alpha)
        
        # Draw labels for high-degree nodes
        high_degree_nodes = [node for node, degree in G.degree() if degree > G.number_of_nodes() / 10]
        nx.draw_networkx_labels(G, pos, {node: node for node in high_degree_nodes}, 
                               font_size=8, font_family='sans-serif')
        
        plt.title('Gene Co-expression Network with Communities')
        plt.axis('off')
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, "gene_co_expression_communities.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Gene co-expression network saved to {output_file}")
        plt.close()
        
        return G, partition
    
    def analyze_spatial_gene_pathways(self):
        """Analyze spatial gene pathways by combining domain and network analysis"""
        print("\nAnalyzing spatial gene pathways...")
        
        # First identify spatial domains
        domains_df = self.identify_spatial_domains()
        
        # Then analyze gene distribution by domain
        gene_domain_results = self.analyze_gene_distribution_by_domain(domains_df)
        
        if gene_domain_results is None:
            print("Could not analyze gene distribution by domain")
            return
        
        gene_domain_df, domain_specific_df = gene_domain_results
        
        # Build gene co-expression network
        network_results = self.build_gene_co_expression_network()
        
        if network_results is None:
            print("Could not build gene co-expression network")
            return
        
        G, partition = network_results
        
        # Combine domain and network information
        if not domain_specific_df.empty:
            # Map communities to domains
            community_domain_map = defaultdict(Counter)
            
            for _, row in domain_specific_df.iterrows():
                gene = row['gene']
                domain = row['primary_domain']
                
                if gene in partition:
                    community = partition[gene]
                    community_domain_map[community][domain] += 1
            
            # Find primary domain for each community
            community_primary_domain = {}
            for community, domain_counts in community_domain_map.items():
                if domain_counts:
                    primary_domain = domain_counts.most_common(1)[0][0]
                    community_primary_domain[community] = primary_domain
            
            # Create a DataFrame with community-domain mapping
            community_domain_data = []
            for community, domain_counts in community_domain_map.items():
                for domain, count in domain_counts.items():
                    community_domain_data.append({
                        'community': community,
                        'domain': domain,
                        'gene_count': count
                    })
            
            community_domain_df = pd.DataFrame(community_domain_data)
            
            if not community_domain_df.empty:
                # Create a heatmap of community-domain relationships
                pivot_df = community_domain_df.pivot_table(
                    index='community', columns='domain', values='gene_count', fill_value=0
                )
                
                plt.figure(figsize=(12, 8))
                sns.heatmap(pivot_df, cmap='YlGnBu', annot=True, fmt='d')
                plt.title('Gene Community to Spatial Domain Mapping')
                plt.xlabel('Spatial Domain')
                plt.ylabel('Gene Community')
                plt.tight_layout()
                
                output_file = os.path.join(self.output_dir, "community_domain_mapping.png")
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"Community-domain mapping saved to {output_file}")
                plt.close()
                
                # Create a summary report
                print("\nSpatial Gene Pathway Summary:")
                print(f"- Identified {len(set(domains_df['domain']))} spatial domains")
                print(f"- Found {len(domain_specific_df)} domain-specific genes")
                print(f"- Detected {len(set(partition.values()))} gene co-expression communities")
                print(f"- Mapped {len(community_primary_domain)} communities to spatial domains")
                
                # For each community, list top genes and associated domain
                for community in sorted(set(partition.values())):
                    if community in community_primary_domain:
                        community_genes = [gene for gene, comm in partition.items() if comm == community]
                        top_genes = sorted(community_genes, key=lambda g: G.degree(g), reverse=True)[:5]
                        
                        print(f"\nCommunity {community} (mapped to domain {community_primary_domain[community]}):")
                        print(f"Top genes: {', '.join(top_genes)}")
                
                return {
                    'domains_df': domains_df,
                    'gene_domain_df': gene_domain_df,
                    'domain_specific_df': domain_specific_df,
                    'network': G,
                    'communities': partition,
                    'community_domain_df': community_domain_df
                }
        
        return None

def main():
    parser = argparse.ArgumentParser(description='Spatial Gene Pathway Analyzer')
    parser.add_argument('--ttl', default="../data/spatial_transcriptomics_advanced.ttl",
                        help='Path to TURTLE file')
    args = parser.parse_args()
    
    analyzer = SpatialGenePathwayAnalyzer(args.ttl)
    results = analyzer.analyze_spatial_gene_pathways()
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main() 