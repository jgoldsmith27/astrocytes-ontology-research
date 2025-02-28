#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF, RDFS
import networkx as nx
from community import community_louvain
from sklearn.metrics.pairwise import euclidean_distances
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict

class SpatialNetworkAnalyzer:
    """
    A class to perform advanced network analysis on spatial gene relationships.
    """
    
    def __init__(self, ttl_file, output_dir="../output/network_analysis"):
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
    
    def build_spatial_proximity_network(self, distance_threshold=100):
        """
        Build a network of spatial points based on proximity.
        
        Parameters:
        -----------
        distance_threshold : float
            Maximum distance for two points to be considered connected
            
        Returns:
        --------
        nx.Graph
            NetworkX graph of spatial points
        """
        print("\nBuilding spatial proximity network...")
        
        # Get spatial points
        points_df = self.get_spatial_points()
        
        if points_df.empty:
            print("No spatial points found")
            return None
        
        # Create a graph
        G = nx.Graph()
        
        # Add nodes (spatial points)
        for _, row in points_df.iterrows():
            G.add_node(row['point'], x=row['x'], y=row['y'])
        
        # Add edges based on proximity
        for i, row1 in points_df.iterrows():
            for j, row2 in points_df.iterrows():
                if i < j:  # Avoid duplicate edges
                    # Calculate Euclidean distance
                    dist = np.sqrt((row1['x'] - row2['x'])**2 + (row1['y'] - row2['y'])**2)
                    
                    # Add edge if distance is below threshold
                    if dist <= distance_threshold:
                        G.add_edge(row1['point'], row2['point'], weight=1.0/dist if dist > 0 else 1.0)
        
        print(f"Created network with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        
        return G
    
    def build_gene_co_expression_network(self, correlation_threshold=0.7):
        """
        Build a network of genes based on spatial co-expression patterns.
        
        Parameters:
        -----------
        correlation_threshold : float
            Minimum correlation coefficient for two genes to be considered connected
            
        Returns:
        --------
        nx.Graph
            NetworkX graph of genes
        """
        print("\nBuilding gene co-expression network...")
        
        # Get genes and their expression points
        genes_df = self.get_genes()
        
        if genes_df.empty:
            print("No gene expression data found")
            return None
        
        # Create a pivot table: genes x points
        pivot_df = genes_df.pivot_table(
            index='gene', 
            columns='point', 
            values='expression_level',
            fill_value=0
        )
        
        # Calculate correlation matrix
        corr_matrix = pivot_df.T.corr(method='spearman')
        
        # Create a graph
        G = nx.Graph()
        
        # Add nodes (genes)
        for gene in pivot_df.index:
            G.add_node(gene)
        
        # Add edges based on correlation
        for i, gene1 in enumerate(corr_matrix.index):
            for j, gene2 in enumerate(corr_matrix.columns):
                if i < j:  # Avoid duplicate edges
                    corr = corr_matrix.loc[gene1, gene2]
                    
                    # Add edge if correlation is above threshold
                    if corr >= correlation_threshold:
                        G.add_edge(gene1, gene2, weight=corr)
        
        print(f"Created gene network with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        
        return G
    
    def analyze_spatial_network(self, G):
        """
        Analyze the spatial proximity network.
        
        Parameters:
        -----------
        G : nx.Graph
            NetworkX graph of spatial points
        """
        if G is None or G.number_of_nodes() == 0:
            print("No network to analyze")
            return
        
        print("\nAnalyzing spatial proximity network...")
        
        # Calculate network metrics
        degree_centrality = nx.degree_centrality(G)
        betweenness_centrality = nx.betweenness_centrality(G)
        closeness_centrality = nx.closeness_centrality(G)
        
        # Identify communities
        communities = community_louvain.best_partition(G)
        
        # Add metrics and community assignments to nodes
        for node in G.nodes():
            G.nodes[node]['degree_centrality'] = degree_centrality[node]
            G.nodes[node]['betweenness_centrality'] = betweenness_centrality[node]
            G.nodes[node]['closeness_centrality'] = closeness_centrality[node]
            G.nodes[node]['community'] = communities[node]
        
        # Count nodes per community
        community_counts = defaultdict(int)
        for node, comm in communities.items():
            community_counts[comm] += 1
        
        print(f"Identified {len(community_counts)} communities")
        
        # Visualize the network
        self.visualize_spatial_network(G)
        
        # Visualize community distribution
        plt.figure(figsize=(10, 6))
        comm_df = pd.DataFrame({
            'community': list(community_counts.keys()),
            'count': list(community_counts.values())
        })
        comm_df = comm_df.sort_values('count', ascending=False)
        
        sns.barplot(x='community', y='count', data=comm_df)
        plt.title('Community Size Distribution')
        plt.xlabel('Community ID')
        plt.ylabel('Number of Points')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, "community_distribution.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Community distribution visualization saved to {output_file}")
        plt.close()
        
        # Analyze centrality metrics
        centrality_df = pd.DataFrame({
            'node': list(G.nodes()),
            'degree_centrality': [G.nodes[n]['degree_centrality'] for n in G.nodes()],
            'betweenness_centrality': [G.nodes[n]['betweenness_centrality'] for n in G.nodes()],
            'closeness_centrality': [G.nodes[n]['closeness_centrality'] for n in G.nodes()],
            'community': [G.nodes[n]['community'] for n in G.nodes()]
        })
        
        # Save centrality metrics
        output_file = os.path.join(self.output_dir, "centrality_metrics.csv")
        centrality_df.to_csv(output_file, index=False)
        print(f"Centrality metrics saved to {output_file}")
        
        # Visualize centrality distributions
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        sns.histplot(centrality_df['degree_centrality'], kde=True, ax=axes[0])
        axes[0].set_title('Degree Centrality Distribution')
        axes[0].set_xlabel('Degree Centrality')
        
        sns.histplot(centrality_df['betweenness_centrality'], kde=True, ax=axes[1])
        axes[1].set_title('Betweenness Centrality Distribution')
        axes[1].set_xlabel('Betweenness Centrality')
        
        sns.histplot(centrality_df['closeness_centrality'], kde=True, ax=axes[2])
        axes[2].set_title('Closeness Centrality Distribution')
        axes[2].set_xlabel('Closeness Centrality')
        
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, "centrality_distributions.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Centrality distributions visualization saved to {output_file}")
        plt.close()
    
    def visualize_spatial_network(self, G):
        """
        Visualize the spatial proximity network.
        
        Parameters:
        -----------
        G : nx.Graph
            NetworkX graph of spatial points
        """
        if G is None or G.number_of_nodes() == 0:
            return
        
        # Extract node positions
        pos = {node: (G.nodes[node]['x'], G.nodes[node]['y']) for node in G.nodes()}
        
        # Extract community assignments
        communities = [G.nodes[node]['community'] for node in G.nodes()]
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Draw edges
        nx.draw_networkx_edges(G, pos, alpha=0.2, width=0.5)
        
        # Draw nodes colored by community
        nx.draw_networkx_nodes(G, pos, node_color=communities, cmap='tab20', 
                              node_size=50, alpha=0.8)
        
        plt.title('Spatial Proximity Network')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.axis('off')
        
        # Add colorbar for communities
        sm = plt.cm.ScalarMappable(cmap='tab20', norm=plt.Normalize(vmin=min(communities), vmax=max(communities)))
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.set_label('Community')
        
        output_file = os.path.join(self.output_dir, "spatial_network.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Spatial network visualization saved to {output_file}")
        plt.close()
    
    def analyze_gene_network(self, G):
        """
        Analyze the gene co-expression network.
        
        Parameters:
        -----------
        G : nx.Graph
            NetworkX graph of genes
        """
        if G is None or G.number_of_nodes() == 0:
            print("No gene network to analyze")
            return
        
        print("\nAnalyzing gene co-expression network...")
        
        # Calculate network metrics
        degree_centrality = nx.degree_centrality(G)
        betweenness_centrality = nx.betweenness_centrality(G)
        eigenvector_centrality = nx.eigenvector_centrality_numpy(G)
        
        # Identify communities
        communities = community_louvain.best_partition(G)
        
        # Add metrics and community assignments to nodes
        for node in G.nodes():
            G.nodes[node]['degree_centrality'] = degree_centrality[node]
            G.nodes[node]['betweenness_centrality'] = betweenness_centrality[node]
            G.nodes[node]['eigenvector_centrality'] = eigenvector_centrality[node]
            G.nodes[node]['community'] = communities[node]
        
        # Count nodes per community
        community_counts = defaultdict(int)
        for node, comm in communities.items():
            community_counts[comm] += 1
        
        print(f"Identified {len(community_counts)} gene communities")
        
        # Visualize the network
        self.visualize_gene_network(G)
        
        # Identify hub genes (high degree centrality)
        hub_threshold = np.percentile(list(degree_centrality.values()), 90)
        hub_genes = [node for node, dc in degree_centrality.items() if dc >= hub_threshold]
        
        print(f"Identified {len(hub_genes)} hub genes")
        print("Top hub genes:")
        for gene in sorted(hub_genes, key=lambda g: degree_centrality[g], reverse=True)[:10]:
            print(f"  {gene}: {degree_centrality[gene]:.4f}")
        
        # Identify bridge genes (high betweenness centrality)
        bridge_threshold = np.percentile(list(betweenness_centrality.values()), 90)
        bridge_genes = [node for node, bc in betweenness_centrality.items() if bc >= bridge_threshold]
        
        print(f"Identified {len(bridge_genes)} bridge genes")
        print("Top bridge genes:")
        for gene in sorted(bridge_genes, key=lambda g: betweenness_centrality[g], reverse=True)[:10]:
            print(f"  {gene}: {betweenness_centrality[g]:.4f}")
        
        # Save gene network metrics
        gene_metrics = pd.DataFrame({
            'gene': list(G.nodes()),
            'degree_centrality': [G.nodes[n]['degree_centrality'] for n in G.nodes()],
            'betweenness_centrality': [G.nodes[n]['betweenness_centrality'] for n in G.nodes()],
            'eigenvector_centrality': [G.nodes[n]['eigenvector_centrality'] for n in G.nodes()],
            'community': [G.nodes[n]['community'] for n in G.nodes()],
            'is_hub': [1 if n in hub_genes else 0 for n in G.nodes()],
            'is_bridge': [1 if n in bridge_genes else 0 for n in G.nodes()]
        })
        
        output_file = os.path.join(self.output_dir, "gene_network_metrics.csv")
        gene_metrics.to_csv(output_file, index=False)
        print(f"Gene network metrics saved to {output_file}")
        
        # Analyze gene communities
        self.analyze_gene_communities(G, communities)
    
    def visualize_gene_network(self, G):
        """
        Visualize the gene co-expression network.
        
        Parameters:
        -----------
        G : nx.Graph
            NetworkX graph of genes
        """
        if G is None or G.number_of_nodes() == 0:
            return
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Extract community assignments
        communities = [G.nodes[node]['community'] for node in G.nodes()]
        
        # Calculate node sizes based on degree centrality
        node_sizes = [300 * G.nodes[node]['degree_centrality'] + 20 for node in G.nodes()]
        
        # Use force-directed layout
        pos = nx.spring_layout(G, seed=42)
        
        # Draw edges with weights as width
        edge_weights = [G[u][v]['weight'] * 2 for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, alpha=0.3, width=edge_weights)
        
        # Draw nodes colored by community
        nx.draw_networkx_nodes(G, pos, node_color=communities, cmap='tab20', 
                              node_size=node_sizes, alpha=0.8)
        
        # Draw labels for hub genes
        hub_threshold = np.percentile([G.nodes[n]['degree_centrality'] for n in G.nodes()], 90)
        hub_genes = {n: n for n in G.nodes() if G.nodes[n]['degree_centrality'] >= hub_threshold}
        nx.draw_networkx_labels(G, pos, labels=hub_genes, font_size=8, font_weight='bold')
        
        plt.title('Gene Co-expression Network')
        plt.axis('off')
        
        # Add colorbar for communities
        sm = plt.cm.ScalarMappable(cmap='tab20', norm=plt.Normalize(vmin=min(communities), vmax=max(communities)))
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.set_label('Community')
        
        output_file = os.path.join(self.output_dir, "gene_network.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Gene network visualization saved to {output_file}")
        plt.close()
    
    def analyze_gene_communities(self, G, communities):
        """
        Analyze gene communities in the co-expression network.
        
        Parameters:
        -----------
        G : nx.Graph
            NetworkX graph of genes
        communities : dict
            Dictionary mapping nodes to community IDs
        """
        # Group genes by community
        community_genes = defaultdict(list)
        for gene, comm in communities.items():
            community_genes[comm].append(gene)
        
        # Analyze each community
        community_data = []
        for comm, genes in community_genes.items():
            # Calculate average centrality metrics
            avg_degree = np.mean([G.nodes[g]['degree_centrality'] for g in genes])
            avg_betweenness = np.mean([G.nodes[g]['betweenness_centrality'] for g in genes])
            avg_eigenvector = np.mean([G.nodes[g]['eigenvector_centrality'] for g in genes])
            
            # Identify top genes in this community
            top_genes = sorted(genes, key=lambda g: G.nodes[g]['degree_centrality'], reverse=True)[:5]
            
            community_data.append({
                'community': comm,
                'size': len(genes),
                'avg_degree_centrality': avg_degree,
                'avg_betweenness_centrality': avg_betweenness,
                'avg_eigenvector_centrality': avg_eigenvector,
                'top_genes': ', '.join(top_genes)
            })
        
        # Convert to DataFrame
        comm_df = pd.DataFrame(community_data)
        
        # Save community data
        output_file = os.path.join(self.output_dir, "gene_communities.csv")
        comm_df.to_csv(output_file, index=False)
        print(f"Gene community data saved to {output_file}")
        
        # Visualize community sizes
        plt.figure(figsize=(12, 6))
        comm_df = comm_df.sort_values('size', ascending=False)
        
        sns.barplot(x='community', y='size', data=comm_df)
        plt.title('Gene Community Size Distribution')
        plt.xlabel('Community ID')
        plt.ylabel('Number of Genes')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, "gene_community_sizes.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Gene community size visualization saved to {output_file}")
        plt.close()
    
    def analyze_spatial_gene_correlations(self):
        """
        Analyze correlations between spatial proximity and gene co-expression.
        """
        print("\nAnalyzing spatial-gene correlations...")
        
        # Get genes and their expression points
        genes_df = self.get_genes()
        
        if genes_df.empty:
            print("No gene expression data found")
            return
        
        # Create a pivot table: genes x points
        pivot_df = genes_df.pivot_table(
            index='gene', 
            columns='point', 
            values='expression_level',
            fill_value=0
        )
        
        # Calculate gene-gene correlation matrix
        gene_corr = pivot_df.T.corr(method='spearman')
        
        # Get spatial points
        points_df = self.get_spatial_points()
        
        if points_df.empty:
            print("No spatial points found")
            return
        
        # Calculate spatial distance matrix
        coords = points_df[['x', 'y']].values
        spatial_dist = euclidean_distances(coords, coords)
        
        # Create a mapping from point URIs to indices
        point_to_idx = {uri: i for i, uri in enumerate(points_df['point'])}
        
        # Calculate correlation between gene expression similarity and spatial proximity
        gene_pairs = []
        for i, gene1 in enumerate(gene_corr.index):
            for j, gene2 in enumerate(gene_corr.columns):
                if i < j:  # Avoid duplicate pairs
                    # Get expression correlation
                    expr_corr = gene_corr.loc[gene1, gene2]
                    
                    # Get points where these genes are expressed
                    gene1_points = genes_df[genes_df['gene'] == gene1]['point'].tolist()
                    gene2_points = genes_df[genes_df['gene'] == gene2]['point'].tolist()
                    
                    # Calculate average spatial distance between expression points
                    if gene1_points and gene2_points:
                        distances = []
                        for p1 in gene1_points:
                            for p2 in gene2_points:
                                if p1 in point_to_idx and p2 in point_to_idx:
                                    idx1 = point_to_idx[p1]
                                    idx2 = point_to_idx[p2]
                                    distances.append(spatial_dist[idx1, idx2])
                        
                        if distances:
                            avg_dist = np.mean(distances)
                            
                            gene_pairs.append({
                                'gene1': gene1,
                                'gene2': gene2,
                                'expression_correlation': expr_corr,
                                'spatial_distance': avg_dist
                            })
        
        # Convert to DataFrame
        pairs_df = pd.DataFrame(gene_pairs)
        
        if pairs_df.empty:
            print("No gene pairs to analyze")
            return
        
        # Calculate correlation between expression similarity and spatial proximity
        corr, p_value = pearsonr(pairs_df['expression_correlation'], -pairs_df['spatial_distance'])
        
        print(f"Correlation between expression similarity and spatial proximity: {corr:.4f} (p={p_value:.4f})")
        
        # Visualize the relationship
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x='spatial_distance', y='expression_correlation', data=pairs_df, alpha=0.5)
        
        # Add regression line
        sns.regplot(x='spatial_distance', y='expression_correlation', data=pairs_df, 
                   scatter=False, line_kws={'color': 'red'})
        
        plt.title(f'Expression Correlation vs. Spatial Distance (r={corr:.4f}, p={p_value:.4f})')
        plt.xlabel('Spatial Distance')
        plt.ylabel('Expression Correlation')
        plt.grid(True, alpha=0.3)
        
        output_file = os.path.join(self.output_dir, "expression_vs_distance.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Expression-distance correlation visualization saved to {output_file}")
        plt.close()
        
        # Save correlation data
        output_file = os.path.join(self.output_dir, "gene_pair_correlations.csv")
        pairs_df.to_csv(output_file, index=False)
        print(f"Gene pair correlation data saved to {output_file}")
    
    def run_analysis(self):
        """
        Run all network analyses.
        """
        # Build and analyze spatial proximity network
        spatial_network = self.build_spatial_proximity_network()
        self.analyze_spatial_network(spatial_network)
        
        # Build and analyze gene co-expression network
        gene_network = self.build_gene_co_expression_network()
        self.analyze_gene_network(gene_network)
        
        # Analyze correlations between spatial proximity and gene co-expression
        self.analyze_spatial_gene_correlations()
        
        print("\nNetwork analysis completed!")

def main():
    parser = argparse.ArgumentParser(description='Spatial Network Analyzer')
    parser.add_argument('--ttl', type=str, default='../data/enhanced_spatial_ontology.ttl',
                        help='Path to the TURTLE file containing the enhanced spatial ontology')
    parser.add_argument('--output', type=str, default='../output/network_analysis',
                        help='Path to the output directory')
    
    args = parser.parse_args()
    
    analyzer = SpatialNetworkAnalyzer(args.ttl, args.output)
    analyzer.run_analysis()
    
    print("\nAnalysis completed successfully!")

if __name__ == "__main__":
    main() 