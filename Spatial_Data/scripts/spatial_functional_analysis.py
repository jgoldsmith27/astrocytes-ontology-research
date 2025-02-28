#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF, RDFS
import requests
import json
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

class SpatialFunctionalAnalyzer:
    """
    A class to perform functional enrichment analysis on genes with specific spatial relationships.
    """
    
    def __init__(self, ttl_file, output_dir="../output/functional_analysis"):
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
        
        # Cache for gene annotations
        self.gene_annotations = {}
    
    def get_all_genes(self):
        """
        Get all genes in the ontology.
        
        Returns:
        --------
        list
            List of gene symbols
        """
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT DISTINCT ?gene
        WHERE {
            ?gene rdf:type :Gene .
        }
        """
        
        results = self.graph.query(query)
        
        genes = []
        for row in results:
            gene_uri = str(row.gene)
            gene_symbol = gene_uri.split('_')[1] if '_' in gene_uri else gene_uri.split('#')[-1]
            genes.append(gene_symbol)
        
        return genes
    
    def get_genes_by_region(self):
        """
        Get genes grouped by spatial region.
        
        Returns:
        --------
        dict
            Dictionary mapping region names to lists of gene symbols
        """
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?regionName
        WHERE {
            ?gene rdf:type :Gene ;
                  :expressedAt ?point .
                  
            ?point :locatedInRegion ?region .
            
            ?region :hasRegionName ?regionName .
        }
        """
        
        results = self.graph.query(query)
        
        region_genes = defaultdict(set)
        for row in results:
            gene_uri = str(row.gene)
            gene_symbol = gene_uri.split('_')[1] if '_' in gene_uri else gene_uri.split('#')[-1]
            region_name = str(row.regionName)
            region_genes[region_name].add(gene_symbol)
        
        # Convert sets to lists
        return {region: list(genes) for region, genes in region_genes.items()}
    
    def get_genes_by_directional_relationship(self):
        """
        Get genes grouped by directional relationship.
        
        Returns:
        --------
        dict
            Dictionary mapping directional relationships to lists of gene pairs
        """
        directional_relationships = [
            'directlyAbove', 'directlyBelow', 'directlyLeftOf', 'directlyRightOf', 'diagonallyRelatedTo'
        ]
        
        relationship_genes = {}
        
        for relationship in directional_relationships:
            query = f"""
            PREFIX : <http://example.org/spatial-transcriptomics#>
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            
            SELECT DISTINCT ?gene1 ?gene2
            WHERE {{
                ?gene1 rdf:type :Gene ;
                       :expressedAt ?point1 .
                       
                ?gene2 rdf:type :Gene ;
                       :expressedAt ?point2 .
                       
                ?point1 :{relationship} ?point2 .
                
                # Ensure we get distinct gene pairs
                FILTER(?gene1 != ?gene2)
            }}
            """
            
            results = self.graph.query(query)
            
            gene_pairs = []
            for row in results:
                gene1_uri = str(row.gene1)
                gene2_uri = str(row.gene2)
                
                gene1_symbol = gene1_uri.split('_')[1] if '_' in gene1_uri else gene1_uri.split('#')[-1]
                gene2_symbol = gene2_uri.split('_')[1] if '_' in gene2_uri else gene2_uri.split('#')[-1]
                
                gene_pairs.append((gene1_symbol, gene2_symbol))
            
            relationship_genes[relationship] = gene_pairs
        
        return relationship_genes
    
    def get_gene_annotations(self, gene_symbols, annotation_type='go'):
        """
        Get functional annotations for a list of genes.
        
        Parameters:
        -----------
        gene_symbols : list
            List of gene symbols to get annotations for
        annotation_type : str
            Type of annotation to retrieve ('go', 'pathway', 'kegg')
            
        Returns:
        --------
        dict
            Dictionary mapping gene symbols to lists of annotations
        """
        # Check cache first
        uncached_genes = [gene for gene in gene_symbols if gene not in self.gene_annotations]
        
        if uncached_genes:
            print(f"Fetching annotations for {len(uncached_genes)} genes...")
            
            # In a real implementation, you would use a proper API to fetch annotations
            # For demonstration purposes, we'll simulate some annotations
            
            # Simulated GO terms
            go_terms = [
                "GO:0006412", "GO:0006457", "GO:0006468", "GO:0006508", "GO:0006511",
                "GO:0006629", "GO:0006810", "GO:0006811", "GO:0006814", "GO:0006886",
                "GO:0007049", "GO:0007155", "GO:0007165", "GO:0007186", "GO:0007268",
                "GO:0008152", "GO:0008380", "GO:0009056", "GO:0009058", "GO:0009987"
            ]
            
            # Simulated pathways
            pathways = [
                "KEGG_CELL_CYCLE", "KEGG_APOPTOSIS", "KEGG_WNT_SIGNALING_PATHWAY",
                "KEGG_NOTCH_SIGNALING_PATHWAY", "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
                "KEGG_TGF_BETA_SIGNALING_PATHWAY", "KEGG_VEGF_SIGNALING_PATHWAY",
                "KEGG_FOCAL_ADHESION", "KEGG_GAP_JUNCTION", "KEGG_TIGHT_JUNCTION",
                "KEGG_ADHERENS_JUNCTION", "KEGG_MAPK_SIGNALING_PATHWAY",
                "KEGG_ERBB_SIGNALING_PATHWAY", "KEGG_CALCIUM_SIGNALING_PATHWAY",
                "KEGG_JAK_STAT_SIGNALING_PATHWAY", "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
            ]
            
            # Assign random annotations to genes
            np.random.seed(42)  # For reproducibility
            for gene in uncached_genes:
                if annotation_type == 'go':
                    # Assign 2-5 random GO terms to each gene
                    num_terms = np.random.randint(2, 6)
                    self.gene_annotations[gene] = np.random.choice(go_terms, size=num_terms, replace=False).tolist()
                elif annotation_type in ['pathway', 'kegg']:
                    # Assign 1-3 random pathways to each gene
                    num_pathways = np.random.randint(1, 4)
                    self.gene_annotations[gene] = np.random.choice(pathways, size=num_pathways, replace=False).tolist()
        
        # Return annotations for requested genes
        return {gene: self.gene_annotations.get(gene, []) for gene in gene_symbols}
    
    def perform_enrichment_analysis(self, gene_set, background_genes=None, annotation_type='go'):
        """
        Perform enrichment analysis on a set of genes.
        
        Parameters:
        -----------
        gene_set : list
            List of gene symbols to analyze
        background_genes : list
            List of background genes to use for comparison (if None, use all genes)
        annotation_type : str
            Type of annotation to use for enrichment ('go', 'pathway', 'kegg')
            
        Returns:
        --------
        pd.DataFrame
            DataFrame containing enrichment results
        """
        if background_genes is None:
            background_genes = self.get_all_genes()
        
        # Get annotations for gene set and background
        gene_set_annotations = self.get_gene_annotations(gene_set, annotation_type)
        background_annotations = self.get_gene_annotations(background_genes, annotation_type)
        
        # Count annotations in gene set
        annotation_counts = defaultdict(int)
        for gene, annotations in gene_set_annotations.items():
            for annotation in annotations:
                annotation_counts[annotation] += 1
        
        # Count annotations in background
        background_counts = defaultdict(int)
        for gene, annotations in background_annotations.items():
            for annotation in annotations:
                background_counts[annotation] += 1
        
        # Perform Fisher's exact test for each annotation
        results = []
        for annotation, count in annotation_counts.items():
            # Contingency table
            # a: genes in set with annotation
            # b: genes in background (not in set) with annotation
            # c: genes in set without annotation
            # d: genes in background (not in set) without annotation
            a = count
            b = background_counts[annotation] - count
            c = len(gene_set) - count
            d = len(background_genes) - len(gene_set) - b
            
            # Fisher's exact test
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
            
            results.append({
                'annotation': annotation,
                'count': count,
                'gene_set_size': len(gene_set),
                'background_count': background_counts[annotation],
                'background_size': len(background_genes),
                'odds_ratio': odds_ratio,
                'p_value': p_value
            })
        
        # Convert to DataFrame
        df = pd.DataFrame(results)
        
        if not df.empty:
            # Multiple testing correction
            df['adjusted_p_value'] = multipletests(df['p_value'], method='fdr_bh')[1]
            
            # Sort by adjusted p-value
            df = df.sort_values('adjusted_p_value')
        
        return df
    
    def analyze_region_enrichment(self):
        """
        Analyze functional enrichment for genes in each spatial region.
        """
        print("\n=== Analyzing Functional Enrichment by Region ===")
        
        # Get genes by region
        region_genes = self.get_genes_by_region()
        
        if not region_genes:
            print("No regions with genes found")
            return
        
        # Get all genes as background
        all_genes = self.get_all_genes()
        
        # Analyze each region
        for region, genes in region_genes.items():
            if len(genes) < 5:
                print(f"Skipping region {region} with only {len(genes)} genes")
                continue
            
            print(f"\nAnalyzing functional enrichment for {len(genes)} genes in {region}")
            
            # Perform GO enrichment
            go_results = self.perform_enrichment_analysis(genes, all_genes, 'go')
            
            if not go_results.empty:
                print(f"Found {len(go_results)} enriched GO terms")
                print(go_results.head(5)[['annotation', 'count', 'odds_ratio', 'adjusted_p_value']])
                
                # Save results
                output_file = os.path.join(self.output_dir, f"go_enrichment_{region}.csv")
                go_results.to_csv(output_file, index=False)
                print(f"GO enrichment results saved to {output_file}")
                
                # Visualize top results
                self.visualize_enrichment_results(go_results, f"GO Enrichment for {region}", f"go_enrichment_{region}")
            else:
                print("No enriched GO terms found")
            
            # Perform pathway enrichment
            pathway_results = self.perform_enrichment_analysis(genes, all_genes, 'pathway')
            
            if not pathway_results.empty:
                print(f"Found {len(pathway_results)} enriched pathways")
                print(pathway_results.head(5)[['annotation', 'count', 'odds_ratio', 'adjusted_p_value']])
                
                # Save results
                output_file = os.path.join(self.output_dir, f"pathway_enrichment_{region}.csv")
                pathway_results.to_csv(output_file, index=False)
                print(f"Pathway enrichment results saved to {output_file}")
                
                # Visualize top results
                self.visualize_enrichment_results(pathway_results, f"Pathway Enrichment for {region}", f"pathway_enrichment_{region}")
            else:
                print("No enriched pathways found")
    
    def analyze_directional_enrichment(self):
        """
        Analyze functional enrichment for genes with directional relationships.
        """
        print("\n=== Analyzing Functional Enrichment by Directional Relationship ===")
        
        # Get genes by directional relationship
        relationship_genes = self.get_genes_by_directional_relationship()
        
        if not relationship_genes:
            print("No directional relationships with genes found")
            return
        
        # Get all genes as background
        all_genes = self.get_all_genes()
        
        # Analyze each relationship
        for relationship, gene_pairs in relationship_genes.items():
            if not gene_pairs:
                print(f"No gene pairs found for {relationship}")
                continue
            
            # Extract unique genes from pairs
            genes = set()
            for gene1, gene2 in gene_pairs:
                genes.add(gene1)
                genes.add(gene2)
            
            genes = list(genes)
            
            if len(genes) < 5:
                print(f"Skipping {relationship} with only {len(genes)} genes")
                continue
            
            print(f"\nAnalyzing functional enrichment for {len(genes)} genes with {relationship} relationship")
            
            # Perform GO enrichment
            go_results = self.perform_enrichment_analysis(genes, all_genes, 'go')
            
            if not go_results.empty:
                print(f"Found {len(go_results)} enriched GO terms")
                print(go_results.head(5)[['annotation', 'count', 'odds_ratio', 'adjusted_p_value']])
                
                # Save results
                output_file = os.path.join(self.output_dir, f"go_enrichment_{relationship}.csv")
                go_results.to_csv(output_file, index=False)
                print(f"GO enrichment results saved to {output_file}")
                
                # Visualize top results
                self.visualize_enrichment_results(go_results, f"GO Enrichment for {relationship}", f"go_enrichment_{relationship}")
            else:
                print("No enriched GO terms found")
            
            # Perform pathway enrichment
            pathway_results = self.perform_enrichment_analysis(genes, all_genes, 'pathway')
            
            if not pathway_results.empty:
                print(f"Found {len(pathway_results)} enriched pathways")
                print(pathway_results.head(5)[['annotation', 'count', 'odds_ratio', 'adjusted_p_value']])
                
                # Save results
                output_file = os.path.join(self.output_dir, f"pathway_enrichment_{relationship}.csv")
                pathway_results.to_csv(output_file, index=False)
                print(f"Pathway enrichment results saved to {output_file}")
                
                # Visualize top results
                self.visualize_enrichment_results(pathway_results, f"Pathway Enrichment for {relationship}", f"pathway_enrichment_{relationship}")
            else:
                print("No enriched pathways found")
    
    def visualize_enrichment_results(self, df, title, filename_prefix, max_terms=10):
        """
        Visualize enrichment results.
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame containing enrichment results
        title : str
            Title for the plot
        filename_prefix : str
            Prefix for the output filename
        max_terms : int
            Maximum number of terms to display
        """
        if df.empty:
            return
        
        # Take top terms
        df_top = df.head(max_terms)
        
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Create bar plot of -log10(adjusted p-value)
        df_top['neg_log_p'] = -np.log10(df_top['adjusted_p_value'])
        
        # Sort by neg_log_p for better visualization
        df_top = df_top.sort_values('neg_log_p')
        
        # Create horizontal bar plot
        sns.barplot(x='neg_log_p', y='annotation', data=df_top)
        
        # Add count labels to bars
        for i, row in enumerate(df_top.itertuples()):
            plt.text(row.neg_log_p + 0.1, i, f"{row.count}/{row.gene_set_size}", 
                    verticalalignment='center', fontsize=10)
        
        # Set labels and title
        plt.xlabel('-log10(adjusted p-value)')
        plt.ylabel('Annotation')
        plt.title(title)
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, f"{filename_prefix}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
        plt.close()
    
    def run_analysis(self):
        """
        Run all functional enrichment analyses.
        """
        # Analyze functional enrichment by region
        self.analyze_region_enrichment()
        
        # Analyze functional enrichment by directional relationship
        self.analyze_directional_enrichment()
        
        print("\nFunctional enrichment analysis completed!")

def main():
    parser = argparse.ArgumentParser(description='Spatial Functional Analyzer')
    parser.add_argument('--ttl', type=str, default='../data/enhanced_spatial_ontology.ttl',
                        help='Path to the TURTLE file containing the enhanced spatial ontology')
    parser.add_argument('--output', type=str, default='../output/functional_analysis',
                        help='Path to the output directory')
    
    args = parser.parse_args()
    
    analyzer = SpatialFunctionalAnalyzer(args.ttl, args.output)
    analyzer.run_analysis()
    
    print("\nAnalysis completed successfully!")

if __name__ == "__main__":
    main() 