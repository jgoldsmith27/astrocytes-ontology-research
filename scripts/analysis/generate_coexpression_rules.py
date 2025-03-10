#!/usr/bin/env python3
"""
Script to analyze single-cell RNA-seq data for co-expression patterns
and generate SPARQL inference rules for cell type identification.

This script reads single-cell data, identifies co-expression patterns for
astrocyte subtypes, and generates SPARQL CONSTRUCT queries that can be used
for inference in spatial data.
"""

import os
import argparse
import numpy as np
import pandas as pd
import anndata as ad
from tqdm import tqdm
import scipy.sparse as sp
from itertools import combinations
import networkx as nx

def analyze_gene_coexpression(adata, cell_types, min_expression=0.2, 
                              min_coexpression_prob=0.6, top_genes_per_type=20):
    """
    Analyze gene co-expression patterns in single-cell data.
    
    Parameters:
    -----------
    adata : AnnData
        Single-cell data object
    cell_types : list
        List of cell types to analyze
    min_expression : float
        Minimum expression level to consider a gene expressed
    min_coexpression_prob : float
        Minimum conditional probability for co-expression rules
    top_genes_per_type : int
        Number of top genes to include for each cell type
    
    Returns:
    --------
    rules_dict : dict
        Dictionary of co-expression rules by cell type
    """
    print("Analyzing gene co-expression patterns...")
    rules_dict = {}
    
    # For each cell type
    for cell_type in tqdm(cell_types):
        # Get cells of this type
        cell_mask = adata.obs['celltype'] == cell_type
        cells = adata[cell_mask]
        
        # Skip if no cells of this type
        if cells.shape[0] == 0:
            print(f"No cells found of type {cell_type}, skipping")
            continue
        
        print(f"\nAnalyzing {cell_type} ({cells.shape[0]} cells)")
        
        # Convert to dense if sparse matrix
        if sp.issparse(cells.X):
            expr_matrix = cells.X.toarray()
        else:
            expr_matrix = cells.X
        
        # Get average expression for each gene
        gene_means = np.mean(expr_matrix, axis=0)
        
        # Get top expressed genes
        top_gene_indices = np.argsort(gene_means)[-top_genes_per_type:]
        top_gene_names = adata.var_names[top_gene_indices].tolist()
        
        print(f"Top {len(top_gene_names)} genes for {cell_type}: {', '.join(top_gene_names[:5])}...")
        
        # Create binary expression matrix (1 if expressed, 0 if not)
        binary_expr = (expr_matrix > min_expression).astype(int)
        
        # Calculate co-expression probabilities
        coexpr_rules = []
        gene_pairs_analyzed = 0
        significant_pairs = 0
        
        # Build co-expression graph
        G = nx.Graph()
        
        # Build co-expression rules for all gene pairs
        for i, gene1 in enumerate(top_gene_indices):
            gene1_name = adata.var_names[gene1]
            G.add_node(gene1_name)
            
            gene1_expressed = binary_expr[:, gene1] == 1
            prob_gene1 = np.mean(gene1_expressed)
            
            if prob_gene1 < 0.1:  # Skip rarely expressed genes
                continue
                
            for j, gene2 in enumerate(top_gene_indices):
                if i >= j:  # Avoid duplicate pairs and self-pairs
                    continue
                    
                gene2_name = adata.var_names[gene2]
                gene_pairs_analyzed += 1
                
                gene2_expressed = binary_expr[:, gene2] == 1
                prob_gene2 = np.mean(gene2_expressed)
                
                if prob_gene2 < 0.1:  # Skip rarely expressed genes
                    continue
                
                # Calculate conditional probabilities
                prob_gene2_given_gene1 = np.mean(gene2_expressed[gene1_expressed]) if np.sum(gene1_expressed) > 0 else 0
                prob_gene1_given_gene2 = np.mean(gene1_expressed[gene2_expressed]) if np.sum(gene2_expressed) > 0 else 0
                
                # Calculate joint probability
                prob_both = np.mean(gene1_expressed & gene2_expressed)
                
                # Add edge with weight based on conditional probabilities
                edge_weight = max(prob_gene2_given_gene1, prob_gene1_given_gene2)
                
                # Only keep strong co-expression relationships
                if edge_weight >= min_coexpression_prob:
                    G.add_edge(gene1_name, gene2_name, weight=edge_weight)
                    significant_pairs += 1
                    
                    # Add rule to list
                    coexpr_rules.append({
                        'gene1': gene1_name,
                        'gene2': gene2_name,
                        'prob_gene1': prob_gene1,
                        'prob_gene2': prob_gene2,
                        'prob_gene2_given_gene1': prob_gene2_given_gene1,
                        'prob_gene1_given_gene2': prob_gene1_given_gene2,
                        'joint_prob': prob_both,
                        'weight': edge_weight
                    })
        
        print(f"Analyzed {gene_pairs_analyzed} gene pairs, found {significant_pairs} significant co-expressions")
        
        # Find cliques in the graph (groups of genes that are all co-expressed)
        cliques = list(nx.find_cliques(G))
        
        # Sort cliques by size (largest first)
        cliques.sort(key=len, reverse=True)
        
        # Only keep cliques of at least 3 genes
        cliques = [c for c in cliques if len(c) >= 3]
        
        # Generate rules from cliques and individual gene pairs
        rules = []
        
        # Add rules from cliques
        for i, clique in enumerate(cliques[:5]):  # Keep top 5 cliques
            if len(clique) >= 3:
                rule = {
                    'type': 'clique',
                    'genes': sorted(clique),
                    'size': len(clique),
                    'confidence': 0.9  # High confidence for cliques
                }
                rules.append(rule)
        
        # Add rules from highest weight gene pairs
        sorted_pairs = sorted(coexpr_rules, key=lambda x: x['weight'], reverse=True)
        for pair in sorted_pairs[:10]:  # Keep top 10 pairs
            rule = {
                'type': 'pair',
                'genes': [pair['gene1'], pair['gene2']],
                'weight': pair['weight'],
                'confidence': pair['weight']
            }
            rules.append(rule)
        
        rules_dict[cell_type] = rules
        
        print(f"Generated {len(rules)} rules for {cell_type}")
    
    return rules_dict

def generate_sparql_rules(rules_dict, output_dir):
    """
    Generate SPARQL CONSTRUCT rules from co-expression patterns.
    
    Parameters:
    -----------
    rules_dict : dict
        Dictionary of co-expression rules by cell type
    output_dir : str
        Directory to save the rules
    """
    print("Generating SPARQL rules...")
    
    os.makedirs(output_dir, exist_ok=True)
    
    all_rules = []
    
    for cell_type, rules in rules_dict.items():
        print(f"Generating rules for {cell_type}...")
        
        # Format the cell type for URIs
        cell_type_uri = cell_type.replace("s", "") if cell_type.endswith("s") else cell_type
        
        # Generate a rule for each co-expression pattern
        for i, rule in enumerate(rules):
            rule_id = f"{cell_type}_{rule['type']}_{i+1}"
            genes = rule['genes']
            confidence = rule['confidence']
            
            if rule['type'] == 'clique':
                sparql_rule = generate_clique_rule(genes, cell_type, 50, None)
            else:  # pair
                sparql_rule = generate_pair_rule(genes[0], genes[1], rule['weight'], cell_type, 50, None)
            
            all_rules.append((rule_id, sparql_rule))
    
    # Write all rules to a file
    with open(os.path.join(output_dir, "coexpression_rules.sparql"), "w") as f:
        for rule_id, rule in all_rules:
            f.write(f"# Rule: {rule_id}\n")
            f.write(rule)
            f.write("\n\n")
    
    print(f"Generated {len(all_rules)} SPARQL rules and saved to {output_dir}")
    return all_rules

def generate_clique_rule(clique, cell_type, distance_threshold, negative_markers=None):
    """Generate a SPARQL CONSTRUCT rule for identifying cells based on a gene clique."""
    # Updated SPARQL template to use cell: prefix
    rule_template = """
# Clique rule for {cell_type} with {n_genes} genes: {genes}
PREFIX cell: <http://example.org/ontology/cell/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

CONSTRUCT {{
  ?cell a cell:SpatialCell ;
        cell:hasCellType cell:{cell_type} ;
        cell:hasX ?centerX ;
        cell:hasY ?centerY ;
        cell:hasRadius ?radius ;
        cell:hasConfidence {confidence} ;
        cell:hasCellID ?cellID .
  
  # Link to constituent points
{point_links}
}}
WHERE {{
{where_clause}

  # Calculate cell center (average of all points)
{center_calc}

  # Calculate maximum distance from center to any point as radius
{radius_calc}

  # Generate a unique ID for the cell
  BIND(CONCAT("cell_", STRUUID()) AS ?cellID)
  
{negative_marker_constraints}
}}
"""

    # Rest of the function implementation...

def generate_pair_rule(gene_i, gene_j, weight, cell_type, distance_threshold, negative_markers=None):
    """Generate a SPARQL CONSTRUCT rule for identifying cells based on a pair of co-expressed genes."""
    # Updated SPARQL template to use cell: prefix
    rule_template = """
# Pair rule for {cell_type} with genes: {gene_i} and {gene_j} (coexpression = {weight:.2f})
PREFIX cell: <http://example.org/ontology/cell/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

CONSTRUCT {{
  ?cell a cell:SpatialCell ;
        cell:hasCellType cell:{cell_type} ;
        cell:hasX ?centerX ;
        cell:hasY ?centerY ;
        cell:hasRadius ?radius ;
        cell:hasConfidence "{weight}"^^xsd:decimal ;
        cell:hasCellID ?cellID .
  
  # Link to constituent points
  ?cell cell:includesPoint ?point0 .
  ?cell cell:includesPoint ?point1 .
}}
WHERE {{
  ?point0 a cell:SpatialDataPoint ;
         cell:expressesGene ?gene0 ;
         cell:hasXCoordinate ?x0 ;
         cell:hasYCoordinate ?y0 .
  ?gene0 cell:hasGeneID "{gene_i}" .
         
  ?point1 a cell:SpatialDataPoint ;
         cell:expressesGene ?gene1 ;
         cell:hasXCoordinate ?x1 ;
         cell:hasYCoordinate ?y1 .
  ?gene1 cell:hasGeneID "{gene_j}" .
  
  # Calculate distance between points
  BIND(SQRT(POW(?x0 - ?x1, 2) + POW(?y0 - ?y1, 2)) AS ?dist)
  FILTER(?dist < {distance_threshold})
  
  # Calculate cell center (midpoint)
  BIND((?x0 + ?x1) / 2 AS ?centerX)
  BIND((?y0 + ?y1) / 2 AS ?centerY)
  
  # Calculate radius (half of distance)
  BIND(?dist / 2 AS ?radius)
  
  # Generate a unique ID for the cell
  BIND(CONCAT("cell_", STRUUID()) AS ?cellID)
  
{negative_marker_constraints}
}}
"""

    # Rest of the function implementation...

def main():
    parser = argparse.ArgumentParser(description="Generate co-expression rules from single-cell data")
    parser.add_argument("--input", "-i", required=True, help="Input h5ad file")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory for SPARQL rules")
    parser.add_argument("--cell-types", "-c", nargs="+", default=["Astrocytes", "Astrocytes1"], 
                        help="Cell types to analyze")
    parser.add_argument("--min-expression", "-e", type=float, default=0.2, 
                        help="Minimum expression level to consider a gene expressed")
    parser.add_argument("--min-coexpression", "-p", type=float, default=0.6, 
                        help="Minimum co-expression probability for rules")
    parser.add_argument("--top-genes", "-g", type=int, default=20, 
                        help="Number of top genes to analyze per cell type")
    args = parser.parse_args()
    
    # Read the single-cell data
    print(f"Reading single-cell data from {args.input}...")
    adata = ad.read_h5ad(args.input)
    
    # Analyze co-expression patterns
    rules_dict = analyze_gene_coexpression(
        adata, 
        args.cell_types, 
        args.min_expression,
        args.min_coexpression,
        args.top_genes
    )
    
    # Generate SPARQL rules
    generate_sparql_rules(rules_dict, args.output_dir)
    
    print("Done!")

if __name__ == "__main__":
    main() 