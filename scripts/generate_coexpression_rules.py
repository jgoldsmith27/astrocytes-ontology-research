#!/usr/bin/env python3
"""
Generate SPARQL rules for identifying spatial cells based on co-expression patterns in single-cell data.
"""

import os
import argparse
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import networkx as nx
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Generate SPARQL rules from single-cell co-expression patterns.")
    parser.add_argument("--input", type=str, required=True, help="Path to the h5ad file containing single-cell data")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to save the generated rules")
    parser.add_argument("--cell-types", nargs='+', default=None, help="Cell types to analyze (if not specified, will analyze all)")
    parser.add_argument("--min-expression", type=float, default=0.2, help="Minimum expression value to consider a gene expressed")
    parser.add_argument("--min-coexpression", type=float, default=0.6, help="Minimum conditional probability for co-expression")
    parser.add_argument("--distance-threshold", type=int, default=50, help="Maximum spatial distance between co-expressed genes")
    parser.add_argument("--top-genes", type=int, default=30, help="Number of top expressed genes to consider for each cell type")
    parser.add_argument("--neg-marker-threshold", type=float, default=0.05, 
                        help="Maximum expression threshold in target type for negative markers")
    parser.add_argument("--min-clique-size", type=int, default=3,
                        help="Minimum number of genes in a clique (default: 3)")
    parser.add_argument("--max-clique-size", type=int, default=None, 
                        help="Maximum number of genes in a clique. If not specified, will find all cliques up to largest possible size.")
    args = parser.parse_args()
    return args

def load_data(input_file):
    """Load single-cell data from h5ad file."""
    logger.info(f"Loading single-cell data from {input_file}")
    adata = sc.read_h5ad(input_file)
    return adata

def get_cell_type_indices(adata, cell_type):
    """Get indices of cells belonging to a specific cell type."""
    return np.where(adata.obs['celltype'] == cell_type)[0]

def get_available_cell_types(adata):
    """Get list of all available cell types in the dataset."""
    return adata.obs['celltype'].unique().tolist()

def get_binary_expression(adata, cell_indices, min_expression):
    """Convert expression matrix to binary (expressed/not expressed) based on threshold."""
    expr_matrix = adata.X[cell_indices]
    binary_matrix = (expr_matrix > min_expression).astype(int)
    return binary_matrix

def calculate_coexpression(binary_matrix, gene_indices, min_coexpression):
    """
    Calculate co-expression probabilities between gene pairs and build a network.
    
    Returns:
        nx.Graph: Network where nodes are genes and edges represent strong co-expression
    """
    n_cells = binary_matrix.shape[0]
    n_genes = len(gene_indices)
    gene_names = [gene_indices[i] for i in range(n_genes)]
    
    # Calculate individual gene expression probabilities
    gene_probs = binary_matrix.mean(axis=0)
    
    # Create co-expression network
    G = nx.Graph()
    
    # Add nodes (genes)
    for i in range(n_genes):
        gene_id = gene_names[i]
        G.add_node(gene_id, probability=gene_probs[i])
    
    # Calculate co-expression and add edges
    for i in range(n_genes):
        gene_i = gene_names[i]
        for j in range(i+1, n_genes):
            gene_j = gene_names[j]
            
            # Calculate joint probability P(A and B)
            joint = np.sum((binary_matrix[:, i] == 1) & (binary_matrix[:, j] == 1)) / n_cells
            
            # Calculate conditional probabilities P(A|B) and P(B|A)
            p_i_given_j = joint / gene_probs[j] if gene_probs[j] > 0 else 0
            p_j_given_i = joint / gene_probs[i] if gene_probs[i] > 0 else 0
            
            # Use maximum of conditional probabilities as the co-expression strength
            coexpr_strength = max(p_i_given_j, p_j_given_i)
            
            # Add edge if co-expression is strong enough
            if coexpr_strength >= min_coexpression:
                G.add_edge(gene_i, gene_j, weight=coexpr_strength, joint_prob=joint)
    
    return G

def identify_negative_markers(adata, target_cell_type, top_genes_per_type=10, 
                             neg_marker_threshold=0.05, min_expr_other=0.5):
    """
    Identify negative markers for a target cell type - genes that are:
    1. Highly expressed in other cell types
    2. Rarely expressed in the target cell type
    
    Args:
        adata: AnnData object with single-cell data
        target_cell_type: The cell type to find negative markers for
        top_genes_per_type: Number of top genes to consider from each other cell type
        neg_marker_threshold: Maximum expression in target type to be a negative marker
        min_expr_other: Minimum expression in other cell types to be a negative marker
    
    Returns:
        list: List of negative marker genes
    """
    target_indices = get_cell_type_indices(adata, target_cell_type)
    target_expr = adata.X[target_indices].mean(axis=0)
    
    negative_markers = []
    
    # For each other cell type, find genes that are highly expressed there but not in target
    for cell_type in adata.obs['celltype'].unique():
        if cell_type == target_cell_type:
            continue
            
        # Get indices for this cell type
        type_indices = get_cell_type_indices(adata, cell_type)
        
        # Calculate mean expression in this cell type
        type_expr = adata.X[type_indices].mean(axis=0)
        
        # Find genes highly expressed in this type but not in target
        for i in range(len(type_expr)):
            if type_expr[i] > min_expr_other and target_expr[i] < neg_marker_threshold:
                gene_name = adata.var_names[i]
                if gene_name not in negative_markers:
                    negative_markers.append(gene_name)
                    
                # Limit to top_genes_per_type per cell type
                if len(negative_markers) >= top_genes_per_type:
                    break
                    
    logger.info(f"Identified {len(negative_markers)} negative markers for {target_cell_type}")
    return negative_markers

def find_cliques(G, min_size=3, max_size=None):
    """
    Find all cliques (fully connected subgraphs) of specified size range.
    
    Parameters:
        G (networkx.Graph): Graph representing gene co-expression relationships
        min_size (int): Minimum clique size to return (default: 3)
        max_size (int, optional): Maximum clique size to return. If None, no upper limit is applied.
    
    Returns:
        list: List of cliques (each clique is a list of node IDs) sorted by size in descending order
    """
    all_cliques = list(nx.find_cliques(G))
    
    if max_size is None:
        # No upper limit, filter only by minimum size
        filtered_cliques = [clique for clique in all_cliques if len(clique) >= min_size]
    else:
        # Filter by both minimum and maximum size
        filtered_cliques = [clique for clique in all_cliques if min_size <= len(clique) <= max_size]
    
    return sorted(filtered_cliques, key=len, reverse=True)

def generate_clique_rule(clique, cell_type, distance_threshold, negative_markers=None):
    """Generate a SPARQL CONSTRUCT rule for identifying cells based on a gene clique."""
    rule_template = """
# Clique rule for {cell_type} with {n_genes} genes: {genes}
PREFIX cell: <http://example.org/ontology/cell/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

CONSTRUCT {{
  ?cell a cell:SpatialCell ;
        cell:hasCellType "{cell_type}" ;
        cell:hasXCoordinate ?centerX ;
        cell:hasYCoordinate ?centerY ;
        cell:hasRadius ?radius ;
        cell:hasConfidence "0.9"^^xsd:decimal ;
        cell:hasCellID ?cellID .
  
  # Link to constituent points
{point_links}
}}
WHERE {{
{points_where}

  # Ensure all points are within distance threshold of each other
{distance_constraints}

  # Calculate cell center (average coordinates)
  BIND(({x_sum}) / {n_genes} AS ?centerX)
  BIND(({y_sum}) / {n_genes} AS ?centerY)
  
  # Calculate radius (half of max distance)
  BIND(MAX({max_distance_options}) / 2 AS ?radius)
  
  # Generate a unique cell ID
  BIND(CONCAT("cell_", STR(?centerX), "_", STR(?centerY), "_", "{cell_type}") AS ?cellID)
  
{negative_marker_constraints}
}}
"""
    
    n_genes = len(clique)
    genes_str = ", ".join(clique)
    
    # Generate WHERE clauses for each point
    points_where = ""
    for i, gene in enumerate(clique):
        points_where += f"""  ?point{i} a cell:SpatialDataPoint ;
        cell:expressesGene ?gene{i} ;
        cell:hasXCoordinate ?x{i} ;
        cell:hasYCoordinate ?y{i} .
  ?gene{i} cell:hasGeneID "{gene}" .
  
"""
    
    # Generate distance constraints between all pairs of points
    distance_constraints = ""
    max_distance_options = []
    for i in range(n_genes):
        for j in range(i+1, n_genes):
            distance_expr = f"SQRT(POW(?x{i} - ?x{j}, 2) + POW(?y{i} - ?y{j}, 2))"
            distance_constraints += f"  BIND({distance_expr} AS ?dist{i}_{j})\n"
            distance_constraints += f"  FILTER(?dist{i}_{j} < {distance_threshold})\n"
            max_distance_options.append(f"?dist{i}_{j}")
    
    # Generate sum expressions for averaging coordinates
    x_sum = " + ".join([f"?x{i}" for i in range(n_genes)])
    y_sum = " + ".join([f"?y{i}" for i in range(n_genes)])
    
    # Generate point links for the CONSTRUCT part
    point_links = ""
    for i in range(n_genes):
        point_links += f"  ?cell cell:includesPoint ?point{i} .\n"
    
    # Generate negative marker constraints if provided
    negative_marker_constraints = ""
    if negative_markers and len(negative_markers) > 0:
        negative_marker_constraints = "  # Ensure no negative markers are present near the cell\n"
        for i, marker in enumerate(negative_markers):
            negative_marker_constraints += f"""  FILTER NOT EXISTS {{
    ?negPoint{i} a cell:SpatialDataPoint ;
               cell:expressesGene ?negGene{i} ;
               cell:hasXCoordinate ?negX{i} ;
               cell:hasYCoordinate ?negY{i} .
    ?negGene{i} cell:hasGeneID "{marker}" .
    BIND(SQRT(POW(?centerX - ?negX{i}, 2) + POW(?centerY - ?negY{i}, 2)) AS ?negDist{i})
    FILTER(?negDist{i} < ?radius * 1.5)
  }}
"""
    
    # Assemble the full rule
    rule = rule_template.format(
        cell_type=cell_type,
        n_genes=n_genes,
        genes=genes_str,
        points_where=points_where,
        distance_constraints=distance_constraints,
        x_sum=x_sum,
        y_sum=y_sum,
        max_distance_options=", ".join(max_distance_options),
        point_links=point_links,
        negative_marker_constraints=negative_marker_constraints
    )
    
    return rule

def generate_pair_rule(gene_i, gene_j, weight, cell_type, distance_threshold, negative_markers=None):
    """Generate a SPARQL CONSTRUCT rule for identifying cells based on a pair of co-expressed genes."""
    rule_template = """
# Pair rule for {cell_type} with genes: {gene1} and {gene2}
# Co-expression weight: {weight:.3f}
PREFIX cell: <http://example.org/ontology/cell/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

CONSTRUCT {{
  ?cell a cell:SpatialCell ;
        cell:hasCellType "{cell_type}" ;
        cell:hasXCoordinate ?centerX ;
        cell:hasYCoordinate ?centerY ;
        cell:hasRadius ?radius ;
        cell:hasConfidence "{confidence}"^^xsd:decimal .
  
  # Link to the two points
  ?cell cell:includesPoint ?point1 .
  ?cell cell:includesPoint ?point2 .
}}
WHERE {{
  # Gene variables
  ?gene1 cell:hasGeneID "{gene1}" .
  ?gene2 cell:hasGeneID "{gene2}" .
  
  # The two spatial data points
  ?point1 a cell:SpatialDataPoint ;
         cell:expressesGene ?gene1 ;
         cell:hasXCoordinate ?x1 ;
         cell:hasYCoordinate ?y1 .
         
  ?point2 a cell:SpatialDataPoint ;
         cell:expressesGene ?gene2 ;
         cell:hasXCoordinate ?x2 ;
         cell:hasYCoordinate ?y2 .
  
  # Calculate distance between points
  BIND(SQRT(POW(?x1 - ?x2, 2) + POW(?y1 - ?y2, 2)) AS ?dist)
  FILTER(?dist < {distance_threshold})
  
  # Calculate cell center (midpoint)
  BIND((?x1 + ?x2) / 2 AS ?centerX)
  BIND((?y1 + ?y2) / 2 AS ?centerY)
  
  # Calculate radius (half of distance)
  BIND(?dist / 2 AS ?radius)
  
  # Generate a unique cell ID
  BIND(CONCAT("cell_", STR(?centerX), "_", STR(?centerY), "_", "{cell_type}") AS ?cellID)
  
{negative_marker_constraints}
}}
"""
    
    # Generate negative marker constraints if provided
    negative_marker_constraints = ""
    if negative_markers and len(negative_markers) > 0:
        negative_marker_constraints = "  # Ensure no negative markers are present near the cell\n"
        for i, marker in enumerate(negative_markers):
            negative_marker_constraints += f"""  FILTER NOT EXISTS {{
    ?negPoint{i} a cell:SpatialDataPoint ;
               cell:expressesGene ?negGene{i} ;
               cell:hasXCoordinate ?negX{i} ;
               cell:hasYCoordinate ?negY{i} .
    ?negGene{i} cell:hasGeneID "{marker}" .
    BIND(SQRT(POW(?centerX - ?negX{i}, 2) + POW(?centerY - ?negY{i}, 2)) AS ?negDist{i})
    FILTER(?negDist{i} < ?radius * 1.5)
  }}
"""
    
    # Assemble the full rule
    rule = rule_template.format(
        cell_type=cell_type,
        gene1=gene_i,
        gene2=gene_j,
        weight=weight,
        distance_threshold=distance_threshold,
        confidence=weight,
        negative_marker_constraints=negative_marker_constraints
    )
    
    return rule

def get_top_genes(adata, cell_indices, n_top):
    """Get the top expressed genes for a specific cell type."""
    # Calculate mean expression for each gene in this cell type
    mean_expr = np.mean(adata.X[cell_indices], axis=0)
    
    # Get indices of top genes
    top_gene_indices = np.argsort(mean_expr.flatten())[-n_top:].tolist()
    
    # Get gene names - convert to strings
    gene_names = [str(adata.var_names[i]) for i in top_gene_indices]
    
    return gene_names, top_gene_indices

def process_cell_type(adata, cell_type, args, output_dir):
    """Process a single cell type and generate rules."""
    logger.info(f"Processing cell type: {cell_type}")
    
    # Get cell indices for this cell type
    cell_indices = get_cell_type_indices(adata, cell_type)
    logger.info(f"Found {len(cell_indices)} cells of type {cell_type}")
    
    # Get top expressed genes for this cell type
    gene_names, gene_indices = get_top_genes(adata, cell_indices, args.top_genes)
    logger.info(f"Top {len(gene_names)} genes: {', '.join(gene_names)}")
    
    # Get binary expression matrix
    binary_matrix = get_binary_expression(adata, cell_indices, args.min_expression)
    
    # Calculate co-expression and build network
    coexpr_network = calculate_coexpression(binary_matrix, gene_indices, args.min_coexpression)
    logger.info(f"Co-expression network has {coexpr_network.number_of_nodes()} nodes and {coexpr_network.number_of_edges()} edges")
    
    # Find cliques in the network for each size from min_size to max_possible
    min_clique_size = args.min_clique_size
    all_cliques = []
    
    # Determine maximum clique size to find
    if args.max_clique_size is not None:
        # Use user-specified maximum
        max_clique_size = args.max_clique_size
        logger.info(f"Using user-specified maximum clique size: {max_clique_size}")
    else:
        # Find the largest possible clique in the network
        largest_cliques = find_cliques(coexpr_network, min_size=min_clique_size)
        max_clique_size = len(largest_cliques[0]) if largest_cliques else min_clique_size
        logger.info(f"Found maximum possible clique size: {max_clique_size}")
    
    # Generate cliques for each size from min_size to max_size
    for size in range(min_clique_size, max_clique_size + 1):
        size_cliques = find_cliques(coexpr_network, min_size=size, max_size=size)
        logger.info(f"Found {len(size_cliques)} cliques of size exactly {size}")
        all_cliques.extend(size_cliques)
    
    logger.info(f"Found a total of {len(all_cliques)} cliques of sizes {min_clique_size} to {max_clique_size}")
    
    # Identify negative markers for this cell type
    negative_markers = identify_negative_markers(adata, cell_type, 
                                              neg_marker_threshold=args.neg_marker_threshold)
    
    # Generate rules
    rules = []
    
    # Generate clique rules for all clique sizes
    for clique in all_cliques:
        rule = generate_clique_rule(clique, cell_type, args.distance_threshold, negative_markers)
        rules.append(rule)
    
    # Find strong edges not in cliques
    used_edges = set()
    for clique in all_cliques:
        for i, gene_i in enumerate(clique):
            for gene_j in clique[i+1:]:
                used_edges.add((gene_i, gene_j))
    
    # Generate pair rules for strong edges not in cliques
    for u, v, data in coexpr_network.edges(data=True):
        if (u, v) not in used_edges and data['weight'] > args.min_coexpression:
            rule = generate_pair_rule(u, v, data['weight'], cell_type, args.distance_threshold, negative_markers)
            rules.append(rule)
    
    # Save rules to file
    rules_file = os.path.join(output_dir, f"{cell_type}_rules.rq")
    with open(rules_file, 'w') as f:
        f.write("\n\n".join(rules))
    
    logger.info(f"Generated {len(rules)} rules for {cell_type}, saved to {rules_file}")
    return len(rules)

def main():
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    adata = load_data(args.input)
    
    # Get list of cell types to process
    all_cell_types = get_available_cell_types(adata)
    logger.info(f"Available cell types: {', '.join(all_cell_types)}")
    
    if args.cell_types:
        cell_types = args.cell_types
        # Validate requested cell types
        for ct in cell_types:
            if ct not in all_cell_types:
                logger.warning(f"Requested cell type '{ct}' not found in dataset")
    else:
        # Process all cell types if none specified
        cell_types = all_cell_types
        logger.info(f"No cell types specified, processing all {len(cell_types)} cell types")
    
    # Process each cell type
    total_rules = 0
    for cell_type in cell_types:
        if cell_type in all_cell_types:
            rules_count = process_cell_type(adata, cell_type, args, args.output_dir)
            total_rules += rules_count
    
    logger.info(f"Generated a total of {total_rules} rules across {len(cell_types)} cell types")
    
    # Create a summary file
    summary_file = os.path.join(args.output_dir, "rules_summary.txt")
    with open(summary_file, 'w') as f:
        f.write(f"# Rules Summary\n\n")
        f.write(f"Generated on: {pd.Timestamp.now()}\n\n")
        f.write(f"## Parameters\n")
        f.write(f"- Input file: {args.input}\n")
        f.write(f"- Cell types: {', '.join(cell_types)}\n")
        f.write(f"- Min expression: {args.min_expression}\n")
        f.write(f"- Min co-expression: {args.min_coexpression}\n")
        f.write(f"- Distance threshold: {args.distance_threshold}\n")
        f.write(f"- Top genes: {args.top_genes}\n")
        f.write(f"- Negative marker threshold: {args.neg_marker_threshold}\n")
        f.write(f"- Minimum clique size: {args.min_clique_size}\n")
        f.write(f"- Maximum clique size: {args.max_clique_size}\n\n")
        f.write(f"## Results\n")
        f.write(f"- Total rules: {total_rules}\n")
        for cell_type in cell_types:
            if cell_type in all_cell_types:
                rules_file = os.path.join(args.output_dir, f"{cell_type}_rules.rq")
                if os.path.exists(rules_file):
                    with open(rules_file, 'r') as rf:
                        rule_count = rf.read().count('CONSTRUCT')
                    f.write(f"- {cell_type}: {rule_count} rules\n")
    
    logger.info(f"Summary saved to {summary_file}")

if __name__ == "__main__":
    main() 