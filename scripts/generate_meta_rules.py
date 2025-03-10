#!/usr/bin/env python3
"""
Generate meta-rules by analyzing relationships between existing co-expression rules.
This script identifies significant associations between rules and creates higher-order
patterns that can be used for more specific cell identification.
"""

import os
import re
import argparse
import logging
import pandas as pd
import numpy as np
import networkx as nx
import scipy.stats
from pathlib import Path
from collections import defaultdict
import rdflib
from rdflib import Graph, Namespace, Literal, URIRef, BNode
from rdflib.namespace import RDF, RDFS, XSD

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Generate meta-rules based on rule co-occurrence.")
    parser.add_argument("--rules-dir", type=str, required=True, help="Directory containing SPARQL rule files")
    parser.add_argument("--single-cell-data", type=str, required=True, help="Path to h5ad file with single-cell data")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to save generated meta-rules")
    parser.add_argument("--min-lift", type=float, default=3.0, 
                        help="Minimum lift value to consider a rule association significant")
    parser.add_argument("--max-p-value", type=float, default=0.01, 
                        help="Maximum p-value for Fisher's exact test to consider association significant")
    parser.add_argument("--spatial-distance", type=int, default=75, 
                        help="Maximum spatial distance between rule patterns to form a meta-rule")
    parser.add_argument("--confidence-boost", type=float, default=0.1, 
                        help="Amount to boost confidence score for meta-rules")
    parser.add_argument("--same-cell-type-only", action="store_true", 
                        help="Only create meta-rules between patterns from the same cell type")
    return parser.parse_args()

def extract_rule_info(rule_file):
    """
    Extract key information from a SPARQL rule file.
    
    Returns:
        list: List of dictionaries containing rule information
    """
    rule_info_list = []
    
    with open(rule_file, 'r') as f:
        content = f.read()
        
    # Split the file into individual rules
    rules = content.split('\n\n')
    
    for rule in rules:
        if not rule.strip() or 'CONSTRUCT' not in rule:
            continue
            
        # Extract rule type (clique or pair)
        if 'Clique rule' in rule:
            rule_type = 'clique'
            # Extract cell type
            cell_type_match = re.search(r'Clique rule for (\w+)', rule)
            cell_type = cell_type_match.group(1) if cell_type_match else "Unknown"
            # Extract genes
            genes_match = re.search(r'with \d+ genes: ([\w, ]+)', rule)
            genes = genes_match.group(1).split(', ') if genes_match else []
            # Default confidence for clique rules
            confidence = 0.9
        elif 'Pair rule' in rule:
            rule_type = 'pair'
            # Extract cell type
            cell_type_match = re.search(r'Pair rule for (\w+)', rule)
            cell_type = cell_type_match.group(1) if cell_type_match else "Unknown"
            # Extract genes
            genes_match = re.search(r'with genes: (\w+) and (\w+)', rule)
            genes = [genes_match.group(1), genes_match.group(2)] if genes_match else []
            # Extract confidence
            conf_match = re.search(r'coexpression = (\d+\.\d+)', rule)
            confidence = float(conf_match.group(1)) if conf_match else 0.6
        else:
            continue
        
        # Create unique rule ID
        rule_id = f"{cell_type}_{rule_type}_{'_'.join(genes)}"
        
        # Store rule information
        rule_info = {
            'rule_id': rule_id,
            'rule_type': rule_type,
            'cell_type': cell_type,
            'genes': genes,
            'confidence': confidence,
            'rule_text': rule
        }
        
        rule_info_list.append(rule_info)
    
    return rule_info_list

def load_all_rules(rules_dir):
    """
    Load all SPARQL rules from the rules directory.
    
    Returns:
        list: List of rule information dictionaries
    """
    all_rules = []
    rule_files = list(Path(rules_dir).glob('*.rq'))
    
    for rule_file in rule_files:
        logger.info(f"Loading rules from {rule_file}")
        rules = extract_rule_info(rule_file)
        all_rules.extend(rules)
    
    logger.info(f"Loaded {len(all_rules)} rules in total")
    return all_rules

def apply_rules_to_single_cell(rules, single_cell_data):
    """
    Apply rules to single-cell data to determine which cells match each rule.
    This allows us to analyze rule co-occurrence patterns.
    
    Returns:
        dict: Dictionary mapping rule IDs to sets of matching cell indices
    """
    import scanpy as sc
    
    logger.info(f"Loading single-cell data from {single_cell_data}")
    adata = sc.read_h5ad(single_cell_data)
    
    # Dictionary to store cells matching each rule
    rule_matches = {}
    
    # Group rules by cell type
    rules_by_cell_type = defaultdict(list)
    for rule in rules:
        rules_by_cell_type[rule['cell_type']].append(rule)
    
    # For each cell type, apply its rules
    for cell_type, cell_type_rules in rules_by_cell_type.items():
        logger.info(f"Analyzing rule matches for {cell_type}")
        
        # Get indices of cells of this type
        cell_indices = np.where(adata.obs['cell_type'] == cell_type)[0]
        if len(cell_indices) == 0:
            logger.warning(f"No cells found of type {cell_type}")
            continue
        
        # Get expression data for these cells
        expr_data = adata.X[cell_indices]
        
        # For each rule, find matching cells
        for rule in cell_type_rules:
            rule_id = rule['rule_id']
            genes = rule['genes']
            
            # Get indices of genes in the expression matrix
            gene_indices = [np.where(adata.var_names == gene)[0][0] for gene in genes if gene in adata.var_names]
            
            if len(gene_indices) != len(genes):
                missing = [gene for gene in genes if gene not in adata.var_names]
                logger.warning(f"Some genes not found in dataset: {missing}")
                continue
            
            # Find cells where all genes in the rule are expressed (above threshold)
            thresh = 0.2  # Same threshold used in rule generation
            matching_cells = set()
            
            for i, cell_idx in enumerate(cell_indices):
                if all(expr_data[i, gene_idx] > thresh for gene_idx in gene_indices):
                    matching_cells.add(cell_idx)
            
            rule_matches[rule_id] = matching_cells
            logger.info(f"Rule {rule_id} matches {len(matching_cells)} cells")
    
    return rule_matches

def calculate_rule_associations(rules, rule_matches, min_lift=3.0, max_p_value=0.01, same_cell_type_only=True):
    """
    Calculate statistical associations between rules based on their co-occurrence in cells.
    
    Returns:
        list: List of significant rule pair associations with statistics
    """
    logger.info("Calculating rule associations")
    
    # Get total number of cells
    all_cell_indices = set()
    for cells in rule_matches.values():
        all_cell_indices.update(cells)
    total_cells = len(all_cell_indices)
    
    # Dictionary to map rule IDs to rule objects
    rule_dict = {rule['rule_id']: rule for rule in rules}
    
    # List to store significant associations
    associations = []
    
    # Calculate associations for each rule pair
    rule_ids = list(rule_matches.keys())
    for i in range(len(rule_ids)):
        for j in range(i+1, len(rule_ids)):
            rule_id1 = rule_ids[i]
            rule_id2 = rule_ids[j]
            
            # Skip if rules are from different cell types and same_cell_type_only is True
            if same_cell_type_only:
                if rule_dict[rule_id1]['cell_type'] != rule_dict[rule_id2]['cell_type']:
                    continue
            
            # Get cells matching each rule
            cells1 = rule_matches[rule_id1]
            cells2 = rule_matches[rule_id2]
            
            # Skip if either rule doesn't match any cells
            if not cells1 or not cells2:
                continue
            
            # Calculate co-occurrence statistics
            n1 = len(cells1)
            n2 = len(cells2)
            n12 = len(cells1.intersection(cells2))
            
            # Skip if no co-occurrence
            if n12 == 0:
                continue
                
            # Expected co-occurrence under independence
            expected = (n1 * n2) / total_cells
            
            # Calculate lift
            lift = n12 / expected if expected > 0 else float('inf')
            
            # Calculate conditional probabilities
            p1_given_2 = n12 / n2
            p2_given_1 = n12 / n1
            
            # Perform Fisher's exact test
            contingency_table = np.array([
                [n12, n1 - n12],
                [n2 - n12, total_cells - n1 - n2 + n12]
            ])
            _, p_value = scipy.stats.fisher_exact(contingency_table)
            
            # Check if association is significant
            if lift >= min_lift and p_value <= max_p_value:
                association = {
                    'rule_id1': rule_id1,
                    'rule_id2': rule_id2,
                    'rule1': rule_dict[rule_id1],
                    'rule2': rule_dict[rule_id2],
                    'n_cells1': n1,
                    'n_cells2': n2,
                    'n_co_occurrence': n12,
                    'lift': lift,
                    'p_value': p_value,
                    'p1_given_2': p1_given_2,
                    'p2_given_1': p2_given_1
                }
                associations.append(association)
    
    # Sort by lift (descending)
    associations.sort(key=lambda x: x['lift'], reverse=True)
    
    logger.info(f"Found {len(associations)} significant rule associations")
    return associations

def generate_meta_rule(rule1, rule2, confidence, spatial_distance=75, cell_type=None):
    """Generate a meta-rule from two co-occurring rules."""
    # Use the more specific cell type if provided, otherwise use the cell type from rule1
    if cell_type is None:
        cell_type = rule1['cell_type']
    
    # Get genes involved in each rule
    genes1 = rule1['genes']
    genes2 = rule2['genes']
    
    # Create a unique ID for the meta-rule
    rule_id = f"meta_rule_{rule1['rule_id']}_{rule2['rule_id']}"
    
    # Generate pattern matching components for each rule
    pattern1 = generate_pattern_matching(1, genes1, pattern_type=rule1['rule_type'])
    pattern2 = generate_pattern_matching(2, genes2, pattern_type=rule2['rule_type'])
    
    # Build point links for SPARQL
    point_links = ""
    for i in range(len(genes1)):
        point_links += f"  ?cell cell:includesPoint ?point1_{i} .\n"
    for i in range(len(genes2)):
        point_links += f"  ?cell cell:includesPoint ?point2_{i} .\n"
    
    # Template for the meta-rule
    meta_rule_template = """
# Meta-rule combining patterns: {rule1_id} and {rule2_id}
# Cell Type: {cell_type}, Combined Confidence: {confidence:.2f}
PREFIX cell: <http://example.org/ontology/cell/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

CONSTRUCT {{
  ?cell a cell:SpatialCell ;
        cell:hasCellType "{cell_type}" ;
        cell:hasXCoordinate ?centerX ;
        cell:hasYCoordinate ?centerY ;
        cell:hasRadius ?combinedRadius ;
        cell:hasConfidence "{confidence}"^^xsd:decimal ;
        cell:hasCellID ?cellID ;
        cell:derivedFromMetaRule "true"^^xsd:boolean .
  
  # Link to the patterns
  ?cell cell:includesPattern ?pattern1 .
  ?cell cell:includesPattern ?pattern2 .
  
  # Link to all included points
{point_links}
}}
WHERE {{
  # Pattern 1: {rule1_id}
  {pattern1}
  
  # Pattern 2: {rule2_id}
  {pattern2}
  
  # Calculate center of the combined cell
  BIND((?centerX1 + ?centerX2) / 2 AS ?centerX)
  BIND((?centerY1 + ?centerY2) / 2 AS ?centerY)
  
  # Check distance between patterns
  BIND(SQRT(POW(?centerX1 - ?centerX2, 2) + POW(?centerY1 - ?centerY2, 2)) AS ?patternDist)
  FILTER(?patternDist < {spatial_distance})
  
  # Calculate combined radius
  BIND(?patternDist / 2 + MAX(?radius1, ?radius2) AS ?combinedRadius)
  
  # Generate a unique ID for the cell
  BIND(CONCAT("meta_", STRUUID()) AS ?cellID)
}}
"""
    
    # Format the template with the specific values
    meta_rule = meta_rule_template.format(
        rule1_id=rule1['rule_id'],
        rule2_id=rule2['rule_id'],
        cell_type=cell_type,
        confidence=confidence,
        spatial_distance=spatial_distance,
        pattern1=pattern1,
        pattern2=pattern2,
        point_links=point_links
    )
    
    return meta_rule

def generate_pattern_matching(pattern_num, genes, pattern_type='clique'):
    """Generate the SPARQL pattern matching for a rule pattern."""
    pattern_where = ""
    
    # Add center variables
    pattern_where += f"  BIND(?centerX{pattern_num} AS ?centerX{pattern_num})\n"
    pattern_where += f"  BIND(?centerY{pattern_num} AS ?centerY{pattern_num})\n"
    pattern_where += f"  BIND(?radius{pattern_num} AS ?radius{pattern_num})\n\n"
    
    # Generate point patterns
    for i, gene in enumerate(genes):
        point_var = f"?point{pattern_num}_{i}"
        gene_var = f"?gene{pattern_num}_{i}"
        x_var = f"?x{pattern_num}_{i}"
        y_var = f"?y{pattern_num}_{i}"
        
        # Add spatial data point pattern
        pattern_where += f"""  {point_var} a cell:SpatialDataPoint ;
        cell:expressesGene {gene_var} ;
        cell:hasXCoordinate {x_var} ;
        cell:hasYCoordinate {y_var} .
  {gene_var} cell:hasGeneID "{gene}" .
  ?pattern{pattern_num}Cell cell:includesPoint {point_var} .
  
"""
    
    # Add pattern cell
    pattern_where += f"?pattern{pattern_num}Cell a cell:SpatialCell .\n    "
    
    # Add pattern-specific constraints based on type
    if pattern_type == 'pair' and len(genes) == 2:
        # For pair rules, we specifically want to find points that match the pattern
        pattern_where += f"""
  ?point{pattern_num}_0 a cell:SpatialDataPoint ;
        cell:expressesGene ?gene{pattern_num}_0 ;
        cell:hasXCoordinate ?x{pattern_num}_0 ;
        cell:hasYCoordinate ?y{pattern_num}_0 .
  ?gene{pattern_num}_0 cell:hasGeneID "{genes[0]}" .
  
  ?point{pattern_num}_1 a cell:SpatialDataPoint ;
        cell:expressesGene ?gene{pattern_num}_1 ;
        cell:hasXCoordinate ?x{pattern_num}_1 ;
        cell:hasYCoordinate ?y{pattern_num}_1 .
  ?gene{pattern_num}_1 cell:hasGeneID "{genes[1]}" .
  
  ?pattern{pattern_num}Cell a cell:SpatialCell ;
        cell:includesPoint ?point{pattern_num}_0 ;
        cell:includesPoint ?point{pattern_num}_1 .
"""
        
    # Add center and radius calculation
    if len(genes) >= 2:
        # Calculate pattern center (average of all points)
        x_vars = [f"?x{pattern_num}_{i}" for i in range(len(genes))]
        y_vars = [f"?y{pattern_num}_{i}" for i in range(len(genes))]
        
        x_sum = " + ".join(x_vars)
        y_sum = " + ".join(y_vars)
        
        pattern_where += f"""
  # Calculate pattern center
  BIND(({x_sum}) / {len(genes)} AS ?centerX{pattern_num})
  BIND(({y_sum}) / {len(genes)} AS ?centerY{pattern_num})
  
  # Calculate pattern radius
"""
        
        # Calculate max distance from center to any point
        for i in range(len(genes)):
            pattern_where += f"""  BIND(SQRT(POW(?centerX{pattern_num} - ?x{pattern_num}_{i}, 2) + 
                     POW(?centerY{pattern_num} - ?y{pattern_num}_{i}, 2)) AS ?dist{pattern_num}_{i})
"""
        
        # Get maximum distance as radius
        dist_vars = [f"?dist{pattern_num}_{i}" for i in range(len(genes))]
        if len(dist_vars) > 1:
            max_dist = f"MAX({', '.join(dist_vars)})"
        else:
            max_dist = dist_vars[0]
        
        pattern_where += f"""  BIND({max_dist} AS ?radius{pattern_num})
"""
    
    return pattern_where

def main():
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load all rules
    all_rules = load_all_rules(args.rules_dir)
    
    # Apply rules to single-cell data to analyze co-occurrence
    rule_matches = apply_rules_to_single_cell(all_rules, args.single_cell_data)
    
    # Calculate rule associations
    associations = calculate_rule_associations(
        all_rules, 
        rule_matches, 
        min_lift=args.min_lift, 
        max_p_value=args.max_p_value,
        same_cell_type_only=args.same_cell_type_only
    )
    
    # Generate meta-rules for significant associations
    meta_rules = []
    for assoc in associations:
        meta_rule = generate_meta_rule(
            assoc['rule1'],
            assoc['rule2'],
            assoc['confidence'],
            spatial_distance=args.spatial_distance,
            cell_type=assoc['rule1']['cell_type']
        )
        meta_rules.append(meta_rule)
    
    # Save meta-rules to files
    meta_rules_file = os.path.join(args.output_dir, "meta_rules.rq")
    with open(meta_rules_file, 'w') as f:
        f.write("\n\n".join(meta_rules))
    
    # Save a summary of rule associations
    summary_file = os.path.join(args.output_dir, "meta_rules_summary.csv")
    summary_data = []
    for assoc in associations:
        summary_data.append({
            'rule_id1': assoc['rule_id1'],
            'rule_id2': assoc['rule_id2'],
            'cell_type': assoc['rule1']['cell_type'],
            'n_cells1': assoc['n_cells1'],
            'n_cells2': assoc['n_cells2'],
            'n_co_occurrence': assoc['n_co_occurrence'],
            'lift': assoc['lift'],
            'p_value': assoc['p_value'],
            'p1_given_2': assoc['p1_given_2'],
            'p2_given_1': assoc['p2_given_1']
        })
    
    pd.DataFrame(summary_data).to_csv(summary_file, index=False)
    
    logger.info(f"Generated {len(meta_rules)} meta-rules, saved to {meta_rules_file}")
    logger.info(f"Rule association summary saved to {summary_file}")
    
    # Generate a visualization of the rule association network
    G = nx.Graph()
    
    # Add nodes (rules)
    for rule in all_rules:
        rule_id = rule['rule_id']
        G.add_node(rule_id, 
                  cell_type=rule['cell_type'], 
                  rule_type=rule['rule_type'],
                  confidence=rule['confidence'],
                  genes=','.join(rule['genes']))
    
    # Add edges (associations)
    for assoc in associations:
        G.add_edge(assoc['rule_id1'], 
                   assoc['rule_id2'], 
                   lift=assoc['lift'], 
                   p_value=assoc['p_value'],
                   co_occurrence=assoc['n_co_occurrence'])
    
    # Save network as GraphML for visualization in external tools
    network_file = os.path.join(args.output_dir, "rule_network.graphml")
    nx.write_graphml(G, network_file)
    
    logger.info(f"Rule network saved to {network_file}")
    
    # Basic reporting
    print(f"\nMeta-Rule Generation Summary:")
    print(f"------------------------------")
    print(f"Found {len(associations)} significant rule associations")
    print(f"Generated {len(meta_rules)} meta-rules")
    print(f"\nTop 5 strongest associations:")
    for i, assoc in enumerate(associations[:5]):
        print(f"{i+1}. {assoc['rule_id1']} + {assoc['rule_id2']} (lift = {assoc['lift']:.2f}, co-occurrence = {assoc['n_co_occurrence']})")
    
    print(f"\nResults saved to:")
    print(f"- Meta-rules: {meta_rules_file}")
    print(f"- Association summary: {summary_file}")
    print(f"- Rule network: {network_file}")

if __name__ == "__main__":
    main() 