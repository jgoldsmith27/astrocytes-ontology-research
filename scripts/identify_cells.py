#!/usr/bin/env python3
"""
Apply the generated co-expression rules to identify whole cells in spatial data.
This script loads SPARQL rules and applies them to spatial data in RDF format.
"""

import os
import sys
import argparse
import logging
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Dict, Optional, Union
import json
from rdflib import Graph, URIRef, Literal, Namespace
from rdflib.plugins.sparql import prepareQuery
from rdflib.namespace import RDF, RDFS, XSD

# Import the conflict resolution system
try:
    from conflict_resolution import resolve_cell_conflicts, ConflictResolutionManager
    CONFLICT_RESOLUTION_AVAILABLE = True
except ImportError:
    CONFLICT_RESOLUTION_AVAILABLE = False
    logging.warning("Enhanced conflict resolution module not available. Using basic conflict resolution.")

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Apply co-expression rules to identify cells")
    parser.add_argument("--spatial-data", required=True, help="Path to spatial data in RDF/Turtle format")
    parser.add_argument("--rules-dir", help="Directory containing SPARQL rule files (paired and clique rules)")
    parser.add_argument("--meta-rules-file", help="Path to meta-rules file (overrides the default rules_dir/meta_rules.rq)")
    parser.add_argument("--output-dir", required=True, help="Directory to save identified cells and visualizations")
    parser.add_argument("--min-confidence", type=float, default=0.7, help="Minimum confidence for cell identification")
    parser.add_argument("--resolve-conflicts", action="store_true", help="Resolve conflicts between overlapping cells")
    parser.add_argument("--overlap-threshold", type=float, default=0.3, help="Overlap threshold for conflict detection")
    parser.add_argument("--prioritize-meta-rules", action="store_true", help="Give meta-rules priority during conflict resolution")
    parser.add_argument("--expression-data", help="Optional JSON file with gene expression data for enhanced conflict resolution")
    parser.add_argument("--enhanced-resolution", action="store_true", help="Use enhanced conflict resolution if available")
    return parser.parse_args()

def load_spatial_data(spatial_data_file: str) -> Graph:
    """Load spatial data from RDF/Turtle file."""
    logger.info(f"Loading spatial data from {spatial_data_file}")
    g = Graph()
    g.parse(spatial_data_file, format="turtle")
    logger.info(f"Loaded {len(g)} triples")
    return g

def load_rules(rules_dir: str, meta_rules_file: Optional[str] = None) -> List[str]:
    """Load SPARQL rules from files."""
    rules = []
    
    # Load rules from directory
    if rules_dir:
        logger.info(f"Loading rules from directory: {rules_dir}")
        rule_files = list(Path(rules_dir).glob("*.rq"))
        for rule_file in rule_files:
            # Skip meta_rules.rq in the directory if meta_rules_file is specified
            if meta_rules_file and rule_file.name == "meta_rules.rq":
                continue
                
            with open(rule_file, "r") as f:
                rule_content = f.read()
                # Split the file by double newlines to get individual rules
                rules.extend([r for r in rule_content.split("\n\n") if r.strip() and "CONSTRUCT" in r])
    
    # Load meta-rules if specified
    if meta_rules_file and os.path.exists(meta_rules_file):
        logger.info(f"Loading meta-rules from: {meta_rules_file}")
        with open(meta_rules_file, "r") as f:
            meta_rules_content = f.read()
            meta_rules = [r for r in meta_rules_content.split("\n\n") if r.strip() and "CONSTRUCT" in r]
            
            # Add meta-rule flag to help with prioritization
            for i, rule in enumerate(meta_rules):
                if "Meta-rule" in rule and "derivedFromMetaRule" not in rule:
                    # Insert the meta-rule flag into the CONSTRUCT block
                    construct_pos = rule.find("CONSTRUCT {")
                    if construct_pos >= 0:
                        insert_pos = rule.find("}", construct_pos)
                        if insert_pos >= 0:
                            meta_rules[i] = (rule[:insert_pos] + 
                                          "\n        cell:derivedFromMetaRule \"true\"^^xsd:boolean ;" +
                                          rule[insert_pos:])
                
            rules.extend(meta_rules)
            
    logger.info(f"Loaded {len(rules)} rules in total")
    return rules

def apply_rules(g: Graph, rules: List[str], min_confidence: float = 0.7) -> pd.DataFrame:
    """
    Apply SPARQL rules to identify cells in the spatial data.
    
    Args:
        g: RDF graph containing spatial data
        rules: List of SPARQL CONSTRUCT queries
        min_confidence: Minimum confidence score to keep an identified cell
        
    Returns:
        DataFrame containing identified cells
    """
    logger.info(f"Applying {len(rules)} rules to identify cells")
    
    # Define namespace
    ASTRO = Namespace("http://example.org/astrocytes#")
    
    # Create a new graph to store the cells
    cells_graph = Graph()
    cells_graph += g  # Start with the spatial data
    
    # Apply each rule and track stats
    rule_counts = {}
    
    for rule_idx, rule in enumerate(rules):
        try:
            # Extract rule name from comments
            rule_name = f"rule_{rule_idx}"
            if "Meta-rule combining patterns" in rule:
                rule_name = "meta_rule"
                rule_type = "meta"
            elif "Clique rule for" in rule:
                rule_type = "clique"
                # Extract cell type
                import re
                m = re.search(r"Clique rule for (\w+)", rule)
                cell_type = m.group(1) if m else "Unknown"
                rule_name = f"{cell_type}_clique_{rule_idx}"
            elif "Pair rule for" in rule:
                rule_type = "pair"
                # Extract cell type
                import re
                m = re.search(r"Pair rule for (\w+)", rule)
                cell_type = m.group(1) if m else "Unknown"
                rule_name = f"{cell_type}_pair_{rule_idx}"
            
            # Apply the rule
            cells_graph.update(rule)
            
            # Count results
            cell_count = len(list(cells_graph.subjects(RDF.type, ASTRO.SpatialCell)))
            rule_counts[rule_name] = rule_counts.get(rule_name, 0) + cell_count
            
            logger.debug(f"Applied rule {rule_idx+1}/{len(rules)}: {rule_name}, total cells: {cell_count}")
            
        except Exception as e:
            logger.error(f"Error applying rule {rule_idx+1}: {str(e)}")
    
    # Extract cells from the graph
    cells = []
    
    for cell in cells_graph.subjects(RDF.type, ASTRO.SpatialCell):
        # Get cell properties
        cell_type = cells_graph.value(cell, ASTRO.hasCellType)
        x = cells_graph.value(cell, ASTRO.hasXCoordinate)
        y = cells_graph.value(cell, ASTRO.hasYCoordinate)
        radius = cells_graph.value(cell, ASTRO.hasRadius)
        confidence = cells_graph.value(cell, ASTRO.hasConfidence)
        cell_id = cells_graph.value(cell, ASTRO.hasCellID)
        
        # Check if derived from meta-rule
        from_meta_rule = cells_graph.value(cell, ASTRO.derivedFromMetaRule)
        is_meta_rule = True if from_meta_rule and from_meta_rule.value else False
        
        # Get rule type
        if is_meta_rule:
            rule_type = "meta"
        else:
            # Try to infer from cell_id if available
            if cell_id:
                if "_clique_" in str(cell_id):
                    rule_type = "clique"
                elif "_pair_" in str(cell_id):
                    rule_type = "pair"
                else:
                    rule_type = "unknown"
            else:
                rule_type = "unknown"
        
        # Get associated genes
        genes = []
        for point in cells_graph.objects(cell, ASTRO.includesPoint):
            for gene_ref in cells_graph.objects(point, ASTRO.expressesGene):
                gene_id = cells_graph.value(gene_ref, ASTRO.hasGeneID)
                if gene_id and str(gene_id) not in genes:
                    genes.append(str(gene_id))
        
        # Skip cells with confidence below threshold
        if not confidence or float(confidence) < min_confidence:
            continue
        
        # Add to cell list
        cells.append({
            "cell_id": str(cell_id) if cell_id else str(cell),
            "cell_type": str(cell_type) if cell_type else "Unknown",
            "x": float(x) if x else 0.0,
            "y": float(y) if y else 0.0,
            "radius": float(radius) if radius else 0.0,
            "confidence": float(confidence) if confidence else 0.0,
            "rule_id": str(cell_id).split("_")[-1] if cell_id else "",
            "rule_type": rule_type,
            "genes": ",".join(genes),
            "derived_from_meta_rule": is_meta_rule
        })
    
    logger.info(f"Identified {len(cells)} cells with confidence >= {min_confidence}")
    
    # Create DataFrame
    df = pd.DataFrame(cells)
    
    # Log rule statistics
    logger.info("Rule application statistics:")
    for rule_name, count in sorted(rule_counts.items(), key=lambda x: x[0]):
        logger.info(f"  {rule_name}: {count} cells")
    
    return df

def basic_conflict_resolution(cells_df: pd.DataFrame, 
                            overlap_threshold: float = 0.3,
                            prioritize_meta_rules: bool = True) -> pd.DataFrame:
    """
    Basic conflict resolution strategy for overlapping cells.
    
    This is a simplified version used as fallback when the enhanced module is not available.
    
    Args:
        cells_df: DataFrame with identified cells
        overlap_threshold: Minimum overlap to consider a conflict
        prioritize_meta_rules: Whether to give meta-rules priority
        
    Returns:
        DataFrame with conflicts resolved
    """
    # If no cells or just one cell, return as-is
    if len(cells_df) <= 1:
        return cells_df
    
    # Create a list to track which cells to keep
    keep_cells = [True] * len(cells_df)
    
    # Process cells with meta-rules first if prioritized
    if prioritize_meta_rules:
        # Sort by whether derived from meta-rule (True first), then by confidence
        sorted_indices = cells_df.sort_values(
            by=["derived_from_meta_rule", "confidence"], 
            ascending=[False, False]
        ).index.tolist()
    else:
        # Sort by confidence only
        sorted_indices = cells_df.sort_values(by="confidence", ascending=False).index.tolist()
    
    # Process cells in order
    for i, idx1 in enumerate(sorted_indices[:-1]):
        # Skip if this cell is already marked for removal
        if not keep_cells[idx1]:
            continue
            
        # Get cell properties
        x1, y1 = cells_df.loc[idx1, ["x", "y"]]
        r1 = cells_df.loc[idx1, "radius"]
        
        # Check for conflicts with remaining cells
        for idx2 in sorted_indices[i+1:]:
            # Skip if this cell is already marked for removal
            if not keep_cells[idx2]:
                continue
                
            # Get cell properties
            x2, y2 = cells_df.loc[idx2, ["x", "y"]]
            r2 = cells_df.loc[idx2, "radius"]
            
            # Calculate distance and overlap
            distance = ((x1 - x2)**2 + (y1 - y2)**2)**0.5
            
            # If centers are far apart, no overlap
            if distance >= r1 + r2:
                continue
            
            # Calculate overlap
            if distance <= abs(r1 - r2):
                # One cell is completely inside the other
                smaller_r = min(r1, r2)
                larger_r = max(r1, r2)
                overlap = (smaller_r / larger_r)**2
            else:
                # Calculate overlap percentage (approximation)
                overlap = 1.0 - (distance / (r1 + r2))
            
            # If overlap exceeds threshold, mark the lower confidence cell for removal
            if overlap >= overlap_threshold:
                keep_cells[idx2] = False
    
    # Apply the filtering
    result_df = cells_df[keep_cells].reset_index(drop=True)
    
    logger.info(f"Basic conflict resolution: kept {len(result_df)} cells out of {len(cells_df)}")
    return result_df

def visualize_cells(cells_df: pd.DataFrame, output_dir: str):
    """Generate visualizations of identified cells."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Create a scatter plot of all cells
    plt.figure(figsize=(10, 10))
    
    # Plot cells by type if available, else by rule type
    if len(cells_df["cell_type"].unique()) > 1:
        for cell_type, group in cells_df.groupby("cell_type"):
            plt.scatter(group["x"], group["y"], label=cell_type, s=30, alpha=0.7)
            
            # Draw circles representing cell boundaries
            for _, cell in group.iterrows():
                circle = plt.Circle((cell["x"], cell["y"]), cell["radius"], 
                                 fill=False, alpha=0.3)
                plt.gca().add_patch(circle)
    else:
        # Distinguish by rule type
        for rule_type, group in cells_df.groupby("rule_type"):
            plt.scatter(group["x"], group["y"], label=rule_type, s=30, alpha=0.7)
            
            # Draw circles representing cell boundaries
            for _, cell in group.iterrows():
                circle = plt.Circle((cell["x"], cell["y"]), cell["radius"], 
                                  fill=False, alpha=0.3)
                plt.gca().add_patch(circle)
    
    plt.legend()
    plt.title(f"Identified Cells (n={len(cells_df)})")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.axis("equal")
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_dir, "identified_cells.png"), dpi=300)
    plt.close()
    
    # Create a heatmap of cell density (using hexbin plot)
    plt.figure(figsize=(10, 10))
    plt.hexbin(cells_df["x"], cells_df["y"], gridsize=50, cmap="viridis")
    plt.colorbar(label="Cell Density")
    plt.title("Cell Density Heatmap")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_dir, "cell_density.png"), dpi=300)
    plt.close()
    
    # Create distribution plots
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    
    # Cell types distribution
    type_counts = cells_df["cell_type"].value_counts()
    type_counts.plot(kind="bar", ax=axs[0, 0])
    axs[0, 0].set_title("Cell Types")
    axs[0, 0].set_ylabel("Count")
    axs[0, 0].tick_params(axis="x", rotation=45)
    
    # Rule types distribution
    rule_counts = cells_df["rule_type"].value_counts()
    rule_counts.plot(kind="bar", ax=axs[0, 1])
    axs[0, 1].set_title("Rule Types")
    axs[0, 1].set_ylabel("Count")
    
    # Confidence distribution
    axs[1, 0].hist(cells_df["confidence"], bins=20)
    axs[1, 0].set_title("Confidence Scores")
    axs[1, 0].set_xlabel("Confidence")
    axs[1, 0].set_ylabel("Count")
    
    # Cell radius distribution
    axs[1, 1].hist(cells_df["radius"], bins=20)
    axs[1, 1].set_title("Cell Radii")
    axs[1, 1].set_xlabel("Radius")
    axs[1, 1].set_ylabel("Count")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cell_statistics.png"), dpi=300)
    plt.close()
    
    logger.info(f"Visualizations saved to {output_dir}")

def main():
    """Main function to run the cell identification pipeline."""
    args = parse_args()
    
    # Create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Load spatial data
    spatial_graph = load_spatial_data(args.spatial_data)
    
    # Load rules
    rules = load_rules(args.rules_dir, args.meta_rules_file)
    
    # Apply rules to identify cells
    cells_df = apply_rules(spatial_graph, rules, args.min_confidence)
    
    # Resolve conflicts if needed
    if args.resolve_conflicts:
        # Use enhanced conflict resolution if available and requested
        if CONFLICT_RESOLUTION_AVAILABLE and args.enhanced_resolution:
            logger.info("Using enhanced conflict resolution")
            
            # Load gene expression data if provided
            cell_type_expression_data = None
            if args.expression_data and os.path.exists(args.expression_data):
                with open(args.expression_data, 'r') as f:
                    cell_type_expression_data = json.load(f)
            
            # Apply enhanced conflict resolution
            cells_df = resolve_cell_conflicts(
                cells_df=cells_df,
                output_dir=args.output_dir,
                overlap_threshold=args.overlap_threshold,
                same_type_overlap_threshold=args.overlap_threshold * 1.5,  # Higher threshold for same type
                cell_type_expression_data=cell_type_expression_data,
                generate_visuals=True
            )
        else:
            logger.info("Using basic conflict resolution")
            cells_df = basic_conflict_resolution(
                cells_df=cells_df,
                overlap_threshold=args.overlap_threshold,
                prioritize_meta_rules=args.prioritize_meta_rules
            )
    
    # Save results
    output_file = os.path.join(args.output_dir, "identified_cells.csv")
    cells_df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(cells_df)} identified cells to {output_file}")
    
    # Create visualizations
    visualize_cells(cells_df, args.output_dir)
    
    # Output summary
    print(f"\nCell Identification Results:")
    print(f"----------------------------")
    print(f"Total cells identified: {len(cells_df)}")
    print(f"Cell types found: {', '.join(cells_df['cell_type'].unique())}")
    print(f"Rule types: {', '.join(cells_df['rule_type'].unique())}")
    print(f"Results saved to {args.output_dir}")

if __name__ == "__main__":
    main() 