#!/usr/bin/env python3
"""
Example script demonstrating how to generate and apply rules for all clique sizes from 3 to N.
This script can be used as a template for testing the variable clique size functionality.
"""

import sys
import os
import logging
import pandas as pd
import numpy as np
import networkx as nx
import tempfile
from pathlib import Path

# Add the scripts directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))

# Import required modules
from generate_coexpression_rules import (
    find_cliques,
    generate_clique_rule,
    generate_pair_rule,
    calculate_coexpression
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def create_test_coexpression_graph():
    """Create a synthetic co-expression graph for testing."""
    G = nx.Graph()
    
    # Define genes
    genes = ["GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A3", "VIM", "SLC1A2", "GJA1", "ATP1A2", "FGFR3"]
    G.add_nodes_from(genes)
    
    # Add edges with varying weights to create different sized cliques
    edges = [
        # Create a complete subgraph (clique) of the first 5 genes
        ("GFAP", "AQP4", 0.9),
        ("GFAP", "S100B", 0.85),
        ("GFAP", "ALDH1L1", 0.75),
        ("GFAP", "SLC1A3", 0.7),
        ("AQP4", "S100B", 0.8),
        ("AQP4", "ALDH1L1", 0.7),
        ("AQP4", "SLC1A3", 0.65),
        ("S100B", "ALDH1L1", 0.75),
        ("S100B", "SLC1A3", 0.7),
        ("ALDH1L1", "SLC1A3", 0.85),
        
        # Add edges to create partial connections to other genes
        ("GFAP", "VIM", 0.8),
        ("AQP4", "VIM", 0.75),
        ("S100B", "VIM", 0.7),
        ("ALDH1L1", "SLC1A2", 0.65),
        ("SLC1A3", "SLC1A2", 0.75),
        ("VIM", "SLC1A2", 0.8),
        
        # Add more edges to create a clique of 3 genes
        ("GJA1", "ATP1A2", 0.85),
        ("GJA1", "FGFR3", 0.8),
        ("ATP1A2", "FGFR3", 0.75)
    ]
    
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)
    
    return G

def main():
    """Generate and test rules for all clique sizes from 3 to N."""
    logger.info("Creating synthetic co-expression graph")
    G = create_test_coexpression_graph()
    
    # Set up a temporary directory for output
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        logger.info(f"Temporary output directory: {output_dir}")
        
        # Find the maximum possible clique size
        largest_cliques = find_cliques(G, min_size=3)
        max_clique_size = len(largest_cliques[0]) if largest_cliques else 3
        logger.info(f"Maximum clique size found: {max_clique_size}")
        
        # Generate rules for each clique size
        all_rules = []
        for size in range(3, max_clique_size + 1):
            size_cliques = find_cliques(G, min_size=size, max_size=size)
            logger.info(f"Found {len(size_cliques)} cliques of size exactly {size}")
            
            # Generate a rule for each clique
            for clique in size_cliques:
                rule_id, rule_content = generate_clique_rule(
                    clique=clique,
                    cell_type="Astrocyte",
                    distance_threshold=50
                )
                
                # Save rule to file
                rule_file = output_dir / f"{rule_id}.rq"
                with open(rule_file, "w") as f:
                    f.write(f"# Rule ID: {rule_id}\n")
                    f.write("# Rule Type: CLIQUE\n")
                    f.write("# Cell Type: Astrocyte\n")
                    f.write(f"# Clique Size: {size}\n")
                    f.write(f"# Genes: {', '.join(clique)}\n\n")
                    f.write(rule_content)
                
                all_rules.append({
                    "rule_id": rule_id,
                    "rule_type": "CLIQUE",
                    "cell_type": "Astrocyte",
                    "clique_size": size,
                    "genes": clique,
                    "file_path": str(rule_file)
                })
        
        # Create rules summary
        rules_df = pd.DataFrame(all_rules)
        rules_summary_file = output_dir / "rules_summary.csv"
        rules_df.to_csv(rules_summary_file, index=False)
        
        # Print summary of generated rules
        logger.info("\n=== Rule Generation Summary ===")
        for size in range(3, max_clique_size + 1):
            size_count = rules_df[rules_df['clique_size'] == size].shape[0]
            logger.info(f"Clique size {size}: {size_count} rules")
        
        # Print the results
        logger.info(f"\nTotal rules generated: {len(all_rules)}")
        logger.info(f"Rules summary saved to: {rules_summary_file}")
        logger.info("\nExample rule IDs:")
        for idx, rule in enumerate(all_rules[:5]):
            logger.info(f"  {idx+1}. {rule['rule_id']} (size {rule['clique_size']})")
        
        if len(all_rules) > 5:
            logger.info(f"  ... and {len(all_rules) - 5} more")

if __name__ == "__main__":
    main() 