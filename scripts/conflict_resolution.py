#!/usr/bin/env python3
"""
Enhanced conflict resolution system for cell identification in spatial transcriptomics data.
This module provides sophisticated resolution of conflicts between overlapping cells identified
by different co-expression rules.
"""

import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple, Optional, Any

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class Cell:
    """Represents an identified cell with all its properties."""
    cell_id: str
    cell_type: str
    x: float
    y: float
    radius: float
    confidence: float
    rule_id: str
    rule_type: str  # 'meta', 'clique', or 'pair'
    genes: List[str]
    negative_markers_missing: bool = True
    derived_from_meta_rule: bool = False
    
    def distance_to(self, other: 'Cell') -> float:
        """Calculate Euclidean distance to another cell."""
        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)
    
    def overlap_percentage(self, other: 'Cell') -> float:
        """Calculate percentage of overlap with another cell."""
        d = self.distance_to(other)
        # If centers are far apart, no overlap
        if d >= self.radius + other.radius:
            return 0.0
        
        # If one cell is completely inside the other
        if d <= abs(self.radius - other.radius):
            smaller_radius = min(self.radius, other.radius)
            larger_radius = max(self.radius, other.radius)
            return (smaller_radius / larger_radius)**2
        
        # Partial overlap - calculate using circle-circle intersection formula
        r1, r2 = self.radius, other.radius
        d2 = d**2
        area1 = np.pi * r1**2
        area2 = np.pi * r2**2
        
        # Calculate intersection area using circle-circle intersection formula
        term1 = r1**2 * np.arccos((d2 + r1**2 - r2**2) / (2 * d * r1))
        term2 = r2**2 * np.arccos((d2 + r2**2 - r1**2) / (2 * d * r2))
        term3 = 0.5 * np.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
        
        intersection_area = term1 + term2 - term3
        smaller_area = min(area1, area2)
        
        return intersection_area / smaller_area

    def get_resolution_score(self, cell_type_hierarchy: Dict[str, int], 
                           gene_specificity: Dict[str, float]) -> float:
        """
        Calculate a comprehensive score for conflict resolution.
        Higher scores indicate stronger evidence for this cell.
        """
        # Start with the base confidence score
        score = self.confidence
        
        # Adjust based on rule type
        rule_type_bonus = {
            'meta': 0.15,  # Meta-rules get highest bonus
            'clique': 0.1,  # Clique rules get medium bonus
            'pair': 0.0     # Pair rules get no bonus
        }
        score += rule_type_bonus.get(self.rule_type, 0.0)
        
        # Adjust based on cell type hierarchy (rarer/more specific types get bonus)
        hierarchy_level = cell_type_hierarchy.get(self.cell_type, 0)
        score += 0.05 * hierarchy_level
        
        # Adjust based on gene specificity (more specific genes give higher scores)
        gene_spec_score = sum(gene_specificity.get(gene, 0.5) for gene in self.genes) / len(self.genes)
        score += 0.1 * gene_spec_score
        
        # Bonus for negative marker validation
        if self.negative_markers_missing:
            score += 0.05
            
        # Bonus for meta-rule derivation
        if self.derived_from_meta_rule:
            score += 0.05
            
        # Bonus based on gene count (more genes = more evidence)
        score += min(0.1, 0.02 * len(self.genes))
        
        return score


class ConflictResolutionManager:
    """Manages the enhanced conflict resolution process for identified cells."""
    
    def __init__(self, 
                 overlap_threshold: float = 0.3,
                 same_type_overlap_threshold: float = 0.7,
                 generate_visuals: bool = True,
                 output_dir: str = None):
        """
        Initialize the conflict resolution manager.
        
        Args:
            overlap_threshold: Minimum overlap percentage to consider a conflict
            same_type_overlap_threshold: Higher threshold for cells of same type
            generate_visuals: Whether to generate visualizations of conflicts
            output_dir: Directory to save visualizations and reports
        """
        self.overlap_threshold = overlap_threshold
        self.same_type_overlap_threshold = same_type_overlap_threshold
        self.generate_visuals = generate_visuals
        self.output_dir = output_dir
        
        # Cell type hierarchy (higher value = more specific/rarer)
        # This would ideally be loaded from a configuration file
        self.cell_type_hierarchy = {
            "Astrocytes": 1,
            "Astrocytes1": 2,
            "Microglia": 2,
            "Neuron": 1,
            "Oligodendrocyte": 1,
        }
        
        # Gene specificity scores (higher = more specific to certain cell types)
        # This would ideally be derived from single-cell data analysis
        self.gene_specificity = {}
        
        # Storage for conflict analysis
        self.cells = []
        self.conflict_groups = []
        self.conflict_stats = {}
        
        # Allow adding to these maps programmatically
        if os.path.exists('configs/cell_type_hierarchy.json'):
            import json
            with open('configs/cell_type_hierarchy.json', 'r') as f:
                self.cell_type_hierarchy.update(json.load(f))
                
        if os.path.exists('configs/gene_specificity.json'):
            import json
            with open('configs/gene_specificity.json', 'r') as f:
                self.gene_specificity.update(json.load(f))
        
    def load_cells(self, cells_df: pd.DataFrame) -> None:
        """
        Load cells from a DataFrame into our internal representation.
        
        Args:
            cells_df: DataFrame containing identified cells with their properties
        """
        self.cells = []
        
        for _, row in cells_df.iterrows():
            # Extract genes from the gene column (comma-separated string)
            genes = row.get('genes', '').split(',') if isinstance(row.get('genes', ''), str) else []
            
            # Create a Cell object
            cell = Cell(
                cell_id=row['cell_id'],
                cell_type=row['cell_type'],
                x=row['x'],
                y=row['y'],
                radius=row['radius'],
                confidence=row['confidence'],
                rule_id=row.get('rule_id', ''),
                rule_type=row.get('rule_type', 'unknown'),
                genes=genes,
                negative_markers_missing=row.get('negative_markers_missing', True),
                derived_from_meta_rule=row.get('derived_from_meta_rule', False)
            )
            
            self.cells.append(cell)
            
        logger.info(f"Loaded {len(self.cells)} cells for conflict resolution")
        
    def calculate_gene_specificity(self, cell_type_gene_expression: Dict[str, Dict[str, float]]) -> None:
        """
        Calculate gene specificity scores based on expression patterns.
        
        Args:
            cell_type_gene_expression: Nested dict with mean expression of each gene in each cell type
        """
        all_types = list(cell_type_gene_expression.keys())
        all_genes = set()
        
        # Get all genes
        for gene_dict in cell_type_gene_expression.values():
            all_genes.update(gene_dict.keys())
            
        # Calculate specificity
        for gene in all_genes:
            # Get expression values for this gene across all cell types
            expressions = [cell_type_gene_expression[ct].get(gene, 0.0) for ct in all_types]
            if sum(expressions) == 0:
                self.gene_specificity[gene] = 0.5  # Default if no expression
                continue
                
            # Normalize
            expressions = [e / sum(expressions) for e in expressions]
            
            # Calculate entropy (lower entropy = more specific)
            from scipy.stats import entropy
            gene_entropy = entropy(expressions)
            
            # Convert to specificity score (1 = very specific, 0 = ubiquitous)
            # Max entropy = log(n) where n is number of cell types
            max_entropy = np.log(len(all_types))
            specificity = 1 - (gene_entropy / max_entropy)
            
            self.gene_specificity[gene] = specificity
    
    def find_conflicts(self) -> None:
        """
        Identify all conflicts between cells based on spatial overlap.
        Groups conflicting cells together for resolution.
        """
        n_cells = len(self.cells)
        # Track which cells are in conflict with each other
        conflict_graph = nx.Graph()
        
        # Add all cells as nodes
        for i, cell in enumerate(self.cells):
            conflict_graph.add_node(i, **cell.__dict__)
        
        # Find conflicts
        for i in range(n_cells):
            for j in range(i+1, n_cells):
                cell_i = self.cells[i]
                cell_j = self.cells[j]
                
                # Calculate overlap
                overlap = cell_i.overlap_percentage(cell_j)
                
                # Determine threshold based on whether cells are same type
                threshold = (self.same_type_overlap_threshold 
                             if cell_i.cell_type == cell_j.cell_type 
                             else self.overlap_threshold)
                
                # If overlap exceeds threshold, add edge to conflict graph
                if overlap >= threshold:
                    conflict_graph.add_edge(i, j, overlap=overlap)
        
        # Find connected components (groups of conflicting cells)
        self.conflict_groups = list(nx.connected_components(conflict_graph))
        
        # Generate conflict statistics
        self.conflict_stats = {
            'total_cells': n_cells,
            'cells_in_conflict': sum(len(group) for group in self.conflict_groups),
            'conflict_groups': len(self.conflict_groups),
            'by_type': defaultdict(int)
        }
        
        # Count conflicts by cell type combinations
        for group in self.conflict_groups:
            types = tuple(sorted(set(self.cells[i].cell_type for i in group)))
            self.conflict_stats['by_type'][types] += 1
            
        logger.info(f"Found {len(self.conflict_groups)} conflict groups involving "
                   f"{self.conflict_stats['cells_in_conflict']} cells")
    
    def resolve_conflicts(self) -> List[Cell]:
        """
        Resolve all identified conflicts and return the final set of cells.
        
        Returns:
            List of cells after conflict resolution
        """
        if not self.conflict_groups:
            logger.warning("No conflicts to resolve. Call find_conflicts() first.")
            return self.cells
        
        # Track which cell indices to keep
        cells_to_keep = set(range(len(self.cells)))
        
        # Process each conflict group
        for group in self.conflict_groups:
            # Skip groups with only one cell (no actual conflict)
            if len(group) <= 1:
                continue
                
            # Get cells in this conflict group
            group_cells = [self.cells[i] for i in group]
            
            # Calculate resolution scores for each cell
            scores = [cell.get_resolution_score(self.cell_type_hierarchy, self.gene_specificity) 
                     for cell in group_cells]
            
            # Find cell with highest score
            best_idx = np.argmax(scores)
            best_cell_idx = list(group)[best_idx]
            
            # Remove other cells in this group
            for idx in group:
                if idx != best_cell_idx:
                    cells_to_keep.discard(idx)
        
        # Create the final list of cells
        resolved_cells = [self.cells[i] for i in sorted(cells_to_keep)]
        
        logger.info(f"Conflict resolution complete. Kept {len(resolved_cells)} cells "
                   f"out of {len(self.cells)} original cells")
        
        # Generate visualization if requested
        if self.generate_visuals and self.output_dir:
            self._visualize_conflicts()
            self._generate_conflict_report()
        
        return resolved_cells
    
    def _visualize_conflicts(self) -> None:
        """Generate visualizations of conflict groups and resolution decisions."""
        if not self.output_dir:
            return
            
        os.makedirs(os.path.join(self.output_dir, 'conflict_visuals'), exist_ok=True)
        
        # Create a scatter plot of all cells, highlighting conflicts
        plt.figure(figsize=(10, 10))
        
        # Plot all cells as circles
        for i, cell in enumerate(self.cells):
            circle = plt.Circle((cell.x, cell.y), cell.radius, fill=False, 
                              color='blue', alpha=0.3)
            plt.gca().add_patch(circle)
            
        # Plot cells in conflict with different color
        conflict_indices = set()
        for group in self.conflict_groups:
            conflict_indices.update(group)
            
        for i in conflict_indices:
            cell = self.cells[i]
            circle = plt.Circle((cell.x, cell.y), cell.radius, fill=False, 
                              color='red', alpha=0.5, linewidth=2)
            plt.gca().add_patch(circle)
        
        plt.axis('equal')
        plt.title('Cell Identifications with Conflicts Highlighted')
        plt.savefig(os.path.join(self.output_dir, 'conflict_visuals', 'all_conflicts.png'))
        plt.close()
        
        # Plot individual conflict groups (up to 20)
        for idx, group in enumerate(self.conflict_groups[:20]):
            if len(group) <= 1:
                continue
                
            plt.figure(figsize=(8, 8))
            
            # Plot cells in this conflict group
            for i in group:
                cell = self.cells[i]
                circle = plt.Circle((cell.x, cell.y), cell.radius, fill=False, 
                                 color='red', alpha=0.5, linewidth=2)
                plt.gca().add_patch(circle)
                plt.text(cell.x, cell.y, f"{cell.cell_type}\n{cell.rule_type}", 
                       ha='center', va='center', fontsize=8)
            
            # Expand limits to show full cells
            max_radius = max(self.cells[i].radius for i in group)
            all_x = [self.cells[i].x for i in group]
            all_y = [self.cells[i].y for i in group]
            plt.xlim(min(all_x) - max_radius*1.2, max(all_x) + max_radius*1.2)
            plt.ylim(min(all_y) - max_radius*1.2, max(all_y) + max_radius*1.2)
            
            plt.axis('equal')
            plt.title(f'Conflict Group {idx+1}: {len(group)} Cells')
            plt.savefig(os.path.join(self.output_dir, 'conflict_visuals', f'conflict_group_{idx+1}.png'))
            plt.close()
    
    def _generate_conflict_report(self) -> None:
        """Generate a detailed report on conflicts and resolution decisions."""
        if not self.output_dir:
            return
            
        report_file = os.path.join(self.output_dir, 'conflict_resolution_report.md')
        
        with open(report_file, 'w') as f:
            f.write("# Conflict Resolution Report\n\n")
            
            # Summary statistics
            f.write("## Summary Statistics\n\n")
            f.write(f"- Total cells: {self.conflict_stats['total_cells']}\n")
            f.write(f"- Cells involved in conflicts: {self.conflict_stats['cells_in_conflict']}\n")
            f.write(f"- Number of conflict groups: {self.conflict_stats['conflict_groups']}\n\n")
            
            # Conflicts by cell type
            f.write("## Conflicts by Cell Type Combination\n\n")
            f.write("| Cell Types | Number of Conflicts |\n")
            f.write("|------------|--------------------|\n")
            for types, count in sorted(self.conflict_stats['by_type'].items(), key=lambda x: -x[1]):
                f.write(f"| {', '.join(types)} | {count} |\n")
            f.write("\n")
            
            # Details of major conflict groups
            f.write("## Major Conflict Groups\n\n")
            for idx, group in enumerate(sorted(self.conflict_groups, key=len, reverse=True)[:10]):
                if len(group) <= 1:
                    continue
                    
                f.write(f"### Group {idx+1}: {len(group)} Cells\n\n")
                f.write("| Cell ID | Cell Type | Rule Type | Confidence | Resolution Score |\n")
                f.write("|---------|-----------|-----------|------------|------------------|\n")
                
                # Get cells and scores
                group_cells = [self.cells[i] for i in group]
                scores = [cell.get_resolution_score(self.cell_type_hierarchy, self.gene_specificity) 
                         for cell in group_cells]
                
                # Find winner
                best_idx = np.argmax(scores)
                
                for i, (cell, score) in enumerate(zip(group_cells, scores)):
                    marker = "✓" if i == best_idx else " "
                    f.write(f"| {cell.cell_id} | {cell.cell_type} | {cell.rule_type} | {cell.confidence:.2f} | {score:.2f} {marker} |\n")
                
                f.write("\n")
            
            # Resolution strategy description
            f.write("## Resolution Strategy\n\n")
            f.write("Conflicts are resolved using a comprehensive scoring system that considers:\n\n")
            f.write("1. Base confidence score from the rule\n")
            f.write("2. Rule type (meta > clique > pair)\n")
            f.write("3. Cell type specificity in the hierarchy\n")
            f.write("4. Gene marker specificity\n")
            f.write("5. Presence/absence of negative markers\n")
            f.write("6. Number of genes involved in the identification\n\n")
            
            f.write("The cell with the highest score in each conflict group is kept, and others are removed.\n")
        
        logger.info(f"Conflict resolution report generated: {report_file}")
    
    def export_resolved_cells(self, resolved_cells: List[Cell], output_file: str) -> None:
        """
        Export the resolved cells to a CSV file.
        
        Args:
            resolved_cells: List of cells after conflict resolution
            output_file: Path to save the CSV file
        """
        # Convert to DataFrame
        data = []
        for cell in resolved_cells:
            data.append({
                'cell_id': cell.cell_id,
                'cell_type': cell.cell_type,
                'x': cell.x,
                'y': cell.y,
                'radius': cell.radius,
                'confidence': cell.confidence,
                'rule_id': cell.rule_id,
                'rule_type': cell.rule_type,
                'genes': ','.join(cell.genes),
                'negative_markers_missing': cell.negative_markers_missing,
                'derived_from_meta_rule': cell.derived_from_meta_rule
            })
        
        df = pd.DataFrame(data)
        df.to_csv(output_file, index=False)
        logger.info(f"Exported {len(resolved_cells)} resolved cells to {output_file}")


def resolve_cell_conflicts(cells_df: pd.DataFrame, 
                          output_dir: str,
                          overlap_threshold: float = 0.3,
                          same_type_overlap_threshold: float = 0.7,
                          cell_type_expression_data: Optional[Dict] = None,
                          generate_visuals: bool = True) -> pd.DataFrame:
    """
    Main function to resolve conflicts between identified cells.
    
    Args:
        cells_df: DataFrame containing identified cells
        output_dir: Directory to save outputs
        overlap_threshold: Minimum overlap percentage to consider a conflict
        same_type_overlap_threshold: Higher threshold for cells of same type
        cell_type_expression_data: Optional dict with gene expression data by cell type
        generate_visuals: Whether to generate visualizations
        
    Returns:
        DataFrame with resolved cells
    """
    # Initialize conflict resolution manager
    manager = ConflictResolutionManager(
        overlap_threshold=overlap_threshold,
        same_type_overlap_threshold=same_type_overlap_threshold,
        generate_visuals=generate_visuals,
        output_dir=output_dir
    )
    
    # Calculate gene specificity if expression data provided
    if cell_type_expression_data:
        manager.calculate_gene_specificity(cell_type_expression_data)
    
    # Load cells
    manager.load_cells(cells_df)
    
    # Find conflicts
    manager.find_conflicts()
    
    # Resolve conflicts
    resolved_cells = manager.resolve_conflicts()
    
    # Export resolved cells
    resolved_file = os.path.join(output_dir, 'resolved_cells.csv')
    manager.export_resolved_cells(resolved_cells, resolved_file)
    
    # Convert back to DataFrame
    data = []
    for cell in resolved_cells:
        data.append({
            'cell_id': cell.cell_id,
            'cell_type': cell.cell_type,
            'x': cell.x,
            'y': cell.y,
            'radius': cell.radius,
            'confidence': cell.confidence,
            'rule_id': cell.rule_id,
            'rule_type': cell.rule_type,
            'genes': ','.join(cell.genes)
        })
    
    return pd.DataFrame(data)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Resolve conflicts between identified cells")
    parser.add_argument("--input", required=True, help="CSV file with identified cells")
    parser.add_argument("--output-dir", required=True, help="Directory to save outputs")
    parser.add_argument("--overlap-threshold", type=float, default=0.3, 
                       help="Minimum overlap to consider a conflict")
    parser.add_argument("--same-type-threshold", type=float, default=0.7,
                       help="Overlap threshold for same cell type")
    parser.add_argument("--expression-data", help="Optional JSON file with gene expression data")
    parser.add_argument("--no-visuals", action="store_true", help="Disable visualization generation")
    
    args = parser.parse_args()
    
    # Load input cells
    cells_df = pd.read_csv(args.input)
    
    # Load expression data if provided
    cell_type_expression_data = None
    if args.expression_data:
        import json
        with open(args.expression_data, 'r') as f:
            cell_type_expression_data = json.load(f)
    
    # Resolve conflicts
    resolved_df = resolve_cell_conflicts(
        cells_df=cells_df,
        output_dir=args.output_dir,
        overlap_threshold=args.overlap_threshold,
        same_type_overlap_threshold=args.same_type_threshold,
        cell_type_expression_data=cell_type_expression_data,
        generate_visuals=not args.no_visuals
    )
    
    # Output success message
    n_original = len(cells_df)
    n_resolved = len(resolved_df)
    logger.info(f"Conflict resolution complete: {n_resolved}/{n_original} cells retained "
               f"({n_original - n_resolved} conflicts resolved)") 