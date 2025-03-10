# Enhanced Conflict Resolution System

This document provides detailed information about the enhanced conflict resolution system used in the Astrocytes Ontology Research Project. This system resolves conflicts between overlapping cell identifications using a sophisticated multi-factor approach.

## Overview

When applying co-expression rules to spatial data, we often identify multiple overlapping cells, which may represent:

1. The same cell identified by different rules
2. Different cells with overlapping boundaries
3. False positives that conflict with more reliable identifications

The enhanced conflict resolution system resolves these conflicts by analyzing various factors beyond simple confidence scores, producing more biologically plausible results.

## Key Features

Our conflict resolution system provides several sophisticated features:

1. **Multi-factor Scoring**: Considers multiple evidence types when resolving conflicts
2. **Precise Overlap Calculation**: Uses exact geometric formulas for circle-circle intersection
3. **Cell Type Awareness**: Applies different thresholds for same vs. different cell types
4. **Rule Type Prioritization**: Gives higher priority to more reliable rule types (e.g., meta > clique > pair)
5. **Gene Specificity Integration**: Incorporates gene marker specificity in resolution decisions
6. **Visualization Tools**: Generates visual representations of conflicts and resolutions
7. **Detailed Reporting**: Produces comprehensive reports explaining resolution decisions
8. **Network-based Conflict Detection**: Uses graph theory to find groups of conflicting cells

## Architecture

The conflict resolution system consists of two main components:

1. **Cell Class**: Represents identified cells with methods for:
   - Calculating distance between cells
   - Determining percentage of overlap
   - Computing comprehensive resolution scores

2. **ConflictResolutionManager Class**: Orchestrates the resolution process:
   - Loads and manages cell data
   - Builds a conflict graph
   - Identifies conflict groups
   - Resolves conflicts using a sophisticated algorithm
   - Generates visualizations and reports

## Conflict Detection

The system detects conflicts using a graph-based approach:

1. Each identified cell is represented as a node in a graph
2. An edge is added between cells if their overlap exceeds a threshold
3. The threshold varies based on whether the cells are of the same type
4. Connected components in this graph represent groups of conflicting cells

The overlap between cells is calculated using precise geometric formulas for circle-circle intersection:

```python
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
```

## Multi-factor Scoring System

The conflict resolution system uses a comprehensive scoring approach that considers:

1. **Base Confidence Score**: The original confidence from the rule
2. **Rule Type Bonus**: Meta-rules get +0.15, clique rules get +0.1, pair rules get +0.0
3. **Cell Type Hierarchy**: More specific cell types get a bonus based on their position in the hierarchy
4. **Gene Specificity**: Cells identified by more specific marker genes receive a bonus
5. **Negative Marker Validation**: Cells with verified absence of negative markers get +0.05
6. **Meta-Rule Derivation**: Cells from meta-rules get an additional +0.05
7. **Gene Count**: More genes involved in identification increases the score (up to +0.1)

This scoring function is implemented as:

```python
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
```

## Resolution Process

The conflict resolution process follows these steps:

1. **Group Identification**: Find all groups of conflicting cells using a graph approach
2. **Score Calculation**: Calculate resolution scores for each cell based on multiple factors
3. **Winner Selection**: Within each conflict group, keep the cell with the highest score
4. **Result Assembly**: Create a final set of non-conflicting cells

## Configuration Options

The system offers several configurable parameters:

| Parameter | Description | Default | When to Adjust |
|---|---|---|---|
| `overlap_threshold` | Minimum overlap to consider a conflict | 0.3 | Higher for more permissive detection, lower for stricter |
| `same_type_overlap_threshold` | Overlap threshold for same cell type | 0.7 | Higher for allowing more same-type cells |
| `generate_visuals` | Whether to generate visualizations | True | Disable for performance in large datasets |
| `cell_type_hierarchy` | Dictionary mapping cell types to hierarchy levels | See code | Adjust for specific research focus |
| `cell_type_expression_data` | Expression data for calculating gene specificity | Optional | Provide for better resolution |

## Gene Specificity Calculation

The system can calculate gene specificity scores to boost confidence in cells identified by more cell-type-specific genes:

```python
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
```

## Visualization and Reporting

The system generates several visualizations to help understand conflict resolution:

1. **Overview visualization**: Shows all cells with conflicts highlighted
2. **Conflict group visualizations**: Detailed views of each conflict group
3. **Detailed report**: Markdown report with statistics and resolution decisions

![Example Conflict Visualization](../data/processed/results/conflict_visuals/all_conflicts.png)
*Example visualization showing cells with conflicts highlighted in red*

## Example: Conflict Resolution

Here's an example of how the system resolves a conflict group with three overlapping cells:

1. **Cell A**:
   - Identified by a clique rule for Astrocytes
   - Confidence: 0.9
   - Genes involved: GPM6A, NRXN1, ADGRB3
   - Resolution score: 1.15

2. **Cell B**:
   - Identified by a pair rule for Astrocytes
   - Confidence: 0.8
   - Genes involved: LSAMP, CNTN1
   - Resolution score: 0.94

3. **Cell C**:
   - Identified by a meta-rule for Astrocytes1
   - Confidence: 0.95
   - Genes involved: LRP1B, NRG3, GPM6A, SLC1A2, CTNNA2
   - Resolution score: 1.28

In this case, Cell C would be kept because it has the highest resolution score, driven by:
- Higher base confidence (0.95)
- Meta-rule bonus (+0.15)
- Cell type hierarchy bonus (Astrocytes1 is more specific, +0.1)
- More genes involved (+0.08)

## Integration with Pipeline

The conflict resolution system is integrated into the pipeline in two ways:

1. **As a standalone module** (`conflict_resolution.py`):
   - Can be run directly on cell identification outputs
   - Offers full configuration and reporting capabilities

2. **Integrated with cell identification** (`identify_cells.py`):
   - Automatically applied during cell identification
   - Controlled via command-line flags

Example usage in the pipeline:

```bash
python scripts/identify_cells.py \
    --spatial-data data/processed/spatial_data.ttl \
    --rules-dir data/processed/rules \
    --meta-rules-file data/processed/meta_rules/meta_rules.rq \
    --output-dir data/processed/results \
    --min-confidence 0.7 \
    --resolve-conflicts \
    --enhanced-resolution \
    --overlap-threshold 0.3
```

## Fallback Mechanism

The system includes a basic conflict resolution fallback when the enhanced module is not available:

```python
def basic_conflict_resolution(cells_df: pd.DataFrame, 
                           overlap_threshold: float = 0.3,
                           prioritize_meta_rules: bool = True) -> pd.DataFrame:
    """
    Basic conflict resolution strategy for overlapping cells.
    
    This is a simplified version used as fallback when the enhanced module is not available.
    """
    # ... (simplified algorithm that prioritizes by rule type and confidence) ...
```

This ensures the pipeline continues to function even if the enhanced resolution module is not available.

## Performance Considerations

For large datasets with many identified cells, the conflict resolution process can be computationally intensive. Consider the following:

1. **Conflict Detection**: O(n²) in the worst case for n cells
2. **Resolution Scoring**: Linear time for each cell in conflict
3. **Visualization**: Can be slow for large numbers of conflict groups

To optimize performance:
- Set `generate_visuals=False` for very large datasets
- Process cell types separately if dealing with many cells
- Consider sampling or pre-filtering obvious conflicts

## Future Enhancements

Planned improvements to the conflict resolution system include:

1. **3D support**: Extending to three-dimensional spatial data
2. **Machine learning integration**: Using ML models to predict resolution decisions
3. **Interactive resolution**: Tools for manual resolution of complex conflicts
4. **Parallel processing**: Optimizations for large-scale datasets
5. **Integration with density-based clustering**: Using cluster information to inform resolution

## References

1. Graph theory for conflict detection:
   - Tarjan, R. E. (1972). "Depth-first search and linear graph algorithms." SIAM Journal on Computing, 1(2), 146-160.

2. Geometric algorithms for overlap calculation:
   - O'Rourke, J. (1998). "Computational Geometry in C." Cambridge University Press.

3. Multi-factor decision making:
   - Saaty, T. L. (1990). "How to make a decision: The analytic hierarchy process." European Journal of Operational Research, 48(1), 9-26. 