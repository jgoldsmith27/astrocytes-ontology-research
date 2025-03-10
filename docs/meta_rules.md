# Meta-Rules: Higher-Order Co-expression Patterns

This document provides detailed documentation on the implementation and usage of meta-rules in the Astrocytes Ontology Research Project.

## Overview

Meta-rules are higher-order co-expression patterns that identify relationships between existing co-expression rules. They provide a more powerful approach to cell identification by requiring multiple co-expression patterns to co-occur, increasing specificity and confidence.

## Key Concepts

### Rule Association

A rule association is a statistically significant relationship between two co-expression rules, indicating they tend to match the same cells in single-cell data. Associations are measured using:

- **Lift**: The ratio of observed co-occurrence to expected co-occurrence under independence
- **Statistical significance**: P-value from Fisher's exact test
- **Conditional probabilities**: P(Rule1|Rule2) and P(Rule2|Rule1)

### Meta-Rule

A meta-rule is a SPARQL CONSTRUCT query that identifies cells in spatial data where two associated co-expression patterns occur in close proximity. Meta-rules have:

- Higher specificity than individual rules
- Increased confidence scores
- Priority during conflict resolution

## Implementation Details

### Workflow

1. **Load Rules**: Extract rule information from all SPARQL rule files
2. **Apply Rules to Single-Cell Data**: Determine which cells match each rule
3. **Calculate Associations**: Identify statistically significant rule relationships
4. **Generate Meta-Rules**: Create SPARQL queries for significant associations
5. **Apply Meta-Rules**: Execute meta-rules alongside regular rules with higher priority

### Applying Rules to Single-Cell Data

A key aspect of our implementation is how we apply co-expression rules to single-cell data. This is done in the `apply_rules_to_single_cell` function in `generate_meta_rules.py`:

```python
def apply_rules_to_single_cell(rules, single_cell_data):
    """
    Apply rules to single-cell data to determine which cells match each rule.
    This allows us to analyze rule co-occurrence patterns.
    
    Returns:
        dict: Dictionary mapping rule IDs to sets of matching cell indices
    """
    import scanpy as sc
    
    # Load single-cell data in h5ad format
    adata = sc.read_h5ad(single_cell_data)
    
    # Dictionary to store cells matching each rule
    rule_matches = {}
    
    # Group rules by cell type
    rules_by_cell_type = defaultdict(list)
    for rule in rules:
        rules_by_cell_type[rule['cell_type']].append(rule)
    
    # For each cell type, apply its rules
    for cell_type, cell_type_rules in rules_by_cell_type.items():
        # Get indices of cells of this type
        cell_indices = np.where(adata.obs['cell_type'] == cell_type)[0]
        if len(cell_indices) == 0:
            continue
        
        # Get expression data for these cells
        expr_data = adata.X[cell_indices]
        
        # For each rule, find matching cells
        for rule in cell_type_rules:
            rule_id = rule['rule_id']
            genes = rule['genes']
            
            # Get indices of genes in the expression matrix
            gene_indices = [np.where(adata.var_names == gene)[0][0] 
                          for gene in genes if gene in adata.var_names]
            
            # Find cells where all genes in the rule are expressed (above threshold)
            thresh = 0.2  # Same threshold used in rule generation
            matching_cells = set()
            
            for i, cell_idx in enumerate(cell_indices):
                if all(expr_data[i, gene_idx] > thresh for gene_idx in gene_indices):
                    matching_cells.add(cell_idx)
            
            rule_matches[rule_id] = matching_cells
    
    return rule_matches
```

Important notes about this process:

1. We **do not use ontology format** for single-cell data analysis. Instead, we:
   - Load the h5ad file directly with scanpy
   - Extract the expression matrix for each cell type
   - Apply the rule logic (gene co-expression) directly to the expression values

2. Rule application logic:
   - For each rule, we identify the genes involved
   - We find cells that express all these genes above the threshold
   - We store the list of matching cells for each rule

3. This approach allows us to:
   - Work with the single-cell data in its native format
   - Efficiently analyze large datasets
   - Generate statistics on rule co-occurrence

### Calculating Rule Associations

After applying rules to single-cell data, we calculate associations between rules:

```python
def calculate_rule_associations(rules, rule_matches, min_lift=3.0, max_p_value=0.01):
    # Get total number of cells
    all_cell_indices = set()
    for cells in rule_matches.values():
        all_cell_indices.update(cells)
    total_cells = len(all_cell_indices)
    
    # List to store significant associations
    associations = []
    
    # Calculate associations for each rule pair
    rule_ids = list(rule_matches.keys())
    for i in range(len(rule_ids)):
        for j in range(i+1, len(rule_ids)):
            rule_id1 = rule_ids[i]
            rule_id2 = rule_ids[j]
            
            # Get cells matching each rule
            cells1 = rule_matches[rule_id1]
            cells2 = rule_matches[rule_id2]
            
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
                    'n_cells1': n1,
                    'n_cells2': n2,
                    'n_co_occurrence': n12,
                    'lift': lift,
                    'p_value': p_value,
                    'p1_given_2': n12 / n2,
                    'p2_given_1': n12 / n1
                }
                associations.append(association)
    
    return associations
```

Key metrics calculated:

- **Co-occurrence count**: Number of cells matching both rules
- **Lift**: Ratio of observed co-occurrence to expected co-occurrence
- **P-value**: Statistical significance of the association
- **Conditional probabilities**: How likely one rule is given the other

### Generating Meta-Rules

For each significant association, we generate a meta-rule that requires both patterns to be present in the spatial data:

```python
def generate_meta_rule(association, spatial_distance=75, confidence_boost=0.1):
    rule1 = association['rule1']
    rule2 = association['rule2']
    cell_type = rule1['cell_type']
    lift = association['lift']
    
    # Calculate confidence score
    base_confidence = max(rule1['confidence'], rule2['confidence'])
    meta_confidence = min(0.99, base_confidence + confidence_boost)
    
    # Generate SPARQL CONSTRUCT query
    meta_rule = """
    # Meta-rule combining patterns: {rule_id1} and {rule_id2}
    # Cell Type: {cell_type}, Association Lift: {lift:.2f}
    
    CONSTRUCT {{
      ?cell a astro:SpatialCell ;
            astro:hasCellType "{cell_type}" ;
            astro:hasXCoordinate ?centerX ;
            astro:hasYCoordinate ?centerY ;
            astro:hasRadius ?combinedRadius ;
            astro:hasConfidence "{confidence}"^^xsd:decimal ;
            astro:derivedFromMetaRule "true"^^xsd:boolean .
    }}
    WHERE {{
      # Find Pattern 1
      {{
        {pattern1_where}
      }}
      
      # Find Pattern 2
      {{
        {pattern2_where}
      }}
      
      # Ensure patterns are close enough
      BIND(SQRT(POW(?centerX1 - ?centerX2, 2) + POW(?centerY1 - ?centerY2, 2)) AS ?patternDist)
      FILTER(?patternDist < {spatial_distance})
      
      # Calculate combined cell
      BIND((?centerX1 + ?centerX2) / 2 AS ?centerX)
      BIND((?centerY1 + ?centerY2) / 2 AS ?centerY)
      BIND(?patternDist / 2 + MAX(?radius1, ?radius2) AS ?combinedRadius)
    }}
    """
    
    # Format with appropriate pattern matching logic
    # ...
    
    return meta_rule
```

Meta-rule structure:

1. **Identification**: Unique identifier and annotation about source rules
2. **Pattern Recognition**: Separate WHERE clauses to find each pattern
3. **Proximity Check**: Ensure patterns are within spatial_distance
4. **Combined Cell Creation**: Calculate center and radius for the new cell

## Key Parameters

| Parameter | Description | Default | When to Adjust |
|-----------|-------------|---------|----------------|
| min_lift | Minimum lift value for significant association | 3.0 | Increase for more stringent associations |
| max_p_value | Maximum p-value for statistical significance | 0.01 | Decrease for higher statistical confidence |
| spatial_distance | Maximum distance between patterns in spatial data | 75 | Adjust based on expected cell sizes |
| confidence_boost | Amount to increase confidence for meta-rules | 0.1 | Adjust to prioritize meta-rules more/less |
| same_cell_type_only | Whether to only associate rules from same cell type | True | Set to False to find cross-type patterns |

## Outputs

The meta-rule generation process produces several outputs:

1. **Meta-rules file**: `meta_rules.rq` with all generated meta-rules
2. **Summary CSV**: `meta_rules_summary.csv` with statistics on all significant associations
3. **Rule network**: `rule_network.graphml` visualizing relationships between rules

### Rule Network Visualization

The rule network provides valuable insights into co-expression pattern relationships:

- **Nodes**: Individual co-expression rules (clique or pair)
- **Edges**: Significant associations between rules
- **Node Attributes**: Cell type, rule type, genes involved
- **Edge Attributes**: Lift, p-value, co-occurrence count

This network can be visualized using tools like Cytoscape or Gephi to:
- Identify clusters of related rules
- Discover higher-order modules of gene co-expression
- Find bridging rules that connect different co-expression patterns

## Integration with Pipeline

Meta-rules are integrated into the pipeline in `run_astrocyte_identification.sh`:

```bash
# Step 3: Generate meta-rules by analyzing rule co-occurrence patterns
echo "Step 3: Generating higher-order co-expression patterns (meta-rules)..."
python scripts/generate_meta_rules.py \
    --rules-dir $RULES_DIR \
    --single-cell-data $SCRNA_DATA \
    --output-dir $META_RULES_DIR \
    --min-lift 3.0 \
    --max-p-value 0.01 \
    --spatial-distance 75 \
    --confidence-boost 0.1 \
    --same-cell-type-only

# Step 4: Apply rules to identify cells in spatial data
echo "Step 4: Applying rules to identify cells in spatial data..."
python scripts/identify_cells.py \
    --spatial-data $SPATIAL_OUTPUT \
    --rules-dir $RULES_DIR \
    --meta-rules-file $META_RULES_DIR/meta_rules.rq \
    --output-dir $RESULTS_DIR \
    --min-confidence 0.7 \
    --resolve-conflicts \
    --prioritize-meta-rules
```

## Example Meta-Rule

Here's an example of a generated meta-rule:

```sparql
# Meta-rule combining patterns: Astrocytes_clique_GPM6A_NRXN1_ADGRB3 and Astrocytes_pair_LSAMP_CNTN1
# Cell Type: Astrocytes, Association Lift: 5.83
# Genes: GPM6A, NRXN1, ADGRB3, LSAMP, CNTN1
PREFIX astro: <http://example.org/astrocytes#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

CONSTRUCT {
  ?cell a astro:SpatialCell ;
        astro:hasCellType "Astrocytes" ;
        astro:hasXCoordinate ?centerX ;
        astro:hasYCoordinate ?centerY ;
        astro:hasRadius ?combinedRadius ;
        astro:hasConfidence "0.97"^^xsd:decimal ;
        astro:hasCellID ?cellID ;
        astro:derivedFromMetaRule "true"^^xsd:boolean .
        
  # Link to constituent patterns
  ?cell astro:includesPattern ?pattern1 .
  ?cell astro:includesPattern ?pattern2 .
  
  # Include all points
  ?cell astro:includesPoint ?point1_0 .
  ?cell astro:includesPoint ?point1_1 .
  ?cell astro:includesPoint ?point1_2 .
  ?cell astro:includesPoint ?point2_0 .
  ?cell astro:includesPoint ?point2_1 .
}
WHERE {
  # Find Pattern 1 (Clique)
  {
    # Points for each gene in the clique
    ?point1_0 a astro:SpatialDataPoint ;
         astro:expressesGene ?gene1_0 ;
         astro:hasXCoordinate ?x1_0 ;
         astro:hasYCoordinate ?y1_0 .
    ?gene1_0 astro:hasGeneID "GPM6A" .
    
    # ... points for NRXN1 and ADGRB3 ...
    
    # Distance constraints between all points
    BIND(SQRT(POW(?x1_0 - ?x1_1, 2) + POW(?y1_0 - ?y1_1, 2)) AS ?dist1_0_1)
    FILTER(?dist1_0_1 < 50)
    # ... other distance constraints ...
    
    # Calculate center and radius
    BIND((?x1_0 + ?x1_1 + ?x1_2) / 3 AS ?centerX1)
    BIND((?y1_0 + ?y1_1 + ?y1_2) / 3 AS ?centerY1)
    BIND(MAX(?dist1_0_1, ?dist1_0_2, ?dist1_1_2) / 2 AS ?radius1)
    
    # Create a temporary cell entity
    ?pattern1Cell a astro:SpatialCell .
    BIND(?pattern1Cell AS ?pattern1)
  }
  
  # Find Pattern 2 (Pair)
  {
    # Points for the pair
    ?point2_0 a astro:SpatialDataPoint ;
         astro:expressesGene ?gene2_0 ;
         astro:hasXCoordinate ?x2_0 ;
         astro:hasYCoordinate ?y2_0 .
    ?gene2_0 astro:hasGeneID "LSAMP" .
    
    ?point2_1 a astro:SpatialDataPoint ;
         astro:expressesGene ?gene2_1 ;
         astro:hasXCoordinate ?x2_1 ;
         astro:hasYCoordinate ?y2_1 .
    ?gene2_1 astro:hasGeneID "CNTN1" .
    
    # Distance constraint
    BIND(SQRT(POW(?x2_0 - ?x2_1, 2) + POW(?y2_0 - ?y2_1, 2)) AS ?dist2)
    FILTER(?dist2 < 50)
    
    # Calculate center and radius
    BIND((?x2_0 + ?x2_1) / 2 AS ?centerX2)
    BIND((?y2_0 + ?y2_1) / 2 AS ?centerY2)
    BIND(?dist2 / 2 AS ?radius2)
    
    # Create a temporary cell entity
    ?pattern2Cell a astro:SpatialCell .
    BIND(?pattern2Cell AS ?pattern2)
  }
  
  # Ensure patterns are close enough
  BIND(SQRT(POW(?centerX1 - ?centerX2, 2) + POW(?centerY1 - ?centerY2, 2)) AS ?patternDist)
  FILTER(?patternDist < 75)
  
  # Calculate combined cell center and radius
  BIND((?centerX1 + ?centerX2) / 2 AS ?centerX)
  BIND((?centerY1 + ?centerY2) / 2 AS ?centerY)
  BIND(?patternDist / 2 + MAX(?radius1, ?radius2) AS ?combinedRadius)
  
  # Generate a unique cell ID
  BIND(CONCAT("meta_cell_", STR(?centerX), "_", STR(?centerY), "_", "Astrocytes") AS ?cellID)
}
```

## Benefits and Limitations

### Benefits

1. **Increased Specificity**: Meta-rules require multiple co-expression patterns, reducing false positives
2. **Biological Insights**: Co-occurring patterns may represent functional modules or biological processes
3. **Cell Type Disambiguation**: Better distinguish between similar cell types that share some markers
4. **Confidence Ranking**: Provides a way to prioritize identifications based on multiple lines of evidence
5. **Discovering Subtypes**: May help identify functional subtypes not explicitly labeled in original annotations

### Limitations

1. **Complexity**: Meta-rules introduce additional computational complexity
2. **Parameter Sensitivity**: Results depend on multiple parameters that may require tuning
3. **Data Requirements**: Requires sufficient single-cell data to detect significant associations
4. **Interpretation Challenges**: Higher-order patterns may be harder to interpret biologically

## Troubleshooting

### Common Issues

1. **No significant associations found**
   - Decrease min_lift or increase max_p_value
   - Check that rules match sufficient cells in the single-cell data
   - Verify the single-cell data format and annotations

2. **Too many associations found**
   - Increase min_lift or decrease max_p_value
   - Consider using same_cell_type_only=True to limit associations

3. **Meta-rules too restrictive**
   - Increase spatial_distance between patterns
   - Decrease min_confidence when applying rules

4. **Memory issues with large datasets**
   - Process cell types in batches
   - Use more efficient data structures for large matrices
   - Consider sampling the single-cell data

## References

1. Agrawal, R., Imielinski, T., & Swami, A. (1993). Mining association rules between sets of items in large databases. ACM SIGMOD Record, 22(2), 207-216.
2. Brin, S., Motwani, R., Ullman, J. D., & Tsur, S. (1997). Dynamic itemset counting and implication rules for market basket data. ACM SIGMOD Record, 26(2), 255-264.
3. Pérez-Suárez, A., Martínez-Trinidad, J. F., Carrasco-Ochoa, J. A., & Medina-Pagola, J. E. (2013). A review of clique-based clustering algorithms. New Frontiers in Applied Data Mining, 30-45.
4. Milo, R., Shen-Orr, S., Itzkovitz, S., Kashtan, N., Chklovskii, D., & Alon, U. (2002). Network motifs: simple building blocks of complex networks. Science, 298(5594), 824-827. 