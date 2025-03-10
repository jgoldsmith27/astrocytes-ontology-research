# Co-expression Rules Guide

This document provides a comprehensive overview of the co-expression rules used in the Astrocytes Ontology Research Project.

## Overview

Co-expression rules are SPARQL CONSTRUCT queries that identify cells in spatial transcriptomics data based on gene co-expression patterns observed in single-cell RNA sequencing data. These rules form the foundation of our cell identification approach.

## Rule Types

### 1. Pair Rules

Pair rules identify potential cells based on the co-expression of two genes that frequently appear together in a specific cell type.

**Structure:**
- Identify spatial points where both genes are expressed
- Ensure the points are within a specified distance threshold
- Create a cell centered at the midpoint between the two points
- Set the cell type based on the source of the co-expression pattern
- Assign a confidence score equal to the co-expression strength

**Example:**
```sparql
# Pair rule for Astrocytes with genes: GPM6A and NRXN1 (coexpression = 0.82)
CONSTRUCT {
  ?cell a astro:SpatialCell ;
        astro:hasCellType "Astrocytes" ;
        astro:hasXCoordinate ?centerX ;
        astro:hasYCoordinate ?centerY ;
        astro:hasRadius ?radius ;
        astro:hasConfidence "0.82"^^xsd:decimal ;
        ...
}
WHERE {
  # Find points expressing each gene and calculate their distance
  ...
}
```

### 2. Clique Rules

Clique rules identify potential cells based on larger groups of co-expressed genes (3 or more) that form a fully connected subgraph in the co-expression network.

**Structure:**
- Identify spatial points where each gene in the clique is expressed
- Ensure all points are within the distance threshold of each other
- Create a cell at the centroid of all the points
- Set the cell type based on the source of the clique
- Assign a high confidence score (0.9) due to the multiple gene evidence

**Example:**
```sparql
# Clique rule for Astrocytes1 with 3 genes: LRP1B, NRG3, GPM6A
CONSTRUCT {
  ?cell a astro:SpatialCell ;
        astro:hasCellType "Astrocytes1" ;
        astro:hasXCoordinate ?centerX ;
        astro:hasYCoordinate ?centerY ;
        astro:hasRadius ?radius ;
        astro:hasConfidence "0.9"^^xsd:decimal ;
        ...
}
WHERE {
  # Find points expressing each gene and ensure proximity
  ...
}
```

### 3. Negative Marker Rules

Negative marker rules enhance the specificity of pair and clique rules by ensuring that negative markers (genes that should NOT be expressed in the cell type) are absent from the identified region.

**Structure:**
- Extend pair or clique rules with FILTER NOT EXISTS clauses
- Check for absence of negative markers within the cell region
- Maintain or boost confidence score when negative markers are confirmed absent

**Example:**
```sparql
# Within a rule's WHERE clause:
FILTER NOT EXISTS {
  ?negPoint a astro:SpatialDataPoint ;
            astro:expressesGene ?negGene ;
            astro:hasXCoordinate ?negX ;
            astro:hasYCoordinate ?negY .
  ?negGene astro:hasGeneID "NEGATIVE_MARKER_GENE" .
  BIND(SQRT(POW(?centerX - ?negX, 2) + POW(?centerY - ?negY, 2)) AS ?negDist)
  FILTER(?negDist < ?radius * 1.5)
}
```

### 4. Meta-Rules (Higher-Order Co-expression Patterns)

Meta-rules represent a significant enhancement to our methodology by identifying relationships between existing co-expression patterns. These rules identify cells with multiple co-expression patterns in proximity, providing higher specificity.

**Structure:**
- Identify regions where two distinct co-expression patterns occur in proximity
- Create a cell that encompasses both patterns
- Assign a higher confidence score than either individual pattern
- Include a flag indicating the cell was identified by a meta-rule

**How Meta-Rules Are Generated:**

1. **Rule Co-occurrence Analysis**: We analyze the single-cell data to identify which rule patterns frequently co-occur in the same cells
   
2. **Statistical Association Measurement**: For each pair of rules, we calculate:
   - **Lift**: The ratio of observed co-occurrence to expected co-occurrence under independence
   - **Fisher's exact test p-value**: Statistical significance of the association
   - **Conditional probabilities**: P(Rule1|Rule2) and P(Rule2|Rule1)

3. **Significant Association Selection**: Only rule pairs with high lift (default: ≥3.0) and low p-value (default: ≤0.01) are used to create meta-rules

4. **Spatial Integration**: Meta-rules require both patterns to be within a specified spatial distance (default: 75 units)

**Example Meta-Rule:**
```sparql
# Meta-rule combining patterns: Astrocytes_clique_GPM6A_NRXN1_ADGRB3 and Astrocytes_pair_LSAMP_CNTN1
# Cell Type: Astrocytes, Association Lift: 5.83
# Genes: GPM6A, NRXN1, ADGRB3, LSAMP, CNTN1
CONSTRUCT {
  ?cell a astro:SpatialCell ;
        astro:hasCellType "Astrocytes" ;
        astro:hasXCoordinate ?centerX ;
        astro:hasYCoordinate ?centerY ;
        astro:hasRadius ?combinedRadius ;
        astro:hasConfidence "0.98"^^xsd:decimal ;
        astro:derivedFromMetaRule "true"^^xsd:boolean ;
        ...
}
WHERE {
  # Find pattern 1 (clique)
  {
    # Clique pattern matching logic
    ...
  }
  
  # Find pattern 2 (pair)
  {
    # Pair pattern matching logic
    ...
  }
  
  # Ensure patterns are spatially related
  BIND(SQRT(POW(?centerX1 - ?centerX2, 2) + POW(?centerY1 - ?centerY2, 2)) AS ?patternDist)
  FILTER(?patternDist < 75)
  
  # Calculate combined cell properties
  ...
}
```

**Benefits of Meta-Rules:**

1. **Increased Specificity**: By requiring multiple co-expression patterns, meta-rules reduce false positives
2. **Biological Relevance**: Co-occurring patterns may represent functional modules or biological processes
3. **Improved Cell Type Disambiguation**: Better distinguish between similar cell types that share some markers
4. **Confidence Ranking**: Cells identified by meta-rules receive higher confidence scores and priority during conflict resolution
5. **Subtypes Discovery**: May help identify functional subtypes not explicitly labeled in original annotations

## Rule Generation Process

1. **Extract Co-expression Patterns**: Analyze single-cell data to identify co-expressed genes within each cell type
2. **Build Co-expression Network**: Create a network where nodes are genes and edges represent strong co-expression
3. **Find Cliques**: Identify fully connected subgraphs in the network
4. **Generate SPARQL Rules**: Convert the patterns into executable SPARQL CONSTRUCT queries
5. **Analyze Rule Co-occurrence**: For meta-rules, analyze which rules tend to match the same cells
6. **Generate Meta-Rules**: Create higher-order rules from statistically significant rule associations

## Rule Application

Rules are applied to the spatial data in RDF format using a SPARQL engine. The application process:

1. Execute each rule against the spatial RDF graph
2. Collect all cells identified by all rules
3. Resolve conflicts when multiple rules identify overlapping cells
4. Filter low-confidence cells and apply proximity validation
5. Export the final set of identified cells

## Customizing Rules

Rules can be customized by adjusting several parameters:

| Parameter | Description | Default | Impact |
|-----------|-------------|---------|--------|
| min_expression | Threshold to consider a gene expressed | 0.2 | Higher = more stringent expression detection |
| min_coexpression | Minimum conditional probability | 0.6 | Higher = more specific co-expression |
| distance_threshold | Maximum spatial distance | 50 | Higher = larger potential cells |
| min_confidence | Minimum confidence to keep a cell | 0.7 | Higher = fewer but more confident cells |
| neg_marker_threshold | Threshold for negative markers | 0.05 | Lower = stricter negative marker filtering |
| min_lift | Minimum lift for meta-rules | 3.0 | Higher = more stringent meta-rule creation |
| max_p_value | Maximum p-value for meta-rules | 0.01 | Lower = more statistically significant meta-rules |
| spatial_distance | Distance between patterns for meta-rules | 75 | Controls how close patterns must be |

## Visualizing Rule Relationships

The meta-rule generation process creates a network visualization (GraphML format) showing relationships between rules:

- **Nodes**: Individual co-expression rules
- **Edges**: Significant associations between rules
- **Edge Weight**: Lift value (strength of association)
- **Node Attributes**: Cell type, rule type, genes involved

This network can be visualized using tools like Cytoscape or Gephi to gain insights into the higher-order structure of co-expression patterns.

## Troubleshooting

### Common Rule Generation Issues

- **Too Few Rules**: Decrease min_coexpression or increase top_genes
- **Too Many Rules**: Increase min_coexpression or decrease top_genes
- **No Meta-Rules**: Decrease min_lift or increase max_p_value
- **Meta-Rules Too Restrictive**: Increase spatial_distance between patterns

### Common Rule Application Issues

- **Too Few Cells**: Lower min_confidence or increase distance_threshold
- **Too Many Cells**: Increase min_confidence or decrease distance_threshold
- **Overlapping Cells**: Refine conflict resolution strategy
- **Unbalanced Cell Types**: Adjust cell type-specific parameters

## References

- The SPARQL 1.1 Query Language: https://www.w3.org/TR/sparql11-query/
- NetworkX Documentation (for clique detection): https://networkx.org/documentation/stable/reference/algorithms/clique.html
- Lift (data mining): https://en.wikipedia.org/wiki/Lift_(data_mining) 