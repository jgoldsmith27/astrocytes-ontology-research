# Methodology for Astrocyte Identification

## Scientific Background

### Single-cell RNA-seq vs. Spatial Transcriptomics

Single-cell RNA sequencing (scRNA-seq) provides high-resolution gene expression data from individual cells but loses spatial context. Spatial transcriptomics preserves spatial information but typically has lower resolution and sensitivity. This project bridges these technologies by transferring cell type information from scRNA-seq to spatial transcriptomics using gene co-expression patterns.

### Astrocyte Cells

Astrocytes are star-shaped glial cells in the central nervous system that perform various functions, including:
- Supporting neurons and synapses
- Regulating blood flow
- Maintaining the blood-brain barrier
- Participating in neurotransmitter recycling

In our dataset, we focus on two types of astrocytes:
1. **Astrocytes (general category)**: 1,287 cells in the single-cell data
2. **Astrocytes1 (subtype)**: 7,116 cells in the single-cell data

## Technical Approach

Our methodology consists of four main steps:

1. Data Conversion
2. Co-expression Analysis
3. SPARQL Rule Generation
4. Rule Application

### 1. Data Conversion

#### Single-cell Data

The single-cell data (h5ad format) is already annotated with cell types. We extract:
- Expression matrices for each cell type
- Gene lists per cell type
- Metadata about cells and genes

#### Spatial Data Conversion

Spatial data is converted to RDF (Resource Description Framework) format following our custom ontology:

```
SpatialData.csv → convert_spatial_to_turtle.py → spatial_data.ttl
```

The conversion process:
1. Reads the CSV file in chunks for memory efficiency
2. Creates RDF entities for:
   - Spatial data points
   - Genes
   - Spatial bins
3. Links these entities with properties like:
   - `expressesGene`
   - `hasXCoordinate`/`hasYCoordinate`
   - `locatedInBin`

### 2. Co-expression Analysis

The co-expression analysis is performed on the single-cell data to identify genes that are reliably co-expressed in specific cell types:

1. **Expression Thresholding**:
   - A gene is considered "expressed" if its normalized expression level exceeds a threshold (default: 0.2)
   - This converts the continuous expression matrix into a binary matrix

2. **Conditional Probability Calculation**:
   - For each gene pair (A, B), we calculate:
     - P(A): Probability that gene A is expressed in the cell population
     - P(B): Probability that gene B is expressed
     - P(B|A): Conditional probability that gene B is expressed given that gene A is expressed
     - P(A|B): Conditional probability that gene A is expressed given that gene B is expressed
     - P(A,B): Joint probability that both genes are expressed

3. **Co-expression Network Construction**:
   - A network is built where:
     - Nodes are genes
     - Edges connect genes with strong co-expression (conditional probability >= threshold, default: 0.6)
     - Edge weights represent the co-expression strength (maximum of the two conditional probabilities)

4. **Clique Detection**:
   - We use NetworkX's clique finding algorithm to identify fully connected subgraphs
   - These represent sets of genes where all genes in the set are co-expressed with each other
   - Cliques are ranked by size (number of genes) and retained if they contain at least 3 genes

5. **Negative Marker Identification** (NEW):
   - For each cell type, we identify negative markers (genes that should NOT be expressed)
   - These are genes that are:
     - Highly expressed in other cell types (top 10% by expression level)
     - Rarely expressed in the target cell type (<5% of cells)
   - Negative markers increase specificity by ensuring identified cells don't express genes characteristic of other cell types

6. **Higher-Order Pattern Detection** (PLANNED):
   - Beyond individual gene co-expression, we analyze how different co-expression patterns and cliques themselves co-occur
   - This creates a "meta-network" where nodes are co-expression rules and edges represent their co-occurrence in cells
   - We can identify significant associations between different co-expression patterns that may represent higher-order biological processes
   - These higher-order patterns provide an additional layer of specificity for cell identification

### 3. SPARQL Rule Generation

Co-expression patterns are converted into SPARQL CONSTRUCT queries:

1. **Clique Rules**:
   - For each clique of co-expressed genes, we generate a rule that:
     - Identifies spatial points where each gene in the clique is expressed
     - Ensures these points are within a proximity threshold (default: 50 units)
     - Creates a cell entity at the centroid of these points
     - Sets the cell type based on the source of the clique (e.g., Astrocyte)
     - Assigns a high confidence score (0.9) due to the multiple gene evidence

2. **Pair Rules**:
   - For strong gene pairs not part of larger cliques, we generate simpler rules that:
     - Identify spatial points where both genes in the pair are expressed
     - Ensure these points are within the proximity threshold
     - Create a cell entity at the midpoint between these points
     - Set the cell type based on the source
     - Assign a confidence score equal to the co-expression strength

3. **Negative Marker Rules** (NEW):
   - Each positive rule (clique or pair) is augmented with negative marker constraints:
   - The rule verifies that identified regions do NOT contain expression of negative markers
   - This is implemented as a FILTER NOT EXISTS clause in the SPARQL query
   - A distance threshold ensures we only exclude negative markers close to the potential cell

4. **Meta-Rules** (PLANNED):
   - For strongly associated co-expression patterns, we generate composite rules that:
     - Require multiple co-expression patterns to be present in proximity to each other
     - Create a cell entity only when both patterns are detected
     - Apply more sophisticated confidence scoring based on the strength of each pattern and their association
     - These rules have the highest specificity but may have lower sensitivity

The mathematical basis for the rule structure is:

- For clique rules with n genes, the distance constraint is:
  ```
  For all pairs (i,j) where i,j ∈ {0,1,...,n-1} and i < j:
    dist(i,j) = sqrt((x_i - x_j)² + (y_i - y_j)²) < threshold
  ```

- Cell center calculation:
  ```
  centerX = mean(x_0, x_1, ..., x_n-1)
  centerY = mean(y_0, y_1, ..., y_n-1)
  ```

- Cell radius calculation:
  ```
  radius = max(dist(0,1), dist(0,2), ..., dist(n-2,n-1)) / 2
  ```

- Negative marker constraint (NEW):
  ```
  FILTER NOT EXISTS {
    ?negPoint a astro:SpatialDataPoint ;
              astro:expressesGene ?negGene ;
              astro:hasXCoordinate ?negX ;
              astro:hasYCoordinate ?negY .
    ?negGene astro:hasGeneID "[NegativeMarkerGene]" .
    BIND(SQRT(POW(?centerX - ?negX, 2) + POW(?centerY - ?negY, 2)) AS ?negDist)
    FILTER(?negDist < ?radius * 1.5)
  }
  ```

- Meta-rule constraint (PLANNED):
  ```
  # First identify cell candidate from pattern 1
  {
    # Pattern 1 matching logic
    ...
  }
  
  # Then verify pattern 2 exists nearby
  {
    # Pattern 2 matching logic
    ...
    
    # Ensure patterns are spatially related
    BIND(SQRT(POW(?centerX1 - ?centerX2, 2) + POW(?centerY1 - ?centerY2, 2)) AS ?patternDist)
    FILTER(?patternDist < ?maxPatternDistance)
  }
  
  # Combine evidence for final cell
  BIND((?centerX1 + ?centerX2) / 2 AS ?centerX)
  BIND((?centerY1 + ?centerY2) / 2 AS ?centerY)
  BIND(MAX(?conf1, ?conf2) + 0.1 AS ?confidence)  # Boost confidence for dual-pattern match
  ```

### 4. Rule Application

The generated SPARQL rules are applied to the spatial data to identify whole cells:

1. **Rule Execution**:
   - Each SPARQL CONSTRUCT query is executed against the spatial RDF graph
   - When a matching pattern is found, new triples are added to create a cell entity
   - Each cell is linked to the spatial points that contributed to its identification

2. **Multi-class Conflict Resolution** (NEW):
   - When multiple cell types are identified in overlapping regions, we resolve conflicts by:
     - Comparing confidence scores and keeping the higher confidence identification
     - Analyzing the density and pattern of marker genes for each competing cell type
     - Considering negative markers to disambiguate cell types
     - Applying a cell type hierarchy when confidence scores are similar (based on biological knowledge)

3. **Cell Validation**:
   - Cells with confidence below a threshold (default: 0.7) are filtered out
   - We apply a proximity filter to remove cells with too many nearby cells (likely false positives)
   - This step ensures that identified cells are distinct and biologically plausible

4. **Visualization and Export**:
   - Identified cells are visualized with:
     - Cell centers as points
     - Cell types distinguished by color
     - Cell radii shown as circles
   - Cells are exported to CSV format with metadata including:
     - Cell type
     - X/Y coordinates
     - Radius
     - Confidence score
     - Contributing genes

## Statistical Considerations

### Minimum Co-expression Probability

The minimum co-expression probability (default: 0.6) determines when two genes are considered to be co-expressed. This threshold:

- Balances sensitivity vs. specificity
- Was determined empirically by testing different values and evaluating the resulting rules
- May need adjustment for different datasets or cell types

Mathematical definition:
```
Two genes A and B are co-expressed if max(P(A|B), P(B|A)) ≥ threshold
```

### Confidence Scores

Confidence scores for identified cells reflect the reliability of the identification:

- Clique rules have high confidence (0.9) because they use multiple genes
- Pair rules have variable confidence based on the actual co-expression strength
- Negative markers boost confidence scores when successfully applied
- Meta-rules receive the highest confidence scores due to multiple layers of evidence
- The final filtering step (min_confidence = 0.7) removes low-confidence cells

### Distance Threshold

The distance threshold (default: 50 units) determines how close gene expressions must be to be considered part of the same cell. This value:

- Represents the maximum expected diameter of astrocyte cells in the tissue
- Is proportional to the resolution of the spatial transcriptomics data
- Can be adjusted based on the known biology of the cell type being identified

### Rule Co-occurrence Significance

For meta-rules, we assess the significance of rule co-occurrence using:

- Conditional probability: P(Rule2 | Rule1) - how often Rule2 is observed given that Rule1 is observed
- Lift: The ratio of the observed co-occurrence to the expected co-occurrence if rules were independent
- Statistical significance: p-values from Fisher's exact test to determine if the co-occurrence is non-random

## Multi-class Approach

### Comprehensive Cell Type Analysis

With our enhanced approach, we now analyze all cell types present in the single-cell dataset, not just astrocytes. This provides several advantages:

1. **Reduced False Positives**: By modeling all cell types, we can better distinguish between similar cell types and reduce misidentifications
   
2. **Comprehensive Tissue Mapping**: Generate a complete spatial map of all cell types in the tissue

3. **Cell Type Interactions**: Enable analysis of spatial relationships between different cell types

4. **Improved Validation**: Each cell type serves as a control for other cell types, improving overall accuracy

### Cell Type Prioritization

When applying multi-class analysis, we prioritize cell types based on:

1. **Abundance**: More abundant cell types in the single-cell data receive higher priority
2. **Distinctiveness**: Cell types with more distinctive marker genes receive higher priority
3. **Biological relevance**: Cell types of particular interest to the research question receive higher priority

### Conflict Resolution Strategy

When multiple cell types are identified in the same region, conflicts are resolved through:

1. **Confidence-based resolution**: Higher confidence identifications are preferred
2. **Marker specificity**: Cell types identified by more specific markers are preferred
3. **Marker count**: Cell types identified by more markers are preferred
4. **Negative marker validation**: Cell types with fewer negative marker violations are preferred
5. **Meta-rule priority**: Identifications from meta-rules are given highest priority due to their increased specificity

## Validation Approach

We validate the identified cells through several methods:

1. **Internal Consistency**:
   - Checking if identified cells have reasonable sizes (radius)
   - Ensuring cells of the same type are not too densely clustered

2. **Spatial Distribution**:
   - Comparing the spatial distribution of identified cells to known brain anatomy
   - Ensuring the relative proportions of cell types are reasonable

3. **Gene Expression Patterns**:
   - Verifying that identified cells contain expected marker genes for their cell type
   - Checking that the co-expression patterns match what is known about each cell type's biology

4. **Cross-validation with Other Cell Types** (NEW):
   - Verifying that the spatial distributions of different cell types match expected tissue architecture
   - Checking for appropriate segregation or co-localization of certain cell types

5. **Higher-Order Pattern Validation** (PLANNED):
   - Validating the biological relevance of identified meta-rules against literature
   - Confirming that cells identified by meta-rules show expected biological characteristics
   - Comparing meta-rule accuracy to single-rule accuracy in cells with known identities

## Parameter Selection and Tuning

The default parameters for the pipeline were selected based on:

1. **Literature Review**:
   - Known marker genes for cell types
   - Typical cell sizes

2. **Empirical Testing**:
   - Testing different thresholds and examining the results
   - Optimizing for both sensitivity (finding most cells) and specificity (avoiding false positives)

3. **Biological Plausibility**:
   - Ensuring the identified cells match known biological constraints
   - Validating against reference data where available

Key parameters that may require tuning for different datasets:

| Parameter | Description | Default | When to Adjust |
|-----------|-------------|---------|----------------|
| min_expression | Threshold to consider a gene expressed | 0.2 | Adjust based on data normalization |
| min_coexpression | Minimum conditional probability for co-expression | 0.6 | Higher for more specific rules, lower for more sensitive rules |
| top_genes | Number of top genes to analyze per cell type | 30 | Increase for more complex cell types |
| distance_threshold | Maximum distance between co-expressed genes | 50 | Adjust based on cell size and data resolution |
| min_confidence | Minimum confidence for identified cells | 0.7 | Higher for more precision, lower for more recall |
| neg_marker_threshold | Expression threshold for negative markers | 0.1 | Adjust based on baseline noise level |
| meta_rule_lift | Minimum lift threshold for rule association | 3.0 | Higher for more stringent meta-rules |
| meta_rule_p_value | Maximum p-value for rule association | 0.01 | Lower for more statistically significant associations |

## Future Improvements

Potential improvements to the methodology include:

1. **Incorporating More Cell Types**:
   - Extending the approach to identify neurons, oligodendrocytes, and other cell types
   - Developing rules for distinguishing closely related cell subtypes

2. **Enhanced Validation**:
   - Implementing cross-validation with known cell markers
   - Developing quantitative metrics for rule quality

3. **Advanced Rule Learning**:
   - Using machine learning to optimize rule parameters
   - Incorporating feedback from domain experts to refine rules

4. **Integration with Other Data Types**:
   - Combining with morphological data
   - Integrating with immunohistochemistry images

5. **Density-Based Clustering** (TODO):
   - Implementing density-based clustering algorithms (e.g., DBSCAN, HDBSCAN) to improve cell boundary detection
   - This would complement our current co-expression approach by:
     - Identifying clusters of gene expressions with irregular shapes
     - Better handling variations in expression density within cells
     - More accurately estimating cell boundaries beyond simple circles
     - Providing an additional layer of validation for cell identification
   - Integration plan:
     - Use co-expression patterns to identify candidate cell regions
     - Apply density-based clustering to refine boundaries
     - Calculate confidence scores based on both approaches
     - Generate more accurate cell visualizations with realistic morphology

6. **Higher-Order Co-expression Patterns** (PLANNED):
   - Analyzing the co-occurrence of co-expression patterns and cliques
   - Implementing meta-rules that incorporate multiple evidence types
   - This enhancement would:
     - Further increase specificity by requiring multiple patterns
     - Better distinguish between closely related cell types
     - Potentially identify functionally distinct cell subtypes
     - Capture complex biological relationships between different gene modules
   - Implementation approach:
     - Construct a network of co-expression rules
     - Calculate statistical significance of rule co-occurrences
     - Generate SPARQL templates for combining rules
     - Develop visualization tools for meta-rule networks

## References

- Cahoy, J.D. et al. (2008). "A transcriptome database for astrocytes, neurons, and oligodendrocytes: a new resource for understanding brain development and function." *Journal of Neuroscience*, 28(1), 264-278.
- Hodge, R.D. et al. (2019). "Conserved cell types with divergent features in human versus mouse cortex." *Nature*, 573(7772), 61-68.
- Zeisel, A. et al. (2018). "Molecular architecture of the mouse nervous system." *Cell*, 174(4), 999-1014.
- Ester, M. et al. (1996). "A density-based algorithm for discovering clusters in large spatial databases with noise." *Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96)*, 226-231.
- McInnes, L. et al. (2017). "HDBSCAN: Hierarchical density based clustering." *Journal of Open Source Software*, 2(11), 205.
- Han, J. et al. (2007). "Frequent pattern mining: current status and future directions." *Data Mining and Knowledge Discovery*, 15(1), 55-86. 