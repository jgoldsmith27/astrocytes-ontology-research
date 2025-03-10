# Examples and Visualizations

This document provides examples of co-expression rules, their application to spatial data, and visualizations of the results to help you understand how the pipeline works in practice.

## Example 1: Basic Co-expression Pattern

### Single-cell Data Analysis

In single-cell RNA-seq data, we observe that genes NRXN1 and GPM6A are highly co-expressed in Astrocytes:

```
Gene pair: NRXN1 and GPM6A
P(NRXN1) = 0.83 (expressed in 83% of Astrocytes)
P(GPM6A) = 0.78 (expressed in 78% of Astrocytes)
P(GPM6A|NRXN1) = 0.91 (91% of cells expressing NRXN1 also express GPM6A)
P(NRXN1|GPM6A) = 0.94 (94% of cells expressing GPM6A also express NRXN1)
Joint probability = 0.75 (75% of Astrocytes express both genes)
```

This strong co-expression relationship (maximum conditional probability = 0.94) leads to the generation of a pair rule.

### Generated SPARQL Rule

```sparql
# Rule: Astrocytes_pair_1

PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX astro: <http://example.org/ontology/astrocyte/>

# Rule for identifying Astrocyte based on co-expression of NRXN1 and GPM6A
CONSTRUCT {
    ?cell a astro:SpatialCell ;
         astro:hasCellType astro:Astrocyte ;
         astro:hasConfidence 0.94 ;
         astro:hasX ?centerX ;
         astro:hasY ?centerY ;
         astro:hasRadius ?radius .
    
    ?cell astro:includesSpatialPoint ?spatialPoint1 .
    ?cell astro:includesSpatialPoint ?spatialPoint2 .
}
WHERE {
    # Gene identifiers
    ?gene1 astro:hasGeneID "NRXN1" .
    ?gene2 astro:hasGeneID "GPM6A" .
    
    # Spatial points for each gene
    ?spatialPoint1 a astro:SpatialDataPoint ;
                   astro:expressesGene ?gene1 ;
                   astro:hasXCoordinate ?x1 ;
                   astro:hasYCoordinate ?y1 .
    
    ?spatialPoint2 a astro:SpatialDataPoint ;
                   astro:expressesGene ?gene2 ;
                   astro:hasXCoordinate ?x2 ;
                   astro:hasYCoordinate ?y2 .
    
    # Distance constraint
    BIND(SQRT(POW(?x1 - ?x2, 2) + POW(?y1 - ?y2, 2)) AS ?dist)
    FILTER(?dist < 50)
    
    # Calculate center coordinates
    BIND((?x1 + ?x2) / 2 AS ?centerX)
    BIND((?y1 + ?y2) / 2 AS ?centerY)
    
    # Set radius
    BIND(?dist / 2 AS ?radius)
    
    # Create a unique cell identifier
    BIND(URI(CONCAT("http://example.org/data/astrocyte/cell_", STRUUID())) AS ?cell)
}
```

### Rule Application

When applied to spatial data, this rule finds instances where NRXN1 and GPM6A are expressed within 50 units of each other:

```
Found: NRXN1 at (x=1235, y=782) and GPM6A at (x=1256, y=798)
Distance: 28.6 units
Created Astrocyte cell at (1245.5, 790) with radius 14.3 and confidence 0.94
```

### Visual Representation

```
Spatial Data Points:                Cell Identification:

    y                                   y
    |                                   |
800 |      •GPM6A                  800 |      •GPM6A
    |                                   |     /   \
    |                                   |    /     \
    |                                   |   •Center •
    |                                   |  /         \
780 |•NRXN1                        780 |•NRXN1       |
    |                                   |             |
    +------------> x                    +------------> x
      1220  1240  1260                    1220  1240  1260
```

## Example 2: Complex Clique Pattern

### Single-cell Data Analysis

In single-cell data, we identify a clique of three highly co-expressed genes in Astrocytes1: LRP1B, NRG3, and GPM6A:

```
Gene Clique: {LRP1B, NRG3, GPM6A}

Pairwise co-expression probabilities:
P(NRG3|LRP1B) = 0.88
P(GPM6A|LRP1B) = 0.85
P(GPM6A|NRG3) = 0.92

All genes are co-expressed with each other above the threshold (0.6)
```

This forms a clique (fully connected subgraph) in the co-expression network, leading to a clique rule.

### Generated SPARQL Rule

```sparql
# Rule: Astrocytes1_clique_1

PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX astro: <http://example.org/ontology/astrocyte/>

# Rule for identifying Astrocyte1 based on clique of 3 co-expressed genes
CONSTRUCT {
    ?cell a astro:SpatialCell ;
         astro:hasCellType astro:Astrocyte1 ;
         astro:hasConfidence 0.9 ;
         astro:hasX ?centerX ;
         astro:hasY ?centerY ;
         astro:hasRadius ?radius .
         
    ?cell astro:includesSpatialPoint ?spatialPoint0 .
    ?cell astro:includesSpatialPoint ?spatialPoint1 .
    ?cell astro:includesSpatialPoint ?spatialPoint2 .
}
WHERE {
    # Gene identifiers
    ?gene0 astro:hasGeneID ?gene0ID .
    ?gene1 astro:hasGeneID ?gene1ID .
    ?gene2 astro:hasGeneID ?gene2ID .
    
    # Gene ID filters
    FILTER(?gene0ID = "LRP1B" && ?gene1ID = "NRG3" && ?gene2ID = "GPM6A")
    
    # Spatial points for each gene
    ?spatialPoint0 a astro:SpatialDataPoint ;
                   astro:expressesGene ?gene0 ;
                   astro:hasXCoordinate ?x0 ;
                   astro:hasYCoordinate ?y0 .
    
    ?spatialPoint1 a astro:SpatialDataPoint ;
                   astro:expressesGene ?gene1 ;
                   astro:hasXCoordinate ?x1 ;
                   astro:hasYCoordinate ?y1 .
    
    ?spatialPoint2 a astro:SpatialDataPoint ;
                   astro:expressesGene ?gene2 ;
                   astro:hasXCoordinate ?x2 ;
                   astro:hasYCoordinate ?y2 .
    
    # Distance constraints
    BIND(SQRT(POW(?x0 - ?x1, 2) + POW(?y0 - ?y1, 2)) AS ?dist0_1)
    FILTER(?dist0_1 < 50)
    BIND(SQRT(POW(?x0 - ?x2, 2) + POW(?y0 - ?y2, 2)) AS ?dist0_2)
    FILTER(?dist0_2 < 50)
    BIND(SQRT(POW(?x1 - ?x2, 2) + POW(?y1 - ?y2, 2)) AS ?dist1_2)
    FILTER(?dist1_2 < 50)
    
    # Calculate center coordinates
    BIND((?x0 + ?x1 + ?x2) / 3 AS ?centerX)
    BIND((?y0 + ?y1 + ?y2) / 3 AS ?centerY)
    
    # Calculate radius
    BIND(MAX(?dist0_1, ?dist0_2, ?dist1_2) / 2 AS ?radius)
    
    # Create a unique cell identifier
    BIND(URI(CONCAT("http://example.org/data/astrocyte/cell_", STRUUID())) AS ?cell)
}
```

### Rule Application

When applied to spatial data, this rule finds instances where all three genes are close to each other:

```
Found: LRP1B at (x=2105, y=1482), NRG3 at (x=2128, y=1465), and GPM6A at (x=2112, y=1490)
Distances: LRP1B-NRG3: 29.1 units, LRP1B-GPM6A: 14.3 units, NRG3-GPM6A: 27.8 units
Created Astrocyte1 cell at (2115, 1479) with radius 14.55 and confidence 0.9
```

### Visual Representation

```
                                  GPM6A
                                    •
                                   / \
                                  /   \
                                 /     \
                                /       \
                               /    •    \
                              /   Center  \
                             /             \
                        LRP1B•-------------•NRG3
```

## Example 3: Multiple Rules Applied

When multiple rules are applied to the same region, cells can be identified with different rules:

### Rules Applied

```
Applied rule Astrocytes_pair_1 (NRXN1-GPM6A): Found 215 cells
Applied rule Astrocytes_clique_1 (NRXN1-GPM6A-SLC1A2): Found 87 cells
Applied rule Astrocytes1_pair_1 (LRP1B-GPM6A): Found 312 cells
Applied rule Astrocytes1_clique_1 (LRP1B-NRG3-GPM6A): Found 128 cells
```

### Validation Results

After applying the confidence threshold and proximity filtering:

```
Filtered 62 cells below confidence threshold (0.7)
Filtered 35 cells with too many nearby cells
Final result: 645 cells (294 Astrocytes, 351 Astrocytes1)
```

### Example Output CSV

The resulting CSV file contains information about each identified cell:

```csv
cell_id,cell_type,x,y,radius,confidence,genes
cell_001,Astrocyte,1245.5,790.0,14.3,0.94,"NRXN1,GPM6A"
cell_002,Astrocyte,1356.2,823.7,13.8,0.92,"NRXN1,GPM6A"
cell_003,Astrocyte1,2115.0,1479.0,14.55,0.9,"LRP1B,NRG3,GPM6A"
cell_004,Astrocyte1,2432.1,1732.6,12.2,0.85,"LRP1B,GPM6A"
...
```

### Visualization

The final visualization shows the spatial distribution of identified cells:

```
      Astrocyte Distribution in Tissue
      
      y
      |
2000  |                                    
      |                           *  *
      |                   *    *  ** *  *
      |              *   ** * ****** * **
      |              ***** ********* ***
1500  |           * ******************  *
      |           **********************
      |         * *********************
      |        ***********************
      |      * *********************    
1000  |     ************************    
      |    *************************
      |    ************************
      |   ************************
      |  ************************     
 500  | ***********************     
      |***********************    
      +------------------------------> x
         500  1000  1500  2000  2500
         
Legend: * Astrocytes    * Astrocytes1
```

## Example 4: Failure Cases

### Insufficient Co-expression

When co-expression probabilities are below the threshold, no rule is generated:

```
Gene pair: NRXN1 and MBP
P(MBP|NRXN1) = 0.12
P(NRXN1|MBP) = 0.08
Max conditional probability (0.12) < threshold (0.6)
No rule generated
```

### Overlapping Cells

When multiple cells are identified in very close proximity, the validation step filters out some of them:

```
Found potential cells:
- Astrocyte at (1245, 790) with radius 14.3
- Astrocyte at (1252, 787) with radius 16.1
- Astrocyte at (1248, 795) with radius 15.2

After proximity filtering:
- Kept: Astrocyte at (1245, 790) with radius 14.3 (highest confidence)
- Filtered: 2 nearby cells (likely duplicates)
```

## Example 5: Comparative Results

By adjusting parameters, we can see how the results change:

### Default Parameters

```
--min-expression 0.2 --min-coexpression 0.6 --min-confidence 0.7 --distance-threshold 50
Results: 645 cells (294 Astrocytes, 351 Astrocytes1)
```

### Higher Precision Parameters

```
--min-expression 0.25 --min-coexpression 0.7 --min-confidence 0.8 --distance-threshold 40
Results: 423 cells (182 Astrocytes, 241 Astrocytes1)
Fewer cells, but higher confidence in their identification
```

### Higher Recall Parameters

```
--min-expression 0.15 --min-coexpression 0.5 --min-confidence 0.6 --distance-threshold 60
Results: 912 cells (421 Astrocytes, 491 Astrocytes1)
More cells identified, but potentially more false positives
```

## Visual Comparison of Results

```
Default             Higher Precision      Higher Recall
Parameters          Parameters            Parameters
*  *  *             *     *               * ** * *
 * ** *              *   *                ********
** * **              * *                 **********
 * * *                *                  **********
```

## Example 6: Rule Generation Process

To illustrate how rules are generated from single-cell data:

### 1. Extract Top Genes

```
Top genes for Astrocytes: NRXN1, RORA, NRG3, GPM6A, SLC1A2, CTNNA2, ...
Top genes for Astrocytes1: LRP1B, NRG3, NRXN1, NPAS3, GPM6A, ...
```

### 2. Calculate Co-expression

```
Calculating co-expression for Astrocytes (45 gene pairs)...
Found 35 significant co-expressions above threshold (0.6)
```

### 3. Build Co-expression Network

```
Gene co-expression network for Astrocytes:
NRXN1 -- (0.94) -- GPM6A -- (0.88) -- SLC1A2
  |                  |
  |                  |
(0.85)            (0.76)
  |                  |
  |                  |
RORA ----------------- CTNNA2
        (0.82)
```

### 4. Find Cliques

```
Found cliques in Astrocytes network:
- Clique 1: {NRXN1, GPM6A, SLC1A2} (size 3)
- Clique 2: {NRXN1, RORA, CTNNA2} (size 3)
- Clique 3: {GPM6A, SLC1A2, CTNNA2} (size 3)
```

### 5. Generate Rules

```
Generated rules for Astrocytes:
- Clique rule for {NRXN1, GPM6A, SLC1A2} with confidence 0.9
- Clique rule for {NRXN1, RORA, CTNNA2} with confidence 0.9
- Pair rule for {NRXN1, GPM6A} with confidence 0.94
- Pair rule for {GPM6A, SLC1A2} with confidence 0.88
...
```

## Example 7: Performance Metrics

By running the pipeline with different parameters, we can measure performance metrics:

### Processing Time

```
Default parameters:
- Data conversion: 25 minutes
- Rule generation: 2 minutes
- Rule application: 15 minutes
Total: 42 minutes

With higher data limit (500,000 rows):
- Data conversion: 118 minutes
- Rule generation: 2 minutes
- Rule application: 42 minutes
Total: 162 minutes
```

### Memory Usage

```
Data conversion: peak of 8.5 GB
Rule generation: peak of 2.1 GB
Rule application: peak of 12.3 GB
```

### Rule Statistics

```
Default parameters:
- Number of rules generated: 22
- Average rule execution time: 40 seconds
- Most productive rule: Astrocytes1_pair_1 (312 cells)
```

## Conclusion

These examples demonstrate how the pipeline processes data from raw single-cell and spatial inputs to identified astrocyte cells. By understanding these examples, you should have a clearer picture of how to interpret the results and how to adjust parameters for your specific research needs.

For a more interactive exploration, consider:
1. Running with sample data and different parameters
2. Visualizing intermediate results (e.g., gene co-expression networks)
3. Comparing outputs with different cell types or datasets 