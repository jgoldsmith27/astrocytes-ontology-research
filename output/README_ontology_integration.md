# Astrocyte Ontology Integration Tool

This tool provides advanced semantic integration between spatial transcriptomics and single-cell RNA sequencing ontologies. It creates explicit semantic relationships between entities in both ontologies, enabling powerful cross-ontology queries and knowledge discovery.

## Features

### 1. Bridge Ontology
- Creates a bridge ontology with properties that connect spatial and single-cell entities
- Defines equivalence classes for astrocyte types across ontologies
- Establishes formal semantic relationships between corresponding entities

### 2. Entity Mapping
- Maps genes between spatial and single-cell ontologies using gene IDs
- Creates explicit connections between astrocyte types in both ontologies
- Establishes `owl:sameAs` relationships for semantic equivalence

### 3. Spatial Inference
- Infers spatial locations for single cells based on cell type
- Assigns single cells to spatial points with matching astrocyte types
- Prioritizes high-probability spatial points for assignment

### 4. Inference Rules
- Implements SPARQL-based inference rules to derive new knowledge
- Propagates gene expression information across ontologies
- Infers cell type information for spatial points based on single-cell data

### 5. Advanced Queries
- Enables powerful cross-ontology queries that span both data modalities
- Supports identification of genes expressed in both modalities
- Facilitates discovery of differentially expressed genes between astrocyte types

## Semantic Relationships

The tool establishes several key semantic relationships:

1. **hasSpatialRepresentation**: Links a cell from single-cell data to its spatial representation
2. **hasSingleCellData**: Links a spatial point to its single-cell data
3. **hasCorrespondingGene**: Links a gene in single-cell data to its corresponding gene in spatial data
4. **hasCorrespondingCellType**: Links a cell type in single-cell data to its corresponding type in spatial data

## Example Queries

The integrated ontology supports sophisticated queries such as:

1. **Cross-Modality Gene Expression**:
   ```sparql
   SELECT ?geneID ?scExprLevel ?stExprLevel
   WHERE {
       # Get gene from single-cell data
       ?scGene rdf:type sc:Gene ;
               sc:hasGeneID ?geneID .
       
       # Get corresponding gene in spatial data
       ?scGene bridge:hasCorrespondingGene ?stGene .
       
       # Get expression in single-cell
       ?scExpr sc:forGene ?scGene ;
               sc:hasExpressionLevel ?scExprLevel .
       
       # Get expression in spatial data
       ?stGene st:expressedAt ?point ;
               st:hasExpressionLevel ?stExprLevel .
   }
   ```

2. **Spatial Regions with Cell Types**:
   ```sparql
   SELECT ?cellType ?x ?y ?geneID ?expressionLevel
   WHERE {
       # Get spatial point with cell type
       ?point st:hasXCoordinate ?x ;
              st:hasYCoordinate ?y ;
              st:hasAstrocyteType ?typeNode .
       
       ?typeNode st:astrocyteType ?cellType ;
                 st:typeProbability ?probability .
       
       # Only consider high probability assignments
       FILTER(?probability > 0.8)
       
       # Get genes expressed at this point
       ?gene st:expressedAt ?point ;
             st:hasExpressionLevel ?expressionLevel ;
             st:hasGeneID ?geneID .
       
       # Only consider high expression
       FILTER(?expressionLevel > 1.0)
   }
   ```

3. **Differentially Expressed Genes**:
   ```sparql
   SELECT ?cellType1 ?cellType2 ?geneID 
          (AVG(?expr1) as ?avgExpr1) (AVG(?expr2) as ?avgExpr2)
          (ABS(AVG(?expr1) - AVG(?expr2)) as ?diffExpr)
   WHERE {
       # Get expression in first cell type
       ?cell1 sc:hasCellType ?cellType1 .
       ?expr1Node sc:belongsToCell ?cell1 ;
                  sc:forGene ?gene ;
                  sc:hasExpressionLevel ?expr1 .
       
       # Get expression in second cell type
       ?cell2 sc:hasCellType ?cellType2 .
       ?expr2Node sc:belongsToCell ?cell2 ;
                  sc:forGene ?gene ;
                  sc:hasExpressionLevel ?expr2 .
       
       # Get gene ID
       ?gene sc:hasGeneID ?geneID .
       
       # Ensure we're comparing different cell types
       FILTER(?cellType1 != ?cellType2)
       
       # Only consider astrocyte types
       FILTER(?cellType1 IN ("Protoplasmic", "Fibrous", "Reactive"))
       FILTER(?cellType2 IN ("Protoplasmic", "Fibrous", "Reactive"))
       
       # Ensure consistent ordering of cell types
       FILTER(STR(?cellType1) < STR(?cellType2))
   }
   GROUP BY ?cellType1 ?cellType2 ?geneID
   HAVING (ABS(AVG(?expr1) - AVG(?expr2)) > 0.5)
   ORDER BY DESC(?diffExpr)
   ```

## Usage

### Command Line Usage

```bash
python ontology_integration.py --spatial <spatial_ontology_file> --single-cell <single_cell_ontology_file> [--output <output_file>]
```

Parameters:
- `--spatial`: Path to the spatial ontology TURTLE file (required)
- `--single-cell`: Path to the single-cell ontology TURTLE file (required)
- `--output`: Path to the output TURTLE file (default: ../output/integrated_ontology.ttl)

### Using the Shell Script

For convenience, you can use the provided shell script:

```bash
./run_ontology_integration.sh
```

This script:
1. Sets up the environment
2. Activates the appropriate virtual environment
3. Checks for required packages
4. Runs the ontology integration tool with the correct parameters

## Requirements

- Python 3.6+
- rdflib
- pandas
- numpy

## Integration with Other Tools

This ontology integration tool is designed to work with:

1. **Spatial Ontology Builder**: Uses the spatial ontology created by the builder
2. **Single-Cell Ontology Builder**: Uses the single-cell ontology created by the builder
3. **Integrated Visualization Tool**: Can use the integrated ontology for enhanced visualizations

## Example Workflow

1. Run the spatial ontology builder to create the spatial ontology
2. Run the single-cell ontology builder to create the single-cell ontology
3. Run this ontology integration tool to create semantic connections
4. Use the integrated ontology for advanced cross-ontology queries
5. Visualize the results using the integrated visualization tool

## Advantages Over Simple Mapping

Unlike simple string-based mapping between datasets, this tool:

1. Creates formal semantic relationships using RDF and OWL
2. Enables reasoning and inference across ontologies
3. Supports complex queries that span both data modalities
4. Provides a standardized way to integrate additional data sources
5. Leverages the full power of semantic web technologies

## Future Enhancements

- Integration with external ontologies (e.g., Gene Ontology, Cell Ontology)
- Support for more complex inference rules
- Implementation of a SPARQL endpoint for web-based querying
- Addition of temporal data for dynamic analysis
- Integration with machine learning for predictive modeling 