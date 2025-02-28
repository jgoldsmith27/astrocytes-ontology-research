# Astrocyte-Focused Ontology Integration Tool

## Overview

The Astrocyte-Focused Ontology Integration Tool is a specialized version of our ontology integration framework that focuses exclusively on astrocytes to simplify the analysis and visualization of astrocyte data across spatial transcriptomics and single-cell RNA sequencing modalities.

## Purpose

This tool addresses the specific challenge of integrating astrocyte data from different experimental modalities by:

1. Creating semantic connections between astrocyte types identified in spatial transcriptomics and single-cell RNA sequencing
2. Enabling cross-modality queries specifically for astrocyte biology
3. Providing specialized visualizations for astrocyte distribution, gene expression, and connectivity
4. Validating astrocyte cell type classifications using industry-standard co-expression analysis

## Features

### Astrocyte-Specific Bridge Ontology

The tool establishes formal semantic relationships between astrocyte entities in spatial and single-cell data:

- **Astrocyte Type Mapping**: Explicitly maps three key astrocyte types (Protoplasmic, Fibrous, and Reactive) between modalities
- **Spatial Representation**: Links single-cell astrocytes to their inferred spatial locations
- **Gene Expression Correlation**: Connects gene expression data across modalities specifically for astrocytes

### Astrocyte-Focused Queries

The integrated ontology supports specialized queries for astrocyte research:

```sparql
# Find genes expressed in astrocytes in both modalities
SELECT ?geneID ?scExprLevel ?stExprLevel ?cellType
WHERE {
    # Get gene from single-cell data
    ?scGene rdf:type sc:Gene ;
            sc:hasGeneID ?geneID .
    
    # Get corresponding gene in spatial data
    ?scGene bridge:hasCorrespondingGene ?stGene .
    
    # Get expression in single-cell for astrocytes only
    ?scExpr sc:forGene ?scGene ;
            sc:hasExpressionLevel ?scExprLevel ;
            sc:belongsToCell ?cell .
    
    ?cell sc:hasCellType ?cellType .
    
    # Only include astrocyte cell types
    FILTER(?cellType IN ("Protoplasmic", "Fibrous", "Reactive"))
    
    # Get expression in spatial data
    ?stGene st:expressedAt ?point ;
            st:hasExpressionLevel ?stExprLevel .
    
    # Get astrocyte type at this point
    ?point st:hasAstrocyteType ?typeNode .
    ?typeNode st:astrocyteType ?astrocyteType .
    
    # Only include points with astrocyte types
    FILTER(?astrocyteType IN ("Protoplasmic", "Fibrous", "Reactive"))
}
```

### Astrocyte-Specific Visualizations

The tool includes specialized visualizations for astrocyte data:

1. **Astrocyte Distribution**: Spatial distribution of different astrocyte types
2. **Brain Region Analysis**: Visualization of astrocytes by brain region
3. **Gene Expression Patterns**: Expression of specific genes across astrocyte types
4. **Astrocyte Connectivity**: Network visualization of astrocyte connectivity
5. **Cross-Modality Integration**: Visual representation of how single-cell and spatial astrocyte data are integrated

### Astrocyte Cell Type Validation

The tool provides comprehensive validation of astrocyte cell type classifications:

1. **Co-expression Network Analysis**: Generates gene co-expression networks for each astrocyte type
2. **Marker Gene Enrichment**: Validates cell types based on expression of canonical marker genes
3. **Differential Expression Analysis**: Identifies genes that distinguish astrocyte types
4. **Literature Comparison**: Compares identified markers with literature-reported markers
5. **Spatial Validation**: Visualizes validation results in spatial context

## Usage

### Command Line Usage

To integrate astrocyte data from spatial and single-cell ontologies:

```bash
python Spatial_Data/scripts/ontology_integration.py \
    --spatial output/enhanced_spatial_ontology.ttl \
    --single-cell Single_Cell/output/single_cell_ontology.ttl \
    --output output/integrated_astrocyte_ontology.ttl
```

To generate astrocyte-specific visualizations:

```bash
python Spatial_Data/scripts/astrocyte_visualization.py \
    --input output/integrated_astrocyte_ontology.ttl \
    --output-dir output/visualizations
```

To validate astrocyte cell type classifications:

```bash
python Spatial_Data/scripts/astrocyte_validation.py \
    --input output/integrated_astrocyte_ontology.ttl \
    --output-dir output/astrocyte_validation
```

### Requirements

- Python 3.8+
- rdflib
- pandas
- numpy
- matplotlib
- seaborn
- scanpy
- anndata
- networkx
- matplotlib-venn

## Advantages of Astrocyte-Focused Integration

By focusing specifically on astrocytes, this tool offers several advantages:

1. **Simplified Analysis**: Reduces complexity by focusing on a single cell type of interest
2. **Specialized Queries**: Enables more targeted questions about astrocyte biology
3. **Optimized Visualizations**: Provides visualizations specifically designed for astrocyte data
4. **Improved Performance**: Reduces computational requirements by limiting the scope of integration
5. **Enhanced Interpretability**: Makes results more interpretable by focusing on a specific biological context
6. **Validation Framework**: Ensures cell type classifications are accurate using industry-standard approaches

## Example Workflow

1. **Prepare Data**: Ensure spatial and single-cell ontologies contain astrocyte annotations
2. **Run Integration**: Execute the integration script to create the astrocyte-focused integrated ontology
3. **Generate Visualizations**: Run the visualization script to create astrocyte-specific visualizations
4. **Validate Cell Types**: Run the validation script to verify astrocyte cell type classifications
5. **Analyze Results**: Examine the visualizations, validation results, and query the integrated ontology for insights into astrocyte biology

## Validation Methodology

The validation framework uses industry-standard approaches to verify astrocyte cell type classifications:

1. **Marker Gene Enrichment**: Validates cell types based on expression of canonical marker genes from literature
2. **Co-expression Analysis**: Identifies modules of co-expressed genes characteristic of each astrocyte type
3. **Differential Expression**: Identifies genes that are significantly differentially expressed between astrocyte types
4. **Spatial Context**: Examines whether cell type classifications are spatially coherent
5. **Literature Comparison**: Compares identified marker genes with established markers from literature

## Future Enhancements

1. **Additional Astrocyte Subtypes**: Expand to include more nuanced astrocyte subtype classifications
2. **Temporal Dynamics**: Add support for analyzing changes in astrocyte populations over time
3. **Disease State Analysis**: Incorporate disease-specific annotations for astrocytes
4. **Integration with External Astrocyte Databases**: Connect to specialized astrocyte databases and resources
5. **Advanced Astrocyte Morphology Analysis**: Add tools for analyzing astrocyte morphology across spatial contexts
6. **Multi-modal Validation**: Extend validation to include additional data modalities like proteomics or ATAC-seq

## Contact

For questions or suggestions regarding the Astrocyte-Focused Ontology Integration Tool, please contact the research team at [your-email@example.com]. 