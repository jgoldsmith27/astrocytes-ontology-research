# Spatial Transcriptomics Ontology

This project implements a Semantic Web ontology for spatial transcriptomics data using RDF and TURTLE format. The ontology models gene expression data with spatial coordinates, allowing for semantic queries and analysis of spatial relationships between genes.

## Dataset

The dataset used in this project is a spatial transcriptomics dataset from a brain tissue sample:
- File: `1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv`
- Size: ~9 million rows
- Contains gene expression data with spatial coordinates (x, y) and expression levels

## Ontology Structure

The ontology models the following key concepts:

- **Gene**: Represents a gene with expression level and other attributes
- **SpatialPoint**: Represents a spatial location with x, y coordinates
- **GeneCluster**: Represents a cluster of genes that are spatially related
- **SpatialRegion**: Represents a region in the spatial domain

Key relationships:
- `expressedAt`: Links a gene to its spatial location
- `locatedNear`: Connects spatially proximate points
- `belongsToCluster`: Associates genes with spatial clusters
- `hasExpressionLevel`: Indicates the expression level of a gene

## Scripts

The project includes several Python scripts for processing and analyzing the data:

1. **spatial_ontology_builder.py**: Converts CSV data to TURTLE format with basic spatial relationships
   - Usage: `python spatial_ontology_builder.py`
   - Parameters:
     - `sample_size`: Number of rows to sample (default: 200)
     - `proximity_threshold`: Distance threshold for spatial relationships (default: 100)
     - `max_neighbors`: Maximum number of neighbors per point (default: 3)

2. **spatial_ontology_advanced.py**: Enhanced ontology builder with clustering and advanced spatial analysis
   - Usage: `python spatial_ontology_advanced.py [--sample N] [--proximity P] [--eps E] [--min_samples M]`
   - Parameters:
     - `--sample`: Number of rows to sample (default: 200)
     - `--proximity`: Distance threshold for spatial relationships (default: 100)
     - `--eps`: DBSCAN clustering epsilon parameter (default: 200)
     - `--min_samples`: DBSCAN minimum samples parameter (default: 3)

3. **query_spatial_ontology.py**: Basic queries for exploring the ontology
   - Usage: `python query_spatial_ontology.py`

4. **advanced_queries.py**: Advanced semantic queries and visualizations
   - Usage: `python advanced_queries.py`
   - Generates various visualizations:
     - Cluster sizes
     - Gene distribution across clusters
     - Spatial network analysis
     - Co-expression networks
     - Spatial expression patterns

## Visualizations

The scripts generate several visualizations to help understand the spatial patterns:

- **spatial_visualization.png**: Basic visualization of gene expression in space
- **gene_clusters.png**: Visualization of gene clusters identified by DBSCAN
- **spatial_network.png**: Network visualization of spatial relationships
- **cluster_sizes.png**: Bar chart of cluster sizes
- **gene_cluster_heatmap.png**: Heatmap of gene distribution across clusters
- **co_expression_network.png**: Network of co-expressed genes
- **expression_pattern_*.png**: Spatial expression patterns for top genes

## Requirements

- Python 3.6+
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn
- networkx
- rdflib

## Usage Example

```bash
# Activate virtual environment
source venv/bin/activate

# Generate basic TURTLE file
python spatial_ontology_builder.py

# Generate advanced TURTLE file with clustering
python spatial_ontology_advanced.py --sample 500 --proximity 100

# Run advanced queries and visualizations
python advanced_queries.py
```

## Future Work

- Development of a web interface for interactive exploration
- Scaling to handle the full dataset with distributed processing
- Integration with single cell dataset ontology for cell detection and segmentation

## References

- RDF: https://www.w3.org/RDF/
- OWL: https://www.w3.org/OWL/
- SPARQL: https://www.w3.org/TR/sparql11-query/
- Spatial Transcriptomics: https://www.nature.com/articles/s41592-019-0548-y
