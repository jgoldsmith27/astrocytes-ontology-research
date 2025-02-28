# Astrocyte Integrated Visualization Tool

This tool provides comprehensive visualization and analysis capabilities for astrocyte data by integrating spatial transcriptomics and single-cell RNA sequencing data. It enables researchers to explore the spatial distribution of different astrocyte types and analyze gene expression patterns within each type.

## Features

### 1. Data Integration
- Combines spatial transcriptomics data with single-cell RNA sequencing data
- Maps astrocyte types between datasets
- Creates an integrated dataset for comprehensive analysis

### 2. Spatial Visualization
- Visualizes the spatial distribution of different astrocyte types (Protoplasmic, Fibrous, Reactive)
- Maps marker gene expression in spatial context
- Provides interactive exploration of spatial patterns

### 3. Gene Expression Analysis
- Analyzes gene expression patterns for each astrocyte type
- Identifies top expressed genes for each type
- Creates heatmaps of marker gene expression

### 4. SPARQL Queries
- Implements semantic queries to analyze gene expression patterns
- Identifies highly expressed genes in specific astrocyte types
- Finds co-expressed genes within each astrocyte type

### 5. Interactive Dashboard
- Provides an interactive web-based dashboard for data exploration
- Allows dynamic selection of visualization types and marker genes
- Displays statistics and insights alongside visualizations

## Output Files

The tool generates several output files:

1. **astrocyte_spatial_distribution.png**: Visualization of the spatial distribution of astrocyte types
2. **gene_expression_heatmap.png**: Heatmap of marker gene expression across astrocyte types
3. **spatial_[GENE]_expression.png**: Spatial distribution of specific marker gene expression
4. **high_expression_genes_by_type.csv**: List of highly expressed genes for each astrocyte type
5. **coexpressed_genes_by_type.csv**: List of co-expressed gene pairs for each astrocyte type

## Usage

### Command Line Usage

```bash
python astrocyte_integrated_visualization.py --spatial <spatial_ontology_file> --single-cell <single_cell_ontology_file> [--output <output_directory>] [--interactive]
```

Parameters:
- `--spatial`: Path to the spatial ontology TURTLE file (required)
- `--single-cell`: Path to the single-cell ontology TURTLE file (required)
- `--output`: Output directory for visualization results (default: ../output/integrated_visualization)
- `--interactive`: Launch the interactive dashboard (optional)

### Using the Shell Script

For convenience, you can use the provided shell script:

```bash
./run_integrated_visualization.sh
```

This script:
1. Sets up the environment
2. Activates the appropriate virtual environment
3. Checks for required packages
4. Runs the visualization tool with the correct parameters

### Interactive Dashboard

The interactive dashboard provides a user-friendly interface for exploring the data. To use it:

1. Run the tool with the `--interactive` flag
2. Open a web browser and navigate to http://127.0.0.1:8050/
3. Use the dropdown menus to select visualization types and marker genes
4. View statistics and insights in the statistics panel

## Requirements

- Python 3.6+
- rdflib
- pandas
- numpy
- matplotlib
- seaborn
- scanpy
- anndata
- plotly
- dash
- dash-bootstrap-components

## Integration with Other Tools

This visualization tool is designed to work with:

1. **Astrocyte Classifier**: Uses the output of the classifier to visualize astrocyte types
2. **Spatial Ontology Builder**: Uses the spatial ontology created by the builder
3. **Single-Cell Ontology Builder**: Uses the single-cell ontology created by the builder

## Example Workflow

1. Run the spatial ontology builder to create the spatial ontology
2. Run the single-cell ontology builder to create the single-cell ontology
3. Run the astrocyte classifier to identify astrocyte types
4. Run this integrated visualization tool to explore the data
5. Use the interactive dashboard to gain insights

## Future Enhancements

- 3D visualization of spatial data
- Temporal analysis of gene expression changes
- Integration with additional data modalities
- Advanced statistical analysis of spatial patterns
- Machine learning-based prediction of astrocyte functions 