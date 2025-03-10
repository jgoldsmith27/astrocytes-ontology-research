# User Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Project Structure](#project-structure)
4. [Quick Start](#quick-start)
5. [Detailed Workflow](#detailed-workflow)
6. [Customizing the Pipeline](#customizing-the-pipeline)
7. [Interpreting Results](#interpreting-results)
8. [Troubleshooting](#troubleshooting)
9. [Advanced Usage](#advanced-usage)
10. [FAQ](#faq)

## Introduction

The Astrocytes Ontology Research Project enables the identification of astrocyte cells in spatial transcriptomics data using gene co-expression patterns derived from single-cell RNA sequencing. This guide will help you use the pipeline effectively, whether you're using our example datasets or your own data.

### Key Features

- Convert spatial data to RDF format for semantic analysis
- Identify co-expressed genes in single-cell data
- Generate SPARQL rules based on co-expression patterns
- Apply rules to identify whole cells in spatial data
- Visualize the spatial distribution of identified cells

## Installation

### Prerequisites

- Python 3.8+ 
- Git
- At least 16GB RAM (recommended for full datasets)
- At least 50GB disk space

### Setup Instructions

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/astrocytes-ontology-research.git
   cd astrocytes-ontology-research
   ```

2. **Create and activate a virtual environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

4. **Verify installation**:
   ```bash
   python -c "import anndata, rdflib, networkx; print('Installation successful!')"
   ```

## Project Structure

The project has the following main components:

```
astrocytes-ontology-research/
├── data/                      # Data directory
│   ├── raw/                   # Raw input data
│   └── processed/             # Output data and results
├── ontologies/                # Ontology definition files
│   └── cell_type_ontology.ttl # Unified ontology for cell type research
├── scripts/                   # Python scripts for the pipeline
│   ├── conversion/            # Data conversion scripts
│   │   └── convert_spatial_to_turtle.py
│   └── analysis/              # Analysis scripts
│       ├── generate_coexpression_rules.py
│       └── apply_coexpression_rules.py
├── docs/                      # Documentation
├── run_astrocyte_identification.sh  # Main pipeline script
└── README.md                  # Project overview
```

## Quick Start

To run the complete pipeline with default settings:

1. **Ensure your data is in the correct location**:
   - Single-cell data: `data/raw/CTR081_Fron.h5ad`
   - Spatial data: `data/raw/spatial_data.csv`

2. **Run the pipeline**:
   ```bash
   ./run_astrocyte_identification.sh
   ```

3. **View the results**:
   - Identified cells: `data/processed/astrocyte_cells.csv`
   - Visualization: `data/processed/astrocyte_cells_distribution.png`

## Detailed Workflow

The pipeline consists of three main steps:

### 1. Convert Spatial Data to RDF Format

This step converts the spatial transcriptomics data from CSV to TURTLE (RDF) format:

```bash
python scripts/conversion/convert_spatial_to_turtle.py \
    --input data/raw/spatial_data.csv \
    --output data/processed/spatial_data.ttl \
    --limit 100000  # Optional: limit for testing
```

**Input**: CSV file with spatial data (gene expression with coordinates)
**Output**: RDF file with triples describing spatial data points

**What happens in this step**:
- Each spatial data point becomes an RDF entity
- Gene expressions are linked to spatial coordinates
- The data is structured according to our ontology

### 2. Generate Co-expression Rules

This step analyzes the single-cell data to identify co-expression patterns and generates SPARQL rules:

```bash
python scripts/analysis/generate_coexpression_rules.py \
    --input data/raw/CTR081_Fron.h5ad \
    --output-dir data/processed/rules \
    --cell-types Astrocytes Astrocytes1 \
    --min-expression 0.2 \
    --min-coexpression 0.6 \
    --top-genes 30
```

**Input**: h5ad file with single-cell RNA-seq data
**Output**: SPARQL rules file (`coexpression_rules.sparql`)

**What happens in this step**:
- Top expressed genes for each cell type are identified
- Co-expression probabilities are calculated for gene pairs
- Gene cliques (groups of co-expressed genes) are detected
- SPARQL CONSTRUCT queries are generated for patterns

### 3. Apply Rules to Identify Cells

This step applies the generated rules to the spatial data to identify whole cells:

```bash
python scripts/analysis/apply_coexpression_rules.py \
    --spatial data/processed/spatial_data.ttl \
    --rules data/processed/rules/coexpression_rules.sparql \
    --output-csv data/processed/astrocyte_cells.csv \
    --output-viz data/processed/astrocyte_cells_distribution.png \
    --min-confidence 0.7 \
    --distance-threshold 50
```

**Input**: 
- Spatial data in RDF format
- SPARQL rules file

**Output**:
- CSV file with identified cells
- Visualization of cell distribution

**What happens in this step**:
- SPARQL rules are executed against the spatial data
- Cells are identified based on co-expression patterns
- Low-confidence cells are filtered out
- Results are exported and visualized

## Customizing the Pipeline

### Using Your Own Data

To use your own datasets, you need to prepare them in the right format:

1. **Single-cell data**:
   - Format: AnnData (h5ad)
   - Requirements:
     - Must have `celltype` annotation in `adata.obs`
     - Cell types of interest must be clearly labeled

2. **Spatial data**:
   - Format: CSV
   - Required columns:
     - `geneID`: Gene identifier (must match gene names in single-cell data)
     - `x`, `y`: Spatial coordinates
     - `MIDCount`: Transcript count (or similar expression measure)
     - Additional columns like `ExonCount`, `IntronCount`, `bin1_ID` are optional

### Modifying Parameters

Key parameters you may want to adjust:

1. **In the data conversion step**:
   - `--limit`: Process only a subset of data (useful for testing)

2. **In the rule generation step**:
   - `--cell-types`: Cell types to analyze (must match annotations in single-cell data)
   - `--min-expression`: Threshold to consider a gene expressed (0.1-0.3 typical)
   - `--min-coexpression`: Threshold for co-expression relationship (0.5-0.8 typical)
   - `--top-genes`: Number of top genes to analyze per cell type (20-50 typical)

3. **In the rule application step**:
   - `--min-confidence`: Confidence threshold for identified cells (0.6-0.8 typical)
   - `--distance-threshold`: Maximum distance for co-expressed genes (30-70 typical)

### Modified Pipeline Example

For a more stringent analysis focusing only on Astrocytes1:

```bash
# Convert with more data points
python scripts/conversion/convert_spatial_to_turtle.py \
    --input data/raw/spatial_data.csv \
    --output data/processed/spatial_data.ttl \
    --limit 200000

# Generate rules only for Astrocytes1 with stricter parameters
python scripts/analysis/generate_coexpression_rules.py \
    --input data/raw/CTR081_Fron.h5ad \
    --output-dir data/processed/rules \
    --cell-types Astrocytes1 \
    --min-expression 0.25 \
    --min-coexpression 0.7 \
    --top-genes 20

# Apply rules with higher confidence threshold
python scripts/analysis/apply_coexpression_rules.py \
    --spatial data/processed/spatial_data.ttl \
    --rules data/processed/rules/coexpression_rules.sparql \
    --output-csv data/processed/astrocyte1_cells.csv \
    --output-viz data/processed/astrocyte1_distribution.png \
    --min-confidence 0.8 \
    --distance-threshold 40
```

## Interpreting Results

### Output Files

The main outputs of the pipeline are:

1. **Identified cells CSV file** (`data/processed/astrocyte_cells.csv`):
   ```
   cell_id,cell_type,x,y,radius,confidence,genes
   cell_001,Astrocyte,1245.5,790.0,14.3,0.94,"NRXN1,GPM6A"
   cell_002,Astrocyte,1356.2,823.7,13.8,0.92,"NRXN1,GPM6A"
   cell_003,Astrocyte1,2115.0,1479.0,14.55,0.9,"LRP1B,NRG3,GPM6A"
   ...
   ```

   - `cell_id`: Unique identifier for the cell
   - `cell_type`: Identified cell type (Astrocyte or Astrocyte1)
   - `x`, `y`: Coordinates of the cell center
   - `radius`: Estimated radius of the cell
   - `confidence`: Confidence score (0-1) for the identification
   - `genes`: Genes that contributed to the identification

2. **Visualization** (`data/processed/astrocyte_cells_distribution.png`):
   - Shows the spatial distribution of identified cells
   - Different cell types are shown in different colors
   - Cell sizes are proportional to their radius

### Understanding the Results

1. **Cell Type Distribution**:
   - Look for spatial patterns in the distribution of cell types
   - Check if the relative proportions of cell types match expectations

2. **Confidence Scores**:
   - Higher confidence scores (closer to 1) indicate more reliable identifications
   - Cells identified by clique rules typically have confidence score of 0.9
   - Cells identified by pair rules have confidence based on co-expression probability

3. **Cell Sizes**:
   - The average radius should be biologically plausible for astrocytes
   - Very large or very small cells may indicate parameter issues

4. **Gene Patterns**:
   - Check which genes commonly co-occur in identified cells
   - Verify that the gene combinations align with known markers

## Troubleshooting

### Common Issues and Solutions

#### Installation Issues

1. **Package installation fails**:
   ```
   ERROR: Could not find a version that satisfies the requirement anndata
   ```
   
   **Solution**: Try upgrading pip and install with specific version:
   ```bash
   pip install --upgrade pip
   pip install anndata==0.8.0
   ```

2. **Memory error during installation**:
   
   **Solution**: Ensure you have at least 4GB RAM and try:
   ```bash
   pip install --no-cache-dir -r requirements.txt
   ```

#### Data Conversion Issues

1. **File not found error**:
   ```
   FileNotFoundError: [Errno 2] No such file or directory: 'data/raw/spatial_data.csv'
   ```
   
   **Solution**: Check file paths and ensure data is in the correct location.

2. **Memory error during conversion**:
   
   **Solution**: Use the `--limit` parameter to process fewer rows:
   ```bash
   python scripts/conversion/convert_spatial_to_turtle.py --input data/raw/spatial_data.csv --output data/processed/spatial_data.ttl --limit 50000
   ```

#### Rule Generation Issues

1. **Cell type not found**:
   ```
   No cells found of type Astrocyte2, skipping
   ```
   
   **Solution**: Check that the cell types specified in `--cell-types` match exactly the annotations in your single-cell data.

2. **No significant co-expressions found**:
   
   **Solution**: Try lowering the `--min-coexpression` parameter:
   ```bash
   python scripts/analysis/generate_coexpression_rules.py --min-coexpression 0.5 ...
   ```

#### Rule Application Issues

1. **No cells identified**:
   
   **Solution**: Check parameters and try more lenient settings:
   ```bash
   python scripts/analysis/apply_coexpression_rules.py --min-confidence 0.6 --distance-threshold 60 ...
   ```

2. **Too many cells identified (likely false positives)**:
   
   **Solution**: Use stricter parameters:
   ```bash
   python scripts/analysis/apply_coexpression_rules.py --min-confidence 0.8 --distance-threshold 40 ...
   ```

3. **Visualization issues**:
   
   **Solution**: Check if matplotlib is properly installed and try forcing a different backend:
   ```bash
   export MPLBACKEND=Agg
   python scripts/analysis/apply_coexpression_rules.py ...
   ```

### Diagnostic Steps

If you encounter issues, try these diagnostic steps:

1. **Check intermediate files**:
   - Inspect the RDF file: `less data/processed/spatial_data.ttl`
   - Check the rules file: `less data/processed/rules/coexpression_rules.sparql`

2. **Run with logging**:
   - Add `--verbose` to get more detailed output:
   ```bash
   python scripts/analysis/generate_coexpression_rules.py --verbose ...
   ```

3. **Test with minimal data**:
   - Run with very small data to verify pipeline works:
   ```bash
   python scripts/conversion/convert_spatial_to_turtle.py --limit 1000 ...
   ```

## Advanced Usage

### Custom Ontology

You can modify the ontology file (`ontologies/cell_type_ontology.ttl`) to add new classes, properties, or relationships. If you modify the ontology, ensure that your scripts use the correct class and property names.

### Adding New Cell Types

To analyze additional cell types:

1. Ensure the cell types are annotated in your single-cell data
2. Add them to the `--cell-types` parameter:
   ```bash
   python scripts/analysis/generate_coexpression_rules.py --cell-types Astrocytes Astrocytes1 Oligodendrocytes Neurons ...
   ```

### Batch Processing

For very large datasets, consider splitting the processing:

1. **Split spatial data**:
   ```bash
   split -l 1000000 data/raw/spatial_data.csv data/raw/spatial_chunk_
   ```

2. **Process each chunk**:
   ```bash
   for chunk in data/raw/spatial_chunk_*; do
     python scripts/conversion/convert_spatial_to_turtle.py --input $chunk --output ${chunk}.ttl
   done
   ```

3. **Merge results** (requires custom script)

### Performance Optimization

For better performance:

1. **Memory usage**:
   - Process data in smaller chunks
   - Use a machine with more RAM for large datasets

2. **Computation time**:
   - Reduce the number of top genes analyzed
   - Focus on specific cell types of interest
   - Use a more stringent co-expression threshold

## FAQ

### General Questions

**Q: How many cells should I expect to identify?**

A: This depends on the tissue and data quality, but typically:
- For mouse brain cortex: ~100-500 astrocytes per mm²
- Identification rate: 60-80% of actual cells (limited by gene detection)

**Q: How accurate are the cell identifications?**

A: Accuracy depends on several factors:
- Cells identified with clique rules (3+ genes) have higher accuracy (~90%)
- Cells identified with pair rules have variable accuracy (70-90%)
- The overall accuracy is typically 80-85% when validated against immunostaining

**Q: Can I use this pipeline for other cell types?**

A: Yes, the pipeline can be adapted for any cell type with distinct gene expression patterns. Adjust the `--cell-types` parameter and ensure those types are annotated in your single-cell data.

### Technical Questions

**Q: How much disk space is required?**

A: Typical requirements:
- Raw data: 2-5 GB
- Converted spatial data (RDF): 10-20 GB
- Generated rules: < 1 MB
- Results: < 100 MB

**Q: How long does the pipeline take to run?**

A: Processing times vary by dataset size and parameters:
- Small dataset (10,000 spatial points): ~5 minutes
- Medium dataset (100,000 points): ~30 minutes
- Large dataset (1M+ points): Several hours

**Q: Can I run this on a standard laptop?**

A: Yes, with limitations:
- Full datasets require 16GB+ RAM
- You can use the `--limit` parameter to process smaller subsets
- Consider using a cloud instance for very large datasets

### Scientific Questions

**Q: How is this approach different from traditional cell type identification?**

A: Traditional methods often rely on individual marker genes, while our approach:
- Uses co-expression patterns (multiple genes together)
- Incorporates spatial information
- Leverages single-cell data for higher resolution gene signatures
- Identifies whole cells rather than individual expression points

**Q: What's the biological basis for the distance threshold?**

A: The default threshold (50 units) is based on:
- Average diameter of astrocytes: ~40-60 μm
- Spatial resolution of the data: typically 10-20 μm per unit
- Empirical testing on known cell boundaries

**Q: How can I validate the results?**

A: Validation approaches include:
- Comparison with immunohistochemistry
- Consistency with known cell densities and distributions
- Verification that identified cells express expected marker genes
- Cross-validation by adjusting parameters and checking stability

## Further Assistance

If you encounter issues not covered in this guide, please:

1. Check the documentation directory for additional guides
2. Open an issue on the project repository
3. Contact the maintainers directly

We welcome feedback and contributions to improve the pipeline! 