# Single Cell Astrocyte Ontology

This project implements a Semantic Web ontology for single-cell RNA sequencing data of astrocytes using RDF and TURTLE format. The ontology models gene expression data at the single-cell level, allowing for semantic queries and analysis of gene expression patterns in astrocytes.

## Dataset

The dataset used in this project is a single-cell RNA sequencing dataset of astrocytes:
- File: `single_cell_data_cleaned.h5ad`
- Contains gene expression data for individual astrocytes
- Additional processed data: `astrocyte_expression.csv`

## Ontology Structure

The ontology models the following key concepts:

- **Gene**: Represents a gene with expression level and other attributes
- **Cell**: Represents an individual cell with its gene expression profile
- **CellCluster**: Represents a cluster of cells with similar expression profiles
- **ExpressionProfile**: Represents the expression profile of a cell

Key relationships:
- `expressedIn`: Links a gene to a cell
- `hasExpressionLevel`: Indicates the expression level of a gene in a cell
- `belongsToCluster`: Associates cells with clusters
- `coExpressedWith`: Connects genes that are co-expressed

## Scripts

The project includes several Python scripts for processing and analyzing the data:

1. **single_cell_ontology_builder.py**: Converts single-cell data to TURTLE format
   - Usage: `python scripts/single_cell_ontology_builder.py`

2. **preprocess_single_cell.py**: Preprocesses the raw single-cell data
   - Usage: `python scripts/preprocess_single_cell.py`

3. **compute_coexpression.py**: Computes co-expression networks from single-cell data
   - Usage: `python scripts/compute_coexpression.py`

4. **generate_h5ad.py**: Generates H5AD files from processed data
   - Usage: `python scripts/generate_h5ad.py`

5. **read_data.py**: Utility functions for reading data
   - Usage: Import in other scripts

6. **main.py**: Main entry point for the application
   - Usage: `python scripts/main.py`

## Testing

The project includes several test scripts:

1. **run_all_tests.py**: Runs all tests
   - Usage: `python scripts/run_all_tests.py`

2. **test_astrocyte_filter.py**: Tests for astrocyte filtering
   - Usage: `python scripts/test_astrocyte_filter.py`

3. **test_integration.py**: Integration tests
   - Usage: `python scripts/test_integration.py`

4. **test_ontology_builder.py**: Tests for ontology building
   - Usage: `python scripts/test_ontology_builder.py`

5. **test_read_data.py**: Tests for data reading functions
   - Usage: `python scripts/test_read_data.py`

## Requirements

- Python 3.6+
- pandas
- numpy
- scanpy
- anndata
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

# Run all tests
python scripts/run_all_tests.py

# Generate ontology
python scripts/single_cell_ontology_builder.py

# Run main application
python scripts/main.py
```

## Future Work

- Integration with spatial transcriptomics data for spatial-temporal analysis
- Enhanced visualization of single-cell clusters
- Development of a web interface for interactive exploration
- Advanced semantic queries for cell type identification

## References

- RDF: https://www.w3.org/RDF/
- OWL: https://www.w3.org/OWL/
- SPARQL: https://www.w3.org/TR/sparql11-query/
- Single-cell RNA sequencing: https://www.nature.com/articles/s41576-018-0088-9
