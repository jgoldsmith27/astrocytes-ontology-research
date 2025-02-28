# Astrocyte Ontology Research: Scripts Organization

This document describes the organization of scripts in the Astrocyte Ontology Research project.

## Overview

The scripts have been reorganized to focus specifically on astrocyte-related functionality and to improve code organization. The scripts are now organized into the following categories:

- **Core**: Essential scripts for data processing, ontology building, and visualization
- **Utils**: Utility functions for data manipulation and analysis
- **Tests**: Test scripts for validating functionality
- **Legacy**: Deprecated scripts kept for reference

## Directory Structure

```
astrocytes-ontology-research/
├── Single_Cell/
│   └── scripts/
│       ├── core/           # Essential single-cell scripts
│       ├── utils/          # Utility functions for single-cell data
│       ├── tests/          # Test scripts for single-cell functionality
│       └── legacy/         # Deprecated single-cell scripts
└── Spatial_Data/
    └── scripts/
        ├── core/           # Essential spatial data scripts
        ├── utils/          # Utility functions for spatial data
        ├── tests/          # Test scripts for spatial data functionality
        └── legacy/         # Deprecated spatial data scripts
```

## Single Cell Scripts

### Core

- **read_data.py**: Functions for reading and filtering single-cell data
- **single_cell_ontology_builder.py**: Tools for building ontologies from single-cell data
- **preprocess_single_cell.py**: Preprocessing functions for single-cell data

### Utils

- **compute_coexpression.py**: Functions for computing gene co-expression networks

### Tests

- **test_astrocyte_filter.py**: Tests for astrocyte filtering functionality
- **test_integration.py**: Tests for integration with other components
- **test_ontology_builder.py**: Tests for the ontology builder
- **test_read_data.py**: Tests for data reading functions

### Legacy

- **generate_h5ad.py**: Legacy script for generating H5AD files
- **main.py**: Legacy main script
- **run_all_tests.py**: Legacy test runner (replaced by the new test framework)

## Spatial Data Scripts

### Core

- **spatial_ontology_builder.py**: Tools for building ontologies from spatial data
- **ontology_integration.py**: Integration of spatial and single-cell ontologies
- **astrocyte_validation.py**: Validation of astrocyte cell type classifications
- **astrocyte_visualization.py**: Visualization of astrocyte data
- **astrocyte_integrated_visualization.py**: Integrated visualization across modalities
- **astrocyte_classifier.py**: Classification of astrocyte cell types

### Utils

- **query_enhanced_ontology.py**: Functions for querying enhanced spatial ontologies
- **query_spatial_ontology.py**: Basic query functions for spatial ontologies
- **enhance_spatial_ontology.py**: Functions for enhancing spatial ontologies with additional information

### Tests

- **test_astrocyte_validation.py**: Tests for astrocyte validation functionality

### Legacy

- **spatial_network_analysis.py**: Legacy network analysis functions
- **temporal_spatial_analysis.py**: Legacy temporal analysis functions
- **spatial_functional_analysis.py**: Legacy functional analysis functions
- **spatial_gene_relationships.py**: Legacy gene relationship analysis
- **spatial_pattern_analyzer.py**: Legacy pattern analysis functions
- **advanced_queries.py**: Legacy advanced query functions
- **spatial_query_builder.py**: Legacy query builder functions
- **spatial_gene_pathway_analyzer.py**: Legacy pathway analysis functions
- **spatial_gene_finder.py**: Legacy gene finder functions
- **spatial_ontology_advanced.py**: Legacy advanced ontology functions
- **astrocyte_spatial_dashboard.py**: Legacy dashboard (replaced by integrated visualization)
