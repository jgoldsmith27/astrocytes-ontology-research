# Astrocyte Ontology Research

## Overview
This project focuses on integrating spatial transcriptomics and single-cell RNA sequencing data to create a comprehensive ontology for astrocyte research. The primary goal is to enhance our understanding of astrocyte heterogeneity and function through multi-modal data integration and validation.

The project has been refactored to specifically focus on astrocyte cell types, with streamlined components for data processing, integration, visualization, and validation.

## Project Context
Astrocytes are a heterogeneous population of glial cells with diverse morphologies, molecular signatures, and functions in the central nervous system. Understanding their diversity requires integrating multiple data modalities, particularly spatial transcriptomics and single-cell RNA sequencing. This project creates a semantic framework to bridge these data types, enabling comprehensive analysis of astrocyte subtypes.

## Components

### 1. Spatial Data Analysis
- **Location**: `/Spatial_Data`
- **Purpose**: Process and analyze spatial transcriptomics data for astrocytes
- **Features**:
  - Spatial ontology creation
  - Astrocyte type classification
  - Spatial gene expression analysis
  - Visualization of spatial patterns

### 2. Single-Cell Analysis
- **Location**: `/Single_Cell`
- **Purpose**: Process and analyze single-cell RNA sequencing data for astrocytes
- **Features**:
  - Single-cell ontology creation
  - Cell clustering and annotation
  - Gene expression profiling
  - Marker gene identification

### 3. Ontology Integration
- **Location**: `/output`
- **Purpose**: Semantically integrate spatial and single-cell data
- **Features**:
  - Bridge ontology creation
  - Cross-modality entity mapping
  - Inference rules for knowledge discovery
  - Advanced SPARQL queries

### 4. Integrated Visualization
- **Location**: `/Spatial_Data/scripts/core/astrocyte_integrated_visualization.py`
- **Purpose**: Visualize integrated data across modalities
- **Features**:
  - Interactive dashboards
  - Spatial distribution visualization
  - Gene expression heatmaps
  - Marker gene spatial mapping

### 5. Astrocyte Validation
- **Location**: `/Spatial_Data/scripts/core/astrocyte_validation.py`
- **Purpose**: Validate astrocyte cell type classifications
- **Features**:
  - Co-expression network analysis
  - Marker gene enrichment
  - Spatial validation visualization
  - Classification accuracy metrics

### 6. Test Suite
- **Location**: `/tests`
- **Purpose**: Verify functionality and prevent regressions
- **Features**:
  - Unit tests for individual components
  - Integration tests for cross-component functionality
  - Regression tests for previously fixed bugs
  - Organized test runner with filtering options

## Script Organization

The codebase has been reorganized to improve maintainability and focus on astrocyte-specific functionality:

### Directory Structure
```
astrocytes-ontology-research/
├── Single_Cell/
│   └── scripts/
│       ├── core/           # Essential single-cell scripts
│       ├── utils/          # Utility functions for single-cell data
│       └── legacy/         # Deprecated single-cell scripts
└── Spatial_Data/
    └── scripts/
        ├── core/           # Essential spatial data scripts
        ├── utils/          # Utility functions for spatial data
        └── legacy/         # Deprecated spatial data scripts
```

For more details on the script organization, see [SCRIPTS_ORGANIZATION.md](SCRIPTS_ORGANIZATION.md).

## Human Protein Atlas Collaboration

This research is being conducted in collaboration with the Human Protein Atlas (HPA) project at Karolinska Institute in Stockholm. The HPA is a Swedish-based program initiated in 2003 with the aim to map all human proteins in cells, tissues, and organs using a combination of various omics technologies, including antibody-based imaging, mass spectrometry, and transcriptomics.

### Contribution to HPA
This project specifically contributes to the HPA Brain Atlas section by:

1. **Enhancing astrocyte classification**: Providing more granular classification of astrocyte subtypes based on integrated spatial and molecular data
2. **Improving spatial mapping**: Creating semantic links between protein expression patterns and spatial locations in brain tissue
3. **Developing new visualization methods**: Offering novel ways to visualize protein expression in different astrocyte subtypes
4. **Creating reusable ontologies**: Building semantic frameworks that can be extended to other cell types in the brain

## Usage

### Prerequisites
- Python 3.8+
- Required packages listed in each component's requirements.txt

### Running the Pipeline

1. **Process spatial data**:
```bash
cd Spatial_Data
python scripts/core/spatial_ontology_builder.py --input data/spatial_data.csv --output output/enhanced_spatial_ontology.ttl
```

2. **Process single-cell data**:
```bash
cd Single_Cell
python scripts/core/single_cell_ontology_builder.py --input data/single_cell_data.csv --output output/single_cell_ontology.ttl
```

3. **Integrate ontologies**:
```bash
./run_ontology_integration.sh
```

4. **Run integrated visualization**:
```bash
cd Spatial_Data
python scripts/core/astrocyte_integrated_visualization.py --spatial ../output/integrated_ontology.ttl --single-cell ../output/integrated_ontology.ttl --interactive
```

5. **Run astrocyte validation**:
```bash
python Spatial_Data/scripts/core/astrocyte_validation.py --input output/integrated_astrocyte_ontology.ttl --output-dir output/astrocyte_validation
```

6. **Run all in one command**:
```bash
./run_astrocyte_integration.sh
```

### Running Tests

Run all tests:
```bash
python tests/run_tests.py
```

Run specific test types:
```bash
python tests/run_tests.py --type unit
python tests/run_tests.py --type integration
python tests/run_tests.py --type regression
```

For more details, see the [Test Suite README](tests/README.md).

## Future Directions

### Planned Enhancements
1. **Improved Gene Mapping**: Enhance the gene mapping between spatial and single-cell data to improve validation accuracy.
2. **Additional Marker Genes**: Incorporate more astrocyte-specific marker genes to better distinguish between subtypes.
3. **Enhanced Validation Algorithms**: Develop more sophisticated validation methods for astrocyte classification.
4. **Interactive Visualization Dashboard**: Create a comprehensive dashboard for exploring astrocyte heterogeneity.
5. **Integration with External Databases**: Connect with additional astrocyte-related databases and resources.

## Contributors
- Research team at Karolinska Institute
- Human Protein Atlas collaborators

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
This work is supported by the Human Protein Atlas project at Karolinska Institute, Stockholm, Sweden. 