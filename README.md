# Astrocytes Ontology Research Project

## Overview

The Astrocytes Ontology Research Project aims to identify various cell types, including astrocytes, in spatial transcriptomics data using co-expression patterns derived from single-cell RNA sequencing. By leveraging both technologies, we can identify whole cells in spatial data where only individual gene expressions are directly measured.

## Datasets

The project works with two main types of data:

1. **Single-Cell RNA Sequencing Data** (scRNA-seq): Contains gene expression profiles for thousands of individual cells, with cell type annotations. This includes:
   - Astrocytes (general category): 1,287 cells
   - Astrocytes1 (subtype): 7,116 cells
   - Neurons (excitatory and inhibitory): 12,583 cells
   - Oligodendrocytes: 3,452 cells
   - Microglia: 1,029 cells
   - Other cell types present in the brain tissue

2. **Spatial Transcriptomics Data**: Contains spatial coordinates of individual gene expressions without cell type information.

## Project Goals

1. Analyze single-cell data to identify co-expression patterns characteristic of different cell types
2. Generate SPARQL inference rules based on these patterns
3. Apply these rules to spatial data to identify whole cells of various types
4. Visualize and validate the results

## Technical Approach

Our approach combines semantic web technologies with Python data analysis:

1. **Co-expression Analysis**: Identify genes that are reliably co-expressed in specific cell types in the single-cell data
2. **Rule Generation**: Convert co-expression patterns into SPARQL CONSTRUCT queries
3. **Cell Identification**: Apply the rules to identify cells in spatial data based on co-localized gene expressions
4. **Validation**: Filter and validate identified cells based on confidence scores and spatial distribution

For more details on why we chose this hybrid approach combining SPARQL/RDF with Python, see our [Architectural Design Document](README_ARCHITECTURE.md).

## Project Structure

```
astrocytes-ontology-research/
├── data/
│   ├── raw/                   # Original data files
│   └── processed/             # Processed data files and results
├── ontologies/
│   └── cell_type_ontology.ttl # Unified ontology for cell types
├── scripts/
│   ├── convert_spatial_to_turtle.py   # Convert spatial data to RDF
│   ├── generate_coexpression_rules.py # Generate rules from scRNA-seq
│   ├── generate_meta_rules.py         # Create meta-rules from rule associations
│   ├── identify_cells.py              # Apply rules to identify cells
│   ├── conflict_resolution.py         # Resolve conflicts between cells
│   └── visualize_cells.py             # Visualize identified cells
├── docs/
│   ├── README_DOCS.md                # Documentation guide and index
│   ├── co_expression_rules.md        # Documentation on rule types and usage
│   ├── methodology.md                # Detailed methodology description
│   ├── meta_rules.md                 # Documentation on higher-order patterns
│   ├── parameters.md                 # Reference for all configurable parameters
│   ├── data_conversion.md            # Documentation on data conversion process
│   ├── rule_application.md           # Documentation on rule application system
│   ├── conflict_resolution.md        # Documentation on conflict resolution system
│   ├── examples.md                   # Examples of rules and applications
│   ├── user_guide.md                 # Step-by-step guide to using the pipeline
│   └── todo_list.md                  # Planned improvements and feature roadmap
├── run_astrocyte_identification.sh   # Main pipeline script
└── README.md                         # This file
```

## Setup

### Prerequisites

- Python 3.8 or higher
- Required Python packages (install with `pip install -r requirements.txt`):
  - rdflib
  - pandas
  - numpy
  - scanpy
  - matplotlib
  - networkx
  - scipy

### Installation

1. Clone this repository:
   ```
   git clone https://github.com/username/astrocytes-ontology-research.git
   cd astrocytes-ontology-research
   ```

2. Create and activate a virtual environment:
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```

## Usage

### Running the Complete Pipeline

To run the complete pipeline with default settings:

```
./run_astrocyte_identification.sh
```

This will:
1. Convert spatial data to RDF format
2. Generate co-expression rules from single-cell data
3. Generate meta-rules by analyzing rule co-occurrence patterns
4. Apply the rules to identify cells in spatial data
5. Resolve conflicts between overlapping cells
6. Generate visualizations and reports

### Running Individual Steps

You can also run individual scripts with custom parameters. For example:

```
python scripts/generate_coexpression_rules.py \
  --input data/raw/CTR081_Fron.h5ad \
  --output-dir data/processed/rules \
  --min-expression 0.2 \
  --min-coexpression 0.6
```

See the [User Guide](docs/user_guide.md) for detailed instructions on running each component.

## Documentation

The `docs/` directory contains comprehensive documentation:

- [Documentation Guide](docs/README_DOCS.md): Overview and index of all documentation
- [Co-expression Rules Guide](docs/co_expression_rules.md): Explains the different types of rules
- [Methodology](docs/methodology.md): Detailed explanation of the scientific approach
- [Meta-Rules](docs/meta_rules.md): Documentation on higher-order co-expression patterns
- [Parameters Reference](docs/parameters.md): Complete list of configurable parameters
- [Data Conversion](docs/data_conversion.md): Documentation on the spatial data conversion process
- [Rule Application](docs/rule_application.md): Documentation on applying rules to identify cells
- [Conflict Resolution](docs/conflict_resolution.md): Documentation on resolving overlapping cells
- [Examples](docs/examples.md): Example rules and their application
- [User Guide](docs/user_guide.md): Step-by-step instructions for using the pipeline
- [Architectural Design](README_ARCHITECTURE.md): Explanation of our hybrid SPARQL/Python approach
- [TODO List](docs/todo_list.md): Planned improvements and features

## Future Improvements

Planned enhancements include:

- ✅ Implement higher-order co-expression patterns (meta-rules)
- ✅ Expand to identify multiple cell types beyond astrocytes
- Implement density-based clustering for improved cell boundary detection
- Enhance conflict resolution between overlapping cell identifications
- Optimize parameters for different cell types
- Add interactive visualization tools
- Explore machine learning approaches for rule generation

See the [TODO List](docs/todo_list.md) for more details.

## Results

After running the pipeline, results are available in `data/processed/`:

- `rules/`: Generated SPARQL rules
- `meta_rules/`: Higher-order co-expression patterns
- `results/identified_cells.csv`: Table of identified cells
- `results/cell_visualization.png`: Visualization of identified cells
- `results/conflict_resolution_report.md`: Detailed report on conflict resolution

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This research builds upon methods from single-cell transcriptomics and spatial genomics, particularly the work on identifying cell types from gene co-expression patterns. 