# Spatial and Single-Cell Integration for Astrocyte Analysis

This directory contains scripts for integrating spatial transcriptomics data with single-cell RNA sequencing data to provide a comprehensive analysis of astrocyte heterogeneity.

## Key Scripts

### `integrate_spatial_single_cell.py`

This script integrates spatial transcriptomics data with single-cell RNA sequencing data to identify and characterize different types of astrocytes.

#### Features

- Loads spatial data from TURTLE ontology files and single-cell data from H5AD files
- Identifies three types of astrocytes (Protoplasmic, Fibrous, and Reactive) in both datasets
- Correlates gene expression patterns between spatial and single-cell data
- Generates comprehensive visualizations of the integrated data
- Exports integrated data for further analysis

#### Usage

```bash
python integrate_spatial_single_cell.py --spatial path/to/enhanced_spatial_ontology.ttl --single-cell path/to/single_cell_data_cleaned.h5ad --output path/to/output/directory
```

#### Command Line Arguments

- `--spatial`: Path to the TURTLE file containing the spatial ontology (required)
- `--single-cell`: Path to the H5AD file containing the single-cell data (required)
- `--output`: Path to the output directory (default: "../output/integrated_analysis")

#### Output Files

- `integrated_gene_data.csv`: Contains correlation data between spatial and single-cell expression for each gene
- `integrated_visualizations.pdf`: PDF containing multiple visualizations of the integrated data

## Astrocyte Types

The scripts identify and analyze three main types of astrocytes:

1. **Protoplasmic Astrocytes**
   - Found mainly in gray matter
   - Marker genes: GFAP, S100B, AQP4, GJA1, SLC1A2, SLC1A3
   - Characterized by: Highly branched processes, high expression of glutamate transporters

2. **Fibrous Astrocytes**
   - Found mainly in white matter
   - Marker genes: GFAP, VIM, ALDH1L1, CD44, CRYAB, HOPX
   - Characterized by: Long, unbranched processes, high expression of intermediate filaments

3. **Reactive Astrocytes**
   - Found in response to injury or disease
   - Marker genes: GFAP, VIM, SERPINA3, C3, CXCL10, STAT3, LCN2
   - Characterized by: Hypertrophy, upregulation of inflammatory genes

## Integration Methodology

The integration approach involves:

1. **Data Loading**: Loading spatial ontology data and single-cell RNA-seq data
2. **Astrocyte Identification**: Identifying astrocyte types in both datasets based on marker gene expression
3. **Common Gene Set**: Finding genes present in both datasets
4. **Correlation Analysis**: Calculating correlation between gene expression patterns in spatial and single-cell data
5. **Visualization**: Generating visualizations to compare astrocyte types across datasets

## Requirements

- Python 3.6+
- pandas
- numpy
- matplotlib
- seaborn
- rdflib
- scanpy
- anndata
- scikit-learn
- scipy

## Future Enhancements

- Integration with temporal data to analyze astrocyte dynamics
- Enhanced clustering to identify additional astrocyte subtypes
- Machine learning-based prediction of spatial location from transcriptional profile
- Integration with functional data to correlate astrocyte types with function

## References

1. Liddelow, S. A., & Barres, B. A. (2017). Reactive astrocytes: production, function, and therapeutic potential. Immunity, 46(6), 957-967.
2. Khakh, B. S., & Deneen, B. (2019). The emerging nature of astrocyte diversity. Annual review of neuroscience, 42, 187-207.
3. Bayraktar, O. A., et al. (2020). Astrocyte layers in the mammalian cerebral cortex revealed by a single-cell in situ transcriptomic map. Nature neuroscience, 23(4), 500-509. 