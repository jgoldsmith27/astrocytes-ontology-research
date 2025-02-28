# Astrocyte Spatial Dashboard

This dashboard provides interactive visualization and analysis of different types of astrocytes in spatial transcriptomics data. It leverages the spatial ontology to identify and characterize three main types of astrocytes: Protoplasmic, Fibrous, and Reactive.

## Features

- **Astrocyte Type Identification**: Automatically identifies three types of astrocytes based on marker gene expression patterns
- **Spatial Visualization**: Displays the spatial distribution of different astrocyte types
- **Marker Gene Analysis**: Visualizes the expression of specific marker genes across the spatial domain
- **Statistical Analysis**: Provides statistics on astrocyte distribution and marker gene expression
- **Data Export**: Exports visualizations and data for further analysis

## Astrocyte Types

The dashboard identifies three main types of astrocytes:

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

## Usage

### Running the Dashboard

```bash
python scripts/astrocyte_spatial_dashboard.py --ttl path/to/enhanced_spatial_ontology.ttl --output path/to/output/directory
```

### Command Line Arguments

- `--ttl`: Path to the TURTLE file containing the spatial ontology (required)
- `--output`: Path to the output directory (default: "../output/astrocyte_analysis")

### Interactive Features

1. **Load Ontology File**: Load a TURTLE file containing spatial transcriptomics data
2. **Show Astrocyte Types**: Visualize the spatial distribution of different astrocyte types
3. **Show Marker Expression**: Visualize the expression of specific marker genes
4. **Export Current View**: Save the current visualization as an image file
5. **Export All Data**: Export all data as CSV files for further analysis

## Output Files

The dashboard generates the following output files:

- **astrocyte_spatial_data.csv**: Contains spatial coordinates and astrocyte type information for each grid cell
- **marker_gene_statistics.csv**: Contains statistics on marker gene expression for each astrocyte type
- **Visualization images**: PNG, JPG, or PDF files of the visualizations

## Integration with Single-Cell Data

This dashboard can be integrated with single-cell RNA sequencing data to provide a more comprehensive view of astrocyte heterogeneity. The integration enables:

1. Validation of astrocyte types identified in spatial data
2. Discovery of additional marker genes for each astrocyte type
3. Correlation of spatial patterns with transcriptional profiles

## Requirements

- Python 3.6+
- pandas
- numpy
- matplotlib
- seaborn
- rdflib
- tkinter
- scikit-learn

## Future Enhancements

- Integration with temporal data to analyze astrocyte dynamics
- Addition of more astrocyte subtypes as they are discovered
- Enhanced visualization with 3D rendering
- Machine learning-based classification of astrocyte types
- Integration with functional data to correlate astrocyte types with function

## References

1. Liddelow, S. A., & Barres, B. A. (2017). Reactive astrocytes: production, function, and therapeutic potential. Immunity, 46(6), 957-967.
2. Khakh, B. S., & Deneen, B. (2019). The emerging nature of astrocyte diversity. Annual review of neuroscience, 42, 187-207.
3. Bayraktar, O. A., et al. (2020). Astrocyte layers in the mammalian cerebral cortex revealed by a single-cell in situ transcriptomic map. Nature neuroscience, 23(4), 500-509. 