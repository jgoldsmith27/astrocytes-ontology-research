# Astrocyte Classification

This directory contains the output of the astrocyte classification process, which identifies and classifies different types of astrocytes in spatial transcriptomics data using machine learning techniques.

## Overview

The astrocyte classifier:
1. Extracts gene expression data from the spatial ontology
2. Creates a point-gene expression matrix
3. Calculates scores for each astrocyte type based on marker gene expression
4. Trains a Random Forest classifier to predict astrocyte types
5. Applies the classifier to all spatial points
6. Enhances the spatial ontology with cell type information
7. Visualizes the spatial distribution of astrocyte types

## Output Files

### Data Files

- **astrocyte_classified_points.csv**: Contains the classification results for all spatial points, including:
  - Spatial coordinates (x, y)
  - Gene expression values
  - Scores for each astrocyte type
  - Dominant astrocyte type (based on scores)
  - Predicted astrocyte type (from the classifier)
  - Prediction probabilities for each astrocyte type

### Model Files

- **astrocyte_classifier.joblib**: The trained Random Forest classifier model, which can be loaded and reused for classifying new data.

### Visualization Files

- **astrocyte_distribution.png**: Visualization of the spatial distribution of different astrocyte types, with points colored by type and sized by classification confidence.

## Astrocyte Types

The classifier identifies three main types of astrocytes:

1. **Protoplasmic Astrocytes** (Blue)
   - Found mainly in gray matter
   - Marker genes: GFAP, S100B, AQP4, GJA1, SLC1A2, SLC1A3
   - Characterized by: Highly branched processes, high expression of glutamate transporters

2. **Fibrous Astrocytes** (Red)
   - Found mainly in white matter
   - Marker genes: GFAP, VIM, ALDH1L1, CD44, CRYAB, HOPX
   - Characterized by: Long, unbranched processes, high expression of intermediate filaments

3. **Reactive Astrocytes** (Green)
   - Found in response to injury or disease
   - Marker genes: GFAP, VIM, SERPINA3, C3, CXCL10, STAT3, LCN2
   - Characterized by: Hypertrophy, upregulation of inflammatory genes

## Classification Methodology

The classification process involves:

1. **Feature Extraction**: Using marker gene expression as features
2. **Score Calculation**: Computing scores for each astrocyte type based on marker gene expression
3. **Initial Classification**: Determining the dominant astrocyte type based on the highest score
4. **Machine Learning**: Training a Random Forest classifier using the dominant types as labels
5. **Prediction**: Applying the trained classifier to all spatial points
6. **Confidence Estimation**: Calculating prediction probabilities for each astrocyte type

## Enhanced Ontology

The classification process enhances the spatial ontology with:

1. **Cell Type Classes**: Adding classes for different astrocyte types
2. **Cell Instances**: Creating instances of cells at each spatial point
3. **Type Information**: Linking cells to their predicted astrocyte type
4. **Confidence Scores**: Adding confidence scores for the classification
5. **Marker Genes**: Linking astrocyte types to their marker genes

## Usage

The enhanced ontology can be used for:

1. **Spatial Queries**: Finding regions with specific astrocyte types
2. **Correlation Analysis**: Correlating astrocyte types with other spatial features
3. **Integration**: Integrating with single-cell data for comprehensive analysis
4. **Visualization**: Creating advanced visualizations of astrocyte distribution

## References

1. Liddelow, S. A., & Barres, B. A. (2017). Reactive astrocytes: production, function, and therapeutic potential. Immunity, 46(6), 957-967.
2. Khakh, B. S., & Deneen, B. (2019). The emerging nature of astrocyte diversity. Annual review of neuroscience, 42, 187-207.
3. Bayraktar, O. A., et al. (2020). Astrocyte layers in the mammalian cerebral cortex revealed by a single-cell in situ transcriptomic map. Nature neuroscience, 23(4), 500-509. 