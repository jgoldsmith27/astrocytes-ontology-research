# Brain Cell-Cell Interaction Analysis

This repository contains code and analysis results for identifying and characterizing potential sites of cell-cell communication in brain tissue using spatial transcriptomics and single-cell RNA-seq data.

## Project Overview

The goal of this project is to identify ligand-receptor interactions in spatial transcriptomic data and predict which cell types are involved in these interactions. We integrate three types of data:

1. **Spatial transcriptomics data**: Gene expression with spatial coordinates
2. **Single-cell RNA-seq data**: Cell type-specific gene expression profiles
3. **Ligand-receptor interaction database**: CellChat database of known interactions

## Key Findings

- Identified 1,195 potential interaction sites between ligands and receptors
- Discovered 9 distinct "hotspots" where multiple interactions co-occur
- Found dominant role of APOE-LRP1 signaling (93% of all detected interactions)
- Characterized cell type-specific communication patterns, with astrocytes and microglia being central players

For detailed findings, see:
- [Spatial Interaction Summary](spatial_interaction_summary.md)
- [Spatial Interaction Conclusions](spatial_interaction_conclusions.md)

## Repository Structure

```
/astrocytes-ontology-research/
├── README.md                         # This file
├── raw/                              # Raw data
│   ├── 1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv  # Spatial data
│   └── CTR081_Fron.h5ad              # Single-cell data
├── cell-chat-data/                   # Ligand-receptor database
│   ├── CellChatDB.human.rda
│   └── PPI.human.rda
├── scripts/                          # Analysis scripts
│   ├── explore_h5ad.py               # Single-cell data exploration
│   ├── explore_h5ad_simplified.py    # Simplified script for cell type analysis
│   ├── analyze_spatial_interactions.py # Spatial proximity analysis
│   └── integrated_spatial_analysis.py  # Combined analysis script
├── results/                          # Analysis outputs
│   ├── cell_type_counts.csv          # Cell type frequencies
│   ├── interaction_sites.csv         # Identified interaction sites
│   ├── interaction_cell_types.csv    # Predicted cell types for each interaction
│   ├── interaction_hotspots.csv      # Clusters of interaction sites
│   └── hotspot_cell_types.csv        # Cell types involved in each hotspot
└── visualizations/                   # Generated figures
    ├── spatial_interactions.png      # Map of interaction sites
    ├── interaction_heatmap.png       # Density of interactions
    ├── cell_interaction_patterns.png # Heatmap of cell-cell interactions
    ├── cell_interaction_network.png  # Circos plot of interactions
    └── interaction_hotspots.png      # Visualization of interaction clusters
```

## Methods

### 1. Spatial Interaction Analysis

We identified potential ligand-receptor interaction sites by analyzing the spatial proximity between ligand and receptor gene expressions in the spatial transcriptomic data. We defined interactions as pairs of ligand and receptor expressions that are within 50 spatial units of each other.

### 2. Cell Type Prediction

For each ligand-receptor pair, we predicted the most likely producing and receiving cell types based on the average expression of the ligand and receptor genes across different cell types in the single-cell RNA-seq data.

### 3. Interaction Hotspot Identification

We used DBSCAN clustering to identify regions in the tissue with high densities of interaction sites, which we call "hotspots". These hotspots likely represent functional microenvironments with specialized cellular activities.

## Key Interactions

| Interaction | Count | Ligand-Producing Cell | Receptor-Expressing Cell | Biological Significance |
|-------------|-------|---------------------------|---------------------------|-------------------------|
| APOE-LRP1 | 1,116 | Astrocytes1 | Mesenchymal cells | Lipid metabolism, amyloid clearance |
| BDNF-NTRK2 | 26 | Excitatory neurons | Astrocytes1 | Neurotrophic support |
| VEGFA-FLT1 | 18 | Mesenchymal cells | Endothelial cells | Angiogenesis |
| CX3CL1-CX3CR1 | 11 | Smooth muscle cells | Microglia | Neuroimmune communication |
| PDGFB-PDGFRB | 10 | Vasculature | Pericytes | Vascular stability |

## Usage

To reproduce the analysis:

1. Clone this repository
2. Set up Python environment with required dependencies:
   ```
   pip install scanpy pandas numpy matplotlib seaborn scipy scikit-learn
   ```
3. Run the analysis scripts in order:
   ```
   python explore_h5ad_simplified.py
   python analyze_spatial_interactions.py
   python integrated_spatial_analysis.py
   ```

## Dependencies

- Python 3.7+
- scanpy
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn

## Future Directions

1. Expand the analysis to include more ligand-receptor pairs
2. Apply the method to disease states (e.g., Alzheimer's disease)
3. Validate key interactions with spatial proteomics methods
4. Develop a more sophisticated model for predicting interaction probabilities

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or collaborations, please open an issue in this repository.

## Acknowledgments

This project uses the CellChat database for ligand-receptor interactions and builds upon methods from spatial transcriptomics and single-cell RNA-seq analysis. 