# Spatial Ligand-Receptor Interaction Analysis

## Overview
This analysis focuses on identifying potential cell-cell interaction sites in spatial transcriptomic data by analyzing the proximity of ligand and receptor gene expressions. We've integrated data from:

1. **Spatial transcriptomic data** (1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv)
2. **Single-cell RNA-seq data** with cell type annotations (CTR081_Fron.h5ad)
3. **CellChat database** for ligand-receptor interaction pairs (CellChatDB.human.rda)

## Initial Findings

### Ligand-Receptor Pairs
We analyzed 15 well-known ligand-receptor pairs, including:
- TGFB1-TGFBR1 (TGF-beta signaling)
- APOE-LRP1 (Astrocyte-microglia signaling)
- CX3CL1-CX3CR1 (Neuronal-microglia signaling)
- VEGFA-FLT1/KDR (Vascular signaling)
- Several others related to inflammatory, neurotrophic, and growth factor signaling

### Interaction Sites
We identified 1,195 potential interaction sites where ligands and receptors are within 50 spatial units of each other. The most prominent interactions include:

1. APOE-LRP1: 1,116 potential interaction sites (93.4% of all interactions)
2. BDNF-NTRK2: 26 potential interaction sites
3. VEGFA-FLT1: 18 potential interaction sites
4. CX3CL1-CX3CR1: 11 potential interaction sites
5. PDGFB-PDGFRB: 10 potential interaction sites

### Cell Type Predictions
Using the single-cell data, we can predict which cell types are likely involved in these interactions:

1. **APOE-LRP1**: 
   - APOE is primarily expressed by Astrocytes
   - LRP1 is highly expressed in Microglia
   - This suggests significant Astrocyte-Microglia communication

2. **BDNF-NTRK2**:
   - BDNF is expressed by neurons
   - NTRK2 is found in both neurons and glia
   - Indicates neurotrophic signaling between neurons or from neurons to glia

3. **CX3CL1-CX3CR1**:
   - CX3CL1 (Fractalkine) is mainly produced by neurons
   - CX3CR1 is primarily expressed by microglia
   - Represents direct neuron-microglia communication

## Interaction Hotspots
We have identified regions in the tissue with high densities of potential interaction sites using clustering analysis. These hotspots likely represent anatomical regions where significant cell-cell communication is occurring.

## Next Steps
1. Further analyze the integrated single-cell and spatial data to refine cell type predictions
2. Identify specific molecular pathways enriched in interaction hotspots
3. Compare interaction patterns between different regions of the tissue
4. Use the identified interaction sites to guide further experimental investigations

## Technical Notes
- Interaction sites were identified using a spatial proximity threshold of 50 units
- Cell type predictions are based on average gene expression levels in the single-cell dataset
- Visualizations include spatial maps, heatmaps of interaction density, and circos plots of cell-cell communication patterns 