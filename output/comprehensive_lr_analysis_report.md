# Comprehensive Ligand-Receptor Analysis Report

Analysis performed on 2025-03-30 16:18:43

## Overview

This report summarizes the analysis of ligand-receptor pairs from the CellChat database in relation to single-cell expression data.

### Gene Coverage

- Total ligands in CellChat: 547
- Total receptors in CellChat: 468
- Ligands found in expression data: 7
- Receptors found in expression data: 8
- Missing ligands: 540
- Missing receptors: 460

### Expression Summary

- Total cell types analyzed: 17
- Expression threshold: 0.1
- Ligand-receptor pairs enabling communication: 8

## Cell Type Expression Profiles

| Cell Type | All Ligands | All Receptors | CellChat Ligands | CellChat Receptors |
|-----------|------------|--------------|------------------|--------------------|
| Astrocytes1 | 3 | 5 | 2 | 5 |
| Oligodendrocytes | 0 | 1 | 0 | 1 |
| Pericytes | 4 | 6 | 3 | 6 |
| Microglia | 3 | 4 | 2 | 4 |
| Astrocytes2 | 4 | 5 | 3 | 5 |
| Endothelial | 3 | 3 | 2 | 3 |
| Neurons | 2 | 3 | 1 | 3 |
| CAMs | 3 | 3 | 2 | 3 |
| Lymphocytes | 2 | 2 | 1 | 2 |
| Mesenchymal | 4 | 3 | 3 | 3 |
| Fibroblasts | 3 | 6 | 2 | 6 |
| SMCs | 3 | 4 | 2 | 4 |
| OPCs | 2 | 5 | 1 | 5 |
| Astrocytes | 3 | 4 | 2 | 4 |
| neurons_ex | 1 | 2 | 1 | 2 |
| Vasculature | 3 | 6 | 2 | 6 |
| neurons_in | 1 | 3 | 1 | 3 |

## Communication Analysis

### Top Communication Pathways

| Pathway | LR Pairs | Total Connections | Avg Connections/Pair | Autocrine | Paracrine |
|---------|----------|-------------------|----------------------|-----------|------------|
| VEGF | 4 | 84 | 21.00 | 4 | 80 |
| PDGF | 1 | 32 | 32.00 | 1 | 31 |
| EGF | 2 | 22 | 11.00 | 2 | 20 |
| CX3C | 1 | 2 | 2.00 | 0 | 2 |

### Top Communication Ligand-Receptor Pairs

| Ligand | Receptor | Pathway | Total Connections | Autocrine | Paracrine |
|--------|----------|---------|-------------------|-----------|------------|
| PDGFB | PDGFRB | PDGF | 32 | 1 | 31 |
| VEGFA | FLT1 | VEGF | 21 | 1 | 20 |
| VEGFA | KDR | VEGF | 21 | 1 | 20 |
| VEGFA | FLT1 | VEGF | 21 | 1 | 20 |
| VEGFA | KDR | VEGF | 21 | 1 | 20 |
| EGF | EGFR | EGF | 11 | 1 | 10 |
| EGF | EGFR | EGF | 11 | 1 | 10 |
| CX3CL1 | CX3CR1 | CX3C | 2 | 0 | 2 |

## Key Findings

### Top Autocrine Signaling Pairs

- **PDGFB-PDGFRB** (PDGF): Used in 1 cell types
- **VEGFA-FLT1** (VEGF): Used in 1 cell types
- **VEGFA-KDR** (VEGF): Used in 1 cell types

### Top Paracrine Signaling Pairs

- **PDGFB-PDGFRB** (PDGF): Enables 31 cell-cell connections
- **VEGFA-FLT1** (VEGF): Enables 20 cell-cell connections
- **VEGFA-KDR** (VEGF): Enables 20 cell-cell connections

## Recommendations for Further Analysis

1. **Validation**: Validate top ligand-receptor pairs with experimental data
2. **Functional Analysis**: Investigate biological functions of top pathways
3. **Network Analysis**: Perform advanced network analysis of the communication patterns
4. **Spatial Context**: Consider spatial relationships between communicating cell types
