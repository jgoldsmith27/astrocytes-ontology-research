# Astrocytes Ontology Research

This project analyzes both spatial and single-cell transcriptomics data to study astrocyte biology, with a focus on cell-cell interactions and ligand-receptor signaling.

## Project Structure

```
.
├── data/                    # All data files
│   ├── spatial/            # Spatial transcriptomics data
│   │   ├── raw/           # Raw spatial data files
│   │   └── processed/     # Processed spatial data files
│   ├── single_cell/       # Single-cell transcriptomics data
│   │   ├── raw/           # Raw single-cell data files
│   │   └── processed/     # Processed single-cell data files
│   └── cellchat/          # CellChat database
│       ├── raw/           # Raw R data files
│       └── processed/     # Processed interaction data
├── src/                    # Source code
│   ├── common/            # Common utilities used by both analyses
│   │   ├── data_utils.py  # General data processing utilities
│   │   └── cellchat_utils.py  # CellChat database utilities
│   ├── spatial/           # Spatial transcriptomics analysis scripts
│   └── single_cell/       # Single-cell transcriptomics analysis scripts
└── notebooks/             # Jupyter notebooks for analysis and visualization
```

## Data Sources

### Spatial Transcriptomics
- File: `1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv`
- Location: `data/spatial/raw/`
- Description: Spatial transcriptomics data of cortical tissue

### Single-cell Transcriptomics
- File: `CTR081_Fron.h5ad`
- Location: `data/single_cell/raw/`
- Description: Single-cell RNA-seq data in AnnData format

### CellChat Database
- Files:
  - `CellChatDB.human.rda`: Main database of ligand-receptor interactions
  - `PPI.human.rda`: Protein-protein interaction data
- Location: `data/cellchat/raw/`
- Description: Database of known ligand-receptor interactions and protein-protein interactions

## Analysis Components

### Common Utilities
Located in `src/common/`
- `data_utils.py`: General data processing utilities for loading, saving, and analyzing data
- `cellchat_utils.py`: Utilities for working with the CellChat database, including:
  - Loading and converting R data files
  - Filtering interactions by cell type
  - Creating interaction networks
  - Processing protein-protein interactions

### Spatial Analysis
Located in `src/spatial/`
- Analysis of spatial transcriptomics data
- Identification of cell-cell interaction sites
- Spatial mapping of ligand-receptor pairs

### Single-cell Analysis
Located in `src/single_cell/`
- Processing of single-cell RNA-seq data
- Cell type identification and characterization
- Expression analysis of ligand-receptor pairs

## Setup and Usage

1. Create and activate virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Dependencies

- Python 3.7+
- pandas
- numpy
- scipy
- scanpy
- anndata
- rpy2
- matplotlib
- seaborn
- scikit-learn

## Future Directions

1. Implement spatial interaction analysis pipeline
2. Develop single-cell analysis pipeline
3. Create visualization tools for interaction networks
4. Integrate spatial and single-cell analyses

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or collaborations, please open an issue in this repository. 