# Processed CellChat Data

This directory contains processed CSV files derived from the CellChat database R data files. The data comes from two source files:

## 1. CellChatDB.human.rda

This R data file contains four separate data frames, which have been exported as:

### interaction.csv (1,939 rows)
Contains ligand-receptor interactions with details including:
- Interaction names and alternative names
- Pathway names
- Ligands and receptors
- Agonists and antagonists
- Co-receptors (activating and inhibiting)
- Evidence sources (e.g., KEGG pathways, PMIDs)
- Annotations (e.g., "Secreted Signaling")

### complex.csv (157 rows)
Lists protein complexes that participate in cell-cell communication:
- Up to 4 subunits per complex (subunit_1 through subunit_4)
- Empty strings ("") indicate unused subunit slots

### cofactor.csv (31 rows)
Contains cofactors that regulate signaling pathways:
- Up to 16 cofactor components per row (cofactor1 through cofactor16)
- Empty strings ("") indicate unused cofactor slots
- Cofactors are molecules that assist in biological processes and help regulate signaling pathways

### gene_info.csv (41,787 rows)
A comprehensive gene database containing:
- Gene symbols
- Full gene names
- EntrezGene IDs
- Ensembl Gene IDs
- MGI (Mouse Genome Informatics) IDs
- Gene group names

## 2. PPI.human.rda

This R data file contains a sparse matrix of protein-protein interactions, exported as:

### ppi.csv (27,702 rows)
Network of protein-protein interactions represented as an edge list:
- protein1: First protein in the interaction
- protein2: Second protein in the interaction
- weight: Interaction weight (all weights are 1 in this dataset)
- Represents interactions between 4,815 different proteins

## Data Processing

These CSV files were created using the R script `src/common/inspect_rdata.R`. The script:
1. Loads the original R data files
2. Extracts each component
3. Converts them to CSV format for easier use in various analysis tools

## File Sizes
- interaction.csv: ~246KB
- complex.csv: ~3.5KB
- cofactor.csv: ~2.0KB
- gene_info.csv: ~3.9MB
- ppi.csv: ~485KB 