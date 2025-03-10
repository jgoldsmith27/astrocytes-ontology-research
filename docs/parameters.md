# Parameters Reference Guide

## Introduction

This document provides a comprehensive reference for all configurable parameters in the Astrocytes Ontology Research Project pipeline. Understanding these parameters will help you customize the analysis for different datasets or research questions.

## Pipeline Parameters

These parameters can be adjusted in the main pipeline script `run_astrocyte_identification.sh`:

| Parameter | Description | Default | Impact |
|-----------|-------------|---------|--------|
| `SPATIAL_DATA` | Path to the spatial transcriptomics CSV file | `data/raw/spatial_data.csv` | Source data for spatial analysis |
| `SCRNA_DATA` | Path to the single-cell RNA-seq h5ad file | `data/raw/CTR081_Fron.h5ad` | Source data for co-expression patterns |
| `ONTOLOGY_FILE` | Path to the ontology TTL file | `ontologies/cell_type_ontology.ttl` | Defines the semantic structure |
| `SPATIAL_OUTPUT` | Output path for converted spatial data | `data/processed/spatial_data.ttl` | Intermediate RDF data |
| `RULES_DIR` | Directory for generated rules | `data/processed/rules` | Location of SPARQL rules |
| `RULES_FILE` | Path to the generated rules file | `$RULES_DIR/coexpression_rules.sparql` | Combined rules file |
| `CELLS_OUTPUT` | Output path for identified cells CSV | `data/processed/astrocyte_cells.csv` | Final results table |
| `VIZ_OUTPUT` | Output path for visualization | `data/processed/astrocyte_cells_distribution.png` | Visualization of results |

## Data Conversion Parameters

These parameters control how spatial data is converted to RDF format (`convert_spatial_to_turtle.py`):

| Parameter | Command Line Flag | Description | Default | Impact |
|-----------|------------------|-------------|---------|--------|
| Input file | `--input`, `-i` | Path to input CSV file | Required | Source spatial data |
| Output file | `--output`, `-o` | Path to output TTL file | Required | Where to save RDF data |
| Limit | `--limit`, `-l` | Limit number of rows to process | None (all rows) | For testing with smaller data |

### Example Usage

```bash
python scripts/conversion/convert_spatial_to_turtle.py \
    --input data/raw/spatial_data.csv \
    --output data/processed/spatial_data.ttl \
    --limit 100000  # Process only 100,000 rows (for testing)
```

## Co-expression Analysis Parameters

These parameters control the co-expression analysis and rule generation (`generate_coexpression_rules.py`):

| Parameter | Command Line Flag | Description | Default | Impact |
|-----------|------------------|-------------|---------|--------|
| Input file | `--input`, `-i` | Path to input h5ad file | Required | Source single-cell data |
| Output directory | `--output-dir`, `-o` | Directory for output files | Required | Where to save generated rules |
| Cell types | `--cell-types`, `-c` | Cell types to analyze | `Astrocytes Astrocytes1` | Target cell types |
| Minimum expression | `--min-expression`, `-e` | Minimum expression threshold | `0.2` | Threshold to consider a gene expressed |
| Minimum co-expression | `--min-coexpression`, `-p` | Minimum co-expression probability | `0.6` | Threshold for co-expression relationship |
| Top genes | `--top-genes`, `-g` | Number of top genes to analyze | `30` | How many genes to include per cell type |

### Impact of Parameters

#### Minimum Expression Threshold (0.2)

This threshold determines when a gene is considered "expressed" in a cell:
- **Lower threshold** (e.g., 0.1): More genes will be considered expressed, increasing sensitivity but potentially including noise
- **Higher threshold** (e.g., 0.3): Stricter definition of expression, reducing noise but might miss relevant genes
- **When to adjust**: Based on your data's normalization method; log-normalized data typically uses lower thresholds than CPM/TPM values

#### Minimum Co-expression Probability (0.6)

This threshold determines when two genes are considered co-expressed:
- **Lower threshold** (e.g., 0.4): More gene pairs will be considered co-expressed, generating more rules but with lower specificity
- **Higher threshold** (e.g., 0.8): Stricter definition of co-expression, resulting in fewer but more reliable rules
- **When to adjust**: If you need higher precision (increase) or higher recall (decrease)

#### Top Genes (30)

This parameter sets how many top genes to analyze for each cell type:
- **Lower value** (e.g., 10-20): Faster analysis focused on the most specific markers
- **Higher value** (e.g., 50-100): More comprehensive analysis that may find additional patterns
- **When to adjust**: Based on how many marker genes are expected for your cell type

### Example Usage

```bash
python scripts/analysis/generate_coexpression_rules.py \
    --input data/raw/CTR081_Fron.h5ad \
    --output-dir data/processed/rules \
    --cell-types Astrocytes Astrocytes1 \
    --min-expression 0.2 \
    --min-coexpression 0.6 \
    --top-genes 30
```

## Rule Application Parameters

These parameters control how rules are applied to identify cells (`apply_coexpression_rules.py`):

| Parameter | Command Line Flag | Description | Default | Impact |
|-----------|------------------|-------------|---------|--------|
| Spatial data | `--spatial`, `-s` | Path to spatial RDF data | Required | Input RDF data |
| Rules file | `--rules`, `-r` | Path to SPARQL rules file | Required | Input rules to apply |
| Output CSV | `--output-csv`, `-o` | Output path for identified cells | Required | Results table |
| Output visualization | `--output-viz`, `-v` | Output path for visualization | Required | Results visualization |
| Minimum confidence | `--min-confidence`, `-c` | Minimum confidence threshold | `0.7` | Threshold for cell filtering |
| Distance threshold | `--distance-threshold`, `-d` | Maximum distance for co-expression | `50` | Maximum distance in spatial units |

### Impact of Parameters

#### Minimum Confidence (0.7)

This threshold filters out low-confidence cell identifications:
- **Lower threshold** (e.g., 0.5): More cells will be included, increasing sensitivity but potentially including false positives
- **Higher threshold** (e.g., 0.8): Stricter filtering, resulting in fewer but more reliable cell identifications
- **When to adjust**: If you need higher precision (increase) or higher recall (decrease)

#### Distance Threshold (50)

This threshold determines how close gene expressions must be to be considered part of the same cell:
- **Lower threshold** (e.g., 30): Stricter spatial constraint, resulting in smaller identified cells
- **Higher threshold** (e.g., 70): More relaxed spatial constraint, allowing larger cell identifications
- **When to adjust**: Based on the expected size of the target cell type and the resolution of your spatial data

### Example Usage

```bash
python scripts/analysis/apply_coexpression_rules.py \
    --spatial data/processed/spatial_data.ttl \
    --rules data/processed/rules/coexpression_rules.sparql \
    --output-csv data/processed/astrocyte_cells.csv \
    --output-viz data/processed/astrocyte_cells_distribution.png \
    --min-confidence 0.7 \
    --distance-threshold 50
```

## Internal Parameters

These parameters are defined within the code and would require code modifications to change:

### In `generate_coexpression_rules.py`

| Parameter | Default | Description | Where Defined | Impact |
|-----------|---------|-------------|--------------|--------|
| Minimum clique size | `3` | Minimum size for gene cliques | Line ~170 | Affects complexity of clique rules |
| Number of cliques kept | `5` | Number of top cliques to keep | Line ~170 | Affects number of clique rules |
| Number of pairs kept | `10` | Number of top pairs to keep | Line ~180 | Affects number of pair rules |
| Minimum expression frequency | `0.1` | Minimum frequency for gene inclusion | Line ~70 | Filters out rarely expressed genes |
| Clique rule confidence | `0.9` | Confidence score for clique rules | Line ~170 | Affects weighting of clique rules |

### In `apply_coexpression_rules.py`

| Parameter | Default | Description | Where Defined | Impact |
|-----------|---------|-------------|--------------|--------|
| Nearby cell threshold | `3` | Maximum nearby cells for validation | Line ~160 | Filters out clustered false positives |
| Visualization point size | `20` | Size of points in visualization | Line ~240 | Affects visualization appearance |
| Visualization alpha | `0.6` | Transparency in visualization | Line ~240 | Affects visualization appearance |

## Parameter Combinations

Certain parameters work together and should be adjusted in tandem:

### For Higher Precision (Fewer False Positives)

```bash
--min-expression 0.25 \
--min-coexpression 0.7 \
--min-confidence 0.8 \
--distance-threshold 40
```

### For Higher Recall (Catch More Cells)

```bash
--min-expression 0.15 \
--min-coexpression 0.5 \
--min-confidence 0.6 \
--distance-threshold 60
```

### For Balanced Approach

```bash
--min-expression 0.2 \
--min-coexpression 0.6 \
--min-confidence 0.7 \
--distance-threshold 50
```

## Troubleshooting Parameter Issues

### Too Few Cells Identified

If your analysis identifies too few cells, try:
1. Lowering the `min-coexpression` parameter
2. Lowering the `min-confidence` parameter
3. Increasing the `distance-threshold` parameter
4. Increasing the `top-genes` parameter

### Too Many Cells (Likely False Positives)

If your analysis identifies too many cells or apparent false positives, try:
1. Increasing the `min-coexpression` parameter
2. Increasing the `min-confidence` parameter
3. Decreasing the `distance-threshold` parameter

### Cells of Wrong Size

If your identified cells are too large or too small, adjust:
1. The `distance-threshold` parameter (primary control of cell size)
2. Check if your spatial data units match your expectations

## Performance Considerations

Parameter choices can significantly affect computational performance:

| Parameter | Performance Impact |
|-----------|-------------------|
| Limit (conversion) | Linear impact on memory and time; use for testing |
| Top genes | Quadratic impact on rule generation time due to pairwise comparisons |
| Min co-expression | Higher values reduce computation by creating fewer rules |
| Distance threshold | Minimal impact on performance |

For large datasets, consider running with a subset first using the `--limit` parameter in the data conversion step. 