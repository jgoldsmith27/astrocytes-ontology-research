# Ligand-Receptor Spatial Analysis

This tool analyzes the spatial relationships between ligand-receptor pairs in spatial transcriptomics data.

## Overview

The `analyze_ligand_receptor_spatial.py` script processes spatial transcriptomics data to:

1. Identify the spatial locations of ligands and their corresponding receptors
2. Calculate a distance/density score between them
3. Visualize their spatial relationships
4. Generate summary statistics and plots

## Spatial Data Context

Understanding spatial distances in the data:

- **Scale**: 1 pixel = 0.5 micrometers
- **Cell size**: Typically 10-50 micrometers (20-100 pixels)
- **RNA location**: RNA transcripts may not be at the exact location of the proteins they code for
- **Interaction distances**: Biologically meaningful interactions typically occur within 20-50 micrometers (40-100 pixels)

## How It Works

### Distance/Density Score

The script calculates a custom score for each ligand-receptor pair that takes into account:

1. **Proximity**: The average distance between a receptor and its nearest ligands
2. **Density**: How clustered the ligands are around each receptor
3. **Relevance**: Whether ligands are within biologically meaningful distance ranges
4. **Abundance**: The number of ligands in proximity to receptors

The formula used is:
```
score = (avg_distance * (1 - density_factor * 0.5) * (1 - count_factor * 0.3)) * distance_penalty
```

Where:
- `avg_distance`: Mean distance to the k-nearest ligands (within threshold)
- `density_factor`: Derived from the coefficient of variation of distances (indicates clustering)
- `count_factor`: Reflects how many ligands are near each receptor
- `distance_penalty`: Applies penalties for distances outside optimal ranges

A **lower score** indicates:
- Shorter average distances between ligands and receptors
- Higher density of ligands around receptors
- Appropriate biological distance ranges
- Higher number of potentially interacting molecules

### Complex Receptors

For complex receptors (made up of multiple genes/components):

1. Identifies locations where all components are within 50 micrometers (100 pixels) of each other
2. Selects the closest pairs when multiple options exist
3. Uses the midpoint between components as the receptor location for visualization
4. Falls back to the first component only if paired locations cannot be found

This approach addresses the biological reality that complex receptors require all components to be in close proximity for functional signaling.

## Usage

```bash
python analyze_ligand_receptor_spatial.py [options]
```

### Options

- `--lr-pairs`: Path to ligand-receptor pairs CSV file (default: 'results/ligand_receptor_pairs.csv')
- `--spatial-data`: Path to spatial data CSV file
- `--num-pairs`: Number of top ligand-receptor pairs to analyze (default: 15)
- `--k-nearest`: Number of nearest ligands to consider for each receptor (default: 50)
- `--output-dir`: Directory to save results (default: 'results/spatial_analysis')

## Output

The script generates:

1. Individual visualization plots for each ligand-receptor pair
2. Density maps/heatmaps showing interaction hotspots
3. A summary CSV file with score results
4. A summary bar plot comparing scores across all analyzed pairs

## Handling Multi-Ligand and Orphaned Receptors

The current implementation:

1. **Orphaned receptors**: Receptors without nearby ligands receive a higher (worse) score due to the distance penalty
2. **Multiple ligands**: Each ligand-receptor pair is analyzed separately, allowing for comparison of scores to determine the most likely interactions
3. **Competitive binding**: Future versions will incorporate competitive binding models where receptors can interact with multiple ligands

## Research Questions

This tool can help explore:

1. **Spatial organization**: How are ligands and receptors spatially organized within tissues?
2. **Signaling hotspots**: Are there regions where specific signaling pathways are concentrated?
3. **Cell-cell communication**: How do spatial arrangements facilitate or inhibit cell-cell communication?
4. **Pathway dominance**: Which ligand-receptor pairs show strongest spatial correlation?
5. **Disease mechanisms**: Do spatial arrangements of signaling molecules differ in disease states?
6. **Regional specialization**: Are certain signaling patterns associated with specific tissue regions?
7. **Complex receptor assembly**: For multi-component receptors, are all parts consistently co-localized?
8. **Competitive interactions**: When multiple ligands can bind the same receptor, which spatial arrangements dominate?

## Future Improvements

1. Advanced spatial statistics using established methods like Ripley's K-function
2. DBSCAN clustering for identifying interaction hotspots
3. Integration with cell type information
4. Statistical testing for significance
5. Interactive visualizations
6. Comparative analysis across samples
7. Competitive binding models for receptors with multiple potential ligands 