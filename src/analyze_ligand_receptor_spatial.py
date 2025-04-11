#!/usr/bin/env python3
"""
Analyze ligand-receptor pairs in spatial data to calculate distance/density scores
and visualize their interactions.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.stats import gaussian_kde
from sklearn.neighbors import NearestNeighbors
import seaborn as sns
from pathlib import Path
import argparse
from tqdm import tqdm
import logging
import networkx as nx
from scipy import stats

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_ligand_receptor_pairs(file_path, num_pairs=15):
    """Load the top N ligand-receptor pairs from CSV file."""
    logger.info(f"Loading top {num_pairs} ligand-receptor pairs from {file_path}")
    lr_pairs = pd.read_csv(file_path)
    return lr_pairs.head(num_pairs)

def load_spatial_data(file_path):
    """Load spatial data from CSV file."""
    logger.info(f"Loading spatial data from {file_path}")
    return pd.read_csv(file_path)

def get_gene_spatial_data(spatial_data, gene_name):
    """Extract spatial coordinates for a specific gene."""
    gene_data = spatial_data[spatial_data['geneID'] == gene_name]
    if len(gene_data) == 0:
        logger.warning(f"No spatial data found for gene: {gene_name}")
        return None
    
    # Return x, y coordinates
    return gene_data[['x', 'y']].values

def calculate_distance_density_score(ligand_coords, receptor_coords, k=50, ligand_name=None, receptor_name=None):
    """
    Calculate an interaction density score based on the product of 
    Ligand KDE and Receptor KDE evaluated at each receptor's location.

    Higher score = higher co-localization density of ligands and receptors.
    
    Args:
        ligand_coords (np.ndarray): Array of (x,y) coordinates for ligand molecules.
        receptor_coords (np.ndarray): Array of (x,y) coordinates for receptor locations (centroids for complex).
        k (int, optional): Not used in this scoring method, kept for compatibility.
        ligand_name (str, optional): Name of the ligand gene (for logging).
        receptor_name (str, optional): Name of the receptor gene (for logging).
    
    Returns:
        dict: Dictionary containing average interaction score and individual 
              receptor scores, or None if calculation fails.
    """
    if ligand_coords is None or receptor_coords is None:
        logger.error(f"Missing ligand or receptor coordinates for {ligand_name}-{receptor_name}.")
        return None
    
    # KDE requires at least 2 points for meaningful calculation
    if len(ligand_coords) < 2:
        logger.warning(f"Cannot calculate KDE score for {ligand_name}-{receptor_name}: requires >= 2 ligand points, found {len(ligand_coords)}.")
        return {"average_score": 0.0, "receptor_scores": {}} 
    if len(receptor_coords) < 2:
        logger.warning(f"Cannot calculate KDE score for {ligand_name}-{receptor_name}: requires >= 2 receptor points/centroids, found {len(receptor_coords)}.")
        return {"average_score": 0.0, "receptor_scores": {}} 
        
    if len(receptor_coords) == 0:
        logger.warning(f"Cannot calculate KDE score for {ligand_name}-{receptor_name}: no receptor points found.")
        return {"average_score": 0.0, "receptor_scores": {}}

    try:
        # 1. Build KDE model from ligand coordinates
        ligand_kde = gaussian_kde(ligand_coords.T)
        
        # 2. Build KDE model from receptor coordinates (centroids)
        receptor_kde = gaussian_kde(receptor_coords.T)
        
        # 3. Evaluate both KDEs at each receptor coordinate
        # receptor_coords.T gives (2, M) where M is number of receptors
        ligand_density_at_receptors = ligand_kde(receptor_coords.T)
        receptor_density_at_receptors = receptor_kde(receptor_coords.T)
        
        # Ensure densities are non-negative
        ligand_density_at_receptors[ligand_density_at_receptors < 0] = 0
        receptor_density_at_receptors[receptor_density_at_receptors < 0] = 0
        
        # 4. Calculate the interaction score (product of densities) for each receptor location
        interaction_scores = ligand_density_at_receptors * receptor_density_at_receptors
        
        # 5. Calculate the average interaction score
        # Use nanmean in case any calculation failed for a point, although KDE evaluation is usually robust
        average_score = np.nanmean(interaction_scores) 
        
        # Store individual scores for potential debugging or detailed analysis
        receptor_scores_detail = { 
            f"receptor_{i}": {
                "interaction_density_score": interaction_scores[i],
                "ligand_density": ligand_density_at_receptors[i],
                "receptor_density": receptor_density_at_receptors[i]
                } 
            for i in range(len(receptor_coords))
        }
        
        # Scale average score for potentially better range/visualization (optional)
        # average_score *= 1e6 # Example scaling
        
        logger.info(f"Calculated Interaction Density score for {ligand_name}-{receptor_name}. Avg score: {average_score:.4g}")
        
        return {
            "average_score": average_score,
            "receptor_scores": receptor_scores_detail # Contains interaction_density_score per receptor
        }

    except Exception as e:
        logger.error(f"Failed to calculate Interaction Density score for {ligand_name}-{receptor_name}: {e}")
        return {"average_score": np.nan, "receptor_scores": {}}

def visualize_interaction(ligand_coords, receptor_clusters, score_data, ligand_name, receptor_name, is_complex, output_dir):
    """
    Visualize the spatial interaction density between ligand and receptor.
    Shows ligand points, receptor locations (centroids for complex), interaction density, and contours.
    Saves the plot to the specified output directory.
    """
    # Extract receptor centroids for overall calculations
    receptor_coords = np.array([cluster['centroid'] for cluster in receptor_clusters]) if receptor_clusters else np.array([])
    
    # Ensure coordinates are valid numpy arrays for calculation
    if ligand_coords is None or not isinstance(ligand_coords, np.ndarray) or ligand_coords.ndim != 2 or ligand_coords.shape[1] != 2:
        logger.warning(f"Invalid or missing ligand coordinates for {ligand_name}-{receptor_name}.")
        return
    if receptor_coords is None or not isinstance(receptor_coords, np.ndarray) or receptor_coords.ndim != 2 or receptor_coords.shape[1] != 2:
         logger.warning(f"Invalid or missing receptor coordinates for {ligand_name}-{receptor_name}.")
         return   
    if len(ligand_coords) == 0 or len(receptor_coords) == 0 or score_data is None or score_data.get('average_score') is None:
        logger.warning(f"Cannot visualize density for {ligand_name}-{receptor_name}: missing data, coordinates, or score.")
        return
    
    avg_score_value = score_data['average_score']
    
    # Get normalized score if available
    normalized_score = score_data.get('normalized_score', None)
    score_display = f'Avg. Interaction Score: {avg_score_value:.4g}'
    if normalized_score is not None:
        score_display += f' (Global Norm: {normalized_score:.3f})'

    # --- REMOVED plot type logic --- 
    title_prefix = "Interaction Density:" # Simplified prefix

    plt.figure(figsize=(12, 10))
    ax_density = plt.gca()
    
    # Define grid for density estimation
    all_points = np.vstack([ligand_coords, receptor_coords])
    x_min, x_max = np.min(all_points[:, 0]), np.max(all_points[:, 0])
    y_min, y_max = np.min(all_points[:, 1]), np.max(all_points[:, 1])
    
    # Add padding
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_padding = x_range * 0.05 if x_range > 0 else 10
    y_padding = y_range * 0.05 if y_range > 0 else 10
    x_min -= x_padding
    x_max += x_padding
    y_min -= y_padding
    y_max += y_padding
    
    # Create meshgrid
    grid_res = 150j # Slightly lower resolution grid for performance
    x_grid, y_grid = np.mgrid[x_min:x_max:grid_res, y_min:y_max:grid_res]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    
    ligand_density = np.zeros(x_grid.shape)
    receptor_density = np.zeros(x_grid.shape)
    interaction_density = np.zeros(x_grid.shape)
    plotted_interaction_density = False

    # --- Calculate Densities (KDE) --- 
    can_calc_ligand_density = False
    if len(ligand_coords) >= 2:
        try:
            ligand_kernel = gaussian_kde(ligand_coords.T)
            ligand_density = np.reshape(ligand_kernel(positions), x_grid.shape)
            can_calc_ligand_density = True
        except Exception as e:
            logger.warning(f"Could not compute ligand density for {ligand_name}: {e}")

    can_calc_receptor_density = False
    if len(receptor_coords) >= 2:
        try:
            receptor_kernel = gaussian_kde(receptor_coords.T)
            receptor_density = np.reshape(receptor_kernel(positions), x_grid.shape)
            can_calc_receptor_density = True
        except Exception as e:
            logger.warning(f"Could not compute receptor density for {receptor_name} centroids: {e}")

    # --- Plot Interaction Density --- 
    if can_calc_ligand_density and can_calc_receptor_density:
        try:
            interaction_density = ligand_density * receptor_density
            max_interaction = np.max(interaction_density)
            
            if max_interaction > 0: 
                interaction_density_normalized = interaction_density / max_interaction 
                im_interaction = ax_density.imshow(np.rot90(interaction_density_normalized), cmap='Greens', 
                                              extent=[x_min, x_max, y_min, y_max], 
                                              alpha=0.7, aspect='auto', vmin=0.01, vmax=1.0, zorder=1)
                plotted_interaction_density = True
                
                density_for_contours = interaction_density_normalized[interaction_density_normalized > 0.01]
                if len(density_for_contours) > 0:
                    density_thresholds = np.percentile(density_for_contours, [60, 80, 95]) 
                    contour_levels = sorted(list(set(d for d in density_thresholds if d > 1e-3)))
                    
                    if len(contour_levels) > 0:
                        contour = ax_density.contour(x_grid, y_grid, interaction_density_normalized, 
                                                levels=contour_levels, colors='purple', alpha=0.8, linewidths=0.8, zorder=2)
                        ax_density.clabel(contour, inline=True, fontsize=6, fmt='%.2f')
                else:
                     logger.info(f"Normalized interaction density too low for contours for {ligand_name}-{receptor_name}")
            else:
                 logger.info(f"Interaction density is zero or negative for {ligand_name}-{receptor_name}, cannot plot heatmap/contours.")
        except Exception as e:
            logger.warning(f"Could not compute/plot interaction density for {ligand_name}-{receptor_name}: {e}")
            
    # --- Plot Raw Points --- 
    point_size = 6
    ligand_alpha = 0.7
    receptor_alpha = 0.8
    
    # Ligands (Blue dots)
    ax_density.scatter(ligand_coords[:, 0], ligand_coords[:, 1], s=point_size, color='blue', alpha=ligand_alpha, label=f'{ligand_name} (Ligand)', marker='.', zorder=3)

    # Receptors (Red dots - centroids for complex)
    plotted_receptor_label = False
    # Limit plotted receptors for clarity if there are too many
    max_receptors_to_plot = 5000 
    indices_to_plot = np.random.choice(len(receptor_clusters), size=min(len(receptor_clusters), max_receptors_to_plot), replace=False) if len(receptor_clusters) > max_receptors_to_plot else range(len(receptor_clusters))

    for i in indices_to_plot:
        cluster = receptor_clusters[i]
        centroid = cluster['centroid']
        receptor_label = f'{receptor_name} (Receptor/Centroid)' if not plotted_receptor_label else ""
        ax_density.scatter(centroid[0], centroid[1], s=point_size+2, color='red', marker='.', label=receptor_label, zorder=4)
        plotted_receptor_label = True

    # --- Final plot settings --- 
    num_ligands_disp = len(ligand_coords)
    num_receptors_disp = len(receptor_coords)
    
    # --- UPDATED TITLE --- 
    plt.title(f'{title_prefix} {ligand_name} ({num_ligands_disp}) - {receptor_name} ({num_receptors_disp})\n' # Title simplified
              f'{score_display} (Higher is Denser Co-localization)')
    ax_density.set_xlabel('X Coordinate (1 pixel = 0.5 micrometers)')
    ax_density.set_ylabel('Y Coordinate (1 pixel = 0.5 micrometers)')
    ax_density.set_xlim(x_min, x_max)
    ax_density.set_ylim(y_min, y_max)
    ax_density.set_aspect('equal', adjustable='box') # Ensure aspect ratio is equal
    ax_density.grid(True, linestyle=':', alpha=0.3)
    
    # Add combined legend
    handles, labels = ax_density.get_legend_handles_labels()
    if plotted_interaction_density: 
        handles.append(plt.Rectangle((0,0),1,1,fc="lightgreen", alpha=0.7))
        labels.append('Interaction Density (Ligand x Receptor KDE)')
        handles.append(plt.Line2D([0], [0], color='purple', lw=1)) 
        labels.append('Interaction Contours')

    # Reduce legend size if too many items
    legend_fontsize = 'small' if len(handles) < 8 else 'x-small'
    # Place legend outside plot area
    ax_density.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.02, 1), 
                      fontsize=legend_fontsize, framealpha=0.8, borderaxespad=0.)

    # Add scale bar
    scale_bar_length = 100
    scale_x_pos = x_min + (x_max - x_min) * 0.05
    scale_y_pos = y_min + (y_max - y_min) * 0.05
    ax_density.plot([scale_x_pos, scale_x_pos + scale_bar_length], 
             [scale_y_pos, scale_y_pos], 'k-', linewidth=3)
    ax_density.text(scale_x_pos + scale_bar_length/2, scale_y_pos - (y_max - y_min) * 0.03, 
             '50 μm', horizontalalignment='center', fontsize=10, color='black')
    
    # Save density figure
    density_output_path = Path(output_dir) / f"{ligand_name}_{receptor_name}_density_map.png"
    plt.savefig(density_output_path, dpi=300, bbox_inches='tight')
    plt.close()
    # --- UPDATED LOG MESSAGE --- 
    logger.info(f"Density map saved to {density_output_path}") # Simplified log

def handle_complex_receptors(receptor_name, is_complex, valid_receptors, spatial_data):
    """
    Handle complex receptors (made up of multiple genes).

    For complex receptors, finds locations where all components are within a threshold
    and returns the centroid and contributing points for each valid cluster.

    Notes on distance units:
    - 1 pixel = 0.5 micrometers in the spatial data
    - For complex receptors, components should be within 50 micrometers (100 pixels)
    """
    if not is_complex:
        # Simple case: return coords directly, format similar to complex for consistency
        coords = get_gene_spatial_data(spatial_data, receptor_name)
        if coords is None:
            return None, 0, {} # No coords, 0 centroids, empty component counts
        # Treat each point as its own 'centroid' for simple receptors
        clusters = [{'centroid': c, 'component_points': [c], 'component_names': [receptor_name]} for c in coords]
        return clusters, len(coords), {receptor_name: len(coords)}

    # --- Complex Receptor Logic ---
    components = valid_receptors.split('_')
    logger.info(f"Processing complex receptor {receptor_name} with components: {components}")

    # Get coordinates for each component
    component_data = {}
    min_component_count = float('inf')
    component_counts = {}
    all_components_found = True
    for component in components:
        coords = get_gene_spatial_data(spatial_data, component)
        if coords is None or len(coords) == 0:
            logger.warning(f"No data found for receptor component: {component} in {receptor_name}")
            all_components_found = False
            # Store 0 count but continue processing if possible
            component_counts[component] = 0
            component_data[component] = {'coords': np.array([])} 
            # If *any* essential component is missing entirely, we cannot form the complex
            # Depending on biological interpretation, could allow complex formation if *some* components found
            # For now, let's require all listed components to be present at least once in the dataset
            # return None, 0, component_counts # Cannot form clusters if a component is entirely missing
        else:
            count = len(coords)
            component_counts[component] = count
            component_data[component] = {'coords': coords}
            min_component_count = min(min_component_count, count)
            
    # If any component was completely missing, report failure
    if not all_components_found:
         logger.warning(f"Cannot form clusters for {receptor_name}: one or more components missing from spatial data.")
         return None, 0, component_counts

    # Threshold in pixels (50 micrometers = 100 pixels)
    proximity_threshold = 100
    
    # Find clusters - anchor search on the component with the fewest points for efficiency
    anchor_component_name = min(component_counts, key=component_counts.get)
    other_component_names = [c for c in components if c != anchor_component_name]
    
    logger.info(f"Anchoring complex search on component: {anchor_component_name} ({component_counts[anchor_component_name]} points)")

    valid_clusters = []
    anchor_coords = component_data[anchor_component_name]['coords']

    for i, anchor_point in enumerate(anchor_coords):
        potential_cluster_points = {anchor_component_name: anchor_point}
        is_valid_cluster = True

        # Check proximity for all other components
        for other_comp_name in other_component_names:
            other_coords = component_data[other_comp_name]['coords']
            distances = np.sqrt(np.sum((other_coords - anchor_point)**2, axis=1))
            
            # Find indices of points within the threshold
            within_threshold_indices = np.where(distances <= proximity_threshold)[0]

            if len(within_threshold_indices) == 0:
                is_valid_cluster = False
                break # This anchor point cannot form a full complex

            # If multiple points are within threshold, choose the closest one
            closest_index = within_threshold_indices[np.argmin(distances[within_threshold_indices])]
            potential_cluster_points[other_comp_name] = other_coords[closest_index]

        if is_valid_cluster:
            # Calculate centroid of this valid cluster
            cluster_points_array = np.array(list(potential_cluster_points.values()))
            centroid = np.mean(cluster_points_array, axis=0)
            
            valid_clusters.append({
                'centroid': centroid,
                'component_points': cluster_points_array, # Points contributing to this centroid
                'component_names': list(potential_cluster_points.keys()) # Names of components in this cluster
            })

    if valid_clusters:
        logger.info(f"Found {len(valid_clusters)} potential complex locations for {receptor_name}")
        # Return all valid clusters found, let visualization handle limits if needed
        return valid_clusters, len(valid_clusters), component_counts
    else:
        logger.warning(f"No valid component clusters found within threshold for {receptor_name}")
        return None, 0, component_counts

def normalize_interaction_scores(scores):
    """
    Normalize interaction scores across different ligand-receptor pairs to make them comparable.
    
    This function applies min-max normalization and different normalization strategies based
    on the distribution of scores.
    
    Args:
        scores (list): List of raw interaction scores (may contain NaN values)
        
    Returns:
        dict: Dictionary with various normalization methods:
            - 'min_max': Min-max normalized scores (0-1 range)
            - 'percentile_rank': Percentile ranks (0-100 range)
            - 'z_score': Z-score normalized values
            - 'robust': Robust scaling using median and IQR
    """
    # Handle case with no valid scores
    valid_scores = [s for s in scores if s is not None and not np.isnan(s)]
    if not valid_scores:
        return {"min_max": [], "percentile_rank": [], "z_score": [], "robust": []}
    
    # Convert to numpy array for faster processing
    score_array = np.array(valid_scores)
    
    # Create mapping of original indices to normalized values
    result = {
        "min_max": [None] * len(scores),
        "percentile_rank": [None] * len(scores),
        "z_score": [None] * len(scores),
        "robust": [None] * len(scores)
    }
    
    # Min-max normalization (0-1 range)
    min_val = np.min(score_array)
    max_val = np.max(score_array)
    
    # Avoid division by zero
    if max_val > min_val:
        min_max_normalized = (score_array - min_val) / (max_val - min_val)
    else:
        # If all scores are identical, set to 1.0 (maximum)
        min_max_normalized = np.ones_like(score_array)
    
    # Calculate percentile ranks (0-100 range)
    percentiles = np.zeros_like(score_array)
    for i, score in enumerate(score_array):
        percentiles[i] = stats.percentileofscore(score_array, score)
    
    # Z-score normalization
    mean = np.mean(score_array)
    std = np.std(score_array)
    if std > 0:
        z_normalized = (score_array - mean) / std
    else:
        z_normalized = np.zeros_like(score_array)
    
    # Robust scaling using median and IQR
    median = np.median(score_array)
    q75, q25 = np.percentile(score_array, [75, 25])
    iqr = q75 - q25
    if iqr > 0:
        robust_normalized = (score_array - median) / iqr
    else:
        robust_normalized = np.zeros_like(score_array)
    
    # Map normalized values back to original indices
    valid_index = 0
    for i, score in enumerate(scores):
        if score is not None and not np.isnan(score):
            result["min_max"][i] = float(min_max_normalized[valid_index])
            result["percentile_rank"][i] = float(percentiles[valid_index])
            result["z_score"][i] = float(z_normalized[valid_index])
            result["robust"][i] = float(robust_normalized[valid_index])
            valid_index += 1
    
    return result

def main():
    parser = argparse.ArgumentParser(description='Analyze ligand-receptor interactions in spatial data')
    parser.add_argument('--lr-pairs', type=str, default='results/ligand_receptor_pairs.csv',
                        help='Path to ligand-receptor pairs CSV file')
    parser.add_argument('--spatial-data', type=str, 
                        default='/Users/jacob/Desktop/astrocytes-ontology-research/data/spatial/raw/1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv',
                        help='Path to spatial data CSV file')
    parser.add_argument('--num-pairs', type=int, default=15,
                        help='Number of top ligand-receptor pairs to analyze')
    parser.add_argument('--k-nearest', type=int, default=50,
                        help='Number of nearest ligands to consider for each receptor')
    parser.add_argument('--output-dir', type=str, default='results/spatial_analysis',
                        help='Directory to save results')
    
    args = parser.parse_args()
    
    # Create the base output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # --- UPDATED: Single density maps directory --- 
    density_maps_dir = output_dir / "density_maps"
    summary_dir = output_dir / "summary"
    
    # Create all directories
    for directory in [density_maps_dir, summary_dir]:
        directory.mkdir(exist_ok=True)
    
    # Load data
    lr_pairs = load_ligand_receptor_pairs(args.lr_pairs, args.num_pairs)
    spatial_data = load_spatial_data(args.spatial_data)
    
    # --- Stage 1: Calculate all raw scores --- 
    raw_results = []
    spatial_data_cache = {} # Cache coordinates to avoid repeated lookups

    logger.info("Stage 1: Calculating raw scores for all pairs...")
    for index, row in tqdm(lr_pairs.iterrows(), total=len(lr_pairs), desc="Calculating Raw Scores"):
        ligand = row['ligand']
        receptor = row['receptor']
        is_complex = row['is_complex_receptor']
        valid_receptors = row['valid_receptors']
        
        # Get ligand coordinates (use cache)
        if ligand not in spatial_data_cache:
            spatial_data_cache[ligand] = get_gene_spatial_data(spatial_data, ligand)
        ligand_coords = spatial_data_cache[ligand]
        num_ligands_found = len(ligand_coords) if ligand_coords is not None else 0

        # Get receptor coordinates/clusters (use cache)
        cache_key_receptor = f"{receptor}_{is_complex}_{valid_receptors}"
        if cache_key_receptor not in spatial_data_cache:
            spatial_data_cache[cache_key_receptor] = handle_complex_receptors(
                receptor, is_complex, valid_receptors, spatial_data
            )
        receptor_clusters, num_centroids, component_counts = spatial_data_cache[cache_key_receptor]

        # Extract centroid coordinates
        receptor_coords = None
        if receptor_clusters is not None and len(receptor_clusters) > 0:
            receptor_coords = np.array([cluster['centroid'] for cluster in receptor_clusters])

        avg_score_value = np.nan
        score_details = {}
        
        if ligand_coords is None or len(ligand_coords) == 0 or receptor_coords is None or len(receptor_coords) == 0 or num_centroids == 0:
            logger.warning(f"Skipping score calculation for {ligand}-{receptor}: insufficient data.")
        else:
            score_data = calculate_distance_density_score(
                ligand_coords, receptor_coords, k=args.k_nearest, ligand_name=ligand, receptor_name=receptor
            )
            if score_data:
                avg_score_value = score_data.get('average_score', np.nan)
                score_details = score_data.get('receptor_scores', {}) # Store detailed scores if needed

        # --- Store component counts as an ordered list --- 
        component_counts_ordered = []
        if is_complex and valid_receptors:
            component_names = valid_receptors.split('_')
            component_counts_ordered = [component_counts.get(comp, 0) for comp in component_names]
        elif not is_complex:
            component_counts_ordered = [component_counts.get(receptor, 0)] # List with single count for simple receptor
        else: # Handle edge cases like complex receptor with empty valid_receptors string
             component_counts_ordered = [0] * len(valid_receptors.split('_')) if valid_receptors else []
             
        # Store raw results
        raw_results.append({
            'index': index, # Keep original index for joining later
            'ligand': ligand,
            'receptor': receptor,
            'is_complex': is_complex,
            'valid_receptors': valid_receptors,
            'avg_score': avg_score_value,
            'num_ligands_found': num_ligands_found,
            'num_receptors_found': num_centroids,
            'receptor_scores_detail': score_details,
            'component_counts_ordered': component_counts_ordered, # Store ordered list
            # Store references needed for visualization
            'ligand_coords_ref': ligand_coords,
            'receptor_clusters_ref': receptor_clusters
        })

    # --- Stage 2: Normalize scores --- 
    logger.info("Stage 2: Normalizing scores across all pairs...")
    raw_scores_list = [res['avg_score'] for res in raw_results]
    normalized_scores = normalize_interaction_scores(raw_scores_list)
    
    # Combine raw results with normalized scores
    results_df = pd.DataFrame(raw_results)
    # Add normalized scores based on the list order (which matches raw_results)
    results_df['score_normalized'] = normalized_scores['min_max']
    results_df['score_percentile'] = normalized_scores['percentile_rank']
    results_df['score_zscore'] = normalized_scores['z_score']
    results_df['score_robust'] = normalized_scores['robust']

    # --- Stage 3: Save results and visualize --- 
    logger.info("Stage 3: Saving results and generating visualizations...")
    if not results_df.empty:
        # --- UPDATED: Select and order columns for CSV, including ordered counts --- 
        base_cols = ['ligand', 'receptor', 'is_complex', 'valid_receptors', 
                     'avg_score', 'score_normalized', 'score_percentile', 'score_zscore', 'score_robust',
                     'num_ligands_found', 'num_receptors_found', 'component_counts_ordered']
        # Ensure only columns present in the DataFrame are selected
        final_cols = [col for col in base_cols if col in results_df.columns] 
        csv_df = results_df[final_cols] # Directly select final columns
        csv_path = summary_dir / 'lr_spatial_analysis_results.csv'
        csv_df.to_csv(csv_path, index=False)
        logger.info(f"Results saved to {csv_path}")

        # Generate density map visualizations (now with normalized scores available)
        for index, row in tqdm(results_df.iterrows(), total=len(results_df), desc="Generating Density Maps"):
            # Visualize only if scoring was successful
            if row['avg_score'] is not None and not np.isnan(row['avg_score']):
                # Prepare score_data dict for visualization function
                viz_score_data = {
                    'average_score': row['avg_score'],
                    'normalized_score': row['score_normalized'], # Pass the normalized score
                    'receptor_scores': row['receptor_scores_detail'] 
                }
                
                # --- CALL VISUALIZE ONCE --- 
                visualize_interaction(
                    row['ligand_coords_ref'],
                    row['receptor_clusters_ref'], 
                    viz_score_data, 
                    row['ligand'], row['receptor'], row['is_complex'],
                    density_maps_dir # Save to single density map dir
                )
            else:
                 logger.warning(f"Skipping visualization for {row['ligand']}-{row['receptor']} due to NaN score.")
    else:
        logger.warning("No results were generated to save to CSV or visualize.")

    # --- Stage 4: Create summary plots --- 
    logger.info("Stage 4: Generating summary score plots...")
    if not results_df.empty:
        results_df['avg_score'] = pd.to_numeric(results_df['avg_score'], errors='coerce')
        results_df['score_normalized'] = pd.to_numeric(results_df['score_normalized'], errors='coerce')
        
        plot_df = results_df.dropna(subset=['avg_score']).copy()
        
        if not plot_df.empty:
            # Create plot for raw scores
            plot_df['plot_id'] = plot_df['ligand'] + '_' + plot_df['receptor'] + '_' + plot_df.index.astype(str)
            plot_df = plot_df.sort_values('avg_score', ascending=False)
            plt.figure(figsize=(15, max(8, len(plot_df) * 0.5)))
            sns.barplot(data=plot_df, x='avg_score', y='ligand', hue='receptor', dodge=True)
            plt.title(f'Top {len(plot_df)} Ligand-Receptor Pairs by Raw Interaction Density Score')
            plt.xlabel('Avg. Interaction Score (Ligand KDE x Receptor KDE - Higher is Better)')
            plt.ylabel('Ligand')
            plt.legend(title='Receptor', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
            plt.tight_layout(rect=[0, 0, 0.85, 1])
            plt.savefig(summary_dir / 'lr_scores_summary_raw.png', dpi=300)
            plt.close()
            
            # Create plot for normalized scores
            plot_df = plot_df.sort_values('score_normalized', ascending=False)
            plt.figure(figsize=(15, max(8, len(plot_df) * 0.5)))
            sns.barplot(data=plot_df, x='score_normalized', y='ligand', hue='receptor', dodge=True)
            plt.title(f'Top {len(plot_df)} Ligand-Receptor Pairs by Normalized Interaction Density Score')
            plt.xlabel('Normalized Interaction Score (0-1 range, Higher is Better)')
            plt.ylabel('Ligand')
            plt.legend(title='Receptor', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
            plt.tight_layout(rect=[0, 0, 0.85, 1])
            plt.savefig(summary_dir / 'lr_scores_summary_normalized.png', dpi=300)
            plt.close()
            
            logger.info(f"Summary score plots saved to {summary_dir}")
        else:
             logger.warning("No valid results to plot for summary scores after cleaning.")
    else:
        logger.warning("No results generated, skipping summary plots.")

if __name__ == "__main__":
    main() 