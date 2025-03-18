#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans, DBSCAN
from sklearn.preprocessing import StandardScaler
from scipy.ndimage import gaussian_filter
import os
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import cKDTree
from sklearn.metrics import silhouette_score
import matplotlib.patches as mpatches

# Set working directory
os.chdir('/Users/jacob/Desktop/astrocytes-ontology-research')

print("Loading spatial data and interaction sites...")
# Load the spatial data
spatial_data = pd.read_csv('raw/1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv')
interaction_sites = pd.read_csv('interaction_sites.csv')

print("\n=== Creating Brain Region and Communication Activity Cluster Map ===")

# 1. Creating a spatial grid for the tissue
print("Creating spatial grid...")
# Determine boundaries of the data
x_min, x_max = spatial_data['x'].min(), spatial_data['x'].max()
y_min, y_max = spatial_data['y'].min(), spatial_data['y'].max()

# Create a grid (adjust resolution as needed)
grid_size = 100  # number of bins in each dimension
x_edges = np.linspace(x_min, x_max, grid_size)
y_edges = np.linspace(y_min, y_max, grid_size)

# 2. Calculate gene expression density for the overall tissue
print("Calculating overall gene expression density...")
gene_density, _, _ = np.histogram2d(
    spatial_data['x'], spatial_data['y'], 
    bins=[x_edges, y_edges]
)
# Smooth the density for better visualization
gene_density_smooth = gaussian_filter(gene_density, sigma=2)

# 3. Calculate interaction density
print("Calculating interaction density...")
interaction_density, _, _ = np.histogram2d(
    interaction_sites['x'], interaction_sites['y'], 
    bins=[x_edges, y_edges]
)
# Smooth the density for better visualization
interaction_density_smooth = gaussian_filter(interaction_density, sigma=2)

# Normalize the densities
gene_density_norm = gene_density_smooth / gene_density_smooth.max()
interaction_density_norm = interaction_density_smooth / interaction_density_smooth.max()

# 4. Create a measure of "communication activity"
# This is the ratio of interaction density to gene expression density
# We add a small constant to avoid division by zero
eps = 1e-10
communication_activity = interaction_density_smooth / (gene_density_smooth + eps)
# Normalize to [0, 1]
communication_activity = (communication_activity - communication_activity.min()) / (communication_activity.max() - communication_activity.min() + eps)

# 5. Use machine learning to identify regions with distinct patterns
# Prepare data for clustering
print("Clustering the spatial data...")
# Get grid cell centers
x_centers = (x_edges[:-1] + x_edges[1:]) / 2
y_centers = (y_edges[:-1] + y_edges[1:]) / 2
X, Y = np.meshgrid(x_centers, y_centers)

# Create feature matrix for clustering
# We include spatial coordinates (X, Y), gene density, and communication activity
features = np.column_stack([
    X.flatten(), 
    Y.flatten(), 
    gene_density_norm.flatten(),
    communication_activity.flatten()
])

# Remove rows with NaN values
valid_mask = ~np.isnan(features).any(axis=1)
features_valid = features[valid_mask]

# Standardize features for clustering
scaler = StandardScaler()
features_scaled = scaler.fit_transform(features_valid)

# 6. Use KMeans to cluster the regions
# Determine optimal number of clusters using silhouette score
max_clusters = 10
silhouette_scores = []

print("Finding optimal number of clusters...")
for n_clusters in range(2, max_clusters + 1):
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(features_scaled)
    
    # Calculate silhouette score if we have more than one cluster
    if n_clusters > 1:
        score = silhouette_score(features_scaled, cluster_labels)
        silhouette_scores.append(score)
        print(f"  {n_clusters} clusters: silhouette score = {score:.3f}")

# Choose number of clusters that maximizes silhouette score
optimal_clusters = np.argmax(silhouette_scores) + 2  # +2 because we started at 2 clusters
print(f"Optimal number of clusters: {optimal_clusters}")

# Fit final model with optimal number of clusters
kmeans = KMeans(n_clusters=optimal_clusters, random_state=42, n_init=10)
cluster_labels = kmeans.fit_predict(features_scaled)

# 7. Create a cluster map
# Prepare a grid for visualization
cluster_grid = np.full((grid_size-1, grid_size-1), np.nan)
grid_indices = np.arange(len(X.flatten()))[valid_mask]
cluster_grid.flat[grid_indices] = cluster_labels

# 8. Analyze the characteristics of each cluster
cluster_stats = []
for cluster_id in range(optimal_clusters):
    mask = cluster_labels == cluster_id
    cluster_stats.append({
        'cluster_id': cluster_id,
        'count': np.sum(mask),
        'avg_x': np.mean(features_valid[mask, 0]),
        'avg_y': np.mean(features_valid[mask, 1]),
        'avg_gene_density': np.mean(features_valid[mask, 2]),
        'avg_communication': np.mean(features_valid[mask, 3])
    })

cluster_stats_df = pd.DataFrame(cluster_stats)
print("\nCluster statistics:")
print(cluster_stats_df.sort_values('avg_communication', ascending=False))

# Save cluster statistics
cluster_stats_df.to_csv('brain_region_clusters.csv', index=False)
print("Cluster statistics saved to brain_region_clusters.csv")

# 9. Calculate interaction counts in each cluster
# Create a mapping from spatial coordinates to cluster ID
cluster_map = {}
for i, (x, y, cluster) in enumerate(zip(
    features_valid[:, 0], 
    features_valid[:, 1], 
    cluster_labels
)):
    cluster_map[(round(x), round(y))] = cluster

# Assign clusters to interaction sites
interaction_sites['cluster'] = -1  # Default to -1 (no cluster)
for i, (x, y) in enumerate(zip(interaction_sites['x'], interaction_sites['y'])):
    # Find the nearest grid point
    x_idx = np.abs(x_centers - x).argmin()
    y_idx = np.abs(y_centers - y).argmin()
    grid_key = (round(x_centers[x_idx]), round(y_centers[y_idx]))
    
    # Assign cluster if the grid point is in a cluster
    if grid_key in cluster_map:
        interaction_sites.loc[i, 'cluster'] = cluster_map[grid_key]

# Count interactions by type in each cluster
interaction_counts = pd.crosstab(
    interaction_sites['cluster'], 
    interaction_sites['pair']
)

# Save interaction counts by cluster
interaction_counts.to_csv('cluster_interaction_counts.csv')
print("Interaction counts by cluster saved to cluster_interaction_counts.csv")

# 10. Visualize the clusters
print("\nCreating visualizations...")
plt.figure(figsize=(20, 16))

# Plot 1: Cluster map with interaction sites
plt.subplot(221)
cluster_cmap = plt.cm.get_cmap('viridis', optimal_clusters)
plt.imshow(cluster_grid.T, origin='lower', extent=[x_min, x_max, y_min, y_max], cmap=cluster_cmap, alpha=0.7)

# Add interaction sites as scatter points
for cluster in range(optimal_clusters):
    cluster_interactions = interaction_sites[interaction_sites['cluster'] == cluster]
    plt.scatter(cluster_interactions['x'], cluster_interactions['y'], 
                s=10, alpha=0.6, label=f'Cluster {cluster}')

plt.colorbar(label='Cluster ID')
plt.title('Brain Region Clusters and Interaction Sites')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')

# Plot 2: Communication activity heatmap
plt.subplot(222)
comm_cmap = plt.cm.YlOrRd
plt.imshow(communication_activity.T, origin='lower', extent=[x_min, x_max, y_min, y_max], 
           cmap=comm_cmap, alpha=0.9)
plt.colorbar(label='Communication Activity (Normalized)')
plt.title('Communication Activity Heatmap')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')

# Plot 3: Cluster map with communication activity overlay
plt.subplot(223)
# Create custom colormap for clusters
cluster_colors = plt.cm.tab10(np.linspace(0, 1, optimal_clusters))
custom_cmap = LinearSegmentedColormap.from_list('custom_cluster', cluster_colors, N=optimal_clusters)

# Plot cluster map
plt.imshow(cluster_grid.T, origin='lower', extent=[x_min, x_max, y_min, y_max], 
           cmap=custom_cmap, alpha=0.7)

# Overlay communication activity contours
comm_levels = np.linspace(0.2, 0.9, 5)
contour = plt.contour(X, Y, communication_activity, levels=comm_levels, 
                      colors='black', alpha=0.7, linewidths=1)
plt.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

plt.colorbar(label='Cluster ID')
plt.title('Brain Regions with Communication Activity Contours')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')

# Plot 4: Interaction type distribution by cluster
plt.subplot(224)

# Get the top interaction types
top_interactions = interaction_sites['pair'].value_counts().head(5).index.tolist()

# Create a matrix of interaction types by cluster
interaction_matrix = pd.DataFrame(index=range(optimal_clusters))
for pair in top_interactions:
    for cluster in range(optimal_clusters):
        count = interaction_sites[(interaction_sites['cluster'] == cluster) & 
                                  (interaction_sites['pair'] == pair)].shape[0]
        interaction_matrix.loc[cluster, pair] = count

# Normalize by cluster size
for cluster in range(optimal_clusters):
    cluster_size = interaction_sites[interaction_sites['cluster'] == cluster].shape[0]
    if cluster_size > 0:
        interaction_matrix.loc[cluster] = interaction_matrix.loc[cluster] / cluster_size * 100

# Plot the heatmap
sns.heatmap(interaction_matrix, annot=True, fmt='.1f', cmap='YlGnBu')
plt.title('Interaction Type Distribution by Cluster (%)')
plt.ylabel('Cluster ID')
plt.xlabel('Interaction Type')

plt.tight_layout()
plt.savefig('brain_region_communication_clusters.png', dpi=300)
print("Visualization saved to brain_region_communication_clusters.png")

# 11. Create a 3D scatter plot showing cluster characteristics
# This helps visualize the relationship between gene density, communication activity, and spatial location
from mpl_toolkits.mplot3d import Axes3D

plt.figure(figsize=(12, 10))
ax = plt.axes(projection='3d')

# Get cluster centers and characteristics
cluster_centers = kmeans.cluster_centers_

# Scatter plot for each cluster
for cluster in range(optimal_clusters):
    mask = cluster_labels == cluster
    points = features_scaled[mask]
    
    # Unscale the features for better interpretation
    points_unscaled = scaler.inverse_transform(points)
    
    # 3D scatter plot with:
    # - X, Y: spatial coordinates
    # - Z: communication activity
    # - Size: gene density
    # - Color: cluster ID
    ax.scatter(
        points_unscaled[:, 0], 
        points_unscaled[:, 1], 
        points_unscaled[:, 3],  # Communication activity
        s=points_unscaled[:, 2]*50 + 10,  # Gene density determines point size
        alpha=0.6,
        c=[cluster_colors[cluster]],
        label=f'Cluster {cluster}'
    )

ax.set_xlabel('X coordinate')
ax.set_ylabel('Y coordinate')
ax.set_zlabel('Communication Activity')
plt.title('3D Visualization of Brain Region Clusters')
plt.legend()
plt.savefig('brain_region_clusters_3d.png', dpi=300)
print("3D visualization saved to brain_region_clusters_3d.png")

# 12. Create a summary of each cluster's characteristics
cluster_summary = []
for cluster in range(optimal_clusters):
    # Get interactions in this cluster
    cluster_inters = interaction_sites[interaction_sites['cluster'] == cluster]
    
    # Top interaction types
    if len(cluster_inters) > 0:
        top_types = cluster_inters['pair'].value_counts().head(3).index.tolist()
        top_types_str = ', '.join(top_types)
    else:
        top_types_str = 'None'
    
    # Cluster location description (quadrant of the tissue)
    center_x = cluster_stats_df.loc[cluster_stats_df['cluster_id'] == cluster, 'avg_x'].values[0]
    center_y = cluster_stats_df.loc[cluster_stats_df['cluster_id'] == cluster, 'avg_y'].values[0]
    
    x_rel = 'east' if center_x > (x_max + x_min)/2 else 'west'
    y_rel = 'north' if center_y > (y_max + y_min)/2 else 'south'
    location = f"{y_rel}-{x_rel}"
    
    # Communication activity level
    comm_activity = cluster_stats_df.loc[cluster_stats_df['cluster_id'] == cluster, 'avg_communication'].values[0]
    if comm_activity > 0.7:
        activity_level = 'Very High'
    elif comm_activity > 0.5:
        activity_level = 'High'
    elif comm_activity > 0.3:
        activity_level = 'Medium'
    else:
        activity_level = 'Low'
    
    # Add to summary
    cluster_summary.append({
        'cluster_id': cluster,
        'location': location,
        'communication_activity': activity_level,
        'interaction_count': len(cluster_inters),
        'top_interaction_types': top_types_str,
        'gene_density': cluster_stats_df.loc[cluster_stats_df['cluster_id'] == cluster, 'avg_gene_density'].values[0]
    })

# Save summary as CSV
cluster_summary_df = pd.DataFrame(cluster_summary)
cluster_summary_df.to_csv('brain_region_cluster_summary.csv', index=False)
print("Cluster summary saved to brain_region_cluster_summary.csv")

print("\nAnalysis complete!")
print(f"Identified {optimal_clusters} distinct brain regions with varying levels of cell communication activity.")
print("Use the generated visualizations and summary files to explore the patterns of cell communication across different brain regions.") 