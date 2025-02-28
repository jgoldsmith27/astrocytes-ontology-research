import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

def convert_csv_to_turtle(csv_file, output_ttl, sample_size=100, proximity_threshold=50, max_neighbors=5):
    """
    Convert spatial transcriptomics data to TURTLE format
    
    Args:
        csv_file: Path to the CSV file
        output_ttl: Path to save the TURTLE file
        sample_size: Number of rows to sample
        proximity_threshold: Distance threshold to define 'locatedNear' relationship
        max_neighbors: Maximum number of neighbors to include for each point
    """
    print(f"Reading data from {csv_file}...")
    
    # Sample the dataset
    if sample_size > 0:
        # Read a random sample of rows
        df = pd.read_csv(csv_file, skiprows=lambda x: x > 0 and np.random.random() > sample_size/9000000)
        # Ensure we don't have too many rows
        if len(df) > sample_size:
            df = df.sample(sample_size, random_state=42)
    else:
        # Read the entire file
        df = pd.read_csv(csv_file)
    
    print(f"Processing {len(df)} records...")
    print(f"Dataset columns: {df.columns.tolist()}")
    
    # Start the TURTLE file with prefixes and ontology definitions
    ttl_lines = [
        "@prefix : <http://example.org/spatial-transcriptomics#> .",
        "@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .",
        "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .",
        "@prefix owl: <http://www.w3.org/2002/07/owl#> .",
        "@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .",
        "",
        "# Define basic ontology classes",
        ":Gene rdf:type owl:Class .",
        ":SpatialPoint rdf:type owl:Class .",
        "",
        "# Define relationships",
        ":expressedAt rdf:type owl:ObjectProperty ;",
        "             rdfs:domain :Gene ;",
        "             rdfs:range :SpatialPoint .",
        "",
        ":hasExpressionLevel rdf:type owl:DatatypeProperty ;",
        "                    rdfs:domain :Gene ;",
        "                    rdfs:range xsd:integer .",
        "",
        ":hasExonCount rdf:type owl:DatatypeProperty ;",
        "              rdfs:domain :Gene ;",
        "              rdfs:range xsd:integer .",
        "",
        ":hasIntronCount rdf:type owl:DatatypeProperty ;",
        "                rdfs:domain :Gene ;",
        "                rdfs:range xsd:integer .",
        "",
        ":locatedNear rdf:type owl:ObjectProperty ;",
        "             rdfs:domain :SpatialPoint ;",
        "             rdfs:range :SpatialPoint .",
        "",
        "# Data Instances"
    ]
    
    # Create a dictionary to store points and their coordinates
    points = {}
    point_to_genes = {}
    
    print("Creating gene and point instances...")
    
    # First pass: create gene and point instances
    for index, row in tqdm(df.iterrows(), total=len(df)):
        gene_id = row['geneID']
        x, y = int(row['x']), int(row['y'])
        expr_level = int(row['MIDCount'])
        exon_count = int(row['ExonCount'])
        intron_count = int(row['IntronCount'])
        
        gene_uri = f":Gene_{gene_id}_{index}"
        point_uri = f":Point_{x}_{y}"
        
        # Store point coordinates for proximity calculation
        points[point_uri] = (x, y)
        
        # Track which genes are expressed at each point
        if point_uri not in point_to_genes:
            point_to_genes[point_uri] = []
        point_to_genes[point_uri].append(gene_uri)
        
        # Add gene instance
        ttl_lines.append(f"{gene_uri} rdf:type :Gene ;")
        ttl_lines.append(f"    :expressedAt {point_uri} ;")
        ttl_lines.append(f"    :hasExpressionLevel \"{expr_level}\"^^xsd:integer ;")
        ttl_lines.append(f"    :hasExonCount \"{exon_count}\"^^xsd:integer ;")
        ttl_lines.append(f"    :hasIntronCount \"{intron_count}\"^^xsd:integer .")
    
    print("Creating spatial point instances...")
    
    # Add point instances
    for point_uri in tqdm(points.keys()):
        ttl_lines.append(f"{point_uri} rdf:type :SpatialPoint .")
    
    # Use KDTree for efficient nearest neighbor search
    if proximity_threshold > 0:
        print("Calculating spatial relationships...")
        
        # Extract coordinates and point URIs
        coords = np.array(list(points.values()))
        point_uris = list(points.keys())
        
        # Build KDTree
        tree = KDTree(coords)
        
        # Query for neighbors within threshold
        for i, point_uri in enumerate(tqdm(point_uris)):
            # Find neighbors within threshold
            neighbors = tree.query_ball_point(coords[i], proximity_threshold)
            
            # Remove self from neighbors
            neighbors = [j for j in neighbors if j != i]
            
            # Limit number of neighbors
            neighbors = neighbors[:max_neighbors]
            
            if neighbors:
                ttl_lines.append(f"{point_uri}")
                for j, neighbor_idx in enumerate(neighbors):
                    neighbor_uri = point_uris[neighbor_idx]
                    if j < len(neighbors) - 1:
                        ttl_lines.append(f"    :locatedNear {neighbor_uri} ;")
                    else:
                        ttl_lines.append(f"    :locatedNear {neighbor_uri} .")
    
    print(f"Writing {len(ttl_lines)} lines to {output_ttl}...")
    
    # Write to file
    with open(output_ttl, 'w') as f:
        f.write('\n'.join(ttl_lines))
    
    print(f"TURTLE file saved: {output_ttl}")
    print(f"Total triples generated: approximately {len(ttl_lines) - 30}")  # Rough estimate
    
    return df

def visualize_sample(df, output_image=None):
    """
    Visualize the spatial distribution of genes in the sample
    """
    plt.figure(figsize=(12, 10))
    
    # Get unique genes for coloring
    unique_genes = df['geneID'].unique()
    gene_to_color = {gene: plt.cm.tab20(i % 20) for i, gene in enumerate(unique_genes)}
    
    # Plot points colored by gene
    for gene in unique_genes:
        gene_df = df[df['geneID'] == gene]
        plt.scatter(gene_df['x'], gene_df['y'], 
                   alpha=0.7, 
                   label=gene if len(unique_genes) <= 20 else None,
                   color=gene_to_color[gene],
                   s=gene_df['MIDCount'] * 5)  # Size by expression level
    
    if len(unique_genes) <= 20:
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.title(f'Spatial Distribution of Genes (Sample of {len(df)} points)')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.tight_layout()
    
    if output_image:
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_image}")
    
    plt.show()

if __name__ == "__main__":
    # File paths
    csv_file = "/Users/jacob/Desktop/Karolinska Research/astrocytes-ontology-research/Spatial_Data/1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv"
    output_dir = "/Users/jacob/Desktop/Karolinska Research/astrocytes-ontology-research/Spatial_Data"
    output_ttl = os.path.join(output_dir, "spatial_transcriptomics_sample.ttl")
    output_image = os.path.join(output_dir, "spatial_visualization.png")
    
    # Parameters
    sample_size = 200  # Number of rows to sample
    proximity_threshold = 100  # Distance threshold for "locatedNear" relationship
    max_neighbors = 3  # Maximum number of neighbors to include for each point
    
    # Convert data to TURTLE format
    df = convert_csv_to_turtle(
        csv_file=csv_file,
        output_ttl=output_ttl,
        sample_size=sample_size,
        proximity_threshold=proximity_threshold,
        max_neighbors=max_neighbors
    )
    
    # Visualize the sample
    visualize_sample(df, output_image) 