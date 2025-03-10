#!/usr/bin/env python3
"""
Script to apply co-expression rules to spatial data for whole-cell identification.

This script reads the SPARQL inference rules generated from single-cell data
and applies them to spatial data to identify whole cells based on co-expression patterns.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdflib import Graph, Namespace, Literal, URIRef, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL
import re
import uuid
from tqdm import tqdm

# Define namespaces
CELL = Namespace("http://example.org/ontology/cell/")
EX = Namespace("http://example.org/data/spatial/")

def load_graphs(spatial_rdf, rules_file):
    """
    Load spatial RDF data and co-expression rules.
    
    Parameters:
    -----------
    spatial_rdf : str
        Path to the spatial RDF file
    rules_file : str
        Path to the co-expression rules file
    
    Returns:
    --------
    g : rdflib.Graph
        The combined RDF graph
    """
    print("Loading RDF graphs...")
    
    # Create a new RDF graph
    g = Graph()
    
    # Load the spatial data
    print(f"Loading spatial data from {spatial_rdf}...")
    g.parse(spatial_rdf, format="turtle")
    
    # Bind namespaces for readability in SPARQL queries
    g.bind("cell", CELL)
    g.bind("ex", EX)
    g.bind("rdf", RDF)
    g.bind("rdfs", RDFS)
    g.bind("xsd", XSD)
    g.bind("owl", OWL)
    
    print(f"Total triples before adding inferences: {len(g)}")
    
    return g

def apply_rules(g, rules_file, distance_threshold=50):
    """
    Apply co-expression rules to identify cells in spatial data.
    
    Parameters:
    -----------
    g : rdflib.Graph
        The RDF graph with spatial data
    rules_file : str
        Path to the co-expression rules file
    distance_threshold : float
        Maximum distance between co-expressed genes to be considered part of the same cell
    
    Returns:
    --------
    g : rdflib.Graph
        The RDF graph with inferred cells
    cells_info : list
        Information about the identified cells
    """
    print(f"Applying co-expression rules from {rules_file}...")
    
    # Read the rules file
    with open(rules_file, 'r') as f:
        rules_content = f.read()
    
    # Extract individual rules using regex
    rule_pattern = r'# Rule: (.+?)\n(PREFIX.+?)(?=# Rule:|$)'
    rules = re.findall(rule_pattern, rules_content, re.DOTALL)
    
    cells_info = []
    
    # Apply each rule
    for rule_id, rule_content in tqdm(rules):
        print(f"Applying rule: {rule_id}")
        
        try:
            # Extract the cell type from the rule ID
            cell_type = rule_id.split('_')[0]
            
            # Execute the CONSTRUCT query
            result = g.update(rule_content)
            
            # Count how many cells were added
            cell_count_query = f"""
            PREFIX cell: <http://example.org/ontology/cell/>
            
            SELECT (COUNT(?cell) as ?count)
            WHERE {{
                ?cell a cell:SpatialCell ;
                     cell:hasCellType cell:{cell_type} .
            }}
            """
            
            count_results = list(g.query(cell_count_query))
            cells_added = int(count_results[0][0]) if count_results else 0
            
            print(f"  Added {cells_added} cells of type {cell_type}")
            
            # Get information about the cells for visualization
            if cells_added > 0:
                cell_info_query = f"""
                PREFIX cell: <http://example.org/ontology/cell/>
                
                SELECT ?cell ?x ?y ?radius ?confidence
                WHERE {{
                    ?cell a cell:SpatialCell ;
                         cell:hasCellType cell:{cell_type} ;
                         cell:hasX ?x ;
                         cell:hasY ?y ;
                         cell:hasRadius ?radius ;
                         cell:hasConfidence ?confidence .
                }}
                """
                
                for row in g.query(cell_info_query):
                    cell_uri, x, y, radius, confidence = row
                    cells_info.append({
                        'cell_uri': str(cell_uri),
                        'cell_type': cell_type,
                        'x': float(x),
                        'y': float(y),
                        'radius': float(radius),
                        'confidence': float(confidence)
                    })
        
        except Exception as e:
            print(f"Error applying rule {rule_id}: {e}")
    
    print(f"Total triples after adding inferences: {len(g)}")
    
    return g, cells_info

def validate_cells(g, cells_info, min_confidence=0.7):
    """
    Validate identified cells based on biological ground truth rules.
    
    Parameters:
    -----------
    g : rdflib.Graph
        The RDF graph with inferred cells
    cells_info : list
        Information about the identified cells
    min_confidence : float
        Minimum confidence threshold
    
    Returns:
    --------
    validated_cells : list
        Information about the validated cells
    """
    print("Validating identified cells...")
    
    validated_cells = []
    
    # Filter by confidence
    high_confidence_cells = [cell for cell in cells_info if cell['confidence'] >= min_confidence]
    print(f"Filtered {len(cells_info) - len(high_confidence_cells)} cells below confidence threshold")
    
    # Implement biological validation rules
    for cell in high_confidence_cells:
        # Get all genes in this cell
        cell_uri = URIRef(cell['cell_uri'])
        
        genes_query = f"""
        PREFIX cell: <http://example.org/ontology/cell/>
        
        SELECT ?geneID
        WHERE {{
            <{cell_uri}> cell:includesSpatialPoint ?point .
            ?point cell:expressesGene ?gene .
            ?gene cell:hasGeneID ?geneID .
        }}
        """
        
        genes = [str(row[0]) for row in g.query(genes_query)]
        
        # Count nearby cells of other types
        nearby_cells_query = f"""
        PREFIX cell: <http://example.org/ontology/cell/>
        
        SELECT (COUNT(?other) as ?count)
        WHERE {{
            ?other a cell:SpatialCell ;
                  cell:hasX ?x ;
                  cell:hasY ?y .
            FILTER(?other != <{cell_uri}>)
            BIND(SQRT(POW(?x - {cell['x']}, 2) + POW(?y - {cell['y']}, 2)) AS ?dist)
            FILTER(?dist < {cell['radius'] * 1.5})
        }}
        """
        
        nearby_count = int(list(g.query(nearby_cells_query))[0][0])
        
        # If cell passes validation, add it to the validated list
        cell['genes'] = genes
        cell['nearby_cells'] = nearby_count
        
        # Skip cells with too many nearby cells (likely false positives)
        if nearby_count <= 3:
            validated_cells.append(cell)
        else:
            print(f"Filtered cell with {nearby_count} nearby cells")
    
    print(f"Validated {len(validated_cells)} of {len(cells_info)} cells")
    return validated_cells

def visualize_cells(cells_info, spatial_points, output_file):
    """
    Visualize the identified cells.
    
    Parameters:
    -----------
    cells_info : list
        Information about the identified cells
    spatial_points : dict
        Dictionary with spatial points data
    output_file : str
        Path to the output visualization file
    """
    print(f"Creating visualization at {output_file}...")
    
    plt.figure(figsize=(12, 10))
    
    # Plot spatial points as small dots
    plt.scatter(
        spatial_points['x'], 
        spatial_points['y'], 
        c='lightgray', 
        alpha=0.3, 
        s=1, 
        label='Spatial Points'
    )
    
    # Define colors for different cell types
    colors = {
        'Astrocyte': 'blue',
        'Astrocyte1': 'red',
    }
    
    # Plot cells by type
    for cell_type in set(cell['cell_type'] for cell in cells_info):
        cell_type_data = [cell for cell in cells_info if cell['cell_type'] == cell_type]
        
        # Plot cell centers
        plt.scatter(
            [cell['x'] for cell in cell_type_data],
            [cell['y'] for cell in cell_type_data],
            c=colors.get(cell_type, 'green'),
            alpha=0.7,
            s=50,
            label=f'{cell_type} ({len(cell_type_data)})'
        )
        
        # Plot cell boundaries as circles
        for cell in cell_type_data:
            circle = plt.Circle(
                (cell['x'], cell['y']), 
                cell['radius'], 
                fill=False, 
                color=colors.get(cell_type, 'green'),
                alpha=0.3
            )
            plt.gca().add_patch(circle)
    
    plt.title("Spatial Distribution of Identified Cells")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.legend()
    plt.grid(alpha=0.3)
    
    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"Visualization saved to {output_file}")

def export_cells_to_csv(cells_info, output_file):
    """
    Export information about identified cells to CSV.
    
    Parameters:
    -----------
    cells_info : list
        Information about the identified cells
    output_file : str
        Path to the output CSV file
    """
    print(f"Exporting cell information to {output_file}...")
    
    cells_df = pd.DataFrame(cells_info)
    cells_df.to_csv(output_file, index=False)
    
    print(f"Exported {len(cells_info)} cells to {output_file}")

def get_cells_by_type(g, cell_type):
    """Get all cells of a specific type."""
    query = """
    PREFIX cell: <http://example.org/ontology/cell/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?cell WHERE {
      ?cell a cell:SpatialCell ;
            cell:hasCellType cell:{cell_type} .
    }
    """
    query = query.format(cell_type=cell_type)
    results = g.query(query)
    return [str(row.cell) for row in results]

def get_cell_info(g, cell_type=None):
    """Get information about all cells, optionally filtered by type."""
    query = """
    PREFIX cell: <http://example.org/ontology/cell/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    
    SELECT ?cell ?cellType ?x ?y ?radius ?confidence WHERE {
      ?cell a cell:SpatialCell ;
            cell:hasCellType ?cellType ;
            cell:hasX ?x ;
            cell:hasY ?y ;
            cell:hasRadius ?radius ;
            cell:hasConfidence ?confidence .
    """
    
    if cell_type:
        query += f"  FILTER(?cellType = \"cell:{cell_type}\")\n"
    
    query += "}"
    
    results = g.query(query)
    return [(str(row.cell), str(row.cellType), 
             float(row.x), float(row.y), 
             float(row.radius), float(row.confidence)) 
            for row in results]

def get_cell_genes(g, cell_uri):
    """Get all genes expressed in a cell."""
    query = """
    PREFIX cell: <http://example.org/ontology/cell/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT DISTINCT ?geneID WHERE {
      <{cell_uri}> cell:includesSpatialPoint ?point .
      ?point cell:expressesGene ?gene .
      ?gene cell:hasGeneID ?geneID .
    }
    """
    query = query.format(cell_uri=cell_uri)
    results = g.query(query)
    return [str(row.geneID) for row in results]

def get_nearby_cells(g, x, y, max_distance=100, cell_type=None):
    """Get cells near a specific location."""
    query = """
    PREFIX cell: <http://example.org/ontology/cell/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    
    SELECT ?other ?x ?y WHERE {
      ?other a cell:SpatialCell ;
             cell:hasX ?x ;
             cell:hasY ?y .
    """
    
    if cell_type:
        query += f"  ?other cell:hasCellType \"cell:{cell_type}\" .\n"
    
    query += f"""
      BIND(SQRT(POW(?x - {x}, 2) + POW(?y - {y}, 2)) AS ?dist)
      FILTER(?dist < {max_distance})
    }}
    ORDER BY ?dist
    """
    
    results = g.query(query)
    return [(str(row.other), float(row.x), float(row.y)) for row in results]

def get_spatial_points(g, x_min, y_min, x_max, y_max):
    """Get spatial points within a bounding box."""
    query = """
    PREFIX cell: <http://example.org/ontology/cell/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    
    SELECT ?point ?x ?y WHERE {
      ?point a cell:SpatialDataPoint ;
             cell:hasXCoordinate ?x ;
             cell:hasYCoordinate ?y .
    
      FILTER(?x >= %d && ?x <= %d && ?y >= %d && ?y <= %d)
    }
    """ % (x_min, x_max, y_min, y_max)
    
    results = g.query(query)
    return [(str(row.point), int(row.x), int(row.y)) for row in results]

def main():
    parser = argparse.ArgumentParser(description="Apply co-expression rules to identify cells in spatial data")
    parser.add_argument("--spatial", "-s", required=True, help="Spatial RDF file")
    parser.add_argument("--rules", "-r", required=True, help="Co-expression rules file")
    parser.add_argument("--output-csv", "-o", required=True, help="Output CSV file for cell information")
    parser.add_argument("--output-viz", "-v", help="Output visualization file")
    parser.add_argument("--min-confidence", "-c", type=float, default=0.7, help="Minimum confidence threshold")
    parser.add_argument("--distance-threshold", "-d", type=float, default=50, help="Maximum distance for co-expression")
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
    if args.output_viz:
        os.makedirs(os.path.dirname(args.output_viz), exist_ok=True)
    
    # Load graphs
    g = load_graphs(args.spatial, args.rules)
    
    # Apply rules
    g, cells_info = apply_rules(g, args.rules, args.distance_threshold)
    
    # Validate cells
    validated_cells = validate_cells(g, cells_info, args.min_confidence)
    
    # Export cell information
    export_cells_to_csv(validated_cells, args.output_csv)
    
    # Get spatial points for visualization
    if args.output_viz:
        spatial_query = """
        PREFIX cell: <http://example.org/ontology/cell/>
        SELECT ?x ?y
        WHERE {
            ?point a cell:SpatialDataPoint ;
                  cell:hasXCoordinate ?x ;
                  cell:hasYCoordinate ?y .
        }
        LIMIT 10000
        """
        
        spatial_points = {
            'x': [float(row[0]) for row in g.query(spatial_query)],
            'y': [float(row[1]) for row in g.query(spatial_query)]
        }
        
        # Visualize cells
        visualize_cells(validated_cells, spatial_points, args.output_viz)
    
    print("Done!")

if __name__ == "__main__":
    main() 