import rdflib
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, XSD, OWL
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

def load_turtle_file(ttl_file):
    """Load a TURTLE file into an RDF graph"""
    print(f"Loading TURTLE file: {ttl_file}")
    g = Graph()
    g.parse(ttl_file, format="turtle")
    print(f"Graph loaded with {len(g)} triples")
    return g

def explore_ontology(g):
    """Print basic information about the ontology"""
    # Define namespaces
    ST = Namespace("http://example.org/spatial-transcriptomics#")
    
    # Count classes
    classes = list(g.subjects(RDF.type, OWL.Class))
    print(f"\nClasses ({len(classes)}):")
    for cls in classes:
        print(f"  - {cls}")
    
    # Count properties
    obj_props = list(g.subjects(RDF.type, OWL.ObjectProperty))
    data_props = list(g.subjects(RDF.type, OWL.DatatypeProperty))
    print(f"\nObject Properties ({len(obj_props)}):")
    for prop in obj_props:
        print(f"  - {prop}")
    print(f"\nDatatype Properties ({len(data_props)}):")
    for prop in data_props:
        print(f"  - {prop}")
    
    # Count instances by type
    for cls in classes:
        instances = list(g.subjects(RDF.type, cls))
        print(f"\nInstances of {cls}: {len(instances)}")
        if len(instances) > 0:
            print(f"  Example: {instances[0]}")

def query_genes(g, limit=10):
    """Query genes and their expression levels"""
    ST = Namespace("http://example.org/spatial-transcriptomics#")
    
    print(f"\nGenes and their expression levels (top {limit}):")
    
    # SPARQL query to get genes and their expression levels
    query = """
    PREFIX : <http://example.org/spatial-transcriptomics#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?gene ?geneID ?point ?x ?y ?expressionLevel
    WHERE {
        ?gene rdf:type :Gene ;
              :expressedAt ?point ;
              :hasExpressionLevel ?expressionLevel .
        
        # Extract gene ID from URI
        BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        
        # Extract coordinates from point URI
        BIND(REPLACE(STR(?point), "^.*Point_([0-9]+)_([0-9]+)$", "$1") AS ?x)
        BIND(REPLACE(STR(?point), "^.*Point_([0-9]+)_([0-9]+)$", "$2") AS ?y)
    }
    LIMIT %d
    """ % limit
    
    results = g.query(query)
    
    for row in results:
        print(f"  Gene: {row.geneID}, Position: ({row.x}, {row.y}), Expression Level: {row.expressionLevel}")
    
    return results

def query_spatial_relationships(g, limit=10):
    """Query spatial points and their nearby points"""
    ST = Namespace("http://example.org/spatial-transcriptomics#")
    
    print(f"\nSpatial relationships (top {limit}):")
    
    # SPARQL query to get points and their nearby points
    query = """
    PREFIX : <http://example.org/spatial-transcriptomics#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?point ?x ?y ?nearPoint ?nearX ?nearY
    WHERE {
        ?point rdf:type :SpatialPoint ;
               :locatedNear ?nearPoint .
        
        # Extract coordinates from point URIs
        BIND(REPLACE(STR(?point), "^.*Point_([0-9]+)_([0-9]+)$", "$1") AS ?x)
        BIND(REPLACE(STR(?point), "^.*Point_([0-9]+)_([0-9]+)$", "$2") AS ?y)
        BIND(REPLACE(STR(?nearPoint), "^.*Point_([0-9]+)_([0-9]+)$", "$1") AS ?nearX)
        BIND(REPLACE(STR(?nearPoint), "^.*Point_([0-9]+)_([0-9]+)$", "$2") AS ?nearY)
    }
    LIMIT %d
    """ % limit
    
    results = g.query(query)
    
    for row in results:
        print(f"  Point ({row.x}, {row.y}) is near Point ({row.nearX}, {row.nearY})")
    
    return results

def visualize_spatial_network(g, max_points=100):
    """Visualize the spatial network of points"""
    ST = Namespace("http://example.org/spatial-transcriptomics#")
    
    # SPARQL query to get points and their nearby points
    query = """
    PREFIX : <http://example.org/spatial-transcriptomics#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?point ?x ?y ?nearPoint ?nearX ?nearY
    WHERE {
        ?point rdf:type :SpatialPoint ;
               :locatedNear ?nearPoint .
        
        # Extract coordinates from point URIs
        BIND(REPLACE(STR(?point), "^.*Point_([0-9]+)_([0-9]+)$", "$1") AS ?x)
        BIND(REPLACE(STR(?point), "^.*Point_([0-9]+)_([0-9]+)$", "$2") AS ?y)
        BIND(REPLACE(STR(?nearPoint), "^.*Point_([0-9]+)_([0-9]+)$", "$1") AS ?nearX)
        BIND(REPLACE(STR(?nearPoint), "^.*Point_([0-9]+)_([0-9]+)$", "$2") AS ?nearY)
    }
    """
    
    results = g.query(query)
    
    # Convert to numeric values
    points = []
    edges = []
    
    for row in results:
        x1, y1 = int(row.x), int(row.y)
        x2, y2 = int(row.nearX), int(row.nearY)
        
        points.append((x1, y1))
        points.append((x2, y2))
        edges.append(((x1, y1), (x2, y2)))
    
    # Limit the number of points to visualize
    unique_points = list(set(points))
    if len(unique_points) > max_points:
        unique_points = unique_points[:max_points]
        edges = [e for e in edges if e[0] in unique_points and e[1] in unique_points]
    
    # Create the plot
    plt.figure(figsize=(12, 10))
    
    # Plot points
    x_coords = [p[0] for p in unique_points]
    y_coords = [p[1] for p in unique_points]
    plt.scatter(x_coords, y_coords, c='blue', alpha=0.7, s=30)
    
    # Plot edges
    for (x1, y1), (x2, y2) in edges:
        plt.plot([x1, x2], [y1, y2], 'k-', alpha=0.2, linewidth=0.5)
    
    plt.title(f'Spatial Network of Gene Expression Points (showing {len(unique_points)} points)')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.tight_layout()
    
    output_image = os.path.join(os.path.dirname(ttl_file), "spatial_network.png")
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"Network visualization saved to {output_image}")
    
    plt.show()

def find_gene_clusters(g, gene_id=None):
    """Find clusters of genes or a specific gene"""
    ST = Namespace("http://example.org/spatial-transcriptomics#")
    
    # If gene_id is provided, filter for that gene
    gene_filter = f"FILTER(CONTAINS(STR(?gene), '{gene_id}'))" if gene_id else ""
    
    # SPARQL query to get genes and their spatial relationships
    query = f"""
    PREFIX : <http://example.org/spatial-transcriptomics#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?gene ?geneID ?point ?x ?y ?expressionLevel ?nearPoint ?nearX ?nearY
    WHERE {{
        ?gene rdf:type :Gene ;
              :expressedAt ?point ;
              :hasExpressionLevel ?expressionLevel .
        
        ?point :locatedNear ?nearPoint .
        
        # Extract gene ID from URI
        BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        
        # Extract coordinates from point URIs
        BIND(REPLACE(STR(?point), "^.*Point_([0-9]+)_([0-9]+)$", "$1") AS ?x)
        BIND(REPLACE(STR(?point), "^.*Point_([0-9]+)_([0-9]+)$", "$2") AS ?y)
        BIND(REPLACE(STR(?nearPoint), "^.*Point_([0-9]+)_([0-9]+)$", "$1") AS ?nearX)
        BIND(REPLACE(STR(?nearPoint), "^.*Point_([0-9]+)_([0-9]+)$", "$2") AS ?nearY)
        
        {gene_filter}
    }}
    """
    
    results = g.query(query)
    
    # Process results
    gene_points = {}
    for row in results:
        gene_id = str(row.geneID)
        if gene_id not in gene_points:
            gene_points[gene_id] = []
        
        gene_points[gene_id].append({
            'x': int(row.x),
            'y': int(row.y),
            'expression': int(row.expressionLevel),
            'near_x': int(row.nearX),
            'near_y': int(row.nearY)
        })
    
    # Print summary
    print(f"\nFound {len(gene_points)} genes with spatial relationships:")
    for gene_id, points in gene_points.items():
        print(f"  Gene {gene_id}: {len(points)} spatial relationships")
    
    return gene_points

if __name__ == "__main__":
    # File paths
    ttl_file = "/Users/jacob/Desktop/Karolinska Research/astrocytes-ontology-research/Spatial_Data/spatial_transcriptomics_sample.ttl"
    
    # Load the TURTLE file
    g = load_turtle_file(ttl_file)
    
    # Explore the ontology
    explore_ontology(g)
    
    # Query genes and their expression levels
    query_genes(g, limit=10)
    
    # Query spatial relationships
    query_spatial_relationships(g, limit=10)
    
    # Find gene clusters
    gene_clusters = find_gene_clusters(g)
    
    # Visualize the spatial network
    visualize_spatial_network(g, max_points=100) 