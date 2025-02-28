import rdflib
from rdflib import Graph, Namespace, URIRef, Literal, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL
import numpy as np
import argparse
import os
import math

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")

class SpatialOntologyEnhancer:
    """A class to enhance spatial relationships in the ontology"""
    
    def __init__(self, ttl_file):
        """Initialize with a TURTLE file"""
        self.ttl_file = ttl_file
        self.graph = self.load_turtle_file(ttl_file)
        self.output_dir = os.path.dirname(ttl_file)
        self.output_file = os.path.join(self.output_dir, "enhanced_spatial_ontology.ttl")
    
    def load_turtle_file(self, ttl_file):
        """Load a TURTLE file into an RDF graph"""
        print(f"Loading TURTLE file: {ttl_file}")
        g = Graph()
        g.parse(ttl_file, format="turtle")
        print(f"Graph loaded with {len(g)} triples")
        return g
    
    def add_directional_relationships(self):
        """Add directional spatial relationships between points"""
        print("\nAdding directional spatial relationships...")
        
        # Define directional properties in the ontology
        self.graph.add((ST.directlyAbove, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.directlyBelow, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.directlyLeftOf, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.directlyRightOf, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.diagonallyRelatedTo, RDF.type, OWL.ObjectProperty))
        
        # Add property descriptions
        self.graph.add((ST.directlyAbove, RDFS.label, Literal("directly above")))
        self.graph.add((ST.directlyBelow, RDFS.label, Literal("directly below")))
        self.graph.add((ST.directlyLeftOf, RDFS.label, Literal("directly left of")))
        self.graph.add((ST.directlyRightOf, RDFS.label, Literal("directly right of")))
        self.graph.add((ST.diagonallyRelatedTo, RDFS.label, Literal("diagonally related to")))
        
        # Define inverse relationships
        self.graph.add((ST.directlyAbove, OWL.inverseOf, ST.directlyBelow))
        self.graph.add((ST.directlyLeftOf, OWL.inverseOf, ST.directlyRightOf))
        
        # Query for all spatial points
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y
        WHERE {
            ?point rdf:type :SpatialPoint ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
        }
        """
        
        results = self.graph.query(query)
        
        # Extract points and coordinates
        points = []
        for row in results:
            points.append({
                'uri': row.point,
                'x': int(row.x),
                'y': int(row.y)
            })
        
        # Define thresholds for directional relationships
        direction_threshold = 50  # Units for considering points as directly above/below/left/right
        
        # Add directional relationships
        relationships_added = 0
        
        for i, point1 in enumerate(points):
            for point2 in points[i+1:]:
                dx = point2['x'] - point1['x']
                dy = point2['y'] - point1['y']
                
                # Check for directional relationships
                if abs(dx) < direction_threshold and dy > direction_threshold:
                    # point2 is directly above point1
                    self.graph.add((point1['uri'], ST.directlyBelow, point2['uri']))
                    self.graph.add((point2['uri'], ST.directlyAbove, point1['uri']))
                    relationships_added += 2
                
                elif abs(dx) < direction_threshold and dy < -direction_threshold:
                    # point2 is directly below point1
                    self.graph.add((point1['uri'], ST.directlyAbove, point2['uri']))
                    self.graph.add((point2['uri'], ST.directlyBelow, point1['uri']))
                    relationships_added += 2
                
                elif abs(dy) < direction_threshold and dx > direction_threshold:
                    # point2 is directly right of point1
                    self.graph.add((point1['uri'], ST.directlyLeftOf, point2['uri']))
                    self.graph.add((point2['uri'], ST.directlyRightOf, point1['uri']))
                    relationships_added += 2
                
                elif abs(dy) < direction_threshold and dx < -direction_threshold:
                    # point2 is directly left of point1
                    self.graph.add((point1['uri'], ST.directlyRightOf, point2['uri']))
                    self.graph.add((point2['uri'], ST.directlyLeftOf, point1['uri']))
                    relationships_added += 2
                
                elif abs(dx) > direction_threshold and abs(dy) > direction_threshold:
                    # points are diagonally related
                    self.graph.add((point1['uri'], ST.diagonallyRelatedTo, point2['uri']))
                    self.graph.add((point2['uri'], ST.diagonallyRelatedTo, point1['uri']))
                    relationships_added += 2
        
        print(f"Added {relationships_added} directional relationships")
    
    def define_spatial_regions(self):
        """Define spatial regions and associate points with them"""
        print("\nDefining spatial regions...")
        
        # Define region class and properties
        self.graph.add((ST.SpatialRegion, RDF.type, OWL.Class))
        self.graph.add((ST.hasRegionName, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasMinX, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasMaxX, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasMinY, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.hasMaxY, RDF.type, OWL.DatatypeProperty))
        self.graph.add((ST.locatedInRegion, RDF.type, OWL.ObjectProperty))
        
        # Add property descriptions
        self.graph.add((ST.SpatialRegion, RDFS.label, Literal("Spatial Region")))
        self.graph.add((ST.hasRegionName, RDFS.label, Literal("has region name")))
        self.graph.add((ST.hasMinX, RDFS.label, Literal("has minimum X coordinate")))
        self.graph.add((ST.hasMaxX, RDFS.label, Literal("has maximum X coordinate")))
        self.graph.add((ST.hasMinY, RDFS.label, Literal("has minimum Y coordinate")))
        self.graph.add((ST.hasMaxY, RDFS.label, Literal("has maximum Y coordinate")))
        self.graph.add((ST.locatedInRegion, RDFS.label, Literal("located in region")))
        
        # Query for all spatial points to determine coordinate ranges
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y
        WHERE {
            ?point rdf:type :SpatialPoint ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
        }
        """
        
        results = self.graph.query(query)
        
        # Extract coordinates
        x_coords = []
        y_coords = []
        points = []
        
        for row in results:
            x = int(row.x)
            y = int(row.y)
            x_coords.append(x)
            y_coords.append(y)
            points.append({
                'uri': row.point,
                'x': x,
                'y': y
            })
        
        if not x_coords or not y_coords:
            print("No spatial points found")
            return
        
        # Determine coordinate ranges
        min_x, max_x = min(x_coords), max(x_coords)
        min_y, max_y = min(y_coords), max(y_coords)
        
        # Define regions (divide the space into a 3x3 grid)
        x_ranges = np.linspace(min_x, max_x, 4)
        y_ranges = np.linspace(min_y, max_y, 4)
        
        regions = []
        for i in range(3):
            for j in range(3):
                region_name = f"Region_{i+1}_{j+1}"
                region_uri = ST[region_name]
                
                # Add region to ontology
                self.graph.add((region_uri, RDF.type, ST.SpatialRegion))
                self.graph.add((region_uri, ST.hasRegionName, Literal(region_name)))
                self.graph.add((region_uri, ST.hasMinX, Literal(x_ranges[i], datatype=XSD.integer)))
                self.graph.add((region_uri, ST.hasMaxX, Literal(x_ranges[i+1], datatype=XSD.integer)))
                self.graph.add((region_uri, ST.hasMinY, Literal(y_ranges[j], datatype=XSD.integer)))
                self.graph.add((region_uri, ST.hasMaxY, Literal(y_ranges[j+1], datatype=XSD.integer)))
                
                regions.append({
                    'uri': region_uri,
                    'name': region_name,
                    'min_x': x_ranges[i],
                    'max_x': x_ranges[i+1],
                    'min_y': y_ranges[j],
                    'max_y': y_ranges[j+1]
                })
        
        print(f"Defined {len(regions)} spatial regions")
        
        # Associate points with regions
        associations_added = 0
        
        for point in points:
            for region in regions:
                if (point['x'] >= region['min_x'] and point['x'] <= region['max_x'] and
                    point['y'] >= region['min_y'] and point['y'] <= region['max_y']):
                    self.graph.add((point['uri'], ST.locatedInRegion, region['uri']))
                    associations_added += 1
                    break
        
        print(f"Associated {associations_added} points with regions")
    
    def add_distance_relationships(self):
        """Add more granular distance-based relationships"""
        print("\nAdding distance-based relationships...")
        
        # Define distance properties
        self.graph.add((ST.veryCloseToPoint, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.moderatelyCloseToPoint, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.farFromPoint, RDF.type, OWL.ObjectProperty))
        
        # Add property descriptions
        self.graph.add((ST.veryCloseToPoint, RDFS.label, Literal("very close to point")))
        self.graph.add((ST.moderatelyCloseToPoint, RDFS.label, Literal("moderately close to point")))
        self.graph.add((ST.farFromPoint, RDFS.label, Literal("far from point")))
        
        # Define distance thresholds
        very_close_threshold = 50
        moderately_close_threshold = 200
        
        # Query for all spatial points
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y
        WHERE {
            ?point rdf:type :SpatialPoint ;
                   :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
        }
        """
        
        results = self.graph.query(query)
        
        # Extract points and coordinates
        points = []
        for row in results:
            points.append({
                'uri': row.point,
                'x': int(row.x),
                'y': int(row.y)
            })
        
        # Add distance relationships
        relationships_added = 0
        
        # Limit the number of points to process to avoid excessive relationships
        max_points = 1000
        if len(points) > max_points:
            print(f"Limiting analysis to {max_points} points to avoid excessive relationships")
            points = points[:max_points]
        
        for i, point1 in enumerate(points):
            for point2 in points[i+1:]:
                # Calculate Euclidean distance
                distance = math.sqrt((point1['x'] - point2['x'])**2 + (point1['y'] - point2['y'])**2)
                
                # Add appropriate relationship based on distance
                if distance <= very_close_threshold:
                    self.graph.add((point1['uri'], ST.veryCloseToPoint, point2['uri']))
                    self.graph.add((point2['uri'], ST.veryCloseToPoint, point1['uri']))
                    relationships_added += 2
                elif distance <= moderately_close_threshold:
                    self.graph.add((point1['uri'], ST.moderatelyCloseToPoint, point2['uri']))
                    self.graph.add((point2['uri'], ST.moderatelyCloseToPoint, point1['uri']))
                    relationships_added += 2
                else:
                    self.graph.add((point1['uri'], ST.farFromPoint, point2['uri']))
                    self.graph.add((point2['uri'], ST.farFromPoint, point1['uri']))
                    relationships_added += 2
        
        print(f"Added {relationships_added} distance-based relationships")
    
    def define_spatial_patterns(self):
        """Define spatial patterns for gene expression"""
        print("\nDefining spatial expression patterns...")
        
        # Define pattern classes
        self.graph.add((ST.SpatialPattern, RDF.type, OWL.Class))
        self.graph.add((ST.ClusteredPattern, RDF.type, OWL.Class))
        self.graph.add((ST.DispersedPattern, RDF.type, OWL.Class))
        self.graph.add((ST.GradientPattern, RDF.type, OWL.Class))
        
        # Add subclass relationships
        self.graph.add((ST.ClusteredPattern, RDFS.subClassOf, ST.SpatialPattern))
        self.graph.add((ST.DispersedPattern, RDFS.subClassOf, ST.SpatialPattern))
        self.graph.add((ST.GradientPattern, RDFS.subClassOf, ST.SpatialPattern))
        
        # Define pattern property
        self.graph.add((ST.exhibitsSpatialPattern, RDF.type, OWL.ObjectProperty))
        self.graph.add((ST.exhibitsSpatialPattern, RDFS.label, Literal("exhibits spatial pattern")))
        self.graph.add((ST.exhibitsSpatialPattern, RDFS.domain, ST.Gene))
        self.graph.add((ST.exhibitsSpatialPattern, RDFS.range, ST.SpatialPattern))
        
        # Query for genes and their expression points
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?geneID ?point ?x ?y
        WHERE {
            ?gene rdf:type :Gene ;
                  :expressedAt ?point .
            
            ?point :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
            
            # Extract gene ID
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        """
        
        results = self.graph.query(query)
        
        # Group expression points by gene
        gene_points = {}
        for row in results:
            gene_uri = row.gene
            gene_id = str(row.geneID)
            x = int(row.x)
            y = int(row.y)
            
            if gene_uri not in gene_points:
                gene_points[gene_uri] = {
                    'id': gene_id,
                    'points': []
                }
            
            gene_points[gene_uri]['points'].append((x, y))
        
        # Analyze spatial patterns for each gene
        patterns_added = 0
        
        for gene_uri, data in gene_points.items():
            points = data['points']
            
            if len(points) < 5:
                continue  # Skip genes with too few expression points
            
            # Calculate mean and standard deviation of distances between points
            distances = []
            for i, p1 in enumerate(points):
                for j, p2 in enumerate(points[i+1:], i+1):
                    dist = math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
                    distances.append(dist)
            
            mean_dist = np.mean(distances) if distances else 0
            std_dist = np.std(distances) if distances else 0
            
            # Calculate spatial distribution metrics
            x_coords = [p[0] for p in points]
            y_coords = [p[1] for p in points]
            
            x_range = max(x_coords) - min(x_coords) if x_coords else 0
            y_range = max(y_coords) - min(y_coords) if y_coords else 0
            
            # Determine pattern based on metrics
            if std_dist < 100 and mean_dist < 200:
                # Points are close together with low variance = clustered
                pattern = ST.ClusteredPattern
            elif std_dist > 300 and mean_dist > 400:
                # Points are far apart with high variance = dispersed
                pattern = ST.DispersedPattern
            elif (x_range > 1000 or y_range > 1000) and std_dist < 200:
                # Points span a large area but with consistent spacing = gradient
                pattern = ST.GradientPattern
            else:
                continue  # Skip if no clear pattern
            
            # Create pattern instance
            pattern_instance = ST[f"Pattern_{data['id']}"]
            self.graph.add((pattern_instance, RDF.type, pattern))
            
            # Associate gene with pattern
            self.graph.add((gene_uri, ST.exhibitsSpatialPattern, pattern_instance))
            patterns_added += 1
        
        print(f"Added {patterns_added} spatial pattern classifications")
    
    def enhance_ontology(self):
        """Enhance the ontology with additional spatial relationships"""
        print("\nEnhancing spatial ontology...")
        
        # Add directional relationships
        self.add_directional_relationships()
        
        # Define spatial regions
        self.define_spatial_regions()
        
        # Add distance relationships
        self.add_distance_relationships()
        
        # Define spatial patterns
        self.define_spatial_patterns()
        
        # Save enhanced ontology
        self.graph.serialize(destination=self.output_file, format="turtle")
        print(f"\nEnhanced ontology saved to {self.output_file}")
        print(f"Total triples in enhanced ontology: {len(self.graph)}")

def main():
    parser = argparse.ArgumentParser(description='Enhance Spatial Ontology')
    parser.add_argument('--ttl', default="../data/spatial_transcriptomics_advanced.ttl",
                        help='Path to TURTLE file')
    args = parser.parse_args()
    
    enhancer = SpatialOntologyEnhancer(args.ttl)
    enhancer.enhance_ontology()

if __name__ == "__main__":
    main() 