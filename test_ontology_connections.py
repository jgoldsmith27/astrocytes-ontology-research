#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import unittest
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, OWL
import pandas as pd
import numpy as np
import argparse

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
SC = Namespace("http://example.org/single-cell#")
CELL = Namespace("http://example.org/cell-ontology#")
GENE = Namespace("http://example.org/gene-ontology#")
BRIDGE = Namespace("http://example.org/ontology-bridge#")

class OntologyConnectionTest(unittest.TestCase):
    """
    Test suite to verify the connections between spatial and single-cell ontologies.
    """
    
    def setUp(self):
        """
        Set up the test environment by creating paths and checking for ontology files.
        """
        # Define paths
        self.project_dir = os.path.dirname(os.path.abspath(__file__))
        self.spatial_dir = os.path.join(self.project_dir, "Spatial_Data")
        self.single_cell_dir = os.path.join(self.project_dir, "Single_Cell")
        self.output_dir = os.path.join(self.project_dir, "output")
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Define ontology file paths
        self.spatial_ontology = os.path.join(self.spatial_dir, "output", "enhanced_spatial_ontology.ttl")
        self.single_cell_ontology = os.path.join(self.single_cell_dir, "output", "single_cell_ontology.ttl")
        self.integrated_ontology = os.path.join(self.output_dir, "integrated_ontology.ttl")
        
        # Check if ontology files exist, if not, create sample ontologies for testing
        self._create_test_ontologies_if_needed()
        
        # Load ontologies
        self.spatial_graph = Graph()
        self.single_cell_graph = Graph()
        self.integrated_graph = Graph()
        
        if os.path.exists(self.spatial_ontology):
            self.spatial_graph.parse(self.spatial_ontology, format="turtle")
            print(f"Loaded spatial ontology with {len(self.spatial_graph)} triples")
        
        if os.path.exists(self.single_cell_ontology):
            self.single_cell_graph.parse(self.single_cell_ontology, format="turtle")
            print(f"Loaded single-cell ontology with {len(self.single_cell_graph)} triples")
        
        if os.path.exists(self.integrated_ontology):
            self.integrated_graph.parse(self.integrated_ontology, format="turtle")
            print(f"Loaded integrated ontology with {len(self.integrated_graph)} triples")
        else:
            # If integrated ontology doesn't exist, create it by running the integration script
            self._run_integration()
            if os.path.exists(self.integrated_ontology):
                self.integrated_graph.parse(self.integrated_ontology, format="turtle")
                print(f"Created and loaded integrated ontology with {len(self.integrated_graph)} triples")
    
    def _create_test_ontologies_if_needed(self):
        """
        Create test ontologies if the real ones don't exist.
        """
        # Create directories if they don't exist
        os.makedirs(os.path.dirname(self.spatial_ontology), exist_ok=True)
        os.makedirs(os.path.dirname(self.single_cell_ontology), exist_ok=True)
        
        # Create spatial ontology if it doesn't exist
        if not os.path.exists(self.spatial_ontology):
            print(f"Creating test spatial ontology at {self.spatial_ontology}")
            self._create_test_spatial_ontology()
        
        # Create single-cell ontology if it doesn't exist
        if not os.path.exists(self.single_cell_ontology):
            print(f"Creating test single-cell ontology at {self.single_cell_ontology}")
            self._create_test_single_cell_ontology()
    
    def _create_test_spatial_ontology(self):
        """
        Create a test spatial ontology with sample data.
        """
        g = Graph()
        
        # Bind namespaces
        g.bind("st", ST)
        g.bind("cell", CELL)
        g.bind("rdf", RDF)
        g.bind("rdfs", RDFS)
        g.bind("owl", OWL)
        
        # Create spatial points
        for i in range(1, 101):
            point_uri = URIRef(ST + f"SpatialPoint_{i}")
            g.add((point_uri, RDF.type, ST.SpatialPoint))
            g.add((point_uri, ST.hasXCoordinate, Literal(np.random.uniform(0, 100))))
            g.add((point_uri, ST.hasYCoordinate, Literal(np.random.uniform(0, 100))))
            
            # Add astrocyte type
            type_node = URIRef(ST + f"AstrocyteType_{i}")
            g.add((point_uri, ST.hasAstrocyteType, type_node))
            
            # Randomly assign astrocyte type
            astrocyte_types = ["Protoplasmic", "Fibrous", "Reactive"]
            astrocyte_type = np.random.choice(astrocyte_types)
            g.add((type_node, ST.astrocyteType, Literal(astrocyte_type)))
            g.add((type_node, ST.typeProbability, Literal(np.random.uniform(0.7, 1.0))))
            
            # Add gene expression
            for j in range(1, 6):  # Add 5 genes per point
                gene_uri = URIRef(ST + f"Gene_{j}")
                g.add((gene_uri, RDF.type, ST.Gene))
                g.add((gene_uri, ST.hasGeneID, Literal(f"GENE{j}")))
                g.add((gene_uri, ST.expressedAt, point_uri))
                g.add((gene_uri, ST.hasExpressionLevel, Literal(np.random.uniform(0, 3.0))))
        
        # Save the ontology
        g.serialize(destination=self.spatial_ontology, format="turtle")
    
    def _create_test_single_cell_ontology(self):
        """
        Create a test single-cell ontology with sample data.
        """
        g = Graph()
        
        # Bind namespaces
        g.bind("sc", SC)
        g.bind("cell", CELL)
        g.bind("rdf", RDF)
        g.bind("rdfs", RDFS)
        g.bind("owl", OWL)
        
        # Create cells
        for i in range(1, 101):
            cell_uri = URIRef(SC + f"Cell_{i}")
            g.add((cell_uri, RDF.type, SC.Cell))
            
            # Randomly assign astrocyte type
            astrocyte_types = ["Protoplasmic", "Fibrous", "Reactive"]
            astrocyte_type = np.random.choice(astrocyte_types)
            g.add((cell_uri, SC.hasCellType, Literal(astrocyte_type)))
            
            # Add gene expression
            for j in range(1, 6):  # Add 5 genes per cell
                gene_uri = URIRef(SC + f"Gene_{j}")
                g.add((gene_uri, RDF.type, SC.Gene))
                g.add((gene_uri, SC.hasGeneID, Literal(f"GENE{j}")))
                
                expr_uri = URIRef(SC + f"Expression_{i}_{j}")
                g.add((expr_uri, SC.belongsToCell, cell_uri))
                g.add((expr_uri, SC.forGene, gene_uri))
                g.add((expr_uri, SC.hasExpressionLevel, Literal(np.random.uniform(0, 3.0))))
        
        # Save the ontology
        g.serialize(destination=self.single_cell_ontology, format="turtle")
    
    def _run_integration(self):
        """
        Run the ontology integration script to create the integrated ontology.
        """
        try:
            # Import the ontology integration module
            sys.path.append(os.path.join(self.spatial_dir, "scripts"))
            from ontology_integration import OntologyIntegration
            
            # Create integration tool
            integration_tool = OntologyIntegration(
                spatial_ttl=self.spatial_ontology,
                single_cell_ttl=self.single_cell_ontology,
                output_ttl=self.integrated_ontology
            )
            
            # Run pipeline
            integration_tool.run_pipeline()
            
            print(f"Integration completed, output saved to {self.integrated_ontology}")
            
        except ImportError:
            print("Could not import ontology_integration module. Make sure it exists in the Spatial_Data/scripts directory.")
            return False
        except Exception as e:
            print(f"Error running integration: {str(e)}")
            return False
        
        return True
    
    def test_spatial_ontology_structure(self):
        """
        Test that the spatial ontology has the expected structure.
        """
        # Check if spatial points exist
        spatial_points_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT (COUNT(?point) as ?count)
        WHERE {
            ?point rdf:type st:SpatialPoint .
        }
        """
        
        results = list(self.spatial_graph.query(spatial_points_query))
        count = int(results[0][0]) if results else 0
        
        self.assertGreater(count, 0, "No spatial points found in the spatial ontology")
        print(f"Found {count} spatial points in the spatial ontology")
        
        # Check if astrocyte types exist
        astrocyte_types_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        
        SELECT DISTINCT ?astrocyteType (COUNT(?point) as ?count)
        WHERE {
            ?point st:hasAstrocyteType ?typeNode .
            ?typeNode st:astrocyteType ?astrocyteType .
        }
        GROUP BY ?astrocyteType
        """
        
        results = list(self.spatial_graph.query(astrocyte_types_query))
        self.assertGreater(len(results), 0, "No astrocyte types found in the spatial ontology")
        
        print("Astrocyte types in spatial ontology:")
        for row in results:
            print(f"  {row[0]}: {row[1]} points")
        
        # Check if genes exist
        genes_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT (COUNT(DISTINCT ?gene) as ?count)
        WHERE {
            ?gene rdf:type st:Gene .
        }
        """
        
        results = list(self.spatial_graph.query(genes_query))
        count = int(results[0][0]) if results else 0
        
        self.assertGreater(count, 0, "No genes found in the spatial ontology")
        print(f"Found {count} genes in the spatial ontology")
    
    def test_single_cell_ontology_structure(self):
        """
        Test that the single-cell ontology has the expected structure.
        """
        # Check if cells exist
        cells_query = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT (COUNT(?cell) as ?count)
        WHERE {
            ?cell rdf:type sc:Cell .
        }
        """
        
        results = list(self.single_cell_graph.query(cells_query))
        count = int(results[0][0]) if results else 0
        
        self.assertGreater(count, 0, "No cells found in the single-cell ontology")
        print(f"Found {count} cells in the single-cell ontology")
        
        # Check if cell types exist
        cell_types_query = """
        PREFIX sc: <http://example.org/single-cell#>
        
        SELECT DISTINCT ?cellType (COUNT(?cell) as ?count)
        WHERE {
            ?cell sc:hasCellType ?cellType .
        }
        GROUP BY ?cellType
        """
        
        results = list(self.single_cell_graph.query(cell_types_query))
        self.assertGreater(len(results), 0, "No cell types found in the single-cell ontology")
        
        print("Cell types in single-cell ontology:")
        for row in results:
            print(f"  {row[0]}: {row[1]} cells")
        
        # Check if genes exist
        genes_query = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT (COUNT(DISTINCT ?gene) as ?count)
        WHERE {
            ?gene rdf:type sc:Gene .
        }
        """
        
        results = list(self.single_cell_graph.query(genes_query))
        count = int(results[0][0]) if results else 0
        
        self.assertGreater(count, 0, "No genes found in the single-cell ontology")
        print(f"Found {count} genes in the single-cell ontology")
    
    def test_integrated_ontology_connections(self):
        """
        Test that the integrated ontology has the expected connections between entities.
        """
        # Check if bridge properties exist
        bridge_properties_query = """
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        
        SELECT ?property
        WHERE {
            ?property rdf:type owl:ObjectProperty .
            FILTER(STRSTARTS(STR(?property), STR(bridge:)))
        }
        """
        
        results = list(self.integrated_graph.query(bridge_properties_query))
        self.assertGreater(len(results), 0, "No bridge properties found in the integrated ontology")
        
        print("Bridge properties in integrated ontology:")
        for row in results:
            print(f"  {row[0]}")
        
        # Check if gene mappings exist
        gene_mappings_query = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT (COUNT(?scGene) as ?count)
        WHERE {
            ?scGene rdf:type sc:Gene .
            ?scGene bridge:hasCorrespondingGene ?stGene .
            ?stGene rdf:type st:Gene .
        }
        """
        
        results = list(self.integrated_graph.query(gene_mappings_query))
        count = int(results[0][0]) if results else 0
        
        self.assertGreater(count, 0, "No gene mappings found in the integrated ontology")
        print(f"Found {count} gene mappings in the integrated ontology")
        
        # Check if cell type mappings exist
        cell_type_mappings_query = """
        PREFIX bridge: <http://example.org/ontology-bridge#>
        
        SELECT ?scType ?stType
        WHERE {
            ?scType bridge:hasCorrespondingCellType ?stType .
        }
        """
        
        results = list(self.integrated_graph.query(cell_type_mappings_query))
        self.assertGreater(len(results), 0, "No cell type mappings found in the integrated ontology")
        
        print("Cell type mappings in integrated ontology:")
        for row in results:
            print(f"  {row[0]} -> {row[1]}")
        
        # Check if spatial inferences exist
        spatial_inferences_query = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT (COUNT(?cell) as ?count)
        WHERE {
            ?cell rdf:type sc:Cell .
            ?cell bridge:hasSpatialRepresentation ?point .
        }
        """
        
        results = list(self.integrated_graph.query(spatial_inferences_query))
        count = int(results[0][0]) if results else 0
        
        self.assertGreater(count, 0, "No spatial inferences found in the integrated ontology")
        print(f"Found {count} spatial inferences in the integrated ontology")
    
    def test_cross_ontology_queries(self):
        """
        Test that cross-ontology queries work as expected.
        """
        # Query 1: Find genes expressed in both modalities
        query1 = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?geneID (COUNT(?scGene) as ?count)
        WHERE {
            ?scGene rdf:type sc:Gene ;
                    sc:hasGeneID ?geneID .
            
            ?scGene bridge:hasCorrespondingGene ?stGene .
        }
        GROUP BY ?geneID
        ORDER BY DESC(?count)
        LIMIT 10
        """
        
        results = list(self.integrated_graph.query(query1))
        self.assertGreater(len(results), 0, "No genes found in both modalities")
        
        print("Genes expressed in both modalities:")
        for row in results:
            print(f"  {row[0]}: {row[1]} occurrences")
        
        # Query 2: Find spatial regions with specific cell types
        query2 = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?cellType (COUNT(?point) as ?count)
        WHERE {
            ?point rdf:type st:SpatialPoint ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?cellType .
            
            # Get cells with this type
            ?cell sc:hasCellType ?cellType ;
                  bridge:hasSpatialRepresentation ?point .
        }
        GROUP BY ?cellType
        ORDER BY DESC(?count)
        """
        
        results = list(self.integrated_graph.query(query2))
        self.assertGreater(len(results), 0, "No spatial regions found with specific cell types")
        
        print("Spatial regions with specific cell types:")
        for row in results:
            print(f"  {row[0]}: {row[1]} regions")
    
    def test_missing_components(self):
        """
        Test to identify missing components in the spatial ontology.
        """
        # Check for missing spatial coordinates
        missing_coords_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point
        WHERE {
            ?point rdf:type st:SpatialPoint .
            FILTER NOT EXISTS { ?point st:hasXCoordinate ?x }
            FILTER NOT EXISTS { ?point st:hasYCoordinate ?y }
        }
        LIMIT 10
        """
        
        results = list(self.spatial_graph.query(missing_coords_query))
        self.assertEqual(len(results), 0, f"Found {len(results)} spatial points without coordinates")
        
        # Check for missing astrocyte types
        missing_types_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point
        WHERE {
            ?point rdf:type st:SpatialPoint .
            FILTER NOT EXISTS { ?point st:hasAstrocyteType ?type }
        }
        LIMIT 10
        """
        
        results = list(self.spatial_graph.query(missing_types_query))
        self.assertEqual(len(results), 0, f"Found {len(results)} spatial points without astrocyte types")
        
        # Check for missing gene expression
        missing_expr_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point
        WHERE {
            ?point rdf:type st:SpatialPoint .
            FILTER NOT EXISTS { ?gene st:expressedAt ?point }
        }
        LIMIT 10
        """
        
        results = list(self.spatial_graph.query(missing_expr_query))
        if len(results) > 0:
            print(f"Warning: Found {len(results)} spatial points without gene expression data")
        
        # Check for missing brain region information
        missing_region_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point
        WHERE {
            ?point rdf:type st:SpatialPoint .
            FILTER NOT EXISTS { ?point st:inBrainRegion ?region }
        }
        LIMIT 10
        """
        
        results = list(self.spatial_graph.query(missing_region_query))
        if len(results) > 0:
            print(f"Warning: Found {len(results)} spatial points without brain region information")
            print("Missing brain region information in the spatial ontology")
        
        # Check for missing cell morphology information
        missing_morphology_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?typeNode
        WHERE {
            ?point st:hasAstrocyteType ?typeNode .
            FILTER NOT EXISTS { ?typeNode st:hasMorphology ?morphology }
        }
        LIMIT 10
        """
        
        results = list(self.spatial_graph.query(missing_morphology_query))
        if len(results) > 0:
            print(f"Warning: Found {len(results)} astrocyte types without morphology information")
            print("Missing morphology information in the spatial ontology")
        
        # Check for missing connectivity information
        missing_connectivity_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point
        WHERE {
            ?point rdf:type st:SpatialPoint .
            FILTER NOT EXISTS { ?point st:connectsTo ?otherPoint }
        }
        LIMIT 10
        """
        
        results = list(self.spatial_graph.query(missing_connectivity_query))
        if len(results) > 0:
            print(f"Warning: Found {len(results)} spatial points without connectivity information")
            print("Missing connectivity information in the spatial ontology")
    
    def test_recommendations(self):
        """
        Test to provide recommendations for improving the ontologies.
        """
        recommendations = []
        
        # Check if brain regions are defined
        brain_regions_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT (COUNT(?region) as ?count)
        WHERE {
            ?region rdf:type st:BrainRegion .
        }
        """
        
        results = list(self.spatial_graph.query(brain_regions_query))
        count = int(results[0][0]) if results else 0
        
        if count == 0:
            recommendations.append("Add brain region definitions to the spatial ontology")
        
        # Check if morphology types are defined
        morphology_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT (COUNT(?morphology) as ?count)
        WHERE {
            ?morphology rdf:type st:Morphology .
        }
        """
        
        results = list(self.spatial_graph.query(morphology_query))
        count = int(results[0][0]) if results else 0
        
        if count == 0:
            recommendations.append("Add astrocyte morphology definitions to the spatial ontology")
        
        # Check if connectivity information is defined
        connectivity_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        
        SELECT (COUNT(?point) as ?count)
        WHERE {
            ?point st:connectsTo ?otherPoint .
        }
        """
        
        results = list(self.spatial_graph.query(connectivity_query))
        count = int(results[0][0]) if results else 0
        
        if count == 0:
            recommendations.append("Add connectivity information between spatial points")
        
        # Check if temporal information is defined
        temporal_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        
        SELECT (COUNT(?point) as ?count)
        WHERE {
            ?point st:hasTimepoint ?timepoint .
        }
        """
        
        results = list(self.spatial_graph.query(temporal_query))
        count = int(results[0][0]) if results else 0
        
        if count == 0:
            recommendations.append("Add temporal information to the spatial ontology")
        
        # Check if functional annotations are defined
        functional_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        
        SELECT (COUNT(?type) as ?count)
        WHERE {
            ?type st:hasFunctionalAnnotation ?function .
        }
        """
        
        results = list(self.spatial_graph.query(functional_query))
        count = int(results[0][0]) if results else 0
        
        if count == 0:
            recommendations.append("Add functional annotations to astrocyte types")
        
        # Print recommendations
        if recommendations:
            print("\nRecommendations for improving the spatial ontology:")
            for i, rec in enumerate(recommendations, 1):
                print(f"{i}. {rec}")
        else:
            print("\nNo recommendations for improving the spatial ontology")

def main():
    """
    Main function to run the ontology connection tests.
    """
    parser = argparse.ArgumentParser(description='Test ontology connections')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    
    args = parser.parse_args()
    
    # Configure unittest
    if args.verbose:
        unittest.main(argv=['first-arg-is-ignored'], verbosity=2)
    else:
        unittest.main(argv=['first-arg-is-ignored'])

if __name__ == "__main__":
    main() 