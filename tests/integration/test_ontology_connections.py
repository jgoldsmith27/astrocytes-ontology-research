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

# Add the Spatial_Data/scripts directory to the path to import the module
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
spatial_scripts_dir = os.path.join(project_root, 'Spatial_Data', 'scripts')
sys.path.insert(0, spatial_scripts_dir)

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
        self.project_dir = project_root
        self.spatial_dir = os.path.join(self.project_dir, "Spatial_Data")
        self.single_cell_dir = os.path.join(self.project_dir, "Single_Cell")
        self.output_dir = os.path.join(self.project_dir, "output")
        self.test_output_dir = os.path.join(self.project_dir, "tests", "output")
        
        # Create output directories if they don't exist
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.test_output_dir, exist_ok=True)
        
        # Define ontology file paths
        self.spatial_ontology = os.path.join(self.test_output_dir, "test_spatial_ontology.ttl")
        self.single_cell_ontology = os.path.join(self.test_output_dir, "test_single_cell_ontology.ttl")
        self.integrated_ontology = os.path.join(self.test_output_dir, "test_integrated_ontology.ttl")
        
        # Create test ontologies
        self._create_test_ontologies()
        
        # Load ontologies
        self.spatial_graph = Graph()
        self.single_cell_graph = Graph()
        self.integrated_graph = Graph()
        
        self.spatial_graph.parse(self.spatial_ontology, format="turtle")
        print(f"Loaded spatial ontology with {len(self.spatial_graph)} triples")
        
        self.single_cell_graph.parse(self.single_cell_ontology, format="turtle")
        print(f"Loaded single-cell ontology with {len(self.single_cell_graph)} triples")
        
        # Run integration
        self._run_integration()
        self.integrated_graph.parse(self.integrated_ontology, format="turtle")
        print(f"Created and loaded integrated ontology with {len(self.integrated_graph)} triples")
    
    def tearDown(self):
        """
        Clean up test files.
        """
        # Remove test ontology files
        for file_path in [self.spatial_ontology, self.single_cell_ontology, self.integrated_ontology]:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"Removed test file: {file_path}")
    
    def _create_test_ontologies(self):
        """
        Create test ontologies for testing.
        """
        # Create spatial ontology
        self._create_test_spatial_ontology()
        
        # Create single-cell ontology
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
        print(f"Created test spatial ontology with {len(g)} triples")
    
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
        print(f"Created test single-cell ontology with {len(g)} triples")
    
    def _run_integration(self):
        """
        Run the ontology integration script to create the integrated ontology.
        """
        try:
            # Import the ontology integration module
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
            self.fail("Could not import ontology_integration module. Make sure it exists in the Spatial_Data/scripts directory.")
        except Exception as e:
            self.fail(f"Error running integration: {str(e)}")
    
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
        
        # Check if astrocyte types exist
        astrocyte_types_query = """
        PREFIX sc: <http://example.org/single-cell#>
        
        SELECT DISTINCT ?cellType (COUNT(?cell) as ?count)
        WHERE {
            ?cell sc:hasCellType ?cellType .
        }
        GROUP BY ?cellType
        """
        
        results = list(self.single_cell_graph.query(astrocyte_types_query))
        self.assertGreater(len(results), 0, "No cell types found in the single-cell ontology")
        
        print("Cell types in single-cell ontology:")
        for row in results:
            print(f"  {row[0]}: {row[1]} cells")
    
    def test_integrated_ontology_connections(self):
        """
        Test that the integrated ontology has connections between spatial and single-cell data.
        """
        # Check if bridge connections exist
        bridge_query = """
        PREFIX bridge: <http://example.org/ontology-bridge#>
        
        SELECT (COUNT(?connection) as ?count)
        WHERE {
            ?connection ?p ?o .
            FILTER(STRSTARTS(STR(?connection), STR(bridge:)))
        }
        """
        
        results = list(self.integrated_graph.query(bridge_query))
        count = int(results[0][0]) if results else 0
        
        self.assertGreater(count, 0, "No bridge connections found in the integrated ontology")
        print(f"Found {count} bridge connections in the integrated ontology")
        
        # Check if spatial points are connected to single-cell data
        connection_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        
        SELECT (COUNT(DISTINCT ?point) as ?count)
        WHERE {
            ?point a st:SpatialPoint .
            ?cell a sc:Cell .
            ?point ?p ?cell .
        }
        """
        
        results = list(self.integrated_graph.query(connection_query))
        count = int(results[0][0]) if results else 0
        
        self.assertGreater(count, 0, "No spatial points connected to single-cell data")
        print(f"Found {count} spatial points connected to single-cell data")

if __name__ == "__main__":
    unittest.main() 