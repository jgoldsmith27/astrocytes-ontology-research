#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import unittest
import tempfile
import shutil
from rdflib import Graph, Namespace, URIRef, Literal, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL

# Add the Spatial_Data/scripts directory to the path to import the module
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
spatial_scripts_dir = os.path.join(project_root, 'Spatial_Data', 'scripts')
sys.path.insert(0, spatial_scripts_dir)

# Import the module
from astrocyte_validation import AstrocyteValidator, CANONICAL_MARKERS

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
SC = Namespace("http://example.org/single-cell#")
BRIDGE = Namespace("http://example.org/ontology-bridge#")

class TestAstrocyteValidator(unittest.TestCase):
    """
    Test cases for the AstrocyteValidator class.
    """
    
    def setUp(self):
        """
        Set up test environment.
        """
        # Create a temporary directory for test outputs
        self.test_dir = tempfile.mkdtemp()
        
        # Create a test ontology
        self.test_ttl = os.path.join(self.test_dir, "test_ontology.ttl")
        self.create_test_ontology()
        
        # Initialize validator
        self.validator = AstrocyteValidator(
            integrated_ttl=self.test_ttl,
            output_dir=os.path.join(self.test_dir, "validation_output")
        )
    
    def tearDown(self):
        """
        Clean up test environment.
        """
        # Remove temporary directory
        shutil.rmtree(self.test_dir)
    
    def create_test_ontology(self):
        """
        Create a test ontology with sample data.
        """
        g = Graph()
        
        # Add namespaces
        g.bind("st", ST)
        g.bind("sc", SC)
        g.bind("bridge", BRIDGE)
        g.bind("rdf", RDF)
        g.bind("rdfs", RDFS)
        g.bind("xsd", XSD)
        g.bind("owl", OWL)
        
        # Create spatial points with astrocyte types
        astrocyte_types = ["Protoplasmic", "Fibrous", "Reactive"]
        
        for i in range(30):
            # Create spatial point
            point_uri = URIRef(ST + f"SpatialPoint_{i}")
            g.add((point_uri, RDF.type, ST.SpatialPoint))
            g.add((point_uri, ST.hasXCoordinate, Literal(i % 10, datatype=XSD.float)))
            g.add((point_uri, ST.hasYCoordinate, Literal(i // 10, datatype=XSD.float)))
            
            # Assign astrocyte type
            type_node = BNode()
            astrocyte_type = astrocyte_types[i % 3]
            g.add((point_uri, ST.hasAstrocyteType, type_node))
            g.add((type_node, ST.astrocyteType, Literal(astrocyte_type)))
            g.add((type_node, ST.typeProbability, Literal(0.8, datatype=XSD.float)))
            
            # Add gene expression
            # Add marker genes for the assigned type with high expression
            for marker in CANONICAL_MARKERS[astrocyte_type][:3]:  # Use first 3 markers
                gene_uri = URIRef(ST + f"Gene_{marker}")
                g.add((gene_uri, RDF.type, ST.Gene))
                g.add((gene_uri, ST.hasGeneID, Literal(marker)))
                g.add((gene_uri, ST.expressedAt, point_uri))
                g.add((gene_uri, ST.hasExpressionLevel, Literal(5.0, datatype=XSD.float)))
                g.add((point_uri, ST.hasGeneExpression, gene_uri))
            
            # Add other genes with lower expression
            for j in range(5):
                gene_uri = URIRef(ST + f"Gene_OTHER{j}")
                g.add((gene_uri, RDF.type, ST.Gene))
                g.add((gene_uri, ST.hasGeneID, Literal(f"OTHER{j}")))
                g.add((gene_uri, ST.expressedAt, point_uri))
                g.add((gene_uri, ST.hasExpressionLevel, Literal(1.0, datatype=XSD.float)))
                g.add((point_uri, ST.hasGeneExpression, gene_uri))
        
        # Save ontology
        g.serialize(destination=self.test_ttl, format="turtle")
        print(f"Created test ontology with {len(g)} triples")
    
    def test_extract_classified_cells(self):
        """
        Test extraction of classified cells.
        """
        self.validator.extract_classified_cells()
        
        # Check that cells were extracted
        self.assertIsNotNone(self.validator.spatial_cells_df)
        self.assertEqual(len(self.validator.spatial_cells_df), 30)
        
        # Check that gene expression matrix was created
        self.assertIsNotNone(self.validator.gene_expression_matrix)
        
        # Check that gene columns include marker genes
        gene_cols = [col for col in self.validator.gene_expression_matrix.columns 
                    if col not in ['point_uri', 'x', 'y', 'astrocyte_type', 'probability']]
        
        # Check that at least some marker genes are present
        markers_present = False
        for cell_type, markers in CANONICAL_MARKERS.items():
            for marker in markers[:3]:  # We added the first 3 markers for each type
                if marker in gene_cols:
                    markers_present = True
                    break
        
        self.assertTrue(markers_present, "No marker genes found in gene expression matrix")
    
    def test_compute_coexpression_networks(self):
        """
        Test computation of co-expression networks.
        """
        self.validator.extract_classified_cells()
        self.validator.compute_coexpression_networks()
        
        # Check that networks were created
        self.assertGreater(len(self.validator.coexpression_networks), 0)
        
        # Check network properties for at least one cell type
        for cell_type, network in self.validator.coexpression_networks.items():
            # Check that network is a DataFrame
            self.assertIsNotNone(network)
            
            # Check that diagonal is 1.0 (self-correlation)
            for i in range(min(3, len(network))):
                self.assertAlmostEqual(network.iloc[i, i], 1.0)
            
            # Only need to check one network
            break
    
    def test_validate_cell_types(self):
        """
        Test validation of cell types.
        """
        self.validator.extract_classified_cells()
        self.validator.validate_cell_types()
        
        # Check that validation results were created
        self.assertIsNotNone(self.validator.validation_results)
        self.assertEqual(len(self.validator.validation_results), 30)
        
        # Check that validation results include required columns
        required_cols = ['cell_uri', 'assigned_type', 'marker_based_type', 'is_valid']
        for col in required_cols:
            self.assertIn(col, self.validator.validation_results.columns)
        
        # Since we designed the test data to have high expression of marker genes
        # for the assigned cell type, most cells should be correctly classified
        valid_count = sum(self.validator.validation_results['is_valid'])
        self.assertGreater(valid_count / len(self.validator.validation_results), 0.7)
    
    def test_full_pipeline(self):
        """
        Test the full validation pipeline.
        """
        results = self.validator.run_validation()
        
        # Check that results include required keys
        required_keys = ['total_cells', 'valid_cells', 'valid_percent', 'results_path']
        for key in required_keys:
            self.assertIn(key, results)
        
        # Check that output files were created
        output_dir = self.validator.output_dir
        expected_files = [
            "validation_results.csv",
            "validation_confusion_matrix.png",
            "marker_expression_by_type.png",
            "spatial_validation_results.png"
        ]
        
        for filename in expected_files:
            filepath = os.path.join(output_dir, filename)
            self.assertTrue(os.path.exists(filepath), f"Expected output file {filepath} not found")

if __name__ == "__main__":
    unittest.main() 