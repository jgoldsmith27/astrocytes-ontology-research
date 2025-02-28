#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import unittest
import tempfile
import shutil
import pandas as pd
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

class TestAstrocyteValidationRegression(unittest.TestCase):
    """
    Regression tests for the AstrocyteValidator class.
    These tests verify that previously fixed bugs remain fixed.
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
    
    def test_missing_point_uri_regression(self):
        """
        Test that the validator handles missing point_uri correctly.
        This was a bug that was fixed previously.
        """
        # Extract cells and validate cell types
        self.validator.extract_classified_cells()
        self.validator.validate_cell_types()
        
        # Create a copy of the data without point_uri
        df_no_uri = self.validator.spatial_cells_df.copy()
        df_no_uri.drop('point_uri', axis=1, inplace=True)
        
        # Save the original and replace with the modified version
        original_df = self.validator.spatial_cells_df
        self.validator.spatial_cells_df = df_no_uri
        
        try:
            # This should not raise an exception
            self.validator.visualize_validation_results()
            
            # Restore the original dataframe
            self.validator.spatial_cells_df = original_df
            
            # Test passed if no exception was raised
            self.assertTrue(True)
        except Exception as e:
            # Restore the original dataframe
            self.validator.spatial_cells_df = original_df
            self.fail(f"visualize_validation_results raised an exception with missing point_uri: {str(e)}")
    
    def test_empty_coexpression_networks_regression(self):
        """
        Test that the validator handles empty co-expression networks correctly.
        This was a bug that was fixed previously.
        """
        # Extract cells
        self.validator.extract_classified_cells()
        
        # Set empty co-expression networks
        self.validator.coexpression_networks = {}
        
        try:
            # This should not raise an exception
            self.validator.visualize_coexpression_networks()
            
            # Test passed if no exception was raised
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"visualize_coexpression_networks raised an exception with empty networks: {str(e)}")
    
    def test_ambiguous_truth_value_regression(self):
        """
        Test that the validator handles ambiguous truth values correctly.
        This was a bug that was fixed previously.
        """
        # Extract cells
        self.validator.extract_classified_cells()
        
        # Create a co-expression network with NaN values
        import numpy as np
        df = pd.DataFrame({
            'gene1': [1.0, np.nan, 0.5],
            'gene2': [np.nan, 1.0, 0.6],
            'gene3': [0.5, 0.6, 1.0]
        }, index=['gene1', 'gene2', 'gene3'])
        
        # Set the co-expression networks
        self.validator.coexpression_networks = {'Protoplasmic': df}
        
        try:
            # This should not raise an exception
            self.validator.visualize_coexpression_networks()
            
            # Test passed if no exception was raised
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"visualize_coexpression_networks raised an exception with NaN values: {str(e)}")
    
    def test_spatial_validation_results_image_creation(self):
        """
        Test that the spatial_validation_results.png file is created.
        This was a bug that was fixed previously.
        """
        # Run the validation
        self.validator.run_validation()
        
        # Check that the file exists
        output_dir = self.validator.output_dir
        filepath = os.path.join(output_dir, "spatial_validation_results.png")
        self.assertTrue(os.path.exists(filepath), f"Expected output file {filepath} not found")

if __name__ == "__main__":
    unittest.main() 