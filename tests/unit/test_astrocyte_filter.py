#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import unittest
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

# Add the Single_Cell/scripts directory to the path to import the module
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
single_cell_scripts_dir = os.path.join(project_root, 'Single_Cell', 'scripts')
sys.path.insert(0, single_cell_scripts_dir)

# Import the module
from read_data import read_data

class TestAstrocyteFilter(unittest.TestCase):
    """
    Test cases for the astrocyte filtering functionality.
    """
    
    def setUp(self):
        """
        Set up test environment.
        """
        # Create test file
        self.test_file = self.create_test_h5ad_with_astrocytes()
    
    def tearDown(self):
        """
        Clean up test environment.
        """
        # Clean up
        if os.path.exists(self.test_file):
            os.remove(self.test_file)
            print(f"🧹 Cleaned up test file")
    
    def create_test_h5ad_with_astrocytes(self):
        """Create a test h5ad file with mock data including astrocyte subtypes"""
        print("🧪 Creating test h5ad file with astrocyte subtypes...")
        
        # Create mock data
        n_cells = 200
        n_genes = 50
        
        # Create expression matrix (sparse)
        X = sparse.csr_matrix(np.random.poisson(1, size=(n_cells, n_genes)))
        
        # Create cell metadata with astrocyte subtypes
        cell_types = [
            'Astrocyte_Type1', 'Astrocyte_Type2', 'Astrocyte_Type3',
            'Neuron', 'Microglia', 'Oligodendrocyte'
        ]
        weights = [0.2, 0.15, 0.15, 0.2, 0.15, 0.15]  # 50% astrocytes
        
        obs = pd.DataFrame({
            'celltype': np.random.choice(cell_types, n_cells, p=weights)
        })
        
        # Create gene metadata
        var = pd.DataFrame(index=[f'gene_{i}' for i in range(n_genes)])
        
        # Create AnnData object
        adata = ad.AnnData(X=X, obs=obs, var=var)
        
        # Save to file
        test_file = "test_astrocytes.h5ad"
        adata.write(test_file)
        print(f"✅ Test file created: {test_file}")
        return test_file
    
    def test_astrocyte_filter(self):
        """Test the read_data function with astrocyte filtering"""
        print("\n🧪 TESTING ASTROCYTE FILTERING")
        
        # Test reading the file
        print("🧪 Testing read_data function with astrocyte filtering...")
        result_df = read_data(self.test_file)
        
        # Verify result is a DataFrame
        self.assertIsInstance(result_df, pd.DataFrame, "Result is not a DataFrame")
        print("✅ Result is a DataFrame")
        
        # Verify only astrocyte types are in the index
        for cell_type in result_df.index:
            self.assertIn('strocyte', cell_type, f"Non-astrocyte type found: {cell_type}")
        print(f"✅ Only astrocyte types present: {list(result_df.index)}")
        
        # Verify gene names are in columns
        adata = sc.read_h5ad(self.test_file)
        for gene in adata.var_names:
            self.assertIn(gene, result_df.columns, f"Gene {gene} not found in result")
        print(f"✅ All genes present: {result_df.shape[1]} genes")
        
        print(f"✅ TEST PASSED: Astrocyte filtering works correctly")

if __name__ == "__main__":
    unittest.main() 