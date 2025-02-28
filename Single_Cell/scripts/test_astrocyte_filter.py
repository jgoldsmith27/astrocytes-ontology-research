import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from read_data import read_data

def create_test_h5ad_with_astrocytes():
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

def test_astrocyte_filter():
    """Test the read_data function with astrocyte filtering"""
    print("\n🧪 TESTING ASTROCYTE FILTERING")
    
    # Create test file
    test_file = create_test_h5ad_with_astrocytes()
    
    try:
        # Test reading the file
        print("🧪 Testing read_data function with astrocyte filtering...")
        result_df = read_data(test_file)
        
        # Verify result is a DataFrame
        assert isinstance(result_df, pd.DataFrame), "Result is not a DataFrame"
        print("✅ Result is a DataFrame")
        
        # Verify only astrocyte types are in the index
        expected_types = ['Astrocyte_Type1', 'Astrocyte_Type2', 'Astrocyte_Type3']
        for cell_type in result_df.index:
            assert 'strocyte' in cell_type, f"Non-astrocyte type found: {cell_type}"
        print(f"✅ Only astrocyte types present: {list(result_df.index)}")
        
        # Verify gene names are in columns
        adata = sc.read_h5ad(test_file)
        assert all(gene in result_df.columns for gene in adata.var_names), "Not all genes in result"
        print(f"✅ All genes present: {result_df.shape[1]} genes")
        
        print(f"✅ TEST PASSED: Astrocyte filtering works correctly")
        return True
        
    except Exception as e:
        print(f"❌ TEST FAILED: {str(e)}")
        return False
    finally:
        # Clean up
        if os.path.exists(test_file):
            os.remove(test_file)
            print(f"🧹 Cleaned up test file")

if __name__ == "__main__":
    test_astrocyte_filter() 