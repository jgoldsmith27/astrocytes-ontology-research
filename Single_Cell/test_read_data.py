import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from read_data import read_data

def create_test_h5ad():
    """Create a test h5ad file with mock data"""
    print("🧪 Creating test h5ad file...")
    
    # Create mock data
    n_cells = 100
    n_genes = 50
    
    # Create expression matrix (sparse)
    X = sparse.csr_matrix(np.random.poisson(1, size=(n_cells, n_genes)))
    
    # Create cell metadata
    obs = pd.DataFrame({
        'celltype': np.random.choice(['Neuron', 'Astrocyte', 'Microglia'], n_cells)
    })
    
    # Create gene metadata
    var = pd.DataFrame(index=[f'gene_{i}' for i in range(n_genes)])
    
    # Create AnnData object
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    # Save to file
    test_file = "test_data.h5ad"
    adata.write(test_file)
    print(f"✅ Test file created: {test_file}")
    return test_file

def test_read_data():
    """Test the read_data function"""
    print("\n🧪 TESTING READ_DATA FUNCTION")
    
    # Create test file
    test_file = create_test_h5ad()
    
    try:
        # Test reading the file
        print("🧪 Testing read_data function...")
        result_df = read_data(test_file)
        
        # Verify result is a DataFrame
        assert isinstance(result_df, pd.DataFrame), "Result is not a DataFrame"
        print("✅ Result is a DataFrame")
        
        # Verify cell types are in the index
        adata = sc.read_h5ad(test_file)
        expected_cell_types = np.unique(adata.obs['celltype'])
        assert all(ct in result_df.index for ct in expected_cell_types), "Not all cell types in result"
        print(f"✅ All cell types present: {list(result_df.index)}")
        
        # Verify gene names are in columns
        assert all(gene in result_df.columns for gene in adata.var_names), "Not all genes in result"
        print(f"✅ All genes present: {result_df.shape[1]} genes")
        
        # Verify values are reasonable
        assert not np.any(np.isnan(result_df.values)), "NaN values in result"
        print("✅ No NaN values in result")
        
        print(f"✅ TEST PASSED: read_data works correctly")
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
    test_read_data() 