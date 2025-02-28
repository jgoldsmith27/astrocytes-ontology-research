import scanpy as sc
import pandas as pd
import numpy as np

def read_data(file_path):
    """Read h5ad file and return average expression for astrocyte subtypes only"""
    print("🔍 Reading h5ad file and filtering for astrocyte subtypes...")
    adata = sc.read_h5ad(file_path)
    
    # Filter for astrocyte cell types only
    astrocyte_types = [ct for ct in adata.obs['celltype'].unique() if 'strocyte' in ct]
    
    if len(astrocyte_types) == 0:
        print("⚠️ No astrocyte types found! Using first three cell types as a fallback.")
        astrocyte_types = list(adata.obs['celltype'].unique()[:3])
    
    print(f"🧠 Found astrocyte types: {astrocyte_types}")
    
    # Calculate average expression per astrocyte type
    avg_df = pd.DataFrame()
    for cell_type in astrocyte_types:
        cell_mask = adata.obs['celltype'] == cell_type
        if sum(cell_mask) > 0:  # Only process if cells exist for this type
            avg_expression = adata[cell_mask].X.mean(axis=0)
            if hasattr(avg_expression, 'A1'):
                avg_df[cell_type] = avg_expression.A1
            else:
                avg_df[cell_type] = avg_expression
    
    avg_df.index = adata.var_names
    return avg_df.T  # Transpose to have cell types as rows