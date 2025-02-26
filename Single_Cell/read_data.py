import anndata as ad
import pandas as pd

def read_data(file_path: str) -> pd.DataFrame:
    """
    Reads an h5ad file, groups data by cell type, and returns a DataFrame 
    with the average gene expression for each cell type.

    Args:
        file_path (str): Path to the h5ad file.

    Returns:
        pd.DataFrame: DataFrame where rows are cell types, columns are genes, 
                      and values are average expression levels for each gene in the cell type.
    """
    # Load h5ad file
    adata = ad.read_h5ad(file_path)

    # Ensure 'celltype' column exists
    if "celltype" not in adata.obs:
        raise ValueError("No 'celltype' column found in metadata. Ensure your dataset has cell type annotations.")

    # Convert data matrix to DataFrame
    df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)

    # Attach cell type labels
    df["celltype"] = adata.obs["celltype"].values

    # Group by cell type and compute the mean expression per gene
    avg_df = df.groupby("celltype").mean()

    return avg_df