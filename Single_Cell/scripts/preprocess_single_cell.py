import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import logging

logger = logging.getLogger("preprocess_single_cell")
logger.setLevel(logging.INFO)

def fix_reserved_index_names(df, new_name):
    """Ensure that the DataFrame's index and any column named '_index' are renamed."""
    if df.index.name is None or df.index.name == "_index":
        df.rename_axis(new_name, inplace=True)
    if "_index" in df.columns:
        df.rename(columns={"_index": new_name + "_col"}, inplace=True)

def preprocess_single_cell(file_path, output_path):
    """Preprocess single-cell data and ensure data integrity for downstream analysis."""
    try:
        logger.info(f"Loading single-cell data from: {file_path}")
        adata = sc.read_h5ad(file_path)

        # Ensure valid structure before proceeding
        if "_index" in adata.obs.columns:
            logger.warning("Renaming '_index' column to 'cell_index' to prevent conflicts.")
            adata.obs.rename(columns={"_index": "cell_index"}, inplace=True)
        
        # Debugging: Check if cells exist before filtering
        print(f"Before filtering: {adata.shape[0]} cells, {adata.shape[1]} genes")

        # Standard Preprocessing
        sc.pp.filter_cells(adata, min_genes=50)  # Lower threshold to prevent data loss
        sc.pp.filter_genes(adata, min_cells=1)  # Keep as many genes as possible

        # Debugging: Check if cells remain after filtering
        print(f"After filtering: {adata.shape[0]} cells, {adata.shape[1]} genes")
        
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Ensure adata.X is numeric (float64) and remove NaN/Inf values
        if sp.issparse(adata.X):
            logger.info("adata.X is sparse. Converting data to float64 and removing NaN/Inf.")
            adata.X.data = np.nan_to_num(adata.X.data.astype(np.float64),
                                          nan=0.0, posinf=0.0, neginf=0.0)
        else:
            logger.info("adata.X is dense. Converting data to float64 and removing NaN/Inf.")
            adata.X = np.nan_to_num(adata.X.astype(np.float64),
                                     nan=0.0, posinf=0.0, neginf=0.0)
        
        # Dimensionality Reduction
        sc.pp.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
        sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
        sc.tl.umap(adata)

        # Ensure adata.obs and adata.var are preserved as copies
        adata.obs = adata.obs.copy()
        adata.var = adata.var.copy()

        # Remove any reserved '_index' columns, if present
        if "_index" in adata.obs.columns:
            print("Removing reserved column '_index' from adata.obs...")
            adata.obs = adata.obs.drop(columns=["_index"], errors='ignore')
        if "_index" in adata.var.columns:
            print("Removing reserved column '_index' from adata.var...")
            adata.var = adata.var.drop(columns=["_index"], errors='ignore')

        # Debugging: Print obs and var before saving
        print("Final adata.obs before saving:", adata.obs.head())
        print("Final adata.var before saving:", adata.var.head())

        # Marker-Based Annotation (Example)
        marker_genes = {
            "Neurons": ["RBFOX3", "MAP2", "SLC17A7"],
            "Astrocytes": ["GFAP", "AQP4"],
            "Microglia": ["ITGAM", "P2RY12"],
            "Oligodendrocytes": ["MBP", "PLP1"],
            "Endothelial": ["CLDN5", "FLT1"]
        }
        adata.obs["Predicted_Cell_Type"] = "Unknown"
        for cell_type, genes in marker_genes.items():
            for gene in genes:
                if gene in adata.var_names:
                    gene_data = adata[:, gene].X
                    if sp.issparse(gene_data):
                        gene_data = gene_data.toarray().ravel()
                    else:
                        gene_data = np.asarray(gene_data, dtype=np.float64).ravel()
                    mask = gene_data > 1
                    adata.obs.loc[mask, "Predicted_Cell_Type"] = cell_type

        # Validate Before Saving
        if adata.obs.empty:
            raise ValueError("adata.obs is missing or empty before saving!")
        if adata.var.empty:
            raise ValueError("adata.var is missing or empty before saving!")

        # Rename the index axes to non-reserved names in the main dataframes.
        fix_reserved_index_names(adata.obs, "cell_id")
        fix_reserved_index_names(adata.var, "gene_id")

        # Also fix raw if it exists.
        if adata.raw is not None:
            # Check and fix raw.obs if it exists.
            if hasattr(adata.raw, 'obs') and isinstance(adata.raw.obs, pd.DataFrame):
                fix_reserved_index_names(adata.raw.obs, "raw_cell_id")
                print("adata.raw.obs index name:", adata.raw.obs.index.name)
            # raw should have var.
            if hasattr(adata.raw, 'var') and isinstance(adata.raw.var, pd.DataFrame):
                fix_reserved_index_names(adata.raw.var, "raw_gene_id")
                print("adata.raw.var index name:", adata.raw.var.index.name)

        print("adata.obs index name:", adata.obs.index.name)
        print("adata.var index name:", adata.var.index.name)

        # Save Processed Data
        adata.write(output_path)
        logger.info(f"Preprocessed data saved successfully to: {output_path}")

    except Exception as e:
        logger.error(f"Error in preprocessing: {e}")

if __name__ == "__main__":
    preprocess_single_cell(
        "../../wetransfer_ctr081_fron-h5ad_2024-02-12_1046/CTR081_Fron.h5ad",
        "../data/single_cell_data_cleaned.h5ad"
    )