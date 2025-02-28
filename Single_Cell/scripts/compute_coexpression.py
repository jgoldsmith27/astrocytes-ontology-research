import numpy as np
import pandas as pd
import scipy.sparse as sp
import logging
from scipy.stats import spearmanr

logger = logging.getLogger("compute_coexpression")
logger.setLevel(logging.INFO)

def compute_coexpression(adata, cell_types, threshold=0.5):
    """
    Computes gene-gene coexpression networks for each specified cell type.
    Saves the coexpression results in an .h5ad file.
    """
    coexpression_data = {}
    
    for cell_type in cell_types:
        logger.info(f"Processing {cell_type}")
        subtype_adata = adata[adata.obs["celltype"] == cell_type]
        
        # Compute average gene expression across all cells of this subtype
        avg_expression = np.mean(subtype_adata.X.toarray(), axis=0)
        gene_names = subtype_adata.var_names
        
        # Compute Spearman correlation matrix
        correlation_matrix = np.zeros((len(gene_names), len(gene_names)))
        for i in range(len(gene_names)):
            for j in range(i+1, len(gene_names)):
                corr, _ = spearmanr(subtype_adata.X[:, i].toarray().ravel(),
                                    subtype_adata.X[:, j].toarray().ravel())
                correlation_matrix[i, j] = corr
                correlation_matrix[j, i] = corr
        
        # Store results
        coexpression_data[cell_type] = pd.DataFrame(correlation_matrix, index=gene_names, columns=gene_names)
    
    return coexpression_data


def save_coexpression_to_h5ad(coexpression_data, output_path="coexpression_results.h5ad"):
    """
    Saves coexpression data into an AnnData object (.h5ad file) for further analysis.
    """
    layers = {}
    for cell_type, df in coexpression_data.items():
        layers[cell_type] = sp.csr_matrix(df.values)  # Store as sparse matrix
    
    coexpression_adata = ad.AnnData(X=np.zeros((1,1)))  # Placeholder matrix
    coexpression_adata.layers = layers
    coexpression_adata.var_names = coexpression_data[next(iter(coexpression_data))].columns
    coexpression_adata.write(output_path)
    logger.info(f"Saved coexpression data to {output_path}")
