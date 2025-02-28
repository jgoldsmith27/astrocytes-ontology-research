import scanpy as sc
import anndata as ad
import logging

logger = logging.getLogger("generate_h5ad")
logger.setLevel(logging.INFO)

def build_h5ad(file_path, output_h5ad):
    """Simplified version that just copies the data"""
    logger.info(f"Reading data from {file_path}")
    adata = sc.read_h5ad(file_path)
    adata.write(output_h5ad)
    logger.info(f"Saved basic data to {output_h5ad}")

if __name__ == "__main__":
    build_h5ad(
        "/Users/jacob/Desktop/Karolinska Research/astrocytes-ontology-research/single_cell_data_cleaned.h5ad",
        "coexpression_results.h5ad"
    )
