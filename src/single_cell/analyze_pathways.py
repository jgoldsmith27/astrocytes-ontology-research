import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Set
import os
import scipy.sparse as sp
import gc

### DOES NOT WORK ### ALL CELLS ARE AUTOCRINE AND WE NEED TO TAKE AVERAGE OR SET THRESHOLD FOR COUNT BECAUSE WE WILL ALWAYS FIND INTERACTIONS AS SINGLE CELL DATA IS NOISY ###

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get the project root directory (two levels up from this script)
PROJECT_ROOT = Path(__file__).parent.parent.parent

def load_interaction_data() -> pd.DataFrame:
    """Load the interaction data from CellChatDB."""
    interaction_path = PROJECT_ROOT / "data/cellchat/processed/interaction.csv"
    if not interaction_path.exists():
        raise FileNotFoundError(f"Interaction data not found at {interaction_path}")
    
    logger.info("Loading interaction data...")
    interactions = pd.read_csv(interaction_path)
    return interactions

def load_single_cell_data(h5ad_path: str) -> sc.AnnData:
    """Load the single-cell data from h5ad file."""
    full_path = PROJECT_ROOT / h5ad_path
    if not full_path.exists():
        raise FileNotFoundError(f"Single-cell data not found at {full_path}")
    
    logger.info(f"Loading single-cell data from {full_path}...")
    adata = sc.read_h5ad(full_path)
    
    # Print available columns in obs
    logger.info("Available columns in single-cell data:")
    for col in adata.obs.columns:
        logger.info(f"- {col}")
    
    return adata

def extract_ligand_receptor_pairs(interactions: pd.DataFrame) -> List[Tuple[str, str]]:
    """Extract unique ligand-receptor pairs from interaction data."""
    logger.info("Extracting ligand-receptor pairs...")
    
    # Create ligand-receptor pairs
    pairs = []
    for _, row in interactions.iterrows():
        ligand = row['ligand']
        # Split receptor string by underscore to handle multiple receptors
        receptors = row['receptor'].split('_')
        for receptor in receptors:
            pairs.append((ligand, receptor))
    
    # Remove duplicates
    unique_pairs = list(set(pairs))
    logger.info(f"Found {len(unique_pairs)} unique ligand-receptor pairs")
    return unique_pairs

def get_gene_expression(adata: sc.AnnData, gene: str, start_idx: int, end_idx: int) -> sp.csr_matrix:
    """Get binary expression vector for a gene for a specific chunk of cells."""
    if gene not in adata.var_names:
        return None
    expr = adata.X[start_idx:end_idx, adata.var_names.get_loc(gene)]
    if sp.issparse(expr):
        return expr > 0
    return sp.csr_matrix(expr > 0)

def analyze_pathway_chunk(
    pathway: str,
    pathway_interactions: pd.DataFrame,
    adata: sc.AnnData,
    cell_type_col: str,
    cell_types: List[str],
    chunk_size: int = 500
) -> List[Dict]:
    """Analyze ligand-receptor pairs for all cell types, processing in chunks."""
    results = []
    
    # Process each ligand-receptor pair separately
    for _, row in pathway_interactions.iterrows():
        ligand = row['ligand']
        receptors = row['receptor'].split('_')
        
        # Skip if ligand is not in data
        if ligand not in adata.var_names:
            continue
        
        # Process cells in chunks
        n_cells = adata.shape[0]
        n_chunks = (n_cells + chunk_size - 1) // chunk_size
        
        for chunk_idx in range(n_chunks):
            chunk_start = chunk_idx * chunk_size
            chunk_end = min(chunk_start + chunk_size, n_cells)
            
            if chunk_idx % 10 == 0:  # Log progress every 10 chunks
                logger.debug(f"Processing chunk {chunk_idx + 1}/{n_chunks} for {ligand} -> {receptors}")
            
            # Get cell type information for this chunk
            chunk_cell_types = adata.obs[cell_type_col].iloc[chunk_start:chunk_end]
            
            # Get ligand expression for chunk
            ligand_idx = adata.var_names.get_loc(ligand)
            ligand_expr = adata.X[chunk_start:chunk_end, ligand_idx]
            if sp.issparse(ligand_expr):
                ligand_expr = ligand_expr > 0
            else:
                ligand_expr = sp.csr_matrix(ligand_expr > 0)
            
            # Process each receptor
            for receptor in receptors:
                if receptor not in adata.var_names:
                    continue
                
                # Get receptor expression for chunk
                receptor_idx = adata.var_names.get_loc(receptor)
                receptor_expr = adata.X[chunk_start:chunk_end, receptor_idx]
                if sp.issparse(receptor_expr):
                    receptor_expr = receptor_expr > 0
                else:
                    receptor_expr = sp.csr_matrix(receptor_expr > 0)
                
                # Find cells expressing both
                expressing_cells = ligand_expr.multiply(receptor_expr)
                
                # Count by cell type
                for source_ct in cell_types:
                    source_mask = (chunk_cell_types == source_ct).values
                    if not np.any(source_mask):
                        continue
                        
                    for target_ct in cell_types:
                        target_mask = (chunk_cell_types == target_ct).values
                        if not np.any(target_mask):
                            continue
                        
                        # Count interactions using sparse matrix operations
                        count = expressing_cells[source_mask & target_mask].sum()
                        
                        if count > 0:
                            results.append({
                                'source_cell_type': source_ct,
                                'target_cell_type': target_ct,
                                'pathway': pathway,
                                'ligand': ligand,
                                'receptor': receptor,
                                'occurrences': int(count),
                                'is_autocrine': source_ct == target_ct
                            })
            
            # Clean up memory
            gc.collect()
    
    return results

def analyze_pathways_by_cell_type(
    adata: sc.AnnData,
    interactions: pd.DataFrame,
    chunk_size: int = 500
) -> pd.DataFrame:
    """Analyze ligand-receptor pairs by cell type, processing in chunks to save memory."""
    logger.info("Analyzing ligand-receptor pairs by cell type...")
    
    # Try different possible cell type column names
    cell_type_cols = ['celltype', 'cell_type', 'CellType', 'Celltype', 'cellType']
    cell_type_col = None
    
    for col in cell_type_cols:
        if col in adata.obs.columns:
            cell_type_col = col
            break
    
    if cell_type_col is None:
        raise KeyError("Could not find cell type column in single-cell data")
    
    # Get unique cell types
    cell_types = adata.obs[cell_type_col].unique()
    logger.info(f"Found {len(cell_types)} cell types")
    
    # Process each ligand-receptor pair
    all_results = []
    total_pairs = len(interactions)
    
    for i, (_, row) in enumerate(interactions.iterrows(), 1):
        ligand = row['ligand']
        receptor_complex = row['receptor']  # Keep the full receptor complex name
        receptors = receptor_complex.split('_')  # Split for checking expression
        pathway = row['pathway_name']
        
        logger.info(f"Processing pair {i}/{total_pairs}: {ligand} -> {receptor_complex} ({pathway})")
        
        # Skip if ligand is not in data
        if ligand not in adata.var_names:
            logger.debug(f"Skipping {ligand} - not found in data")
            continue
            
        # Skip if any receptor in the complex is not in data
        if not all(r in adata.var_names for r in receptors):
            missing = [r for r in receptors if r not in adata.var_names]
            logger.debug(f"Skipping {receptor_complex} - missing receptors: {missing}")
            continue
        
        # Process cells in chunks
        n_cells = adata.shape[0]
        n_chunks = (n_cells + chunk_size - 1) // chunk_size
        
        for chunk_idx in range(n_chunks):
            chunk_start = chunk_idx * chunk_size
            chunk_end = min(chunk_start + chunk_size, n_cells)
            
            if chunk_idx % 10 == 0:  # Log progress every 10 chunks
                logger.debug(f"Processing chunk {chunk_idx + 1}/{n_chunks} for {ligand} -> {receptor_complex}")
            
            # Get cell type information for this chunk
            chunk_cell_types = adata.obs[cell_type_col].iloc[chunk_start:chunk_end]
            
            # Get ligand expression for chunk
            ligand_idx = adata.var_names.get_loc(ligand)
            ligand_expr = adata.X[chunk_start:chunk_end, ligand_idx]
            if sp.issparse(ligand_expr):
                ligand_expr = ligand_expr > 0
            else:
                ligand_expr = sp.csr_matrix(ligand_expr > 0)
            
            # Get combined receptor expression (all receptors must be expressed)
            receptor_exprs = []
            for receptor in receptors:
                receptor_idx = adata.var_names.get_loc(receptor)
                receptor_expr = adata.X[chunk_start:chunk_end, receptor_idx]
                if sp.issparse(receptor_expr):
                    receptor_expr = receptor_expr > 0
                else:
                    receptor_expr = sp.csr_matrix(receptor_expr > 0)
                receptor_exprs.append(receptor_expr)
            
            # Combine receptor expressions (all must be expressed)
            combined_receptor_expr = receptor_exprs[0]
            for expr in receptor_exprs[1:]:
                combined_receptor_expr = combined_receptor_expr.multiply(expr)
            
            # Find cells expressing both ligand and all receptors
            expressing_cells = ligand_expr.multiply(combined_receptor_expr)
            
            # Count by cell type
            for source_ct in cell_types:
                source_mask = (chunk_cell_types == source_ct).values
                if not np.any(source_mask):
                    continue
                    
                for target_ct in cell_types:
                    target_mask = (chunk_cell_types == target_ct).values
                    if not np.any(target_mask):
                        continue
                    
                    # Count interactions using sparse matrix operations
                    count = expressing_cells[source_mask & target_mask].sum()
                    
                    if count > 0:
                        all_results.append({
                            'source_cell_type': source_ct,
                            'target_cell_type': target_ct,
                            'pathway': pathway,
                            'ligand': ligand,
                            'receptor_complex': receptor_complex,  # Full complex name
                            'receptors': receptors,  # Array of individual receptors
                            'n_receptors': len(receptors),  # Number of receptors in complex
                            'occurrences': int(count),
                            'is_autocrine': source_ct == target_ct
                        })
            
            # Clean up memory
            gc.collect()
        
        # Save intermediate results every 100 pairs
        if i % 100 == 0:
            temp_df = pd.DataFrame(all_results)
            if len(temp_df) > 0:
                output_path = PROJECT_ROOT / "data/single_cell/processed/pathway_analysis_temp.csv"
                temp_df.to_csv(output_path, index=False)
                logger.info(f"Saved intermediate results after {i} pairs")
    
    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)
    
    if len(results_df) > 0:
        # Sort by source cell type, target cell type, and occurrences
        results_df = results_df.sort_values(
            ['source_cell_type', 'target_cell_type', 'occurrences'],
            ascending=[True, True, False]
        )
    
    return results_df

def main():
    """Main function to run the pathway analysis."""
    try:
        # Load data
        interactions = load_interaction_data()
        adata = load_single_cell_data("data/single_cell/raw/CTR081_Fron.h5ad")
        
        # Analyze pathways
        results = analyze_pathways_by_cell_type(adata, interactions, chunk_size=500)
        
        # Save results
        output_path = PROJECT_ROOT / "data/single_cell/processed/pathway_analysis.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")
        
        # Print summary
        if len(results) > 0:
            logger.info("\nAnalysis Summary:")
            logger.info(f"Total pathways analyzed: {len(results['pathway'].unique())}")
            logger.info(f"Total cell types: {len(results['source_cell_type'].unique())}")
            logger.info(f"Total interactions: {len(results)}")
        else:
            logger.warning("No significant interactions found")
        
    except Exception as e:
        logger.error(f"Error in pathway analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main() 