import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Set
import scipy.sparse as sp
import gc

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent

def load_interaction_data() -> pd.DataFrame:
    """Load the interaction data from CellChatDB."""
    interaction_path = PROJECT_ROOT / "data/cellchat/processed/interaction.csv"
    if not interaction_path.exists():
        raise FileNotFoundError(f"Interaction data not found at {interaction_path}")
    
    logger.info("Loading interaction data...")
    interactions = pd.read_csv(interaction_path)
    
    # Validate required columns
    required_columns = ['ligand', 'receptor']
    missing_columns = [col for col in required_columns if col not in interactions.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in interaction data: {missing_columns}")
    
    # Validate no empty values in key columns
    if interactions['ligand'].isna().any() or interactions['receptor'].isna().any():
        raise ValueError("Found empty values in ligand or receptor columns")
    
    logger.info(f"Loaded {len(interactions)} interactions")
    return interactions

def load_spatial_data(h5ad_path: str) -> pd.DataFrame:
    """Load the spatial data from CSV file."""
    full_path = PROJECT_ROOT / h5ad_path
    if not full_path.exists():
        raise FileNotFoundError(f"Spatial data not found at {full_path}")
    
    logger.info(f"Loading spatial data from {full_path}...")
    # Read the CSV file
    df = pd.read_csv(full_path)
    
    # Validate required columns
    required_columns = ['geneID', 'bin1_ID']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in spatial data: {missing_columns}")
    
    # Validate no empty values in key columns
    if df['geneID'].isna().any() or df['bin1_ID'].isna().any():
        raise ValueError("Found empty values in geneID or bin1_ID columns")
    
    # Print available columns
    logger.info("Available columns in spatial data:")
    for col in df.columns:
        logger.info(f"- {col}")
    
    logger.info(f"Loaded {len(df)} rows of spatial data")
    return df

def extract_unique_genes(interactions: pd.DataFrame) -> Dict[str, Set[str]]:
    """Extract unique ligands and receptors from interaction data."""
    logger.info("\nExtracting unique ligands and receptors from CellChatDB interactions...")
    
    # Get unique ligands
    ligands = set(interactions['ligand'].str.upper())
    logger.info(f"Total unique ligands in CellChatDB: {len(ligands)}")
    logger.info("Sample of ligands (first 10):")
    for ligand in sorted(list(ligands))[:10]:
        logger.info(f"- {ligand}")
    
    # Get unique receptors (split complexes)
    all_receptors = interactions['receptor'].str.upper()
    receptors = set()
    complex_receptors = set()
    
    for receptor in all_receptors:
        if '_' in receptor:
            complex_receptors.add(receptor)
            # Split complex into individual genes
            receptors.update(receptor.split('_'))
        else:
            receptors.add(receptor)
    
    logger.info(f"\nTotal unique receptor genes in CellChatDB: {len(receptors)}")
    logger.info("Sample of receptor genes (first 10):")
    for receptor in sorted(list(receptors))[:10]:
        logger.info(f"- {receptor}")
    
    logger.info(f"\nTotal receptor complexes: {len(complex_receptors)}")
    logger.info("Sample of receptor complexes (first 10):")
    for complex_r in sorted(list(complex_receptors))[:10]:
        logger.info(f"- {complex_r}")
    
    # Validate no empty sets
    if not ligands or not receptors:
        raise ValueError("No ligands or receptors found in interaction data")
    
    return {
        'ligands': ligands,
        'receptors': receptors
    }

def analyze_gene_expression(
    df: pd.DataFrame,
    genes: Set[str],
    gene_type: str,
    chunk_size: int = 1000
) -> pd.DataFrame:
    """Analyze presence of genes in spatial spots."""
    logger.info(f"\nAnalyzing {gene_type} presence in spatial spots...")
    logger.info(f"Total {gene_type} to analyze: {len(genes)}")
    
    # Validate input
    if not genes:
        raise ValueError(f"No {gene_type} provided for analysis")
    if df.empty:
        raise ValueError("Empty spatial data provided")
    
    results = []
    
    # Convert gene names to uppercase for case-insensitive matching
    genes_upper = {g.upper() for g in genes}
    df['geneID_upper'] = df['geneID'].str.upper()
    
    # Find valid genes
    valid_genes = set(df['geneID_upper'].unique()) & genes_upper
    missing_genes = genes_upper - valid_genes
    
    logger.info(f"\nGene matching summary for {gene_type}:")
    logger.info(f"- Total {gene_type} in CellChatDB: {len(genes)}")
    logger.info(f"- Found in spatial data: {len(valid_genes)}")
    logger.info(f"- Missing in spatial data: {len(missing_genes)}")
    
    if missing_genes:
        logger.info(f"\nFirst 20 missing {gene_type}:")
        for gene in sorted(list(missing_genes))[:20]:
            logger.info(f"- {gene}")
    
    if len(valid_genes) == 0:
        logger.warning(f"No {gene_type} found in the dataset!")
        return pd.DataFrame()
    
    # Get total number of unique spots
    total_spots = len(df['bin1_ID'].unique())
    logger.info(f"\nTotal number of spatial spots in dataset: {total_spots}")
    
    # Validate total spots
    if total_spots == 0:
        raise ValueError("No spatial spots found in dataset")
    
    # Process genes in chunks
    valid_genes = list(valid_genes)
    n_genes = len(valid_genes)
    n_chunks = (n_genes + chunk_size - 1) // chunk_size
    
    logger.info(f"\nProcessing {n_genes} {gene_type} in {n_chunks} chunks...")
    
    for chunk_idx in range(n_chunks):
        chunk_start = chunk_idx * chunk_size
        chunk_end = min(chunk_start + chunk_size, n_genes)
        chunk_genes = valid_genes[chunk_start:chunk_end]
        
        logger.info(f"\nProcessing chunk {chunk_idx + 1}/{n_chunks} ({chunk_start+1}-{chunk_end} of {n_genes} {gene_type})")
        
        for i, gene in enumerate(chunk_genes, 1):
            # Get spots where this gene is present
            gene_spots = df[df['geneID_upper'] == gene]['bin1_ID'].unique()
            n_spots = len(gene_spots)
            percent_spots = (n_spots / total_spots) * 100
            
            # Validate spot counts
            if n_spots > total_spots:
                raise ValueError(f"Gene {gene} has more spots ({n_spots}) than total spots ({total_spots})")
            
            # Find original case of gene name
            orig_gene = df[df['geneID_upper'] == gene]['geneID'].iloc[0]
            
            # Log progress for each gene
            logger.info(f"Processing {gene_type} {i}/{len(chunk_genes)}: {orig_gene} - {n_spots} spots ({percent_spots:.2f}%)")
            
            results.append({
                'gene': orig_gene,  # Use original case
                'spots_present': int(n_spots),
                'percent_spots': float(percent_spots)
            })
        
        # Clean up memory
        gc.collect()
    
    # Convert to DataFrame and sort
    if results:
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('spots_present', ascending=False)
        
        # Validate results
        if len(results_df) != len(valid_genes):
            raise ValueError(f"Number of results ({len(results_df)}) doesn't match number of valid genes ({len(valid_genes)})")
        
        logger.info(f"\nTop 10 most prevalent {gene_type} (by number of spots):")
        for _, row in results_df.head(10).iterrows():
            logger.info(f"- {row['gene']}: {row['spots_present']} spots ({row['percent_spots']:.2f}%)")
    else:
        results_df = pd.DataFrame()
    
    return results_df

def main():
    """Main function to run the spatial expression analysis."""
    try:
        # Load data
        interactions = load_interaction_data()
        df = load_spatial_data("data/spatial/raw/1996-081_GFM_SS200000954BR_A2_tissue_cleaned_cortex_crop.csv")
        
        # Extract unique genes
        genes = extract_unique_genes(interactions)
        
        # Analyze ligand presence
        logger.info("\n=== Starting Ligand Analysis ===")
        ligand_results = analyze_gene_expression(df, genes['ligands'], "ligands")
        
        # Analyze receptor presence
        logger.info("\n=== Starting Receptor Analysis ===")
        receptor_results = analyze_gene_expression(df, genes['receptors'], "receptors")
        
        # Save results
        output_dir = PROJECT_ROOT / "data/spatial/processed"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        ligand_path = output_dir / "ligand_spots.csv"
        receptor_path = output_dir / "receptor_spots.csv"
        
        # Validate before saving
        if ligand_results.empty or receptor_results.empty:
            raise ValueError("Empty results detected - analysis failed")
        
        ligand_results.to_csv(ligand_path, index=False)
        receptor_results.to_csv(receptor_path, index=False)
        
        logger.info(f"\nResults saved to:")
        logger.info(f"- Ligands: {ligand_path}")
        logger.info(f"- Receptors: {receptor_path}")
        
        # Print final summary
        logger.info("\nFinal Analysis Summary:")
        logger.info(f"Total ligands analyzed: {len(ligand_results)}")
        logger.info(f"Total receptors analyzed: {len(receptor_results)}")
        
    except Exception as e:
        logger.error(f"Error in spatial expression analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main() 