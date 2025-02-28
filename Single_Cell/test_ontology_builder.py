import os
import pandas as pd
import numpy as np
from single_cell_ontology_builder import GeneExpressionOntologyBuilder
from rdflib import Graph

def create_test_dataframe():
    """Create a test DataFrame with mock expression data"""
    print("🧪 Creating test expression data...")
    
    # Create mock data
    cell_types = ['Neuron', 'Astrocyte', 'Microglia']
    genes = [f'gene_{i}' for i in range(10)]
    
    # Random expression values
    data = np.random.rand(len(cell_types), len(genes))
    
    # Create DataFrame
    df = pd.DataFrame(data, index=cell_types, columns=genes)
    print(f"✅ Test DataFrame created with shape {df.shape}")
    return df

def test_ontology_builder():
    """Test the GeneExpressionOntologyBuilder"""
    print("\n🧪 TESTING ONTOLOGY BUILDER")
    
    # Create test data
    test_df = create_test_dataframe()
    test_file = "test_ontology.ttl" 