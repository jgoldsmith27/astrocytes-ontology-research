import pytest
import pandas as pd
import networkx as nx
import os
import json
from pathlib import Path

@pytest.fixture
def sample_gene_expression_data():
    """Fixture providing a small sample of gene expression data in pandas DataFrame format."""
    return pd.DataFrame({
        'cell_id': ['cell1', 'cell2', 'cell3', 'cell4', 'cell5'],
        'GFAP': [1.2, 0.1, 2.3, 0.4, 1.8],
        'AQP4': [1.5, 0.2, 1.9, 0.3, 2.1],
        'S100B': [0.9, 0.3, 1.4, 0.1, 1.2],
        'ALDH1L1': [1.1, 0.1, 1.3, 0.2, 0.9],
        'SLC1A3': [0.8, 0.2, 1.2, 0.3, 1.0],
        'cluster': [1, 2, 1, 2, 1]
    })

@pytest.fixture
def sample_spatial_data():
    """Fixture providing spatial transcriptomics data with coordinates."""
    return pd.DataFrame({
        'spot_id': ['spot1', 'spot2', 'spot3', 'spot4', 'spot5'],
        'x': [1.0, 10.0, 20.0, 30.0, 40.0],
        'y': [1.0, 5.0, 15.0, 25.0, 35.0],
        'GFAP': [1.2, 0.1, 2.3, 0.4, 1.8],
        'AQP4': [1.5, 0.2, 1.9, 0.3, 2.1],
        'S100B': [0.9, 0.3, 1.4, 0.1, 1.2],
        'ALDH1L1': [1.1, 0.1, 1.3, 0.2, 0.9],
        'SLC1A3': [0.8, 0.2, 1.2, 0.3, 1.0],
    })

@pytest.fixture
def sample_coexpression_graph():
    """Fixture providing a small gene co-expression graph."""
    G = nx.Graph()
    
    # Add genes as nodes
    genes = ['GFAP', 'AQP4', 'S100B', 'ALDH1L1', 'SLC1A3']
    G.add_nodes_from(genes)
    
    # Add co-expression edges with weights
    edges = [
        ('GFAP', 'AQP4', 0.8),
        ('GFAP', 'S100B', 0.7),
        ('AQP4', 'ALDH1L1', 0.6),
        ('S100B', 'SLC1A3', 0.5),
        ('ALDH1L1', 'SLC1A3', 0.4)
    ]
    
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)
    
    return G

@pytest.fixture
def sample_rule_data():
    """Fixture providing sample rule data in JSON format."""
    return [
        {
            "rule_id": "rule1",
            "rule_type": "CLIQUE",
            "genes": ["GFAP", "AQP4", "S100B"],
            "threshold": 0.5
        },
        {
            "rule_id": "rule2",
            "rule_type": "PAIR",
            "gene1": "ALDH1L1",
            "gene2": "SLC1A3",
            "threshold": 0.4
        }
    ]

@pytest.fixture
def sample_rule_sparql():
    """Fixture providing sample SPARQL rule content."""
    return """
    CONSTRUCT {
        ?cell a astro:Astrocyte ;
            astro:identifiedBy "GFAP_AQP4_S100B_rule" .
    }
    WHERE {
        ?cell astro:expressesGene ?gene1 .
        ?cell astro:expressesGene ?gene2 .
        ?cell astro:expressesGene ?gene3 .
        FILTER(?gene1 = astro:GFAP && ?gene2 = astro:AQP4 && ?gene3 = astro:S100B)
    }
    """

@pytest.fixture
def test_output_dir(tmp_path):
    """Create a temporary directory for test outputs."""
    output_dir = tmp_path / "test_outputs"
    output_dir.mkdir()
    return output_dir 