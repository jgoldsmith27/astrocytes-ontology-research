import pytest
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
import sys
import os

# Add the scripts directory to the path so we can import the modules
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))

# Import the module under test with a try/except to provide better error messages
try:
    from generate_coexpression_rules import (
        get_binary_expression, 
        calculate_coexpression, 
        find_cliques,
        generate_clique_rule,
        generate_pair_rule
    )
except ImportError as e:
    pytest.skip(f"Could not import generate_coexpression_rules module: {e}")

class TestBinaryExpression:
    def test_get_binary_expression(self):
        """Test conversion of expression matrix to binary based on threshold."""
        # Mock adata.X as a numpy array
        mock_expr_matrix = np.array([
            [0.1, 0.5, 0.3],
            [0.6, 0.2, 0.8],
            [0.4, 0.7, 0.1]
        ])
        
        # Mock a simple AnnData-like object with X attribute
        class MockAnnData:
            def __init__(self, X):
                self.X = X
        
        adata = MockAnnData(mock_expr_matrix)
        cell_indices = [0, 1, 2]  # All cells
        min_expression = 0.4
        
        # Expected binary matrix based on min_expression=0.4
        expected = np.array([
            [0, 1, 0],
            [1, 0, 1],
            [1, 1, 0]
        ])
        
        result = get_binary_expression(adata, cell_indices, min_expression)
        np.testing.assert_array_equal(result, expected)

class TestCoexpressionCalculation:
    def test_calculate_coexpression(self):
        """Test calculation of coexpression probabilities between genes."""
        # Create a binary expression matrix:
        # Genes A, B, C for cells 1-5
        binary_matrix = np.array([
            [1, 1, 0],  # Cell 1: A and B expressed
            [1, 0, 0],  # Cell 2: A expressed
            [1, 1, 1],  # Cell 3: A, B, and C expressed
            [0, 1, 1],  # Cell 4: B and C expressed
            [1, 1, 0]   # Cell 5: A and B expressed
        ])
        
        # Consider all genes
        gene_indices = [0, 1, 2]  # genes A, B, C
        min_coexpression = 0.5
        
        # Calculate the coexpression graph
        G = calculate_coexpression(binary_matrix, gene_indices, min_coexpression)
        
        # Verify the graph structure
        assert isinstance(G, nx.Graph)
        
        # Check nodes (should be 0, 1, 2 for genes A, B, C)
        assert set(G.nodes()) == {0, 1, 2}
        
        # Check edges and weights
        # P(B|A) = 3/4 = 0.75 (A is expressed in 4 cells, and both A and B in 3 cells)
        # P(A|B) = 3/4 = 0.75 (B is expressed in 4 cells, and both A and B in 3 cells)
        # P(C|A) = 1/4 = 0.25 (A is expressed in 4 cells, and both A and C in 1 cell)
        # P(A|C) = 1/2 = 0.5 (C is expressed in 2 cells, and both A and C in 1 cell)
        # P(C|B) = 2/4 = 0.5 (B is expressed in 4 cells, and both B and C in 2 cells)
        # P(B|C) = 2/2 = 1.0 (C is expressed in 2 cells, and both B and C in 2 cells)
        
        # With min_coexpression=0.5, we expect edges (0,1), (1,2), (2,1)
        expected_edges = [(0, 1), (1, 2)]  # Undirected graph merges (1,2) and (2,1)
        for u, v in expected_edges:
            assert G.has_edge(u, v)
        
        # Verify edge weights
        assert G[0][1]['weight'] == 0.75  # max(P(B|A), P(A|B)) = max(0.75, 0.75) = 0.75
        assert G[1][2]['weight'] == 1.0   # max(P(C|B), P(B|C)) = max(0.5, 1.0) = 1.0

class TestCliqueDetection:
    def test_find_cliques(self):
        """Test finding cliques in a coexpression graph."""
        # Create a simple graph with multiple cliques
        G = nx.Graph()
        
        # Add nodes for genes
        genes = list(range(6))  # 6 genes labeled 0-5
        G.add_nodes_from(genes)
        
        # Add edges to create:
        # - a clique of genes 0, 1, 2
        # - a clique of genes 2, 3, 4
        # - a clique of genes 0, 1, 2, 5
        edges = [
            (0, 1, 0.8),
            (0, 2, 0.7),
            (0, 5, 0.9),
            (1, 2, 0.9),
            (1, 5, 0.8),
            (2, 3, 0.8),
            (2, 4, 0.7),
            (2, 5, 0.7),
            (3, 4, 0.6)
        ]
        
        for u, v, w in edges:
            G.add_edge(u, v, weight=w)
        
        # Test finding cliques of size 3
        cliques_size_3 = find_cliques(G, min_size=3, max_size=3)
        expected_cliques_size_3 = [{0, 1, 2}, {2, 3, 4}, {0, 1, 5}, {0, 2, 5}, {1, 2, 5}]
        
        # Convert found cliques to sets for comparison
        found_cliques_as_sets = [set(clique) for clique in cliques_size_3]
        assert len(cliques_size_3) == len(expected_cliques_size_3)
        assert set(map(frozenset, found_cliques_as_sets)) == set(map(frozenset, expected_cliques_size_3))
        
        # Test finding cliques of size 4
        cliques_size_4 = find_cliques(G, min_size=4, max_size=4)
        expected_cliques_size_4 = [{0, 1, 2, 5}]
        
        found_cliques_as_sets = [set(clique) for clique in cliques_size_4]
        assert len(cliques_size_4) == len(expected_cliques_size_4)
        assert set(map(frozenset, found_cliques_as_sets)) == set(map(frozenset, expected_cliques_size_4))
        
        # Test finding all cliques from size 3 to N (in this case N=4)
        all_cliques = find_cliques(G, min_size=3)
        
        # Should include both size 3 and size 4 cliques
        assert len(all_cliques) == len(expected_cliques_size_3) + len(expected_cliques_size_4)
        
        # Test that each clique size is represented in the results
        clique_sizes = [len(clique) for clique in all_cliques]
        assert 3 in clique_sizes
        assert 4 in clique_sizes
    
    def test_find_cliques_with_variable_sizes(self):
        """Test finding cliques with variable sizes from min_size to max_size."""
        # Create a more complex graph with cliques of varying sizes
        G = nx.Graph()
        
        # Add 8 nodes to potentially form cliques of up to size 8
        genes = list(range(8))
        G.add_nodes_from(genes)
        
        # Create a complete graph of size 8 (every node connected to every other node)
        for i in range(8):
            for j in range(i+1, 8):
                G.add_edge(i, j, weight=0.7 + (i*j % 3) * 0.1)  # Varying weights
        
        # Find cliques with different size ranges
        # Find all cliques of size exactly 3
        cliques_size_3 = find_cliques(G, min_size=3, max_size=3)
        # Find all cliques of sizes 3 through 5
        cliques_size_3_to_5 = find_cliques(G, min_size=3, max_size=5)
        # Find all cliques of sizes 6 and above
        cliques_size_6_and_up = find_cliques(G, min_size=6)
        
        # In a complete graph of size n, the number of k-cliques is binomial(n,k)
        # For n=8, k=3, we have binomial(8,3) = 56 cliques of size 3
        assert len(cliques_size_3) == 56
        
        # Check some expected properties
        clique_sizes_3_to_5 = [len(clique) for clique in cliques_size_3_to_5]
        assert set(clique_sizes_3_to_5) == {3, 4, 5}
        
        clique_sizes_6_and_up = [len(clique) for clique in cliques_size_6_and_up]
        assert set(clique_sizes_6_and_up) == {6, 7, 8}
        
        # Test that the max_size parameter is respected
        for clique in cliques_size_3_to_5:
            assert len(clique) <= 5
            
        # Test that the min_size parameter is respected
        for clique in cliques_size_6_and_up:
            assert len(clique) >= 6
            
        # Test finding all possible cliques from size 3 up
        all_cliques = find_cliques(G, min_size=3)
        all_clique_sizes = [len(clique) for clique in all_cliques]
        assert set(all_clique_sizes) == {3, 4, 5, 6, 7, 8}

class TestRuleGeneration:
    def test_generate_clique_rule(self):
        """Test generation of SPARQL rule for a gene clique."""
        # Test a clique of genes of size 3
        clique_size_3 = ["GFAP", "AQP4", "S100B"]
        cell_type = "Astrocyte"
        distance_threshold = 50
        
        # Generate the clique rule
        rule_id, rule_content = generate_clique_rule(clique_size_3, cell_type, distance_threshold)
        
        # Check the rule ID format
        assert "GFAP" in rule_id
        assert "AQP4" in rule_id
        assert "S100B" in rule_id
        assert "CLIQUE" in rule_id
        assert cell_type in rule_id
        
        # Check rule content
        assert "CONSTRUCT" in rule_content
        assert "astro:Astrocyte" in rule_content
        assert all(gene in rule_content for gene in clique_size_3)
        assert str(distance_threshold) in rule_content
        
        # Test a larger clique (size 5)
        clique_size_5 = ["GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A3"]
        
        # Generate the clique rule for larger clique
        rule_id, rule_content = generate_clique_rule(clique_size_5, cell_type, distance_threshold)
        
        # Check rule ID includes all genes
        for gene in clique_size_5:
            assert gene in rule_id
        
        # Check SPARQL content has all gene references
        for gene in clique_size_5:
            assert gene in rule_content
        
        # Check that the rule has the correct number of gene clauses
        gene_clauses_count = rule_content.count("expressesGene")
        assert gene_clauses_count >= len(clique_size_5)
        
    def test_generate_pair_rule(self):
        """Test generation of SPARQL rule for a gene pair."""
        # Test a pair of genes
        gene_i = "ALDH1L1"
        gene_j = "SLC1A3"
        weight = 0.75
        cell_type = "Astrocyte"
        distance_threshold = 50
        
        # Generate the pair rule
        rule_id, rule_content = generate_pair_rule(gene_i, gene_j, weight, cell_type, distance_threshold)
        
        # Check the rule ID format
        assert gene_i in rule_id
        assert gene_j in rule_id
        assert "PAIR" in rule_id
        assert cell_type in rule_id
        
        # Check rule content
        assert "CONSTRUCT" in rule_content
        assert "astro:Astrocyte" in rule_content
        assert gene_i in rule_content
        assert gene_j in rule_content
        assert str(distance_threshold) in rule_content
        
    def test_generate_rules_with_negative_markers(self):
        """Test generation of rules with negative markers."""
        # Test a clique of genes with negative markers
        clique = ["GFAP", "AQP4", "S100B"]
        cell_type = "Astrocyte"
        distance_threshold = 50
        negative_markers = ["CD45", "PDGFRA"]
        
        # Generate the clique rule with negative markers
        rule_id, rule_content = generate_clique_rule(clique, cell_type, distance_threshold, negative_markers)
        
        # Check for negative marker filtering in the rule
        assert all(marker in rule_content for marker in negative_markers)
        assert "FILTER NOT EXISTS" in rule_content or "MINUS" in rule_content 