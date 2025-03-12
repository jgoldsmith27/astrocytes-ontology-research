import pytest
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os
from typing import Dict, List

# Add the scripts directory to the path so we can import the modules
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))

# Import the module under test with a try/except to provide better error messages
try:
    from conflict_resolution import (
        Cell,
        ConflictResolutionManager,
        resolve_cell_conflicts
    )
except ImportError as e:
    pytest.skip(f"Could not import conflict_resolution module: {e}")

class TestCellClass:
    def test_distance_calculation(self):
        """Test calculation of distance between cells."""
        cell1 = Cell(
            cell_id="cell1",
            cell_type="Astrocyte",
            x=10.0,
            y=10.0,
            radius=5.0,
            confidence=0.8,
            rule_id="GFAP_AQP4_PAIR",
            rule_type="pair",
            genes=["GFAP", "AQP4"]
        )
        
        cell2 = Cell(
            cell_id="cell2",
            cell_type="Astrocyte",
            x=20.0,
            y=10.0,
            radius=5.0,
            confidence=0.7,
            rule_id="S100B_ALDH1L1_PAIR",
            rule_type="pair",
            genes=["S100B", "ALDH1L1"]
        )
        
        # Distance should be 10.0 (cells are 10 units apart on x-axis)
        assert cell1.distance_to(cell2) == 10.0
        
        # Distance should be symmetrical
        assert cell2.distance_to(cell1) == 10.0
    
    def test_overlap_percentage(self):
        """Test calculation of overlap percentage between cells."""
        cell1 = Cell(
            cell_id="cell1",
            cell_type="Astrocyte",
            x=10.0,
            y=10.0,
            radius=5.0,
            confidence=0.8,
            rule_id="GFAP_AQP4_PAIR",
            rule_type="pair",
            genes=["GFAP", "AQP4"]
        )
        
        # Test no overlap (cells far apart)
        cell2 = Cell(
            cell_id="cell2",
            cell_type="Astrocyte",
            x=20.0,
            y=10.0,
            radius=4.0,
            confidence=0.7,
            rule_id="S100B_ALDH1L1_PAIR",
            rule_type="pair",
            genes=["S100B", "ALDH1L1"]
        )
        assert cell1.overlap_percentage(cell2) == 0.0
        
        # Test partial overlap
        cell3 = Cell(
            cell_id="cell3",
            cell_type="Astrocyte",
            x=14.0,
            y=10.0,
            radius=5.0,
            confidence=0.7,
            rule_id="S100B_ALDH1L1_PAIR",
            rule_type="pair",
            genes=["S100B", "ALDH1L1"]
        )
        # Distance is 4, radii are both 5, so there should be significant overlap
        overlap = cell1.overlap_percentage(cell3)
        assert overlap > 0.0 and overlap < 1.0
        
        # Test complete overlap (one cell inside another)
        cell4 = Cell(
            cell_id="cell4",
            cell_type="Astrocyte",
            x=10.0,
            y=10.0,
            radius=3.0,
            confidence=0.7,
            rule_id="S100B_ALDH1L1_PAIR",
            rule_type="pair",
            genes=["S100B", "ALDH1L1"]
        )
        # Cell4 is completely inside cell1, so smaller_r²/larger_r² = (3/5)² = 0.36
        assert np.isclose(cell4.overlap_percentage(cell1), (3/5)**2, rtol=1e-3)

    def test_resolution_score(self):
        """Test calculation of conflict resolution score."""
        # Create cell with various properties
        cell = Cell(
            cell_id="cell1",
            cell_type="Astrocyte",
            x=10.0,
            y=10.0,
            radius=5.0,
            confidence=0.8,
            rule_id="GFAP_AQP4_S100B_CLIQUE",
            rule_type="clique",
            genes=["GFAP", "AQP4", "S100B"],
            negative_markers_missing=True,
            derived_from_meta_rule=False
        )
        
        # Create mock hierarchies and specificities
        cell_type_hierarchy = {
            "Astrocyte": 2,
            "Neuron": 1,
            "Microglia": 3
        }
        
        gene_specificity = {
            "GFAP": 0.9,
            "AQP4": 0.8,
            "S100B": 0.7,
            "ALDH1L1": 0.6
        }
        
        # Calculate score
        score = cell.get_resolution_score(cell_type_hierarchy, gene_specificity)
        
        # Score should be positive
        assert score > 0
        
        # Score should include base confidence (0.8)
        assert score > 0.8
        
        # Create another cell with different properties for comparison
        cell2 = Cell(
            cell_id="cell2",
            cell_type="Astrocyte",
            x=10.0,
            y=10.0,
            radius=5.0,
            confidence=0.7,  # Lower confidence
            rule_id="GFAP_AQP4_PAIR",
            rule_type="pair",  # Lower scoring rule type
            genes=["GFAP", "AQP4"],  # Fewer genes
            negative_markers_missing=False,  # No negative marker bonus
            derived_from_meta_rule=False
        )
        
        score2 = cell2.get_resolution_score(cell_type_hierarchy, gene_specificity)
        
        # First cell should score higher
        assert score > score2

class TestConflictResolutionManager:
    def setup_method(self):
        """Setup test data for conflict resolution tests."""
        # Create test cells dataframe
        self.cells_df = pd.DataFrame({
            'cell_id': ['cell1', 'cell2', 'cell3', 'cell4'],
            'cell_type': ['Astrocyte', 'Astrocyte', 'Neuron', 'Microglia'],
            'x': [10.0, 12.0, 30.0, 50.0],
            'y': [10.0, 11.0, 30.0, 50.0],
            'radius': [5.0, 5.0, 4.0, 3.0],
            'confidence': [0.8, 0.7, 0.9, 0.85],
            'rule_id': ['GFAP_AQP4_CLIQUE', 'S100B_ALDH1L1_PAIR', 'MAP2_NEUN_CLIQUE', 'CD45_CD68_PAIR'],
            'rule_type': ['clique', 'pair', 'clique', 'pair'],
            'genes': [['GFAP', 'AQP4'], ['S100B', 'ALDH1L1'], ['MAP2', 'NEUN'], ['CD45', 'CD68']],
            'negative_markers_missing': [True, True, True, True],
            'derived_from_meta_rule': [False, False, False, False]
        })
        
        # Create a temporary directory for output
        self.test_output_dir = "tests/unit/test_output"
        os.makedirs(self.test_output_dir, exist_ok=True)
    
    def test_conflict_detection(self):
        """Test detection of conflicts between cells."""
        # Initialize resolution manager
        manager = ConflictResolutionManager(
            overlap_threshold=0.3,
            same_type_overlap_threshold=0.7,
            generate_visuals=False,
            output_dir=self.test_output_dir
        )
        
        # Load cells
        manager.load_cells(self.cells_df)
        
        # Find conflicts
        manager.find_conflicts()
        
        # Check that a conflict was detected between cell1 and cell2 (they overlap)
        conflict_pairs = [(c1.cell_id, c2.cell_id) for c1, c2 in manager.conflicts]
        assert ('cell1', 'cell2') in conflict_pairs or ('cell2', 'cell1') in conflict_pairs
        
        # No conflict should be detected between distant cells
        assert ('cell1', 'cell3') not in conflict_pairs
        assert ('cell2', 'cell3') not in conflict_pairs
    
    def test_conflict_resolution(self):
        """Test resolution of conflicts between cells."""
        # Initialize resolution manager
        manager = ConflictResolutionManager(
            overlap_threshold=0.3,
            same_type_overlap_threshold=0.7,
            generate_visuals=False,
            output_dir=self.test_output_dir
        )
        
        # Load cells
        manager.load_cells(self.cells_df)
        
        # Set up gene specificity info
        cell_type_gene_expression = {
            'Astrocyte': {'GFAP': 0.9, 'AQP4': 0.8, 'S100B': 0.7, 'ALDH1L1': 0.6},
            'Neuron': {'MAP2': 0.9, 'NEUN': 0.8},
            'Microglia': {'CD45': 0.9, 'CD68': 0.8}
        }
        manager.calculate_gene_specificity(cell_type_gene_expression)
        
        # Find and resolve conflicts
        manager.find_conflicts()
        resolved_cells = manager.resolve_conflicts()
        
        # Check that we got a valid resolution
        assert isinstance(resolved_cells, list)
        assert all(isinstance(cell, Cell) for cell in resolved_cells)
        
        # Either cell1 or cell2 should be removed (but not both), as they overlap
        cell_ids = [cell.cell_id for cell in resolved_cells]
        assert not ('cell1' in cell_ids and 'cell2' in cell_ids)
        assert 'cell1' in cell_ids or 'cell2' in cell_ids
        
        # cell3 and cell4 should be kept (no conflicts)
        assert 'cell3' in cell_ids
        assert 'cell4' in cell_ids

class TestEndToEndResolution:
    def test_resolve_cell_conflicts_function(self, tmp_path):
        """Test the end-to-end conflict resolution function."""
        # Create test cells dataframe
        cells_df = pd.DataFrame({
            'cell_id': ['cell1', 'cell2', 'cell3', 'cell4'],
            'cell_type': ['Astrocyte', 'Astrocyte', 'Neuron', 'Microglia'],
            'x': [10.0, 12.0, 30.0, 50.0],
            'y': [10.0, 11.0, 30.0, 50.0],
            'radius': [5.0, 5.0, 4.0, 3.0],
            'confidence': [0.8, 0.7, 0.9, 0.85],
            'rule_id': ['GFAP_AQP4_CLIQUE', 'S100B_ALDH1L1_PAIR', 'MAP2_NEUN_CLIQUE', 'CD45_CD68_PAIR'],
            'rule_type': ['clique', 'pair', 'clique', 'pair'],
            'genes': [['GFAP', 'AQP4'], ['S100B', 'ALDH1L1'], ['MAP2', 'NEUN'], ['CD45', 'CD68']],
            'negative_markers_missing': [True, True, True, True],
            'derived_from_meta_rule': [False, False, False, False]
        })
        
        # Create a temporary directory for output
        output_dir = str(tmp_path / "conflicts_output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Create mock cell type expression data
        cell_type_expression_data = {
            'Astrocyte': {'GFAP': 0.9, 'AQP4': 0.8, 'S100B': 0.7, 'ALDH1L1': 0.6},
            'Neuron': {'MAP2': 0.9, 'NEUN': 0.8},
            'Microglia': {'CD45': 0.9, 'CD68': 0.8}
        }
        
        # Run the conflict resolution
        resolved_df = resolve_cell_conflicts(
            cells_df=cells_df,
            output_dir=output_dir,
            overlap_threshold=0.3,
            same_type_overlap_threshold=0.7,
            cell_type_expression_data=cell_type_expression_data,
            generate_visuals=False
        )
        
        # Check that the output is a DataFrame
        assert isinstance(resolved_df, pd.DataFrame)
        
        # Check that there are fewer rows in the result (conflicts resolved)
        assert len(resolved_df) < len(cells_df)
        
        # Check that output files were created
        assert any(Path(output_dir).glob("resolved_cells*.csv")) 