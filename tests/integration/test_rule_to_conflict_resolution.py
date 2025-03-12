import pytest
import os
import pandas as pd
import numpy as np
import tempfile
from pathlib import Path
import sys
import shutil
import json

# Add the scripts directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))

# Import components with try/except to provide better error messages
try:
    from generate_coexpression_rules import generate_clique_rule, generate_pair_rule
    from conflict_resolution import Cell, resolve_cell_conflicts
except ImportError as e:
    pytest.skip(f"Could not import required modules: {e}")

class TestRulesToCellIdentification:
    """Integration test for the rule generation to conflict resolution pipeline."""
    
    def setup_method(self):
        """Set up the test environment with necessary data."""
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()
        self.rules_dir = os.path.join(self.test_dir, "rules")
        self.output_dir = os.path.join(self.test_dir, "output")
        os.makedirs(self.rules_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Define test genes and cell types
        self.genes = {
            "Astrocyte": ["GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A3"],
            "Neuron": ["MAP2", "NEUN", "SYP"],
            "Microglia": ["CD45", "CD68", "IBA1"]
        }
        
        # Create some sample rules
        self._create_sample_rules()
        
        # Create sample cell data that would be identified using these rules
        self._create_sample_identified_cells()
        
    def teardown_method(self):
        """Clean up temporary test directory."""
        shutil.rmtree(self.test_dir)
        
    def _create_sample_rules(self):
        """Create sample SPARQL rule files with different clique sizes."""
        # Create clique rule for astrocytes (size 3)
        rule_id, rule_content = generate_clique_rule(
            clique=["GFAP", "AQP4", "S100B"],
            cell_type="Astrocyte",
            distance_threshold=50
        )
        with open(os.path.join(self.rules_dir, f"{rule_id}.rq"), "w") as f:
            f.write(f"# Rule ID: {rule_id}\n")
            f.write("# Rule Type: CLIQUE\n")
            f.write("# Cell Type: Astrocyte\n")
            f.write("# Confidence Score: 0.85\n")
            f.write("# Genes: GFAP, AQP4, S100B\n\n")
            f.write(rule_content)
        
        # Create clique rule for astrocytes (size 4)
        rule_id, rule_content = generate_clique_rule(
            clique=["GFAP", "AQP4", "S100B", "ALDH1L1"],
            cell_type="Astrocyte",
            distance_threshold=50
        )
        with open(os.path.join(self.rules_dir, f"{rule_id}.rq"), "w") as f:
            f.write(f"# Rule ID: {rule_id}\n")
            f.write("# Rule Type: CLIQUE\n")
            f.write("# Cell Type: Astrocyte\n")
            f.write("# Confidence Score: 0.9\n")  # Higher confidence for larger clique
            f.write("# Genes: GFAP, AQP4, S100B, ALDH1L1\n\n")
            f.write(rule_content)
            
        # Create pair rule for astrocytes
        rule_id, rule_content = generate_pair_rule(
            gene_i="ALDH1L1",
            gene_j="SLC1A3",
            weight=0.75,
            cell_type="Astrocyte",
            distance_threshold=50
        )
        with open(os.path.join(self.rules_dir, f"{rule_id}.rq"), "w") as f:
            f.write(f"# Rule ID: {rule_id}\n")
            f.write("# Rule Type: PAIR\n")
            f.write("# Cell Type: Astrocyte\n")
            f.write("# Confidence Score: 0.75\n")
            f.write("# Genes: ALDH1L1, SLC1A3\n\n")
            f.write(rule_content)
            
        # Create rule for neurons
        rule_id, rule_content = generate_clique_rule(
            clique=["MAP2", "NEUN"],
            cell_type="Neuron",
            distance_threshold=50
        )
        with open(os.path.join(self.rules_dir, f"{rule_id}.rq"), "w") as f:
            f.write(f"# Rule ID: {rule_id}\n")
            f.write("# Rule Type: CLIQUE\n")
            f.write("# Cell Type: Neuron\n")
            f.write("# Confidence Score: 0.9\n")
            f.write("# Genes: MAP2, NEUN\n\n")
            f.write(rule_content)
    
    def _create_sample_identified_cells(self):
        """Create sample data for cells identified by the rules."""
        # Create a DataFrame representing cells identified by our rules
        # with some overlapping cells to test conflict resolution
        self.cells_df = pd.DataFrame({
            'cell_id': ['cell1', 'cell2', 'cell3', 'cell4', 'cell5', 'cell6'],
            'cell_type': ['Astrocyte', 'Astrocyte', 'Astrocyte', 'Astrocyte', 'Neuron', 'Neuron'],
            'x': [10.0, 12.0, 30.0, 25.0, 50.0, 100.0],
            'y': [10.0, 11.0, 30.0, 28.0, 50.0, 100.0],
            'radius': [5.0, 5.0, 4.0, 5.0, 3.0, 3.5],
            'confidence': [0.85, 0.75, 0.8, 0.9, 0.9, 0.85],
            'rule_id': [
                'GFAP_AQP4_S100B_CLIQUE_Astrocyte',  # Size 3 clique
                'ALDH1L1_SLC1A3_PAIR_Astrocyte',     # Pair
                'GFAP_AQP4_S100B_CLIQUE_Astrocyte',  # Size 3 clique
                'GFAP_AQP4_S100B_ALDH1L1_CLIQUE_Astrocyte',  # Size 4 clique
                'MAP2_NEUN_CLIQUE_Neuron',
                'MAP2_NEUN_CLIQUE_Neuron'
            ],
            'rule_type': ['clique', 'pair', 'clique', 'clique', 'clique', 'clique'],
            'genes': [
                ['GFAP', 'AQP4', 'S100B'],           # Size 3 clique
                ['ALDH1L1', 'SLC1A3'],               # Pair
                ['GFAP', 'AQP4', 'S100B'],           # Size 3 clique
                ['GFAP', 'AQP4', 'S100B', 'ALDH1L1'], # Size 4 clique
                ['MAP2', 'NEUN'],
                ['MAP2', 'NEUN']
            ],
            'negative_markers_missing': [True, True, True, True, True, True],
            'derived_from_meta_rule': [False, False, False, False, False, False]
        })
        
        # Save the cell data
        self.cells_df.to_csv(os.path.join(self.test_dir, "identified_cells.csv"), index=False)
    
    def test_end_to_end_pipeline(self):
        """Test the full pipeline from rule generation to conflict resolution."""
        # 1. Verify that we have the expected rule files
        rule_files = list(Path(self.rules_dir).glob("*.rq"))
        assert len(rule_files) == 4
        
        # 2. Check rule content
        for rule_file in rule_files:
            with open(rule_file, "r") as f:
                content = f.read()
                assert "CONSTRUCT" in content
                assert "WHERE" in content
                assert "astro:expressesGene" in content
        
        # 3. Create mock cell type expression data for specificity calculation
        cell_type_expression_data = {
            'Astrocyte': {'GFAP': 0.9, 'AQP4': 0.8, 'S100B': 0.7, 'ALDH1L1': 0.6, 'SLC1A3': 0.5},
            'Neuron': {'MAP2': 0.9, 'NEUN': 0.8, 'SYP': 0.7},
            'Microglia': {'CD45': 0.9, 'CD68': 0.8, 'IBA1': 0.7}
        }
        
        # 4. Resolve conflicts between identified cells
        resolved_df = resolve_cell_conflicts(
            cells_df=self.cells_df,
            output_dir=self.output_dir,
            overlap_threshold=0.3,
            same_type_overlap_threshold=0.7,
            cell_type_expression_data=cell_type_expression_data,
            generate_visuals=False
        )
        
        # 5. Verify conflict resolution results
        assert isinstance(resolved_df, pd.DataFrame)
        
        # Cell1 and cell2 should have a conflict (overlap) - one should be removed
        # Should not have both cell1 and cell2
        cell_ids = resolved_df['cell_id'].tolist()
        assert not ('cell1' in cell_ids and 'cell2' in cell_ids)
        
        # Should have either cell1 or cell2 (not both removed)
        assert 'cell1' in cell_ids or 'cell2' in cell_ids
        
        # Cell3 should remain (no conflict)
        assert 'cell3' in cell_ids
        
        # Cell4 (size 4 clique) and cell3 (size 3 clique) might overlap, but size 4 should be prioritized
        if 'cell3' in cell_ids and 'cell4' in cell_ids:
            # If both are kept, they didn't overlap enough to cause a conflict
            pass
        elif 'cell4' in cell_ids and 'cell3' not in cell_ids:
            # Correct: larger clique (cell4) was prioritized over smaller clique (cell3)
            pass
        elif 'cell3' in cell_ids and 'cell4' not in cell_ids:
            # This is unexpected - the larger clique should be prioritized
            pytest.fail("Conflict resolution unexpectedly favored smaller clique over larger clique")
        
        # Cell5 and cell6 (neurons) should remain (no conflict)
        assert 'cell5' in cell_ids
        assert 'cell6' in cell_ids
        
        # 6. Check that output files were created
        assert any(Path(self.output_dir).glob("resolved_cells*.csv"))
        
        # 7. Verify that the resolution favored the clique rule over the pair rule
        if 'cell1' in cell_ids and 'cell2' not in cell_ids:
            # Cell1 (clique rule) was kept, cell2 (pair rule) was removed
            pass
        elif 'cell2' in cell_ids and 'cell1' not in cell_ids:
            # This is unexpected based on our scoring system 
            # (clique rules should score higher than pair rules)
            pytest.fail("Conflict resolution unexpectedly favored pair rule over clique rule")
            
    def test_rule_type_prioritization(self):
        """Test that conflict resolution correctly prioritizes rule types."""
        # Create mock cells with different rule types but same cell type
        meta_rule_cell = {
            'cell_id': 'meta_cell',
            'cell_type': 'Astrocyte',
            'x': 10.0,
            'y': 10.0,
            'radius': 5.0,
            'confidence': 0.8,
            'rule_id': 'META_RULE_Astrocyte',
            'rule_type': 'meta',
            'genes': ['GFAP', 'AQP4', 'S100B', 'ALDH1L1'],
            'negative_markers_missing': True,
            'derived_from_meta_rule': True
        }
        
        clique_rule_cell = {
            'cell_id': 'clique_cell',
            'cell_type': 'Astrocyte',
            'x': 12.0,
            'y': 11.0,
            'radius': 5.0,
            'confidence': 0.8,  # Same confidence
            'rule_id': 'GFAP_AQP4_S100B_CLIQUE_Astrocyte',
            'rule_type': 'clique',
            'genes': ['GFAP', 'AQP4', 'S100B'],
            'negative_markers_missing': True,
            'derived_from_meta_rule': False
        }
        
        pair_rule_cell = {
            'cell_id': 'pair_cell',
            'cell_type': 'Astrocyte',
            'x': 11.0,
            'y': 12.0,
            'radius': 5.0,
            'confidence': 0.8,  # Same confidence
            'rule_id': 'GFAP_AQP4_PAIR_Astrocyte',
            'rule_type': 'pair',
            'genes': ['GFAP', 'AQP4'],
            'negative_markers_missing': True,
            'derived_from_meta_rule': False
        }
        
        # Create a DataFrame with all three cells that overlap
        test_cells_df = pd.DataFrame([meta_rule_cell, clique_rule_cell, pair_rule_cell])
        
        # Create mock cell type expression data
        cell_type_expression_data = {
            'Astrocyte': {'GFAP': 0.9, 'AQP4': 0.8, 'S100B': 0.7, 'ALDH1L1': 0.6}
        }
        
        # Resolve conflicts
        resolved_df = resolve_cell_conflicts(
            cells_df=test_cells_df,
            output_dir=self.output_dir,
            overlap_threshold=0.1,  # Low threshold to ensure conflict detection
            same_type_overlap_threshold=0.1,
            cell_type_expression_data=cell_type_expression_data,
            generate_visuals=False
        )
        
        # Meta rule should win over both clique and pair
        assert 'meta_cell' in resolved_df['cell_id'].tolist()
        
        # Create another test case: clique vs pair
        test_cells_df = pd.DataFrame([clique_rule_cell, pair_rule_cell])
        
        # Resolve conflicts
        resolved_df = resolve_cell_conflicts(
            cells_df=test_cells_df,
            output_dir=self.output_dir,
            overlap_threshold=0.1,  # Low threshold to ensure conflict detection
            same_type_overlap_threshold=0.1,
            cell_type_expression_data=cell_type_expression_data,
            generate_visuals=False
        )
        
        # Clique rule should win over pair
        assert 'clique_cell' in resolved_df['cell_id'].tolist() 