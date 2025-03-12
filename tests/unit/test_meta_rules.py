import pytest
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
import sys
import os
import tempfile

# Add the scripts directory to the path so we can import the modules
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))

# Import the module under test with a try/except to provide better error messages
try:
    from generate_meta_rules import (
        extract_rule_info,
        calculate_rule_associations,
        generate_meta_rule
    )
except ImportError as e:
    pytest.skip(f"Could not import generate_meta_rules module: {e}")

class TestRuleInfoExtraction:
    def test_extract_rule_info(self, tmp_path):
        """Test extracting rule information from SPARQL files."""
        # Create a temporary rule file
        rule_file = tmp_path / "test_rule.rq"
        rule_content = """
        # Rule ID: GFAP_AQP4_S100B_CLIQUE_Astrocyte
        # Rule Type: CLIQUE
        # Cell Type: Astrocyte
        # Confidence Score: 0.85
        # Genes: GFAP, AQP4, S100B
        
        PREFIX astro: <http://example.org/astrocytes/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        CONSTRUCT {
            ?cell a astro:Astrocyte ;
                astro:identifiedBy "GFAP_AQP4_S100B_rule" ;
                astro:confidence "0.85"^^xsd:float .
        }
        WHERE {
            ?cell astro:expressesGene astro:GFAP .
            ?cell astro:expressesGene astro:AQP4 .
            ?cell astro:expressesGene astro:S100B .
        }
        """
        rule_file.write_text(rule_content)
        
        # Extract rule info
        rule_info_list = extract_rule_info(str(rule_file))
        
        # Check that rule info was correctly extracted
        assert len(rule_info_list) == 1
        rule_info = rule_info_list[0]
        
        assert rule_info["rule_id"] == "GFAP_AQP4_S100B_CLIQUE_Astrocyte"
        assert rule_info["rule_type"] == "CLIQUE"
        assert rule_info["cell_type"] == "Astrocyte"
        assert rule_info["confidence"] == 0.85
        assert set(rule_info["genes"]) == {"GFAP", "AQP4", "S100B"}
        assert rule_info["file_path"] == str(rule_file)

class TestRuleAssociation:
    def test_calculate_rule_associations(self):
        """Test calculation of significant associations between rules."""
        # Create mock rule data
        rules = [
            {
                "rule_id": "GFAP_AQP4_PAIR_Astrocyte",
                "rule_type": "PAIR",
                "cell_type": "Astrocyte",
                "genes": ["GFAP", "AQP4"],
                "confidence": 0.8
            },
            {
                "rule_id": "S100B_ALDH1L1_PAIR_Astrocyte",
                "rule_type": "PAIR",
                "cell_type": "Astrocyte",
                "genes": ["S100B", "ALDH1L1"],
                "confidence": 0.7
            },
            {
                "rule_id": "SLC1A3_GFAP_PAIR_Astrocyte",
                "rule_type": "PAIR",
                "cell_type": "Astrocyte",
                "genes": ["SLC1A3", "GFAP"],
                "confidence": 0.75
            },
            {
                "rule_id": "CD45_CD68_PAIR_Microglia",
                "rule_type": "PAIR",
                "cell_type": "Microglia",
                "genes": ["CD45", "CD68"],
                "confidence": 0.9
            }
        ]
        
        # Create mock rule matches (cells x rules matrix)
        # Cells matching each rule (1=match, 0=no match)
        # Strongly associate rule 0 and rule 1, but not others
        rule_matches = pd.DataFrame({
            "GFAP_AQP4_PAIR_Astrocyte": [1, 1, 1, 0, 0],
            "S100B_ALDH1L1_PAIR_Astrocyte": [1, 1, 0, 0, 0],
            "SLC1A3_GFAP_PAIR_Astrocyte": [0, 0, 1, 1, 0],
            "CD45_CD68_PAIR_Microglia": [0, 0, 0, 0, 1]
        })
        
        # Calculate rule associations
        associations = calculate_rule_associations(
            rules, 
            rule_matches, 
            min_lift=1.5, 
            max_p_value=0.05,
            same_cell_type_only=True
        )
        
        # Check that associations were found
        assert len(associations) > 0
        
        # The strongest association should be between rules 0 and 1
        # Check if first association is between the expected rules
        first_assoc = associations[0]
        rule1_id = first_assoc["rule1_id"]
        rule2_id = first_assoc["rule2_id"]
        
        expected_pair = {
            "GFAP_AQP4_PAIR_Astrocyte", 
            "S100B_ALDH1L1_PAIR_Astrocyte"
        }
        assert {rule1_id, rule2_id} == expected_pair
        
        # Verify that no associations were found between different cell types
        for assoc in associations:
            rule1_idx = next(i for i, r in enumerate(rules) if r["rule_id"] == assoc["rule1_id"])
            rule2_idx = next(i for i, r in enumerate(rules) if r["rule_id"] == assoc["rule2_id"])
            assert rules[rule1_idx]["cell_type"] == rules[rule2_idx]["cell_type"]

class TestMetaRuleGeneration:
    def test_generate_meta_rule(self):
        """Test generation of meta-rules from associated rules."""
        # Create two simple rules
        rule1 = {
            "rule_id": "GFAP_AQP4_PAIR_Astrocyte",
            "rule_type": "PAIR",
            "cell_type": "Astrocyte",
            "genes": ["GFAP", "AQP4"],
            "confidence": 0.8,
            "sparql": """
            CONSTRUCT {
                ?cell a astro:Astrocyte ;
                    astro:identifiedBy "GFAP_AQP4_rule" .
            }
            WHERE {
                ?cell astro:expressesGene astro:GFAP .
                ?cell astro:expressesGene astro:AQP4 .
            }
            """
        }
        
        rule2 = {
            "rule_id": "S100B_ALDH1L1_PAIR_Astrocyte",
            "rule_type": "PAIR",
            "cell_type": "Astrocyte",
            "genes": ["S100B", "ALDH1L1"],
            "confidence": 0.7,
            "sparql": """
            CONSTRUCT {
                ?cell a astro:Astrocyte ;
                    astro:identifiedBy "S100B_ALDH1L1_rule" .
            }
            WHERE {
                ?cell astro:expressesGene astro:S100B .
                ?cell astro:expressesGene astro:ALDH1L1 .
            }
            """
        }
        
        confidence = 0.9
        spatial_distance = 75
        
        # Generate meta-rule
        rule_id, rule_sparql = generate_meta_rule(rule1, rule2, confidence, spatial_distance)
        
        # Check meta-rule properties
        assert "META" in rule_id
        assert "Astrocyte" in rule_id
        assert all(gene in rule_id for gene in ["GFAP", "AQP4", "S100B", "ALDH1L1"])
        
        # Check SPARQL content
        assert "CONSTRUCT" in rule_sparql
        assert "astro:Astrocyte" in rule_sparql
        assert all(gene in rule_sparql for gene in ["GFAP", "AQP4", "S100B", "ALDH1L1"])
        assert str(spatial_distance) in rule_sparql
        
        # Check confidence value
        assert str(confidence) in rule_sparql 