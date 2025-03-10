#!/bin/bash
# Script to run the pipeline for astrocyte identification in spatial data
# This pipeline identifies whole cells (Astrocytes and Astrocytes1) in spatial data
# using co-expression patterns derived from single-cell RNA-seq data.

set -e  # Exit on error

# Activate the virtual environment
source venv/bin/activate

# Define paths
SPATIAL_DATA="data/raw/spatial_data.csv"
SCRNA_DATA="data/raw/CTR081_Fron.h5ad"
ONTOLOGY_FILE="ontologies/cell_type_ontology.ttl"
SPATIAL_OUTPUT="data/processed/spatial_data.ttl"
RULES_DIR="data/processed/rules"
META_RULES_DIR="data/processed/meta_rules"
RESULTS_DIR="data/processed/results"

# Create output directories if they don't exist
mkdir -p data/processed
mkdir -p $RULES_DIR
mkdir -p $META_RULES_DIR
mkdir -p $RESULTS_DIR

echo "==========================================="
echo "Astrocyte Identification Pipeline"
echo "==========================================="
echo "Input single-cell data: $SCRNA_DATA"
echo "Input spatial data: $SPATIAL_DATA"
echo "Output directory: data/processed"
echo "Ontology file: $ONTOLOGY_FILE"
echo

# Step 1: Convert spatial data to RDF/Turtle format
echo "Step 1: Converting spatial data to RDF..."
python scripts/convert_spatial_to_turtle.py \
    --input $SPATIAL_DATA \
    --output $SPATIAL_OUTPUT \
    --ontology $ONTOLOGY_FILE
echo "Conversion complete. Output saved to $SPATIAL_OUTPUT"
echo

# Step 2: Generate co-expression rules for all cell types
echo "Step 2: Generating co-expression rules from single-cell data..."
python scripts/generate_coexpression_rules.py \
    --input $SCRNA_DATA \
    --output-dir $RULES_DIR \
    --min-expression 0.2 \
    --min-coexpression 0.6 \
    --top-genes 30 \
    --neg-marker-threshold 0.05
echo "Rule generation complete. Rules saved to $RULES_DIR"
echo

# Step 3: Generate meta-rules by analyzing rule co-occurrence patterns
echo "Step 3: Generating higher-order co-expression patterns (meta-rules)..."
python scripts/generate_meta_rules.py \
    --rules-dir $RULES_DIR \
    --single-cell-data $SCRNA_DATA \
    --output-dir $META_RULES_DIR \
    --min-lift 3.0 \
    --max-p-value 0.01 \
    --spatial-distance 75 \
    --confidence-boost 0.1 \
    --same-cell-type-only
echo "Meta-rule generation complete. Meta-rules saved to $META_RULES_DIR"
echo

# Step 4: Apply rules to identify cells in spatial data
echo "Step 4: Applying rules to identify cells in spatial data..."
python scripts/identify_cells.py \
    --spatial-data $SPATIAL_OUTPUT \
    --rules-dir $RULES_DIR \
    --meta-rules-file $META_RULES_DIR/meta_rules.rq \
    --output-dir $RESULTS_DIR \
    --min-confidence 0.7 \
    --resolve-conflicts \
    --prioritize-meta-rules \
    --enhanced-resolution \
    --overlap-threshold 0.3 \
    --expression-data $META_RULES_DIR/gene_expression_data.json
echo "Cell identification complete. Results saved to $RESULTS_DIR"
echo

# Step 5: Visualize results
echo "Step 5: Generating visualization of identified cells..."
python scripts/visualize_cells.py \
    --input $RESULTS_DIR/identified_cells.csv \
    --output $RESULTS_DIR/cell_visualization.png \
    --plot-by-type \
    --highlight-meta-rule-cells
echo "Visualization complete. Output saved to $RESULTS_DIR/cell_visualization.png"
echo

echo "Pipeline completed successfully!"
echo "Results:"
echo "- Identified cells: $RESULTS_DIR/identified_cells.csv"
echo "- Cell visualization: $RESULTS_DIR/cell_visualization.png"
echo "- Rules summary: $RULES_DIR/rules_summary.txt"
echo "- Meta-rules summary: $META_RULES_DIR/meta_rules_summary.csv"
echo "- Rule network visualization: $META_RULES_DIR/rule_network.graphml"
echo "- Conflict resolution report: $RESULTS_DIR/conflict_resolution_report.md"

# Display summary statistics
echo
echo "Summary Statistics:"
CELL_COUNT=$(wc -l < $RESULTS_DIR/identified_cells.csv)
CELL_COUNT=$((CELL_COUNT - 1))  # Subtract header row
echo "- Total cells identified: $CELL_COUNT"

# Count cells by type
echo "- Cells by type:"
tail -n +2 $RESULTS_DIR/identified_cells.csv | cut -d, -f2 | sort | uniq -c | while read -r count type; do
    echo "  - $type: $count"
done

# Count cells by rule type
echo "- Cells by rule type:"
META_RULE_COUNT=$(grep "meta" $RESULTS_DIR/identified_cells.csv | wc -l)
CLIQUE_RULE_COUNT=$(grep "clique" $RESULTS_DIR/identified_cells.csv | grep -v "meta" | wc -l)
PAIR_RULE_COUNT=$(grep "pair" $RESULTS_DIR/identified_cells.csv | grep -v "meta" | wc -l)
echo "  - Meta-rules: $META_RULE_COUNT"
echo "  - Clique rules: $CLIQUE_RULE_COUNT"
echo "  - Pair rules: $PAIR_RULE_COUNT" 