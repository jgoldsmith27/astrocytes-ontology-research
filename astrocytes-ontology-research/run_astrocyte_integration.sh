#!/bin/bash
# Astrocyte-Focused Ontology Integration and Visualization Script
# This script automates the process of integrating astrocyte data from spatial and single-cell ontologies
# and generating visualizations from the integrated ontology.

# Activate virtual environment if it exists
if [ -d "venv" ]; then
    echo "Activating virtual environment..."
    source venv/bin/activate
fi

# Check for required packages
echo "Checking for required packages..."
python -c "import rdflib, pandas, numpy, matplotlib, seaborn, scanpy, anndata, networkx, matplotlib_venn" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip install rdflib pandas numpy matplotlib seaborn scanpy anndata networkx matplotlib-venn
fi

# Create output directories
echo "Creating output directories..."
mkdir -p output/visualizations
mkdir -p output/astrocyte_validation

# Step 1: Run the ontology integration
echo "Step 1: Running astrocyte-focused ontology integration..."
python Spatial_Data/scripts/core/ontology_integration.py \
    --spatial output/enhanced_spatial_ontology.ttl \
    --single-cell Single_Cell/output/single_cell_ontology.ttl \
    --output output/integrated_astrocyte_ontology.ttl

# Check if integration was successful
if [ $? -ne 0 ]; then
    echo "Error: Ontology integration failed."
    exit 1
fi

# Step 2: Run the visualization
echo "Step 2: Generating astrocyte-specific visualizations..."
python Spatial_Data/scripts/core/astrocyte_visualization.py \
    --input output/integrated_astrocyte_ontology.ttl \
    --output-dir output/visualizations

# Check if visualization was successful
if [ $? -ne 0 ]; then
    echo "Error: Visualization generation failed."
    exit 1
fi

# Step 3: Run the validation
echo "Step 3: Validating astrocyte cell type classifications..."
python Spatial_Data/scripts/core/astrocyte_validation.py \
    --input output/integrated_astrocyte_ontology.ttl \
    --output-dir output/astrocyte_validation

# Check if validation was successful
if [ $? -ne 0 ]; then
    echo "Error: Astrocyte validation failed."
    exit 1
fi

echo "Astrocyte integration, visualization, and validation completed successfully!"
echo "Integrated ontology saved to: output/integrated_astrocyte_ontology.ttl"
echo "Visualizations saved to: output/visualizations/"
echo "Validation results saved to: output/astrocyte_validation/"

# List generated visualizations
echo "Generated visualizations:"
ls -l output/visualizations/

# List validation results
echo "Validation results:"
ls -l output/astrocyte_validation/

echo "To view the README for the astrocyte-focused integration, see: output/README_astrocyte_integration.md"
