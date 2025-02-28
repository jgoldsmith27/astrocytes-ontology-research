#!/bin/bash

# Run Ontology Integration Tool for Astrocyte Analysis
# This script runs the ontology integration tool that creates semantic connections
# between spatial and single-cell ontologies for comprehensive analysis

# Set up environment
echo "Setting up environment..."
SPATIAL_DIR="astrocytes-ontology-research/Spatial_Data"
SINGLE_CELL_DIR="astrocytes-ontology-research/Single_Cell"
OUTPUT_DIR="astrocytes-ontology-research/output"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Activate virtual environment
if [ -d "$SPATIAL_DIR/venv" ]; then
    echo "Activating Spatial Data virtual environment..."
    source $SPATIAL_DIR/venv/bin/activate
elif [ -d "$SINGLE_CELL_DIR/venv" ]; then
    echo "Activating Single Cell virtual environment..."
    source $SINGLE_CELL_DIR/venv/bin/activate
else
    echo "No virtual environment found. Please create one first."
    exit 1
fi

# Check for required Python packages
echo "Checking required packages..."
python -c "import rdflib" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip install rdflib
fi

# Define paths to ontology files
SPATIAL_ONTOLOGY="$SPATIAL_DIR/output/enhanced_spatial_ontology.ttl"
SINGLE_CELL_ONTOLOGY="$SINGLE_CELL_DIR/output/single_cell_ontology.ttl"
INTEGRATED_ONTOLOGY="$OUTPUT_DIR/integrated_ontology.ttl"

# Check if ontology files exist
if [ ! -f "$SPATIAL_ONTOLOGY" ]; then
    echo "Spatial ontology file not found at $SPATIAL_ONTOLOGY"
    echo "Please run the spatial ontology builder first."
    exit 1
fi

if [ ! -f "$SINGLE_CELL_ONTOLOGY" ]; then
    echo "Single-cell ontology file not found at $SINGLE_CELL_ONTOLOGY"
    echo "Please run the single-cell ontology builder first."
    exit 1
fi

# Run the ontology integration tool
echo "Running ontology integration tool..."
python $SPATIAL_DIR/scripts/ontology_integration.py \
    --spatial "$SPATIAL_ONTOLOGY" \
    --single-cell "$SINGLE_CELL_ONTOLOGY" \
    --output "$INTEGRATED_ONTOLOGY"

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Ontology integration completed successfully!"
    echo "Integrated ontology saved to $INTEGRATED_ONTOLOGY"
    echo ""
    echo "You can now use this integrated ontology for advanced cross-ontology queries."
    echo "Example usage:"
    echo "  - Find genes expressed in both modalities"
    echo "  - Identify spatial regions with specific cell types"
    echo "  - Discover differentially expressed genes between astrocyte types"
else
    echo "Error running ontology integration tool."
    exit 1
fi 