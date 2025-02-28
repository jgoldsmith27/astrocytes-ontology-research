#!/bin/bash

# Run Integrated Visualization Tool for Astrocyte Analysis
# This script runs the integrated visualization tool that combines
# spatial and single-cell data for comprehensive astrocyte analysis

# Set up environment
echo "Setting up environment..."
SPATIAL_DIR="astrocytes-ontology-research/Spatial_Data"
SINGLE_CELL_DIR="astrocytes-ontology-research/Single_Cell"
OUTPUT_DIR="astrocytes-ontology-research/output/integrated_visualization"

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
python -c "import dash, plotly, dash_bootstrap_components" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip install dash plotly dash-bootstrap-components
fi

# Define paths to ontology files
SPATIAL_ONTOLOGY="$SPATIAL_DIR/output/enhanced_spatial_ontology.ttl"
SINGLE_CELL_ONTOLOGY="$SINGLE_CELL_DIR/output/single_cell_ontology.ttl"

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

# Run the integrated visualization tool
echo "Running integrated visualization tool..."
python $SPATIAL_DIR/scripts/astrocyte_integrated_visualization.py \
    --spatial "$SPATIAL_ONTOLOGY" \
    --single-cell "$SINGLE_CELL_ONTOLOGY" \
    --output "$OUTPUT_DIR" \
    --interactive

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Visualization completed successfully!"
    echo "Results saved to $OUTPUT_DIR"
    echo "Interactive dashboard is running at http://127.0.0.1:8050/"
    echo "Press Ctrl+C to stop the dashboard when finished."
else
    echo "Error running visualization tool."
    exit 1
fi 