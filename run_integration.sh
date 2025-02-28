#!/bin/bash

# Run Spatial and Single-Cell Integration
# This script activates the virtual environment and runs the integration between spatial and single-cell data

# Set the base directory to the script's location
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Define paths
SPATIAL_DIR="$BASE_DIR/Spatial_Data"
SINGLE_CELL_DIR="$BASE_DIR/Single_Cell"
SPATIAL_TTL="$SPATIAL_DIR/data/enhanced_spatial_ontology.ttl"
SINGLE_CELL_H5AD="$SINGLE_CELL_DIR/data/single_cell_data_cleaned.h5ad"
OUTPUT_DIR="$BASE_DIR/output/integrated_analysis"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if virtual environment exists
if [ -d "$SPATIAL_DIR/venv" ]; then
    VENV_DIR="$SPATIAL_DIR/venv"
elif [ -d "$SINGLE_CELL_DIR/venv" ]; then
    VENV_DIR="$SINGLE_CELL_DIR/venv"
else
    echo "Error: Virtual environment not found in either Spatial_Data or Single_Cell directories."
    exit 1
fi

# Activate virtual environment
echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

# Check if required packages are installed
echo "Checking required packages..."
python -c "import pandas, numpy, matplotlib, seaborn, rdflib, scanpy, anndata, sklearn, scipy" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip install pandas numpy matplotlib seaborn rdflib scanpy anndata scikit-learn scipy
fi

# Run the integration
echo "Running Spatial and Single-Cell Integration..."
python "$BASE_DIR/scripts/integrate_spatial_single_cell.py" --spatial "$SPATIAL_TTL" --single-cell "$SINGLE_CELL_H5AD" --output "$OUTPUT_DIR"

# Deactivate virtual environment
deactivate

echo "Integration complete. Results saved to $OUTPUT_DIR" 