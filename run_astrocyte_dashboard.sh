#!/bin/bash

# Run Astrocyte Spatial Dashboard
# This script activates the virtual environment and runs the astrocyte spatial dashboard

# Set the base directory to the script's location
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Define paths
SPATIAL_DIR="$BASE_DIR/Spatial_Data"
SINGLE_CELL_DIR="$BASE_DIR/Single_Cell"
SPATIAL_TTL="$SPATIAL_DIR/data/enhanced_spatial_ontology.ttl"
SINGLE_CELL_H5AD="$SINGLE_CELL_DIR/data/single_cell_data_cleaned.h5ad"
OUTPUT_DIR="$BASE_DIR/output/astrocyte_analysis"

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
python -c "import pandas, numpy, matplotlib, seaborn, rdflib, tkinter, sklearn" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip install pandas numpy matplotlib seaborn rdflib scikit-learn
fi

# Run the dashboard
echo "Running Astrocyte Spatial Dashboard..."
python "$SPATIAL_DIR/scripts/astrocyte_spatial_dashboard.py" --ttl "$SPATIAL_TTL" --output "$OUTPUT_DIR"

# Deactivate virtual environment
deactivate

echo "Dashboard execution complete." 