#!/bin/bash

# Run Astrocyte Classifier
# This script activates the virtual environment and runs the astrocyte classifier
# to classify astrocytes, apply the classifier to all spatial points, and enhance
# the spatial ontology with cell type information.

# Set the base directory to the script's location
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Define paths
SPATIAL_DIR="$BASE_DIR/Spatial_Data"
SINGLE_CELL_DIR="$BASE_DIR/Single_Cell"
SPATIAL_TTL="$SPATIAL_DIR/data/enhanced_spatial_ontology.ttl"
OUTPUT_TTL="$SPATIAL_DIR/data/astrocyte_enhanced_spatial_ontology.ttl"
OUTPUT_DIR="$BASE_DIR/output/astrocyte_classification"

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
python -c "import pandas, numpy, matplotlib, seaborn, rdflib, sklearn, joblib" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip install pandas numpy matplotlib seaborn rdflib scikit-learn joblib
fi

# Run the classifier
echo "Running Astrocyte Classifier..."
python "$SPATIAL_DIR/scripts/astrocyte_classifier.py" --input "$SPATIAL_TTL" --output "$OUTPUT_TTL" --output-dir "$OUTPUT_DIR"

# Check if the classifier ran successfully
if [ $? -eq 0 ]; then
    echo "Astrocyte classification completed successfully!"
    echo "Enhanced ontology saved to: $OUTPUT_TTL"
    echo "Classification results saved to: $OUTPUT_DIR"
else
    echo "Error: Astrocyte classification failed."
    exit 1
fi

# Deactivate virtual environment
deactivate

echo "Classification process complete." 