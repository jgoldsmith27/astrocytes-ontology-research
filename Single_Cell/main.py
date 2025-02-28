from read_data import read_data
from single_cell_ontology_builder import GeneExpressionOntologyBuilder
import pandas as pd
import os

def main():
    # Define file paths
    file_path = "/Users/jacob/Desktop/Karolinska Research/astrocytes-ontology-research/single_cell_data_cleaned.h5ad"
    output_ttl = "astrocyte_ontology.ttl"
    output_csv = "astrocyte_expression.csv"
    
    try:
        # Step 1: Read the data and compute average expression for astrocyte types
        print("🔍 Reading data and filtering for astrocyte subtypes...")
        avg_df = read_data(file_path)
        
        # Save the DataFrame
        avg_df.to_csv(output_csv, float_format='%.3f')
        print(f"💾 Saved astrocyte expression data to {output_csv} ({avg_df.shape[0]} astrocyte types, {avg_df.shape[1]} genes)")
        
        # Step 2: Build the ontology
        print("🧠 Building astrocyte ontology...")
        builder = GeneExpressionOntologyBuilder()
        builder.create_classes()
        builder.add_data(avg_df)
        builder.save_to_file(output_ttl)
        print(f"✅ Saved astrocyte ontology to {output_ttl}")
        
    except Exception as e:
        print(f"❌ Error: {str(e)}")

if __name__ == "__main__":
    main()
