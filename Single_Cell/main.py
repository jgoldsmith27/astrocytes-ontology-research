from optimized_pipeline.Single_Cell.read_data import read_data
from single_cell_ontology_builder import GeneExpressionOntologyBuilder

def main():
    # Define file paths
    file_path = "/Users/jacob/Desktop/Karolinska Research/optimized_pipeline/single_cell_data_cleaned.h5ad"
    output_path = "output_ontology.ttl"

    # Read the data and compute average expression per cell type
    avg_df = read_data(file_path)

    # Build the ontology
    builder = GeneExpressionOntologyBuilder()
    builder.create_classes()
    builder.add_data(avg_df)
    builder.save_to_file(output_path)

if __name__ == "__main__":
    main()