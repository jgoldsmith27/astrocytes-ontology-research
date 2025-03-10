#!/usr/bin/env python3
"""
Script to convert spatial transcriptomics data to TURTLE (RDF) format.

This script reads the spatial transcriptomics dataset and converts it into
TURTLE format using the cell type ontology.
"""

import os
import pandas as pd
import numpy as np
from rdflib import Graph, Namespace, Literal, URIRef, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL
import argparse
from tqdm import tqdm

# Define namespaces
CELL = Namespace("http://example.org/ontology/cell/")
EX = Namespace("http://example.org/data/spatial/")

def create_spatial_rdf(input_file, output_file, limit=None):
    """
    Convert spatial transcriptomics data to RDF using the cell type ontology.
    
    Parameters:
    -----------
    input_file : str
        Path to the input CSV file containing spatial data
    output_file : str
        Path to the output TURTLE file
    limit : int, optional
        Limit the number of rows to process (for testing)
    """
    print(f"Reading spatial data from {input_file}...")
    
    # Create a new RDF graph
    g = Graph()
    
    # Bind namespaces for readability
    g.bind("cell", CELL)
    g.bind("ex", EX)
    g.bind("rdf", RDF)
    g.bind("rdfs", RDFS)
    g.bind("xsd", XSD)
    g.bind("owl", OWL)
    
    # Read the spatial data
    # Using chunksize for memory efficiency with large files
    chunk_size = 100000
    reader = pd.read_csv(input_file, chunksize=chunk_size)
    
    # Track unique genes and bins to avoid duplicates
    unique_genes = set()
    unique_bins = set()
    
    # Process data in chunks
    row_count = 0
    for chunk_idx, chunk in enumerate(reader):
        print(f"Processing chunk {chunk_idx+1}...")
        
        if limit and row_count >= limit:
            break
        
        # Process each row in the chunk
        for idx, row in tqdm(chunk.iterrows(), total=len(chunk)):
            if limit and row_count >= limit:
                break
                
            # Create URIs for the entities
            data_point_uri = EX[f"datapoint_{row['Unnamed: 0']}"]
            gene_id = str(row['geneID'])
            gene_uri = EX[f"gene_{gene_id}"]
            bin_id = int(row['bin1_ID'])
            bin_uri = EX[f"bin_{bin_id}"]
            
            # Add the data point
            g.add((data_point_uri, RDF.type, CELL.SpatialDataPoint))
            g.add((data_point_uri, CELL.hasXCoordinate, Literal(int(row['x']), datatype=XSD.integer)))
            g.add((data_point_uri, CELL.hasYCoordinate, Literal(int(row['y']), datatype=XSD.integer)))
            g.add((data_point_uri, CELL.hasMIDCount, Literal(int(row['MIDCount']), datatype=XSD.integer)))
            g.add((data_point_uri, CELL.hasExonCount, Literal(int(row['ExonCount']), datatype=XSD.integer)))
            g.add((data_point_uri, CELL.hasIntronCount, Literal(int(row['IntronCount']), datatype=XSD.integer)))
            
            # Add gene information (only if not already added)
            if gene_id not in unique_genes:
                g.add((gene_uri, RDF.type, CELL.Gene))
                g.add((gene_uri, CELL.hasGeneID, Literal(gene_id, datatype=XSD.string)))
                unique_genes.add(gene_id)
            
            # Link data point to gene
            g.add((data_point_uri, CELL.expressesGene, gene_uri))
            
            # Add bin information (only if not already added)
            if bin_id not in unique_bins:
                g.add((bin_uri, RDF.type, CELL.SpatialBin))
                g.add((bin_uri, CELL.hasBinID, Literal(bin_id, datatype=XSD.integer)))
                unique_bins.add(bin_id)
            
            # Link data point to bin
            g.add((data_point_uri, CELL.locatedInBin, bin_uri))
            
            row_count += 1
    
    print(f"Processed {row_count} data points, {len(unique_genes)} unique genes, and {len(unique_bins)} unique bins.")
    
    # Write the graph to a TURTLE file
    print(f"Writing RDF to {output_file}...")
    g.serialize(destination=output_file, format="turtle")
    print(f"RDF graph written to {output_file}")
    print(f"Total triples: {len(g)}")

def main():
    parser = argparse.ArgumentParser(description="Convert spatial transcriptomics data to TURTLE format")
    parser.add_argument("--input", "-i", required=True, help="Input CSV file")
    parser.add_argument("--output", "-o", required=True, help="Output TURTLE file")
    parser.add_argument("--limit", "-l", type=int, help="Limit number of rows (for testing)")
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Convert data to RDF
    create_spatial_rdf(args.input, args.output, args.limit)

if __name__ == "__main__":
    main() 