#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef, Literal, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL
import argparse
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import joblib

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
CELL = Namespace("http://example.org/cell-ontology#")

# Define astrocyte marker genes for different types
# Based on literature and single-cell data
ASTROCYTE_MARKERS = {
    'Protoplasmic': ['GFAP', 'S100B', 'AQP4', 'GJA1', 'SLC1A2', 'SLC1A3'],
    'Fibrous': ['GFAP', 'VIM', 'ALDH1L1', 'CD44', 'CRYAB', 'HOPX'],
    'Reactive': ['GFAP', 'VIM', 'SERPINA3', 'C3', 'CXCL10', 'STAT3', 'LCN2']
}

class AstrocyteClassifier:
    """
    A class for classifying astrocytes in spatial transcriptomics data and
    enhancing the spatial ontology with cell type information.
    """
    
    def __init__(self, input_ttl, output_ttl, output_dir="../output/astrocyte_classification"):
        """
        Initialize the classifier with input and output file paths.
        
        Parameters:
        -----------
        input_ttl : str
            Path to the input TURTLE file containing the spatial ontology
        output_ttl : str
            Path to the output TURTLE file where the enhanced ontology will be saved
        output_dir : str
            Path to the directory where analysis results will be saved
        """
        self.input_ttl = input_ttl
        self.output_ttl = output_ttl
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize data structures
        self.graph = None
        self.gene_coordinates = {}
        self.spatial_points = {}
        self.point_gene_matrix = None
        self.classifier = None
        self.astrocyte_data = None
        
        # Load data
        self.load_data()
    
    def load_data(self):
        """
        Load the spatial ontology data from a TURTLE file.
        """
        print(f"Loading TURTLE file: {self.input_ttl}")
        self.graph = Graph()
        self.graph.parse(self.input_ttl, format="turtle")
        print(f"Graph loaded with {len(self.graph)} triples")
        
        # Extract gene coordinates and spatial points
        self.extract_data_from_ontology()
    
    def extract_data_from_ontology(self):
        """
        Extract gene coordinates and spatial points from the ontology.
        """
        # Extract gene coordinates
        query_genes = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?geneID ?point ?x ?y ?expressionLevel
        WHERE {
            ?gene rdf:type :Gene ;
                  :expressedAt ?point ;
                  :hasExpressionLevel ?expressionLevel .
            
            ?point :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        """
        
        results = self.graph.query(query_genes)
        
        # Store gene coordinates
        gene_coords = {}
        for row in results:
            gene_id = str(row.geneID)
            gene_uri = row.gene
            point_uri = row.point
            x = float(row.x)
            y = float(row.y)
            expr = float(row.expressionLevel)
            
            if gene_id not in gene_coords:
                gene_coords[gene_id] = []
            
            gene_coords[gene_id].append((gene_uri, point_uri, x, y, expr))
            
            # Also store by spatial point
            if point_uri not in self.spatial_points:
                self.spatial_points[point_uri] = {
                    'x': x,
                    'y': y,
                    'genes': {}
                }
            
            self.spatial_points[point_uri]['genes'][gene_id] = {
                'gene_uri': gene_uri,
                'expression': expr
            }
        
        self.gene_coordinates = gene_coords
        
        print(f"Extracted {len(gene_coords)} genes and {len(self.spatial_points)} spatial points")
    
    def create_point_gene_matrix(self):
        """
        Create a matrix of spatial points and gene expression values.
        """
        print("Creating point-gene expression matrix...")
        
        # Get all unique gene IDs
        all_genes = sorted(self.gene_coordinates.keys())
        
        # Create a DataFrame with spatial points as rows and genes as columns
        data = []
        for point_uri, point_data in self.spatial_points.items():
            row = {
                'point_uri': point_uri,
                'x': point_data['x'],
                'y': point_data['y']
            }
            
            # Add gene expression values
            for gene_id in all_genes:
                if gene_id in point_data['genes']:
                    row[gene_id] = point_data['genes'][gene_id]['expression']
                else:
                    row[gene_id] = 0
            
            data.append(row)
        
        # Convert to DataFrame
        self.point_gene_matrix = pd.DataFrame(data)
        
        print(f"Created point-gene matrix with {len(self.point_gene_matrix)} points and {len(all_genes)} genes")
    
    def calculate_astrocyte_scores(self):
        """
        Calculate scores for each astrocyte type based on marker gene expression.
        """
        print("Calculating astrocyte type scores...")
        
        if self.point_gene_matrix is None:
            self.create_point_gene_matrix()
        
        # Calculate scores for each astrocyte type
        for cell_type, markers in ASTROCYTE_MARKERS.items():
            # Initialize score column
            self.point_gene_matrix[f'{cell_type}_score'] = 0
            
            # Sum expression of marker genes (if present in data)
            for marker in markers:
                if marker in self.point_gene_matrix.columns:
                    self.point_gene_matrix[f'{cell_type}_score'] += self.point_gene_matrix[marker].fillna(0)
            
            # Normalize by number of markers found
            markers_found = sum(1 for marker in markers if marker in self.point_gene_matrix.columns)
            if markers_found > 0:
                self.point_gene_matrix[f'{cell_type}_score'] /= markers_found
        
        # Determine dominant astrocyte type for each point
        type_columns = [f'{cell_type}_score' for cell_type in ASTROCYTE_MARKERS.keys()]
        self.point_gene_matrix['dominant_type'] = self.point_gene_matrix[type_columns].idxmax(axis=1)
        self.point_gene_matrix['dominant_type'] = self.point_gene_matrix['dominant_type'].str.replace('_score', '')
        
        # Store results
        self.astrocyte_data = self.point_gene_matrix.copy()
        
        print(f"Classified {len(self.astrocyte_data)} spatial points into astrocyte types")
        for cell_type in ASTROCYTE_MARKERS.keys():
            count = sum(self.astrocyte_data['dominant_type'] == cell_type)
            print(f"  - {cell_type}: {count} points")
    
    def train_classifier(self):
        """
        Train a machine learning classifier for astrocyte types.
        """
        print("Training astrocyte classifier...")
        
        if self.astrocyte_data is None:
            self.calculate_astrocyte_scores()
        
        # Get all marker genes
        all_markers = sorted(set(gene for genes in ASTROCYTE_MARKERS.values() for gene in genes))
        
        # Filter for marker genes that are present in the data
        present_markers = [marker for marker in all_markers if marker in self.astrocyte_data.columns]
        
        print(f"Using {len(present_markers)} marker genes for classification: {', '.join(present_markers)}")
        
        # Prepare features and target
        X = self.astrocyte_data[present_markers].values
        y = self.astrocyte_data['dominant_type'].values
        
        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        
        # Standardize features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Train a Random Forest classifier
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(X_train_scaled, y_train)
        
        # Evaluate the classifier
        y_pred = clf.predict(X_test_scaled)
        
        print("\nClassification Report:")
        print(classification_report(y_test, y_pred))
        
        print("\nConfusion Matrix:")
        print(confusion_matrix(y_test, y_pred))
        
        # Store the classifier and scaler
        self.classifier = {
            'model': clf,
            'scaler': scaler,
            'features': present_markers
        }
        
        # Save the classifier
        model_file = os.path.join(self.output_dir, "astrocyte_classifier.joblib")
        joblib.dump(self.classifier, model_file)
        print(f"Classifier saved to {model_file}")
    
    def apply_classifier_to_all_points(self):
        """
        Apply the trained classifier to all spatial points.
        """
        print("Applying classifier to all spatial points...")
        
        if self.classifier is None:
            self.train_classifier()
        
        # Get features for all points
        features = self.classifier['features']
        X = self.astrocyte_data[features].values
        
        # Scale features
        X_scaled = self.classifier['scaler'].transform(X)
        
        # Predict astrocyte types
        y_pred = self.classifier['model'].predict(X_scaled)
        
        # Get prediction probabilities
        y_prob = self.classifier['model'].predict_proba(X_scaled)
        
        # Add predictions to the data
        self.astrocyte_data['predicted_type'] = y_pred
        
        # Add prediction probabilities
        for i, cell_type in enumerate(self.classifier['model'].classes_):
            self.astrocyte_data[f'{cell_type}_probability'] = y_prob[:, i]
        
        # Compare with dominant type
        match_count = sum(self.astrocyte_data['dominant_type'] == self.astrocyte_data['predicted_type'])
        match_percentage = match_count / len(self.astrocyte_data) * 100
        
        print(f"Classifier predictions match dominant type for {match_count} points ({match_percentage:.2f}%)")
        
        # Save the classified data
        data_file = os.path.join(self.output_dir, "astrocyte_classified_points.csv")
        self.astrocyte_data.to_csv(data_file, index=False)
        print(f"Classified data saved to {data_file}")
    
    def enhance_ontology(self):
        """
        Enhance the spatial ontology with astrocyte type information.
        """
        print("Enhancing spatial ontology with astrocyte type information...")
        
        if 'predicted_type' not in self.astrocyte_data.columns:
            self.apply_classifier_to_all_points()
        
        # Add cell type namespace to the graph
        self.graph.bind("cell", CELL)
        
        # Add cell type classes to the ontology
        for cell_type in ASTROCYTE_MARKERS.keys():
            cell_type_uri = CELL[cell_type + "Astrocyte"]
            
            # Add cell type class
            self.graph.add((cell_type_uri, RDF.type, OWL.Class))
            self.graph.add((cell_type_uri, RDFS.subClassOf, CELL.Astrocyte))
            self.graph.add((cell_type_uri, RDFS.label, Literal(f"{cell_type} Astrocyte", datatype=XSD.string)))
            
            # Add description
            if cell_type == "Protoplasmic":
                description = "Astrocytes found mainly in gray matter with highly branched processes and high expression of glutamate transporters."
            elif cell_type == "Fibrous":
                description = "Astrocytes found mainly in white matter with long, unbranched processes and high expression of intermediate filaments."
            elif cell_type == "Reactive":
                description = "Astrocytes found in response to injury or disease, characterized by hypertrophy and upregulation of inflammatory genes."
            else:
                description = f"{cell_type} type of astrocyte."
            
            self.graph.add((cell_type_uri, RDFS.comment, Literal(description, datatype=XSD.string)))
            
            # Add marker genes
            for marker in ASTROCYTE_MARKERS[cell_type]:
                marker_uri = URIRef(f"http://example.org/gene-ontology#{marker}")
                self.graph.add((cell_type_uri, CELL.hasMarkerGene, marker_uri))
                self.graph.add((marker_uri, RDF.type, CELL.MarkerGene))
                self.graph.add((marker_uri, RDFS.label, Literal(marker, datatype=XSD.string)))
        
        # Add cell type information to spatial points
        for _, row in self.astrocyte_data.iterrows():
            point_uri = row['point_uri']
            cell_type = row['predicted_type']
            probability = row[f'{cell_type}_probability']
            
            # Create a cell instance
            cell_instance_uri = BNode()
            
            # Add cell type information
            self.graph.add((cell_instance_uri, RDF.type, CELL[cell_type + "Astrocyte"]))
            self.graph.add((cell_instance_uri, CELL.hasConfidence, Literal(probability, datatype=XSD.float)))
            
            # Link cell to spatial point
            self.graph.add((point_uri, CELL.containsCell, cell_instance_uri))
            self.graph.add((cell_instance_uri, CELL.locatedAt, point_uri))
            
            # Add scores for each cell type
            for ct in ASTROCYTE_MARKERS.keys():
                score = row[f'{ct}_score']
                self.graph.add((cell_instance_uri, CELL[f"has{ct}Score"], Literal(score, datatype=XSD.float)))
        
        # Add object properties
        self.graph.add((CELL.containsCell, RDF.type, OWL.ObjectProperty))
        self.graph.add((CELL.containsCell, RDFS.domain, ST.SpatialPoint))
        self.graph.add((CELL.containsCell, RDFS.range, CELL.Cell))
        
        self.graph.add((CELL.locatedAt, RDF.type, OWL.ObjectProperty))
        self.graph.add((CELL.locatedAt, RDFS.domain, CELL.Cell))
        self.graph.add((CELL.locatedAt, RDFS.range, ST.SpatialPoint))
        
        self.graph.add((CELL.hasMarkerGene, RDF.type, OWL.ObjectProperty))
        self.graph.add((CELL.hasMarkerGene, RDFS.domain, OWL.Class))
        self.graph.add((CELL.hasMarkerGene, RDFS.range, CELL.MarkerGene))
        
        # Add data properties
        self.graph.add((CELL.hasConfidence, RDF.type, OWL.DatatypeProperty))
        self.graph.add((CELL.hasConfidence, RDFS.domain, CELL.Cell))
        self.graph.add((CELL.hasConfidence, RDFS.range, XSD.float))
        
        for cell_type in ASTROCYTE_MARKERS.keys():
            score_property = CELL[f"has{cell_type}Score"]
            self.graph.add((score_property, RDF.type, OWL.DatatypeProperty))
            self.graph.add((score_property, RDFS.domain, CELL.Cell))
            self.graph.add((score_property, RDFS.range, XSD.float))
        
        # Save the enhanced ontology
        self.graph.serialize(destination=self.output_ttl, format="turtle")
        print(f"Enhanced ontology saved to {self.output_ttl}")
    
    def visualize_astrocyte_distribution(self):
        """
        Visualize the spatial distribution of different astrocyte types.
        """
        print("Visualizing astrocyte distribution...")
        
        if 'predicted_type' not in self.astrocyte_data.columns:
            self.apply_classifier_to_all_points()
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Define colors for each astrocyte type
        colors = {
            'Protoplasmic': 'blue',
            'Fibrous': 'red',
            'Reactive': 'green'
        }
        
        # Plot each astrocyte type
        for cell_type, color in colors.items():
            # Filter data for this cell type
            type_data = self.astrocyte_data[self.astrocyte_data['predicted_type'] == cell_type]
            
            if not type_data.empty:
                # Plot points
                ax.scatter(
                    type_data['x'], 
                    type_data['y'],
                    c=color,
                    s=type_data[f'{cell_type}_probability'] * 100,  # Size based on probability
                    alpha=0.7,
                    label=f'{cell_type} Astrocytes'
                )
        
        # Add legend and labels
        ax.legend(fontsize=12)
        ax.set_xlabel('X Coordinate', fontsize=14)
        ax.set_ylabel('Y Coordinate', fontsize=14)
        ax.set_title('Spatial Distribution of Astrocyte Types', fontsize=16)
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Save the figure
        fig_file = os.path.join(self.output_dir, "astrocyte_distribution.png")
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        print(f"Visualization saved to {fig_file}")
    
    def run_pipeline(self):
        """
        Run the full pipeline: classify astrocytes, apply to all points, enhance ontology.
        """
        self.create_point_gene_matrix()
        self.calculate_astrocyte_scores()
        self.train_classifier()
        self.apply_classifier_to_all_points()
        self.enhance_ontology()
        self.visualize_astrocyte_distribution()
        
        print("Pipeline completed successfully!")

def main():
    """
    Main function to run the astrocyte classifier.
    """
    parser = argparse.ArgumentParser(description="Classify astrocytes in spatial transcriptomics data and enhance the ontology")
    parser.add_argument("--input", type=str, required=True, help="Path to the input TURTLE file containing the spatial ontology")
    parser.add_argument("--output", type=str, required=True, help="Path to the output TURTLE file where the enhanced ontology will be saved")
    parser.add_argument("--output-dir", type=str, default="../output/astrocyte_classification", help="Path to the output directory")
    
    args = parser.parse_args()
    
    classifier = AstrocyteClassifier(
        input_ttl=args.input,
        output_ttl=args.output,
        output_dir=args.output_dir
    )
    
    classifier.run_pipeline()

if __name__ == "__main__":
    main() 