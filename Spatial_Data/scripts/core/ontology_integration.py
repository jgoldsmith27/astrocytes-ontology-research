#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
from rdflib import Graph, Namespace, URIRef, Literal, BNode
from rdflib.namespace import RDF, RDFS, XSD, OWL
import argparse

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")
SC = Namespace("http://example.org/single-cell#")
CELL = Namespace("http://example.org/cell-ontology#")
GENE = Namespace("http://example.org/gene-ontology#")
BRIDGE = Namespace("http://example.org/ontology-bridge#")

class OntologyIntegration:
    """
    A class for semantically integrating spatial and single-cell ontologies
    to enable advanced cross-ontology queries and knowledge discovery.
    """
    
    def __init__(self, spatial_ttl, single_cell_ttl, output_ttl="../output/integrated_ontology.ttl"):
        """
        Initialize the integration tool with input and output file paths.
        
        Parameters:
        -----------
        spatial_ttl : str
            Path to the TURTLE file containing the spatial ontology
        single_cell_ttl : str
            Path to the TURTLE file containing the single-cell ontology
        output_ttl : str
            Path to the output TURTLE file where the integrated ontology will be saved
        """
        self.spatial_ttl = spatial_ttl
        self.single_cell_ttl = single_cell_ttl
        self.output_ttl = output_ttl
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(self.output_ttl), exist_ok=True)
        
        # Initialize graphs
        self.spatial_graph = None
        self.single_cell_graph = None
        self.integrated_graph = None
        
        # Load data
        self.load_data()
    
    def load_data(self):
        """
        Load the spatial and single-cell ontology data from TURTLE files.
        """
        print(f"Loading spatial TURTLE file: {self.spatial_ttl}")
        self.spatial_graph = Graph()
        self.spatial_graph.parse(self.spatial_ttl, format="turtle")
        print(f"Spatial graph loaded with {len(self.spatial_graph)} triples")
        
        print(f"Loading single-cell TURTLE file: {self.single_cell_ttl}")
        self.single_cell_graph = Graph()
        self.single_cell_graph.parse(self.single_cell_ttl, format="turtle")
        print(f"Single-cell graph loaded with {len(self.single_cell_graph)} triples")
        
        # Initialize integrated graph
        self.integrated_graph = Graph()
        
        # Add all triples from both graphs
        for triple in self.spatial_graph:
            self.integrated_graph.add(triple)
        
        for triple in self.single_cell_graph:
            self.integrated_graph.add(triple)
        
        print(f"Integrated graph initialized with {len(self.integrated_graph)} triples")
    
    def define_bridge_ontology(self):
        """
        Define a bridge ontology with properties and classes that connect
        the spatial and single-cell ontologies.
        """
        print("Defining bridge ontology...")
        
        # Define bridge properties
        self.integrated_graph.add((BRIDGE.hasSpatialRepresentation, RDF.type, OWL.ObjectProperty))
        self.integrated_graph.add((BRIDGE.hasSpatialRepresentation, RDFS.domain, SC.Cell))
        self.integrated_graph.add((BRIDGE.hasSpatialRepresentation, RDFS.range, ST.SpatialPoint))
        self.integrated_graph.add((BRIDGE.hasSpatialRepresentation, RDFS.comment, Literal("Links a cell from single-cell data to its spatial representation")))
        
        self.integrated_graph.add((BRIDGE.hasSingleCellData, RDF.type, OWL.ObjectProperty))
        self.integrated_graph.add((BRIDGE.hasSingleCellData, RDFS.domain, ST.SpatialPoint))
        self.integrated_graph.add((BRIDGE.hasSingleCellData, RDFS.range, SC.Cell))
        self.integrated_graph.add((BRIDGE.hasSingleCellData, RDFS.comment, Literal("Links a spatial point to its single-cell data")))
        
        self.integrated_graph.add((BRIDGE.hasCorrespondingGene, RDF.type, OWL.ObjectProperty))
        self.integrated_graph.add((BRIDGE.hasCorrespondingGene, RDFS.domain, SC.Gene))
        self.integrated_graph.add((BRIDGE.hasCorrespondingGene, RDFS.range, ST.Gene))
        self.integrated_graph.add((BRIDGE.hasCorrespondingGene, RDFS.comment, Literal("Links a gene in single-cell data to its corresponding gene in spatial data")))
        
        self.integrated_graph.add((BRIDGE.hasCorrespondingCellType, RDF.type, OWL.ObjectProperty))
        self.integrated_graph.add((BRIDGE.hasCorrespondingCellType, RDFS.domain, SC.CellType))
        self.integrated_graph.add((BRIDGE.hasCorrespondingCellType, RDFS.range, ST.AstrocyteType))
        self.integrated_graph.add((BRIDGE.hasCorrespondingCellType, RDFS.comment, Literal("Links a cell type in single-cell data to its corresponding type in spatial data")))
        
        # Define equivalence classes for astrocyte types
        self.integrated_graph.add((BRIDGE.ProtoplasticAstrocyte, RDF.type, OWL.Class))
        self.integrated_graph.add((BRIDGE.ProtoplasticAstrocyte, RDFS.subClassOf, CELL.Astrocyte))
        self.integrated_graph.add((BRIDGE.ProtoplasticAstrocyte, OWL.equivalentClass, URIRef(ST + "Protoplasmic")))
        self.integrated_graph.add((BRIDGE.ProtoplasticAstrocyte, OWL.equivalentClass, URIRef(SC + "Protoplasmic")))
        
        self.integrated_graph.add((BRIDGE.FibrousAstrocyte, RDF.type, OWL.Class))
        self.integrated_graph.add((BRIDGE.FibrousAstrocyte, RDFS.subClassOf, CELL.Astrocyte))
        self.integrated_graph.add((BRIDGE.FibrousAstrocyte, OWL.equivalentClass, URIRef(ST + "Fibrous")))
        self.integrated_graph.add((BRIDGE.FibrousAstrocyte, OWL.equivalentClass, URIRef(SC + "Fibrous")))
        
        self.integrated_graph.add((BRIDGE.ReactiveAstrocyte, RDF.type, OWL.Class))
        self.integrated_graph.add((BRIDGE.ReactiveAstrocyte, RDFS.subClassOf, CELL.Astrocyte))
        self.integrated_graph.add((BRIDGE.ReactiveAstrocyte, OWL.equivalentClass, URIRef(ST + "Reactive")))
        self.integrated_graph.add((BRIDGE.ReactiveAstrocyte, OWL.equivalentClass, URIRef(SC + "Reactive")))
        
        print("Bridge ontology defined")
    
    def map_genes_between_ontologies(self):
        """
        Create explicit mappings between genes in the spatial and single-cell ontologies.
        """
        print("Mapping genes between ontologies...")
        
        # Get all genes from spatial ontology
        spatial_genes_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?geneID
        WHERE {
            ?gene rdf:type st:Gene .
            ?gene st:hasGeneID ?geneID .
        }
        """
        
        spatial_genes = {}
        for row in self.spatial_graph.query(spatial_genes_query):
            gene_uri = row.gene
            gene_id = str(row.geneID)
            spatial_genes[gene_id] = gene_uri
        
        # Get all genes from single-cell ontology
        sc_genes_query = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?geneID
        WHERE {
            ?gene rdf:type sc:Gene .
            ?gene sc:hasGeneID ?geneID .
        }
        """
        
        sc_genes = {}
        for row in self.single_cell_graph.query(sc_genes_query):
            gene_uri = row.gene
            gene_id = str(row.geneID)
            sc_genes[gene_id] = gene_uri
        
        # Create mappings for genes that exist in both ontologies
        mapped_count = 0
        for gene_id, sc_gene_uri in sc_genes.items():
            if gene_id in spatial_genes:
                spatial_gene_uri = spatial_genes[gene_id]
                
                # Add correspondence relationship
                self.integrated_graph.add((sc_gene_uri, BRIDGE.hasCorrespondingGene, spatial_gene_uri))
                self.integrated_graph.add((spatial_gene_uri, BRIDGE.hasCorrespondingGene, sc_gene_uri))
                
                # Add owl:sameAs relationship for semantic equivalence
                self.integrated_graph.add((sc_gene_uri, OWL.sameAs, spatial_gene_uri))
                
                mapped_count += 1
        
        print(f"Mapped {mapped_count} genes between ontologies")
    
    def map_cell_types(self):
        """
        Create explicit mappings between cell types in the spatial and single-cell ontologies.
        """
        print("Mapping cell types between ontologies...")
        
        # Define mappings between cell types
        cell_type_mappings = {
            "Protoplasmic": (URIRef(ST + "Protoplasmic"), URIRef(SC + "Protoplasmic")),
            "Fibrous": (URIRef(ST + "Fibrous"), URIRef(SC + "Fibrous")),
            "Reactive": (URIRef(ST + "Reactive"), URIRef(SC + "Reactive"))
        }
        
        # Add mappings to the integrated graph
        for cell_type, (spatial_type, sc_type) in cell_type_mappings.items():
            # Add correspondence relationship
            self.integrated_graph.add((sc_type, BRIDGE.hasCorrespondingCellType, spatial_type))
            self.integrated_graph.add((spatial_type, BRIDGE.hasCorrespondingCellType, sc_type))
            
            # Add owl:sameAs relationship for semantic equivalence
            self.integrated_graph.add((sc_type, OWL.sameAs, spatial_type))
        
        print(f"Mapped {len(cell_type_mappings)} cell types between ontologies")
    
    def infer_spatial_locations_for_cells(self):
        """
        Infer spatial locations for single cells based on cell type and gene expression.
        Focus only on astrocytes to simplify the integration.
        """
        print("Inferring spatial locations for single cells...")
        
        # Get all single cells with astrocyte cell types
        sc_cells_query = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?cell ?cellType
        WHERE {
            ?cell rdf:type sc:Cell ;
                  sc:hasCellType ?cellType .
            
            # Only include astrocyte cell types
            FILTER(?cellType IN ("Protoplasmic", "Fibrous", "Reactive"))
        }
        """
        
        # Get all spatial points with their types
        spatial_points_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y ?astrocyteType ?probability
        WHERE {
            ?point rdf:type st:SpatialPoint ;
                   st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?astrocyteType ;
                      st:typeProbability ?probability .
            
            # Only include astrocyte types
            FILTER(?astrocyteType IN ("Protoplasmic", "Fibrous", "Reactive"))
        }
        """
        
        # Execute queries
        sc_cells = list(self.integrated_graph.query(sc_cells_query))
        spatial_points = list(self.integrated_graph.query(spatial_points_query))
        
        # Group spatial points by astrocyte type
        spatial_points_by_type = {}
        for row in spatial_points:
            astrocyte_type = str(row.astrocyteType)
            if astrocyte_type not in spatial_points_by_type:
                spatial_points_by_type[astrocyte_type] = []
            
            spatial_points_by_type[astrocyte_type].append({
                'point': row.point,
                'x': float(row.x),
                'y': float(row.y),
                'probability': float(row.probability)
            })
        
        # Sort spatial points by probability (highest first)
        for cell_type in spatial_points_by_type:
            spatial_points_by_type[cell_type].sort(key=lambda p: p['probability'], reverse=True)
        
        # Map single cells to spatial points based on cell type
        mapped_count = 0
        for row in sc_cells:
            cell = row.cell
            cell_type = str(row.cellType)
            
            if cell_type in spatial_points_by_type and spatial_points_by_type[cell_type]:
                # Get the highest probability spatial point for this cell type
                # that hasn't been assigned yet (or the first one if all assigned)
                spatial_point = spatial_points_by_type[cell_type][0]
                
                # Add relationships between cell and spatial point
                self.integrated_graph.add((cell, BRIDGE.hasSpatialRepresentation, spatial_point['point']))
                self.integrated_graph.add((spatial_point['point'], BRIDGE.hasSingleCellData, cell))
                
                # Move this point to the end of the list to prioritize unmapped points
                spatial_points_by_type[cell_type].append(spatial_points_by_type[cell_type].pop(0))
                
                mapped_count += 1
        
        print(f"Inferred spatial locations for {mapped_count} single cells")
    
    def add_inference_rules(self):
        """
        Add inference rules to the integrated ontology to enable advanced reasoning.
        """
        print("Adding inference rules...")
        
        # Rule 1: If a gene is expressed in a cell and that cell has a spatial representation,
        # then the gene is expressed at that spatial location
        rule1 = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        INSERT {
            ?gene st:expressedAt ?point .
            ?point st:hasGeneExpression ?gene .
        }
        WHERE {
            ?expr sc:belongsToCell ?cell ;
                  sc:forGene ?gene .
            
            ?cell bridge:hasSpatialRepresentation ?point .
        }
        """
        
        # Rule 2: If a cell has a certain type and has a spatial representation,
        # then that spatial point has the same cell type
        rule2 = """
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        INSERT {
            ?point st:hasCellType ?cellType .
        }
        WHERE {
            ?cell sc:hasCellType ?cellType ;
                  bridge:hasSpatialRepresentation ?point .
        }
        """
        
        # Execute inference rules
        self.integrated_graph.update(rule1)
        self.integrated_graph.update(rule2)
        
        print(f"Added inference rules, graph now has {len(self.integrated_graph)} triples")
    
    def save_integrated_ontology(self):
        """
        Save the integrated ontology to a TURTLE file.
        """
        print(f"Saving integrated ontology to {self.output_ttl}...")
        
        # Add prefixes for readability
        self.integrated_graph.bind("st", ST)
        self.integrated_graph.bind("sc", SC)
        self.integrated_graph.bind("cell", CELL)
        self.integrated_graph.bind("gene", GENE)
        self.integrated_graph.bind("bridge", BRIDGE)
        self.integrated_graph.bind("owl", OWL)
        self.integrated_graph.bind("rdf", RDF)
        self.integrated_graph.bind("rdfs", RDFS)
        self.integrated_graph.bind("xsd", XSD)
        
        # Save to file
        self.integrated_graph.serialize(destination=self.output_ttl, format="turtle")
        
        print(f"Integrated ontology saved with {len(self.integrated_graph)} triples")
    
    def run_example_queries(self):
        """
        Run example queries that demonstrate the power of the integrated ontology.
        Focus only on astrocytes to simplify the analysis.
        """
        print("Running example queries on the integrated ontology...")
        
        # Query 1: Find genes that are expressed in astrocytes in both single-cell and spatial data
        query1 = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?geneID ?scExprLevel ?stExprLevel ?cellType
        WHERE {
            # Get gene from single-cell data
            ?scGene rdf:type sc:Gene ;
                    sc:hasGeneID ?geneID .
            
            # Get corresponding gene in spatial data
            ?scGene bridge:hasCorrespondingGene ?stGene .
            
            # Get expression in single-cell for astrocytes only
            ?scExpr sc:forGene ?scGene ;
                    sc:hasExpressionLevel ?scExprLevel ;
                    sc:belongsToCell ?cell .
            
            ?cell sc:hasCellType ?cellType .
            
            # Only include astrocyte cell types
            FILTER(?cellType IN ("Protoplasmic", "Fibrous", "Reactive"))
            
            # Get expression in spatial data
            ?stGene st:expressedAt ?point ;
                    st:hasExpressionLevel ?stExprLevel .
            
            # Get astrocyte type at this point
            ?point st:hasAstrocyteType ?typeNode .
            ?typeNode st:astrocyteType ?astrocyteType .
            
            # Only include points with astrocyte types
            FILTER(?astrocyteType IN ("Protoplasmic", "Fibrous", "Reactive"))
        }
        LIMIT 10
        """
        
        # Query 2: Find spatial regions where specific astrocyte types are predominantly found
        # and the genes that are highly expressed in those regions
        query2 = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?cellType ?x ?y ?geneID ?expressionLevel
        WHERE {
            # Get spatial point with astrocyte type
            ?point st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?cellType ;
                      st:typeProbability ?probability .
            
            # Only include astrocyte types
            FILTER(?cellType IN ("Protoplasmic", "Fibrous", "Reactive"))
            
            # Only consider high probability assignments
            FILTER(?probability > 0.8)
            
            # Get genes expressed at this point
            ?gene st:expressedAt ?point ;
                  st:hasExpressionLevel ?expressionLevel ;
                  st:hasGeneID ?geneID .
            
            # Only consider high expression
            FILTER(?expressionLevel > 1.0)
        }
        ORDER BY ?cellType ?x ?y
        LIMIT 20
        """
        
        # Query 3: Find genes that are differentially expressed between astrocyte types
        # in both spatial and single-cell data
        query3 = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX sc: <http://example.org/single-cell#>
        PREFIX bridge: <http://example.org/ontology-bridge#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?cellType1 ?cellType2 ?geneID 
               (AVG(?expr1) as ?avgExpr1) (AVG(?expr2) as ?avgExpr2)
               (ABS(AVG(?expr1) - AVG(?expr2)) as ?diffExpr)
        WHERE {
            # Get expression in first cell type
            ?cell1 sc:hasCellType ?cellType1 .
            ?expr1Node sc:belongsToCell ?cell1 ;
                       sc:forGene ?gene ;
                       sc:hasExpressionLevel ?expr1 .
            
            # Get expression in second cell type
            ?cell2 sc:hasCellType ?cellType2 .
            ?expr2Node sc:belongsToCell ?cell2 ;
                       sc:forGene ?gene ;
                       sc:hasExpressionLevel ?expr2 .
            
            # Get gene ID
            ?gene sc:hasGeneID ?geneID .
            
            # Ensure we're comparing different cell types
            FILTER(?cellType1 != ?cellType2)
            
            # Only consider astrocyte types
            FILTER(?cellType1 IN ("Protoplasmic", "Fibrous", "Reactive"))
            FILTER(?cellType2 IN ("Protoplasmic", "Fibrous", "Reactive"))
            
            # Ensure consistent ordering of cell types
            FILTER(STR(?cellType1) < STR(?cellType2))
        }
        GROUP BY ?cellType1 ?cellType2 ?geneID
        HAVING (ABS(AVG(?expr1) - AVG(?expr2)) > 0.5)
        ORDER BY DESC(?diffExpr)
        LIMIT 20
        """
        
        # Execute queries
        print("\nQuery 1: Genes expressed in astrocytes in both modalities")
        results1 = self.integrated_graph.query(query1)
        for row in results1:
            # Convert Literal objects to appropriate Python types before formatting
            gene_id = str(row.geneID)
            sc_expr = float(row.scExprLevel) if row.scExprLevel else 0.0
            st_expr = float(row.stExprLevel) if row.stExprLevel else 0.0
            cell_type = str(row.cellType) if hasattr(row, 'cellType') else "Unknown"
            print(f"Gene: {gene_id}, Cell Type: {cell_type}, SC Expr: {sc_expr:.2f}, ST Expr: {st_expr:.2f}")
        
        print("\nQuery 2: Spatial regions with cell types and highly expressed genes")
        results2 = self.integrated_graph.query(query2)
        for row in results2:
            # Convert Literal objects to appropriate Python types before formatting
            cell_type = str(row.cellType)
            x_coord = float(row.x) if row.x else 0.0
            y_coord = float(row.y) if row.y else 0.0
            gene_id = str(row.geneID)
            expr_level = float(row.expressionLevel) if row.expressionLevel else 0.0
            print(f"Cell Type: {cell_type}, Position: ({x_coord:.1f}, {y_coord:.1f}), Gene: {gene_id}, Expr: {expr_level:.2f}")
        
        print("\nQuery 3: Differentially expressed genes between astrocyte types")
        results3 = self.integrated_graph.query(query3)
        for row in results3:
            # Convert Literal objects to appropriate Python types before formatting
            cell_type1 = str(row.cellType1)
            cell_type2 = str(row.cellType2)
            gene_id = str(row.geneID)
            avg_expr1 = float(row.avgExpr1) if row.avgExpr1 else 0.0
            avg_expr2 = float(row.avgExpr2) if row.avgExpr2 else 0.0
            diff_expr = float(row.diffExpr) if row.diffExpr else 0.0
            print(f"Gene: {gene_id}, {cell_type1} vs {cell_type2}, Expr: {avg_expr1:.2f} vs {avg_expr2:.2f}, Diff: {diff_expr:.2f}")
    
    def run_pipeline(self):
        """
        Run the full integration pipeline.
        """
        print("Running ontology integration pipeline...")
        
        # Define bridge ontology
        self.define_bridge_ontology()
        
        # Map entities between ontologies
        self.map_genes_between_ontologies()
        self.map_cell_types()
        
        # Infer relationships
        self.infer_spatial_locations_for_cells()
        self.add_inference_rules()
        
        # Save integrated ontology
        self.save_integrated_ontology()
        
        # Run example queries
        self.run_example_queries()
        
        print("Integration pipeline completed")

def main():
    """
    Main function to run the ontology integration tool.
    """
    parser = argparse.ArgumentParser(description='Semantic integration of spatial and single-cell ontologies')
    parser.add_argument('--spatial', required=True, help='Path to spatial ontology TURTLE file')
    parser.add_argument('--single-cell', required=True, help='Path to single-cell ontology TURTLE file')
    parser.add_argument('--output', default='../output/integrated_ontology.ttl', help='Output TURTLE file')
    
    args = parser.parse_args()
    
    # Create integration tool
    integration_tool = OntologyIntegration(
        spatial_ttl=args.spatial,
        single_cell_ttl=args.single_cell,
        output_ttl=args.output
    )
    
    # Run pipeline
    integration_tool.run_pipeline()

if __name__ == "__main__":
    main() 