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
CELL = Namespace("http://example.org/cell-ontology#")
GENE = Namespace("http://example.org/gene-ontology#")
UBERON = Namespace("http://purl.obolibrary.org/obo/UBERON_")
GO = Namespace("http://purl.obolibrary.org/obo/GO_")

class SpatialOntologyEnhancer:
    """
    A class for enhancing the spatial ontology with additional information
    such as brain regions, morphology, connectivity, and functional annotations.
    """
    
    def __init__(self, input_ttl, output_ttl="../output/enhanced_spatial_ontology.ttl"):
        """
        Initialize the enhancer with input and output file paths.
        
        Parameters:
        -----------
        input_ttl : str
            Path to the input TURTLE file containing the spatial ontology
        output_ttl : str
            Path to the output TURTLE file where the enhanced ontology will be saved
        """
        self.input_ttl = input_ttl
        self.output_ttl = output_ttl
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(self.output_ttl), exist_ok=True)
        
        # Initialize graph
        self.graph = None
        
        # Load data
        self.load_data()
    
    def load_data(self):
        """
        Load the spatial ontology data from TURTLE file.
        """
        print(f"Loading spatial TURTLE file: {self.input_ttl}")
        self.graph = Graph()
        self.graph.parse(self.input_ttl, format="turtle")
        print(f"Spatial graph loaded with {len(self.graph)} triples")
    
    def add_brain_regions(self):
        """
        Add brain region information to the spatial ontology.
        """
        print("Adding brain region information...")
        
        # Define brain regions
        brain_regions = {
            "Cortex": {
                "description": "The cerebral cortex, the outermost layer of the brain",
                "subregions": ["PrefrontalCortex", "MotorCortex", "SensoryCortex", "VisualCortex"]
            },
            "Hippocampus": {
                "description": "A brain structure embedded deep in the temporal lobe",
                "subregions": ["CA1", "CA2", "CA3", "DentateGyrus"]
            },
            "Cerebellum": {
                "description": "A structure at the back of the brain involved in motor control",
                "subregions": ["CerebellarCortex", "DeepCerebellarNuclei"]
            },
            "Thalamus": {
                "description": "A midline symmetrical structure in the brain that relays sensory and motor signals",
                "subregions": ["MedialDorsalNucleus", "VentralPosteriorNucleus"]
            },
            "Hypothalamus": {
                "description": "A region of the brain that controls many body functions",
                "subregions": ["SupraopticNucleus", "ParaventricularNucleus"]
            }
        }
        
        # Add brain regions to the ontology
        for region_name, region_info in brain_regions.items():
            region_uri = URIRef(ST + region_name)
            self.graph.add((region_uri, RDF.type, ST.BrainRegion))
            self.graph.add((region_uri, RDFS.label, Literal(region_name)))
            self.graph.add((region_uri, RDFS.comment, Literal(region_info["description"])))
            
            # Link to UBERON ontology if applicable
            if region_name == "Cortex":
                self.graph.add((region_uri, OWL.sameAs, UBERON["0000956"]))  # cerebral cortex
            elif region_name == "Hippocampus":
                self.graph.add((region_uri, OWL.sameAs, UBERON["0002421"]))  # hippocampus
            elif region_name == "Cerebellum":
                self.graph.add((region_uri, OWL.sameAs, UBERON["0002037"]))  # cerebellum
            elif region_name == "Thalamus":
                self.graph.add((region_uri, OWL.sameAs, UBERON["0001897"]))  # thalamus
            elif region_name == "Hypothalamus":
                self.graph.add((region_uri, OWL.sameAs, UBERON["0001898"]))  # hypothalamus
            
            # Add subregions
            for subregion in region_info["subregions"]:
                subregion_uri = URIRef(ST + subregion)
                self.graph.add((subregion_uri, RDF.type, ST.BrainRegion))
                self.graph.add((subregion_uri, RDFS.label, Literal(subregion)))
                self.graph.add((subregion_uri, ST.isPartOf, region_uri))
        
        # Assign brain regions to spatial points
        # Get all spatial points
        spatial_points_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y
        WHERE {
            ?point rdf:type st:SpatialPoint ;
                   st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y .
        }
        """
        
        results = list(self.graph.query(spatial_points_query))
        
        # Assign brain regions based on spatial coordinates
        # This is a simplified approach; in a real scenario, this would be based on actual anatomical data
        for row in results:
            point_uri = row.point
            x = float(row.x)
            y = float(row.y)
            
            # Assign brain region based on x, y coordinates (simplified example)
            if x < 25 and y < 25:
                region_uri = URIRef(ST + "PrefrontalCortex")
                self.graph.add((point_uri, ST.inBrainRegion, region_uri))
            elif x < 25 and y >= 25:
                region_uri = URIRef(ST + "MotorCortex")
                self.graph.add((point_uri, ST.inBrainRegion, region_uri))
            elif x >= 25 and x < 50 and y < 50:
                region_uri = URIRef(ST + "Hippocampus")
                self.graph.add((point_uri, ST.inBrainRegion, region_uri))
            elif x >= 50 and x < 75 and y < 75:
                region_uri = URIRef(ST + "Cerebellum")
                self.graph.add((point_uri, ST.inBrainRegion, region_uri))
            else:
                region_uri = URIRef(ST + "Thalamus")
                self.graph.add((point_uri, ST.inBrainRegion, region_uri))
        
        print("Brain region information added")
    
    def add_morphology_information(self):
        """
        Add morphology information to astrocyte types.
        """
        print("Adding morphology information...")
        
        # Define morphology types
        morphology_types = {
            "Protoplasmic": {
                "description": "Highly branched with short, thick processes",
                "features": {
                    "process_length": "short",
                    "process_thickness": "thick",
                    "branching_density": "high",
                    "domain_radius": "40-60 μm"
                }
            },
            "Fibrous": {
                "description": "Long, thin, unbranched processes",
                "features": {
                    "process_length": "long",
                    "process_thickness": "thin",
                    "branching_density": "low",
                    "domain_radius": "100-200 μm"
                }
            },
            "Reactive": {
                "description": "Hypertrophic with increased GFAP expression",
                "features": {
                    "process_length": "variable",
                    "process_thickness": "thick",
                    "branching_density": "medium",
                    "domain_radius": "variable"
                }
            }
        }
        
        # Add morphology types to the ontology
        for morph_name, morph_info in morphology_types.items():
            morph_uri = URIRef(ST + f"Morphology_{morph_name}")
            self.graph.add((morph_uri, RDF.type, ST.Morphology))
            self.graph.add((morph_uri, RDFS.label, Literal(f"{morph_name} Morphology")))
            self.graph.add((morph_uri, RDFS.comment, Literal(morph_info["description"])))
            
            # Add features
            for feature_name, feature_value in morph_info["features"].items():
                feature_uri = URIRef(ST + feature_name)
                self.graph.add((feature_uri, RDF.type, ST.MorphologicalFeature))
                self.graph.add((feature_uri, RDFS.label, Literal(feature_name)))
                self.graph.add((morph_uri, ST.hasFeature, feature_uri))
                self.graph.add((morph_uri, feature_uri, Literal(feature_value)))
        
        # Assign morphology to astrocyte types
        astrocyte_types_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        
        SELECT ?typeNode ?astrocyteType
        WHERE {
            ?point st:hasAstrocyteType ?typeNode .
            ?typeNode st:astrocyteType ?astrocyteType .
        }
        """
        
        results = list(self.graph.query(astrocyte_types_query))
        
        for row in results:
            type_node = row.typeNode
            astrocyte_type = str(row.astrocyteType)
            
            if astrocyte_type in morphology_types:
                morph_uri = URIRef(ST + f"Morphology_{astrocyte_type}")
                self.graph.add((type_node, ST.hasMorphology, morph_uri))
        
        print("Morphology information added")
    
    def add_connectivity_information(self):
        """
        Add connectivity information between spatial points.
        """
        print("Adding connectivity information...")
        
        # Get all spatial points
        spatial_points_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point ?x ?y ?astrocyteType
        WHERE {
            ?point rdf:type st:SpatialPoint ;
                   st:hasXCoordinate ?x ;
                   st:hasYCoordinate ?y ;
                   st:hasAstrocyteType ?typeNode .
            
            ?typeNode st:astrocyteType ?astrocyteType .
        }
        """
        
        results = list(self.graph.query(spatial_points_query))
        
        # Convert to list of dictionaries for easier processing
        points = []
        for row in results:
            points.append({
                'uri': row.point,
                'x': float(row.x),
                'y': float(row.y),
                'type': str(row.astrocyteType)
            })
        
        # Add connectivity based on proximity and astrocyte type
        # This is a simplified approach; in a real scenario, this would be based on actual connectivity data
        for i, point1 in enumerate(points):
            # Find nearest neighbors
            neighbors = []
            for j, point2 in enumerate(points):
                if i != j:
                    # Calculate Euclidean distance
                    distance = np.sqrt((point1['x'] - point2['x'])**2 + (point1['y'] - point2['y'])**2)
                    neighbors.append((j, distance))
            
            # Sort by distance
            neighbors.sort(key=lambda x: x[1])
            
            # Connect to 3 nearest neighbors
            for j, distance in neighbors[:3]:
                point2 = points[j]
                
                # Add connectivity relationship
                self.graph.add((point1['uri'], ST.connectsTo, point2['uri']))
                
                # Add connection strength based on distance (inverse relationship)
                strength = 1.0 / (1.0 + distance)
                conn_node = BNode()
                self.graph.add((point1['uri'], ST.hasConnection, conn_node))
                self.graph.add((conn_node, ST.toPoint, point2['uri']))
                self.graph.add((conn_node, ST.connectionStrength, Literal(strength)))
                
                # Add connection type based on astrocyte types
                if point1['type'] == point2['type']:
                    self.graph.add((conn_node, ST.connectionType, Literal("homotypic")))
                else:
                    self.graph.add((conn_node, ST.connectionType, Literal("heterotypic")))
        
        print("Connectivity information added")
    
    def add_functional_annotations(self):
        """
        Add functional annotations to astrocyte types.
        """
        print("Adding functional annotations...")
        
        # Define functional annotations for each astrocyte type
        functional_annotations = {
            "Protoplasmic": [
                {
                    "function": "GlutamateUptake",
                    "description": "Uptake of glutamate from synaptic cleft",
                    "go_term": "0015813",  # L-glutamate transport
                    "genes": ["SLC1A2", "SLC1A3"]
                },
                {
                    "function": "PotassiumBuffering",
                    "description": "Regulation of extracellular potassium levels",
                    "go_term": "0006813",  # potassium ion transport
                    "genes": ["KCNJ10", "AQP4"]
                },
                {
                    "function": "SynapticModulation",
                    "description": "Modulation of synaptic transmission",
                    "go_term": "0050804",  # modulation of chemical synaptic transmission
                    "genes": ["GLUL", "GJA1"]
                }
            ],
            "Fibrous": [
                {
                    "function": "StructuralSupport",
                    "description": "Structural support for axons",
                    "go_term": "0045103",  # intermediate filament-based process
                    "genes": ["GFAP", "VIM"]
                },
                {
                    "function": "MyelinMaintenance",
                    "description": "Support for myelin maintenance",
                    "go_term": "0032288",  # myelin assembly
                    "genes": ["CRYAB", "HOPX"]
                },
                {
                    "function": "BloodBrainBarrierMaintenance",
                    "description": "Maintenance of blood-brain barrier integrity",
                    "go_term": "0043085",  # positive regulation of catalytic activity
                    "genes": ["CD44", "ALDH1L1"]
                }
            ],
            "Reactive": [
                {
                    "function": "InflammatoryResponse",
                    "description": "Response to inflammation and injury",
                    "go_term": "0006954",  # inflammatory response
                    "genes": ["SERPINA3", "C3"]
                },
                {
                    "function": "CytokineProduction",
                    "description": "Production of cytokines",
                    "go_term": "0001816",  # cytokine production
                    "genes": ["CXCL10", "LCN2"]
                },
                {
                    "function": "GlialScarFormation",
                    "description": "Formation of glial scar",
                    "go_term": "0007050",  # cell cycle arrest
                    "genes": ["STAT3", "GFAP"]
                }
            ]
        }
        
        # Add functional annotations to the ontology
        for astro_type, functions in functional_annotations.items():
            for func in functions:
                # Create function node
                func_uri = URIRef(ST + func["function"])
                self.graph.add((func_uri, RDF.type, ST.Function))
                self.graph.add((func_uri, RDFS.label, Literal(func["function"])))
                self.graph.add((func_uri, RDFS.comment, Literal(func["description"])))
                
                # Link to GO term if available
                if "go_term" in func:
                    self.graph.add((func_uri, OWL.sameAs, GO[func["go_term"]]))
                
                # Add genes associated with this function
                for gene in func["genes"]:
                    gene_uri = URIRef(ST + gene)
                    self.graph.add((func_uri, ST.associatedWithGene, gene_uri))
        
        # Assign functions to astrocyte types
        astrocyte_types_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        
        SELECT ?typeNode ?astrocyteType
        WHERE {
            ?point st:hasAstrocyteType ?typeNode .
            ?typeNode st:astrocyteType ?astrocyteType .
        }
        """
        
        results = list(self.graph.query(astrocyte_types_query))
        
        for row in results:
            type_node = row.typeNode
            astrocyte_type = str(row.astrocyteType)
            
            if astrocyte_type in functional_annotations:
                for func in functional_annotations[astrocyte_type]:
                    func_uri = URIRef(ST + func["function"])
                    self.graph.add((type_node, ST.hasFunctionalAnnotation, func_uri))
        
        print("Functional annotations added")
    
    def add_temporal_information(self):
        """
        Add temporal information to the spatial ontology.
        """
        print("Adding temporal information...")
        
        # Define timepoints
        timepoints = ["Baseline", "Day1", "Day3", "Day7", "Day14", "Day28"]
        
        # Add timepoints to the ontology
        for tp in timepoints:
            tp_uri = URIRef(ST + tp)
            self.graph.add((tp_uri, RDF.type, ST.Timepoint))
            self.graph.add((tp_uri, RDFS.label, Literal(tp)))
        
        # Assign timepoints to spatial points
        # In this example, we'll assign all points to Baseline
        # In a real scenario, this would be based on actual temporal data
        spatial_points_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?point
        WHERE {
            ?point rdf:type st:SpatialPoint .
        }
        """
        
        results = list(self.graph.query(spatial_points_query))
        
        for row in results:
            point_uri = row.point
            tp_uri = URIRef(ST + "Baseline")
            self.graph.add((point_uri, ST.hasTimepoint, tp_uri))
        
        print("Temporal information added")
    
    def add_human_protein_atlas_links(self):
        """
        Add links to the Human Protein Atlas for genes in the ontology.
        """
        print("Adding Human Protein Atlas links...")
        
        # Define HPA namespace
        HPA = Namespace("https://www.proteinatlas.org/")
        
        # Get all genes
        genes_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene ?geneID
        WHERE {
            ?gene rdf:type st:Gene ;
                  st:hasGeneID ?geneID .
        }
        """
        
        results = list(self.graph.query(genes_query))
        
        # Add HPA links for common astrocyte marker genes
        hpa_links = {
            "GFAP": "ENSG00000131095-GFAP",
            "S100B": "ENSG00000160307-S100B",
            "AQP4": "ENSG00000171885-AQP4",
            "GJA1": "ENSG00000152661-GJA1",
            "SLC1A2": "ENSG00000110436-SLC1A2",
            "SLC1A3": "ENSG00000079215-SLC1A3",
            "VIM": "ENSG00000026025-VIM",
            "ALDH1L1": "ENSG00000144908-ALDH1L1",
            "CD44": "ENSG00000026508-CD44",
            "CRYAB": "ENSG00000109846-CRYAB",
            "HOPX": "ENSG00000171476-HOPX",
            "SERPINA3": "ENSG00000196136-SERPINA3",
            "C3": "ENSG00000125730-C3",
            "CXCL10": "ENSG00000169245-CXCL10",
            "STAT3": "ENSG00000168610-STAT3",
            "LCN2": "ENSG00000148346-LCN2"
        }
        
        for row in results:
            gene_uri = row.gene
            gene_id = str(row.geneID)
            
            if gene_id in hpa_links:
                hpa_uri = URIRef(HPA + hpa_links[gene_id])
                self.graph.add((gene_uri, ST.hasHPAEntry, hpa_uri))
                self.graph.add((hpa_uri, RDF.type, ST.HPAEntry))
                self.graph.add((hpa_uri, RDFS.label, Literal(f"HPA Entry for {gene_id}")))
        
        # Add HPA Brain Atlas link
        hpa_brain_atlas = URIRef(HPA + "humanproteome/brain")
        self.graph.add((hpa_brain_atlas, RDF.type, ST.HPABrainAtlas))
        self.graph.add((hpa_brain_atlas, RDFS.label, Literal("Human Protein Atlas Brain Atlas")))
        
        # Link astrocyte types to HPA Brain Atlas
        astrocyte_types_query = """
        PREFIX st: <http://example.org/spatial-transcriptomics#>
        
        SELECT DISTINCT ?astrocyteType
        WHERE {
            ?typeNode st:astrocyteType ?astrocyteType .
        }
        """
        
        results = list(self.graph.query(astrocyte_types_query))
        
        for row in results:
            astrocyte_type = str(row.astrocyteType)
            self.graph.add((URIRef(ST + astrocyte_type), ST.referencedInHPA, hpa_brain_atlas))
        
        print("Human Protein Atlas links added")
    
    def enhance_ontology(self):
        """
        Enhance the spatial ontology with additional information.
        """
        print("Enhancing spatial ontology...")
        
        # Add brain regions
        self.add_brain_regions()
        
        # Add morphology information
        self.add_morphology_information()
        
        # Add connectivity information
        self.add_connectivity_information()
        
        # Add functional annotations
        self.add_functional_annotations()
        
        # Add temporal information
        self.add_temporal_information()
        
        # Add Human Protein Atlas links
        self.add_human_protein_atlas_links()
        
        # Save enhanced ontology
        self.save_enhanced_ontology()
        
        print("Spatial ontology enhancement completed")
    
    def save_enhanced_ontology(self):
        """
        Save the enhanced ontology to a TURTLE file.
        """
        print(f"Saving enhanced ontology to {self.output_ttl}...")
        
        # Add prefixes for readability
        self.graph.bind("st", ST)
        self.graph.bind("cell", CELL)
        self.graph.bind("gene", GENE)
        self.graph.bind("uberon", UBERON)
        self.graph.bind("go", GO)
        self.graph.bind("owl", OWL)
        self.graph.bind("rdf", RDF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind("xsd", XSD)
        
        # Save to file
        self.graph.serialize(destination=self.output_ttl, format="turtle")
        
        print(f"Enhanced ontology saved with {len(self.graph)} triples")

def main():
    """
    Main function to run the spatial ontology enhancer.
    """
    parser = argparse.ArgumentParser(description='Enhance spatial ontology with additional information')
    parser.add_argument('--input', required=True, help='Path to input spatial ontology TURTLE file')
    parser.add_argument('--output', default='../output/enhanced_spatial_ontology.ttl', help='Output TURTLE file')
    
    args = parser.parse_args()
    
    # Create enhancer
    enhancer = SpatialOntologyEnhancer(
        input_ttl=args.input,
        output_ttl=args.output
    )
    
    # Enhance ontology
    enhancer.enhance_ontology()

if __name__ == "__main__":
    main() 