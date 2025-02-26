from abc import ABC, abstractmethod
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, XSD
import pandas as pd

class OntologyBuilder(ABC):
    """
    Abstract base class for creating an ontology using the Template Design Pattern.
    """
    def __init__(self):
        # Initialize the RDF graph and namespace
        self.graph = Graph()
        self.EX = Namespace("http://example.org/ontology#")
        self.graph.bind("ex", self.EX)

    @abstractmethod
    def create_classes(self):
        """ Define ontology classes (to be implemented by subclasses). """
        pass

    @abstractmethod
    def add_data(self, df):
        """ Add data to the ontology (to be implemented by subclasses). """
        pass

    def save_to_file(self, output_path: str):
        """ Serialize the ontology to a TURTLE file. """
        self.graph.serialize(output_path, format="turtle")
        print(f"Ontology saved to {output_path}")

class GeneExpressionOntologyBuilder(OntologyBuilder):
    """
    Ontology Builder for Gene Expression Data Averaged by Cell Type.
    """
    def create_classes(self):
        """ Define ontology classes and properties for cell types and genes. """
        OWL = Namespace("http://www.w3.org/2002/07/owl#")
        
        # Define ontology classes
        self.graph.add((self.EX.CellType, RDF.type, self.EX.Class))
        self.graph.add((self.EX.Gene, RDF.type, self.EX.Class))
        self.graph.add((self.EX.GeneExpression, RDF.type, self.EX.Class))

        # Define ontology properties
        self.graph.add((self.EX.hasGene, RDF.type, self.EX.Property))
        self.graph.add((self.EX.belongsToCellType, RDF.type, self.EX.Property))
        self.graph.add((self.EX.hasAverageExpression, RDF.type, self.EX.Property))
        
        # Define inverse properties
        self.graph.add((self.EX.isExpressedIn, RDF.type, OWL.InverseFunctionalProperty))
        self.graph.add((self.EX.containsExpression, RDF.type, OWL.InverseFunctionalProperty))

        # Declare inverse relationships
        self.graph.add((self.EX.hasGene, OWL.inverseOf, self.EX.isExpressedIn))
        self.graph.add((self.EX.belongsToCellType, OWL.inverseOf, self.EX.containsExpression))

    def add_data(self, df):
        """
        Add average gene expression data to the ontology, filtering only Astrocytes.

        Args:
            df (pd.DataFrame): DataFrame where rows are cell types, columns are genes, 
                               and values are average expression levels.
        """
        # Filter only Astrocyte cell types (case-insensitive)
        astrocyte_df = df[df.index.str.contains("astrocyte", case=False, na=False)]
        
        if astrocyte_df.empty:
            print("Warning: No Astrocytes found in dataset!")
            return

        for celltype, row in astrocyte_df.iterrows():
            # Create a URI for the cell type
            celltype_uri = URIRef(self.EX[celltype.replace(" ", "_")])
            self.graph.add((celltype_uri, RDF.type, self.EX.CellType))

            for gene, avg_expression in row.items():
                if avg_expression > 0:  # Only include expressed genes
                    gene_uri = URIRef(self.EX[gene])
                    self.graph.add((gene_uri, RDF.type, self.EX.Gene))

                    # Create a unique instance for the gene expression in this cell type
                    gene_expression_uri = URIRef(self.EX[f"{gene}_in_{celltype.replace(' ', '_')}"])
                    self.graph.add((gene_expression_uri, RDF.type, self.EX.GeneExpression))

                    # Link gene expression instance to cell type and gene
                    self.graph.add((gene_expression_uri, self.EX.hasGene, gene_uri))
                    self.graph.add((gene_expression_uri, self.EX.belongsToCellType, celltype_uri))
                    self.graph.add((gene_expression_uri, self.EX.hasAverageExpression, Literal(avg_expression, datatype=XSD.float)))
