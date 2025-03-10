# Data Conversion: Spatial Data to RDF/Turtle

This document provides detailed information about the data conversion system used in the Astrocytes Ontology Research Project. This module transforms spatial transcriptomics data into RDF/Turtle format for semantic analysis.

## Overview

The data conversion module is a critical first step in our pipeline, transforming raw spatial transcriptomics data into a semantic RDF representation. This conversion:

1. Preserves spatial coordinates of each gene expression event
2. Establishes semantic relationships between genes, expressions, and spatial points
3. Creates a standardized format for rule application
4. Integrates with existing biological ontologies

Converting to RDF/Turtle format enables powerful semantic queries via SPARQL, allowing complex pattern matching that would be difficult to express using conventional programming methods.

## Key Features

Our data conversion system provides several important features:

1. **Semantic Representation**: Maps raw data to meaningful RDF entities and relationships
2. **Ontology Integration**: Links data to established biological ontologies
3. **Spatial Preservation**: Maintains accurate spatial coordinates in the semantic model
4. **Expression Thresholding**: Filters expression data based on confidence scores
5. **Data Validation**: Verifies data integrity during conversion
6. **Efficient Storage**: Optimizes RDF representation for query performance
7. **Provenance Tracking**: Records metadata about the conversion process

## RDF Data Model

The data conversion system maps spatial transcriptomics data to this RDF model:

```turtle
# Namespaces
@prefix astro: <http://example.org/astrocytes/> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

# Gene definitions
astro:gene_GPM6A a astro:Gene ;
    rdfs:label "GPM6A" ;
    astro:ensemblID "ENSG00000150625" .

# Spatial data points
astro:point_123 a astro:SpatialDataPoint ;
    astro:hasXCoordinate "1234.56"^^xsd:decimal ;
    astro:hasYCoordinate "789.12"^^xsd:decimal ;
    astro:expressesGene astro:gene_GPM6A ;
    astro:expressionIntensity "2.45"^^xsd:decimal .
```

Key entity types include:
- **Gene**: Represents a gene that can be expressed
- **SpatialDataPoint**: Represents a point in space where a gene is expressed
- **SpatialCell**: (Created later by rule application) Represents an identified cell

## Implementation

The conversion process is implemented in `scripts/conversion/convert_spatial_to_turtle.py`:

```python
class SpatialDataConverter:
    def __init__(self, input_file, output_file, min_expression=0.1):
        """
        Initialize the converter.
        
        Args:
            input_file: Path to input spatial data file (CSV, H5AD, or other formats)
            output_file: Path where Turtle output will be written
            min_expression: Minimum expression threshold to include a data point
        """
        self.input_file = input_file
        self.output_file = output_file
        self.min_expression = min_expression
        self.gene_counter = 0
        self.point_counter = 0
        self.graph = rdflib.Graph()
        
        # Add namespaces
        self.astro = rdflib.Namespace("http://example.org/astrocytes/")
        self.graph.bind("astro", self.astro)
        self.graph.bind("rdfs", rdflib.RDFS)
        self.graph.bind("xsd", rdflib.XSD)
        
    def load_data(self):
        """Load the spatial data from the input file."""
        if self.input_file.endswith('.h5ad'):
            return self._load_h5ad()
        elif self.input_file.endswith('.csv'):
            return self._load_csv()
        else:
            raise ValueError(f"Unsupported file format: {self.input_file}")
            
    def _load_h5ad(self):
        """Load data from an H5AD file."""
        import scanpy as sc
        data = sc.read_h5ad(self.input_file)
        return self._extract_spatial_data(data)
    
    def _load_csv(self):
        """Load data from a CSV file."""
        import pandas as pd
        data = pd.read_csv(self.input_file)
        return data
        
    def _extract_spatial_data(self, data):
        """Extract spatial coordinates and gene expressions from the data."""
        # Implementation depends on input format
        pass
        
    def create_gene_entities(self, gene_names):
        """Create RDF entities for each gene."""
        gene_entities = {}
        
        for gene_name in gene_names:
            gene_uri = self.astro[f"gene_{gene_name}"]
            gene_entities[gene_name] = gene_uri
            
            self.graph.add((gene_uri, rdflib.RDF.type, self.astro.Gene))
            self.graph.add((gene_uri, rdflib.RDFS.label, rdflib.Literal(gene_name)))
            
        return gene_entities
        
    def create_spatial_points(self, spatial_data, gene_entities):
        """Create RDF entities for spatial data points."""
        for _, row in spatial_data.iterrows():
            # Skip low expression values
            if row['expression'] < self.min_expression:
                continue
                
            # Create point entity
            point_id = f"point_{self.point_counter}"
            self.point_counter += 1
            point_uri = self.astro[point_id]
            
            # Add point properties
            self.graph.add((point_uri, rdflib.RDF.type, self.astro.SpatialDataPoint))
            self.graph.add((point_uri, self.astro.hasXCoordinate, 
                           rdflib.Literal(row['x'], datatype=rdflib.XSD.decimal)))
            self.graph.add((point_uri, self.astro.hasYCoordinate, 
                           rdflib.Literal(row['y'], datatype=rdflib.XSD.decimal)))
            
            # Link to gene and add expression intensity
            gene_uri = gene_entities[row['gene']]
            self.graph.add((point_uri, self.astro.expressesGene, gene_uri))
            self.graph.add((point_uri, self.astro.expressionIntensity, 
                           rdflib.Literal(row['expression'], datatype=rdflib.XSD.decimal)))
    
    def convert(self):
        """Convert spatial data to RDF/Turtle format."""
        # Load data
        spatial_data = self.load_data()
        
        # Get unique gene names
        gene_names = spatial_data['gene'].unique()
        
        # Create gene entities
        gene_entities = self.create_gene_entities(gene_names)
        
        # Create spatial point entities
        self.create_spatial_points(spatial_data, gene_entities)
        
        # Add metadata
        self._add_metadata()
        
        # Write to output file
        self.graph.serialize(destination=self.output_file, format="turtle")
        
        return {
            'genes': len(gene_names),
            'points': self.point_counter
        }
        
    def _add_metadata(self):
        """Add metadata about the conversion process."""
        import datetime
        
        metadata_uri = self.astro.conversionMetadata
        self.graph.add((metadata_uri, rdflib.RDF.type, self.astro.ConversionMetadata))
        self.graph.add((metadata_uri, self.astro.conversionDate, 
                       rdflib.Literal(datetime.datetime.now().isoformat())))
        self.graph.add((metadata_uri, self.astro.inputFile, 
                       rdflib.Literal(os.path.basename(self.input_file))))
        self.graph.add((metadata_uri, self.astro.expressionThreshold, 
                       rdflib.Literal(self.min_expression, datatype=rdflib.XSD.decimal)))
```

## Input File Formats

The data conversion module supports multiple input formats:

### H5AD Format (AnnData)

The H5AD format is common for single-cell and spatial transcriptomics data:

```python
def _load_h5ad(self):
    """Load data from an H5AD file."""
    import scanpy as sc
    data = sc.read_h5ad(self.input_file)
    
    # Extract spatial coordinates
    if 'spatial' in data.obsm:
        coordinates = data.obsm['spatial']
    else:
        raise ValueError("H5AD file does not contain spatial coordinates")
        
    # Extract gene expression
    genes = data.var_names.tolist()
    spatial_data = []
    
    for i in range(data.n_obs):
        for j, gene in enumerate(genes):
            expression = data.X[i, j]
            if expression > self.min_expression:
                spatial_data.append({
                    'x': coordinates[i, 0],
                    'y': coordinates[i, 1],
                    'gene': gene,
                    'expression': expression
                })
                
    return pd.DataFrame(spatial_data)
```

### CSV Format

The CSV format should contain columns for coordinates and gene expression:

```python
def _load_csv(self):
    """Load data from a CSV file."""
    import pandas as pd
    
    # Validate CSV format
    data = pd.read_csv(self.input_file)
    required_columns = ['x', 'y', 'gene', 'expression']
    for col in required_columns:
        if col not in data.columns:
            raise ValueError(f"CSV file missing required column: {col}")
            
    return data
```

## Command-Line Interface

The module provides a comprehensive command-line interface:

```python
def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Convert spatial data to RDF/Turtle format.')
    
    parser.add_argument('--input', required=True,
                      help='Path to input spatial data file')
    parser.add_argument('--output', required=True,
                      help='Path for output Turtle file')
    parser.add_argument('--min-expression', type=float, default=0.1,
                      help='Minimum expression threshold (default: 0.1)')
    parser.add_argument('--gene-mapping', default=None,
                      help='Path to gene ID mapping file (optional)')
    parser.add_argument('--add-ensembl-ids', action='store_true',
                      help='Add Ensembl IDs to gene entities')
    parser.add_argument('--add-gene-names', action='store_true',
                      help='Add gene names in addition to symbols')
    parser.add_argument('--spatial-units', default='pixels',
                      help='Units of spatial coordinates (default: pixels)')
    
    return parser.parse_args()
```

This interface provides control over:
- Input and output file paths
- Expression thresholds for filtering
- Gene ID mapping options
- Additional gene metadata
- Spatial coordinate units

## Example Usage

Basic usage for data conversion:

```bash
python scripts/conversion/convert_spatial_to_turtle.py \
    --input data/raw/spatial_data.h5ad \
    --output data/processed/spatial_data.ttl \
    --min-expression 0.1
```

Advanced usage with gene mapping:

```bash
python scripts/conversion/convert_spatial_to_turtle.py \
    --input data/raw/spatial_data.h5ad \
    --output data/processed/spatial_data.ttl \
    --min-expression 0.1 \
    --gene-mapping data/raw/gene_id_mapping.csv \
    --add-ensembl-ids \
    --spatial-units micrometers
```

## Gene ID Mapping

To enhance interoperability, the system can add external gene identifiers:

```python
def add_gene_identifiers(self, gene_entities, mapping_file):
    """Add external identifiers (e.g., Ensembl IDs) to gene entities."""
    import pandas as pd
    
    # Load gene mapping
    mapping = pd.read_csv(mapping_file)
    mapping_dict = dict(zip(mapping['symbol'], mapping['ensembl_id']))
    
    # Add identifiers to gene entities
    for gene_name, gene_uri in gene_entities.items():
        if gene_name in mapping_dict:
            ensembl_id = mapping_dict[gene_name]
            self.graph.add((gene_uri, self.astro.ensemblID, rdflib.Literal(ensembl_id)))
```

## Performance Optimization

For large datasets, the system includes performance optimizations:

```python
def optimize_for_large_datasets(self):
    """Apply performance optimizations for large datasets."""
    # Use bulk loading instead of individual triples
    from rdflib.plugins.parsers.ntriples import NTriplesParser
    
    # Process in chunks
    chunk_size = 100000
    for i in range(0, len(self.spatial_data), chunk_size):
        chunk = self.spatial_data[i:i+chunk_size]
        self.create_spatial_points(chunk, self.gene_entities)
        
    # Use optimized serialization
    self.graph.serialize(destination=self.output_file, format="turtle", encoding="utf-8", 
                        compression=None)
```

## Output Format Details

The output Turtle file follows this structure:

```turtle
@prefix astro: <http://example.org/astrocytes/> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

# Metadata
astro:conversionMetadata a astro:ConversionMetadata ;
    astro:conversionDate "2023-10-15T14:30:22.451234" ;
    astro:inputFile "spatial_data.h5ad" ;
    astro:expressionThreshold "0.1"^^xsd:decimal .

# Genes
astro:gene_GPM6A a astro:Gene ;
    rdfs:label "GPM6A" ;
    astro:ensemblID "ENSG00000150625" .

astro:gene_NRG3 a astro:Gene ;
    rdfs:label "NRG3" ;
    astro:ensemblID "ENSG00000185737" .

# Spatial data points
astro:point_0 a astro:SpatialDataPoint ;
    astro:hasXCoordinate "1234.56"^^xsd:decimal ;
    astro:hasYCoordinate "789.12"^^xsd:decimal ;
    astro:expressesGene astro:gene_GPM6A ;
    astro:expressionIntensity "2.45"^^xsd:decimal .

astro:point_1 a astro:SpatialDataPoint ;
    astro:hasXCoordinate "1240.34"^^xsd:decimal ;
    astro:hasYCoordinate "795.67"^^xsd:decimal ;
    astro:expressesGene astro:gene_NRG3 ;
    astro:expressionIntensity "1.98"^^xsd:decimal .
```

## Integration with Other Tools

The RDF/Turtle format enables integration with external tools:

1. **SPARQL Endpoints**: The data can be loaded into a SPARQL endpoint
2. **RDF Visualization Tools**: Tools like Graphviz can visualize the RDF graph
3. **Ontology Editors**: Protégé can be used to explore the data structure
4. **Other RDF Formats**: The data can be converted to other RDF formats (JSON-LD, RDF/XML)
5. **Triple Stores**: For large datasets, the data can be loaded into a triple store

## Validation

The system includes validation to ensure data integrity:

```python
def validate_output(self):
    """Validate the generated RDF/Turtle file."""
    # Check if file was created
    if not os.path.exists(self.output_file):
        return False, "Output file was not created"
        
    # Try loading the file back to verify it's valid Turtle
    try:
        test_graph = rdflib.Graph()
        test_graph.parse(self.output_file, format="turtle")
        
        # Check if we have the expected entities
        genes_count = len(list(test_graph.subjects(rdflib.RDF.type, self.astro.Gene)))
        points_count = len(list(test_graph.subjects(rdflib.RDF.type, self.astro.SpatialDataPoint)))
        
        if genes_count == 0:
            return False, "No gene entities found in output"
        
        if points_count == 0:
            return False, "No spatial data points found in output"
            
        return True, f"Valid RDF with {genes_count} genes and {points_count} data points"
        
    except Exception as e:
        return False, f"Invalid Turtle format: {str(e)}"
```

## Troubleshooting

### Common Issues

| Issue | Possible Solution |
|---|---|
| Memory errors with large files | Use the `optimize_for_large_datasets` method |
| Missing gene IDs | Provide a gene mapping file with `--gene-mapping` |
| Invalid coordinates | Check the input format and coordinate system |
| Slow conversion | Increase `--min-expression` to filter more points |
| Empty output | Verify input file format and content |

### Logging and Diagnostics

The module includes detailed logging:

```python
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Example log messages
logger.info(f"Processing input file: {args.input}")
logger.info(f"Expression threshold: {args.min_expression}")
logger.info(f"Created {stats['genes']} gene entities and {stats['points']} data points")
```

## Future Enhancements

Planned improvements to the data conversion system include:

1. **3D Data Support**: Extension to three-dimensional spatial data
2. **Batch Processing**: Improved handling of very large datasets
3. **Additional Metadata**: Support for experimental and sample metadata
4. **Alternative Formats**: Support for additional input and output formats
5. **Performance Improvements**: Optimizations for very large spatial datasets
6. **Ontology Integration**: Enhanced integration with external biological ontologies

## References

1. RDF and Turtle format:
   - W3C. (2014). "RDF 1.1 Turtle: Terse RDF Triple Language." W3C Recommendation.

2. Spatial transcriptomics data formats:
   - Deng, Y., et al. (2021). "Spatial transcriptomics for the construction of a high-resolution cell atlas." Journal of Genetics and Genomics, 48(10), 878-886.

3. Gene identifiers:
   - Braschi, B., et al. (2019). "Genenames.org: the HGNC and VGNC resources in 2019." Nucleic Acids Research, 47(D1), D786-D792. 