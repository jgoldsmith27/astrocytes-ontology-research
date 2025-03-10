# Rule Application and Cell Identification

This document provides detailed information about the rule application and cell identification system used in the Astrocytes Ontology Research Project. This module applies SPARQL rules to spatial transcriptomics data to identify cells based on gene co-expression patterns.

## Overview

The rule application module is the core component that connects gene co-expression patterns to spatial transcriptomics data. It:

1. Loads spatial data in RDF format
2. Applies SPARQL rules to identify potential cells
3. Integrates with meta-rules for higher-order pattern matching
4. Resolves conflicts between overlapping cell identifications
5. Generates comprehensive output files and visualizations

This module represents the culmination of the pipeline, where co-expression patterns derived from single-cell data are used to identify whole cells in spatial data.

## Key Features

Our rule application system provides several sophisticated features:

1. **Multi-Rule Integration**: Applies both standard co-expression rules and meta-rules
2. **Dynamic Cell Creation**: Constructs cell entities with appropriate properties
3. **Confidence Scoring**: Assigns confidence scores to each identified cell
4. **Spatial Awareness**: Ensures proper spatial relationships between gene expressions
5. **Conflict Resolution**: Implements sophisticated resolution of overlapping cells
6. **Visualization**: Generates informative visualizations of identified cells
7. **Comprehensive Reporting**: Produces detailed reports on identification results

## Architecture

The rule application system consists of several main components:

1. **RDF Data Manager**: Handles loading and querying spatial data in RDF format
2. **Rule Processor**: Loads and applies SPARQL rules to identify cells
3. **Conflict Resolver**: Resolves conflicts between overlapping cells
4. **Results Manager**: Handles output generation and visualization

## RDF Data Manager

The RDF Data Manager handles all interactions with the RDF representation of spatial data:

```python
class RDFDataManager:
    def __init__(self, ttl_file):
        """Initialize with path to Turtle file containing spatial data."""
        self.graph = rdflib.Graph()
        self.graph.parse(ttl_file, format="turtle")
        self.namespace_manager = self.graph.namespace_manager
        
    def get_namespaces(self):
        """Return all namespaces used in the graph."""
        return dict(self.graph.namespaces())
        
    def execute_query(self, query_string):
        """Execute a SPARQL query on the graph."""
        return self.graph.query(query_string)
        
    def execute_construct(self, query_string):
        """Execute a SPARQL CONSTRUCT query and return the new triples."""
        result_graph = self.graph.query(query_string).graph
        return result_graph
```

This component provides a clean interface for:
- Loading spatial data from Turtle files
- Executing SPARQL queries to identify patterns
- Managing RDF namespaces for proper URIs
- Constructing new graph entities (cells)

## Rule Processor

The Rule Processor is responsible for loading and applying SPARQL rules:

```python
class RuleProcessor:
    def __init__(self, rdf_manager, rule_dir):
        """Initialize with RDF manager and directory containing rule files."""
        self.rdf_manager = rdf_manager
        self.rule_dir = rule_dir
        self.rules = []
        self.load_rules()
        
    def load_rules(self):
        """Load all SPARQL rule files from the rules directory."""
        rule_files = glob.glob(os.path.join(self.rule_dir, "*.rq"))
        for rule_file in rule_files:
            with open(rule_file, 'r') as f:
                rule_content = f.read()
                rule_id = os.path.basename(rule_file).replace('.rq', '')
                self.rules.append({
                    'id': rule_id,
                    'content': rule_content,
                    'type': 'standard'
                })
    
    def load_meta_rules(self, meta_rule_file):
        """Load meta-rules from the specified file."""
        if not os.path.exists(meta_rule_file):
            return
            
        with open(meta_rule_file, 'r') as f:
            content = f.read()
            
        # Split content by rule delimiter
        rule_delimiter = "# Meta-rule combining patterns"
        meta_rules_content = content.split(rule_delimiter)
        
        # Skip the first empty part
        if meta_rules_content[0].strip() == '':
            meta_rules_content = meta_rules_content[1:]
            
        for i, rule_content in enumerate(meta_rules_content):
            if rule_content.strip():
                self.rules.append({
                    'id': f"meta_rule_{i+1}",
                    'content': rule_delimiter + rule_content,
                    'type': 'meta'
                })
    
    def apply_rule(self, rule):
        """Apply a SPARQL rule to identify cells."""
        result_graph = self.rdf_manager.execute_construct(rule['content'])
        return result_graph
        
    def apply_all_rules(self):
        """Apply all loaded rules to identify cells."""
        all_cells = rdflib.Graph()
        
        for rule in self.rules:
            result = self.apply_rule(rule)
            if result and len(result) > 0:
                all_cells += result
                
        return all_cells
```

Key functions:
- Loading standard SPARQL rules from files
- Loading meta-rules from the meta-rule file
- Applying rules to the spatial data
- Aggregating results from all rules

## Rule Application Process

The rule application process follows these steps:

1. **Load Spatial Data**: Convert spatial transcriptomics data to RDF format
2. **Load Rules**: Read SPARQL CONSTRUCT queries from rule files
3. **Apply Rules**: Execute each rule against the spatial data
4. **Collect Cells**: Gather all identified cells from rule application
5. **Extract Properties**: Extract properties of each identified cell
6. **Remove Duplicates**: Ensure each unique cell is only represented once
7. **Filter Low Confidence**: Remove cells below confidence threshold
8. **Resolve Conflicts**: Handle overlapping cell identifications

### SPARQL Rule Application

The core of the rule application process is the execution of SPARQL CONSTRUCT queries. Here's an example of how a rule is executed:

```python
def apply_rule_to_spatial_data(rule_content, spatial_graph):
    """Apply a SPARQL rule to spatial data."""
    # Create a query object
    query = prepareQuery(rule_content, initNs=namespaces)
    
    # Execute the CONSTRUCT query
    result_graph = spatial_graph.query(query).graph
    
    # Return the constructed cells
    return result_graph
```

### Cell Extraction and Processing

After rules are applied, we extract the identified cells from the result graph:

```python
def extract_cells_from_graph(cells_graph):
    """Extract cell information from RDF graph."""
    cell_query = """
    SELECT DISTINCT ?cell ?cellType ?x ?y ?radius ?confidence ?ruleType
    WHERE {
        ?cell a astro:SpatialCell ;
              astro:hasCellType ?cellType ;
              astro:hasXCoordinate ?x ;
              astro:hasYCoordinate ?y ;
              astro:hasRadius ?radius ;
              astro:hasConfidence ?confidence .
        OPTIONAL { ?cell astro:derivedFromMetaRule ?derivedFromMeta . }
        BIND(IF(BOUND(?derivedFromMeta) && ?derivedFromMeta = true, "meta", "standard") AS ?ruleType)
    }
    """
    
    cells_result = cells_graph.query(cell_query)
    
    cells = []
    for row in cells_result:
        cell = {
            'uri': str(row.cell),
            'cell_type': str(row.cellType),
            'x': float(row.x),
            'y': float(row.y),
            'radius': float(row.radius),
            'confidence': float(row.confidence),
            'rule_type': str(row.ruleType)
        }
        cells.append(cell)
        
    return cells
```

## Visualization and Output

The module generates several output formats:

1. **CSV File**: Tabular data with all cell information
2. **Visualization**: Plot showing all identified cells
3. **RDF Export**: Complete RDF graph with all identified cells
4. **Conflict Report**: Detailed report on conflict resolution

Example of the visualization generation:

```python
def visualize_cells(cells_df, output_file):
    """Generate visualization of identified cells."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Define colors for different cell types
    cell_type_colors = {
        'Astrocytes': 'blue',
        'Astrocytes1': 'red',
        'Unknown': 'gray'
    }
    
    # Plot each cell as a circle
    for _, cell in cells_df.iterrows():
        cell_type = cell['cell_type']
        color = cell_type_colors.get(cell_type, 'gray')
        alpha = min(1.0, 0.5 + cell['confidence'] * 0.5)  # Higher confidence = more opaque
        
        circle = Circle((cell['x'], cell['y']), cell['radius'], 
                        fill=True, alpha=alpha, color=color, 
                        linewidth=1, edgecolor='black')
        ax.add_patch(circle)
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                   label=cell_type, markersize=10)
        for cell_type, color in cell_type_colors.items()
        if cell_type in cells_df['cell_type'].values
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Set axis labels and title
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_title('Identified Cells')
    
    # Equal aspect ratio and limits
    ax.set_aspect('equal')
    x_min, x_max = cells_df['x'].min() - 100, cells_df['x'].max() + 100
    y_min, y_max = cells_df['y'].min() - 100, cells_df['y'].max() + 100
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
```

## Integration with Enhanced Conflict Resolution

The rule application module integrates with the enhanced conflict resolution system:

```python
def process_with_enhanced_resolution(cells_df, overlap_threshold=0.3):
    """Process cells using the enhanced conflict resolution system."""
    from conflict_resolution import ConflictResolutionManager
    
    # Create a conflict resolution manager
    resolver = ConflictResolutionManager(
        cells_df=cells_df,
        overlap_threshold=overlap_threshold,
        generate_visuals=True
    )
    
    # Resolve conflicts
    resolved_cells = resolver.resolve_conflicts()
    
    # Generate report
    resolver.generate_report()
    
    return resolved_cells
```

## Parallelization for Large Datasets

For large datasets, the module implements parallel rule application:

```python
def apply_rules_parallel(rules, spatial_graph, num_processes=4):
    """Apply rules in parallel for large datasets."""
    from multiprocessing import Pool
    
    # Function to apply a single rule
    def apply_single_rule(rule):
        result = apply_rule_to_spatial_data(rule['content'], spatial_graph)
        return rule['id'], result
    
    # Apply rules in parallel
    with Pool(processes=num_processes) as pool:
        results = pool.map(apply_single_rule, rules)
    
    # Combine results
    all_cells = rdflib.Graph()
    for rule_id, result in results:
        if result and len(result) > 0:
            all_cells += result
    
    return all_cells
```

## Command-Line Interface

The module provides a comprehensive command-line interface:

```python
def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Apply SPARQL rules to identify cells in spatial data.')
    
    parser.add_argument('--spatial-data', required=True,
                      help='Path to Turtle file containing spatial data')
    parser.add_argument('--rules-dir', required=True,
                      help='Directory containing SPARQL rule files')
    parser.add_argument('--output-dir', required=True,
                      help='Directory for output files')
    parser.add_argument('--meta-rules-file', default=None,
                      help='Path to file containing meta-rules')
    parser.add_argument('--min-confidence', type=float, default=0.7,
                      help='Minimum confidence threshold (default: 0.7)')
    parser.add_argument('--resolve-conflicts', action='store_true',
                      help='Apply conflict resolution')
    parser.add_argument('--enhanced-resolution', action='store_true',
                      help='Use enhanced conflict resolution')
    parser.add_argument('--overlap-threshold', type=float, default=0.3,
                      help='Overlap threshold for conflict detection (default: 0.3)')
    parser.add_argument('--parallel', action='store_true',
                      help='Use parallel processing for rule application')
    parser.add_argument('--num-processes', type=int, default=4,
                      help='Number of parallel processes to use (default: 4)')
    
    return parser.parse_args()
```

This interface provides fine-grained control over:
- Input data and rule locations
- Output formatting and destinations
- Confidence thresholds for filtering
- Conflict resolution settings
- Performance optimization options

## Example Usage

Basic usage of rule application:

```bash
python scripts/identify_cells.py \
    --spatial-data data/processed/spatial_data.ttl \
    --rules-dir data/processed/rules \
    --output-dir data/processed/results \
    --min-confidence 0.7
```

Advanced usage with meta-rules and enhanced conflict resolution:

```bash
python scripts/identify_cells.py \
    --spatial-data data/processed/spatial_data.ttl \
    --rules-dir data/processed/rules \
    --meta-rules-file data/processed/meta_rules/meta_rules.rq \
    --output-dir data/processed/results \
    --min-confidence 0.7 \
    --resolve-conflicts \
    --enhanced-resolution \
    --overlap-threshold 0.3
```

## Output Files

The rule application module generates several output files:

| File | Description | Format |
|---|---|---|
| `identified_cells.csv` | Table of all identified cells | CSV |
| `cell_visualization.png` | Visual representation of cells | PNG |
| `cells.ttl` | RDF representation of identified cells | Turtle |
| `conflict_resolution_report.md` | Detailed report on conflict resolution | Markdown |
| `conflict_visuals/` | Directory with conflict visualizations | PNG |

Example output format for the CSV file:

```
cell_id,cell_type,x,y,radius,confidence,rule_type,genes
cell_1,Astrocytes,1230.5,840.2,45.8,0.95,meta,"NRXN1,GPM6A,ADGRB3"
cell_2,Astrocytes1,980.7,730.5,38.2,0.85,clique,"LRP1B,NRG3,GPM6A"
cell_3,Astrocytes,1450.2,920.8,42.5,0.72,pair,"LSAMP,CNTN1"
```

## Performance Considerations

For efficient rule application, consider the following:

1. **RDF Store Selection**: The backend RDF store can significantly impact performance
2. **Parallelization**: Use the `--parallel` flag for large datasets
3. **Rule Optimization**: More specific rules perform better than generic ones
4. **Confidence Filtering**: Use appropriate confidence thresholds to reduce post-processing
5. **Selective Meta-Rule Application**: Only apply meta-rules when needed

## Use Cases

### Regular Rule Application

For basic identification of cells:

```python
# Simple rule application
rdf_manager = RDFDataManager(args.spatial_data)
rule_processor = RuleProcessor(rdf_manager, args.rules_dir)
cells_graph = rule_processor.apply_all_rules()
cells = extract_cells_from_graph(cells_graph)
```

### Filtering by Cell Type

For focusing on specific cell types:

```python
# Filter for specific cell types
filtered_cells = [cell for cell in cells if cell['cell_type'] in target_cell_types]
```

### Interactive Analysis

For exploring results interactively:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load results
cells_df = pd.read_csv('data/processed/results/identified_cells.csv')

# Count by cell type
cell_type_counts = cells_df['cell_type'].value_counts()
print(cell_type_counts)

# Plot distributions
plt.figure(figsize=(10, 6))
cells_df.groupby('cell_type')['confidence'].plot.kde()
plt.xlabel('Confidence Score')
plt.ylabel('Density')
plt.title('Confidence Score Distribution by Cell Type')
plt.legend()
plt.show()
```

## Troubleshooting

### Common Issues

| Issue | Possible Solution |
|---|---|
| No cells identified | Check rule files and spatial data format |
| Low confidence scores | Adjust rule parameters or check data quality |
| Overlapping cells | Enable enhanced conflict resolution |
| Performance issues | Use parallelization or optimize rules |
| Memory errors | Process cell types separately |

### Diagnostic Messages

The module provides detailed logging:

```python
import logging

def setup_logging(log_file=None):
    """Set up logging configuration."""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format,
                       filename=log_file)
    
    # Add console handler
    if not log_file:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(logging.Formatter(log_format))
        logging.getLogger('').addHandler(console)
```

Important log messages to watch for:

- "Loaded X rules from directory Y"
- "Applied rule Z, identified N cells"
- "Found N conflicts between cells"
- "Final result: X cells after conflict resolution"

## References

1. RDF and SPARQL:
   - W3C. (2013). "SPARQL 1.1 Query Language." W3C Recommendation.

2. Spatial data processing:
   - Deng, Y., et al. (2021). "Spatial transcriptomics for the construction of a high-resolution cell atlas." Journal of Genetics and Genomics, 48(10), 878-886.

3. Cell identification algorithms:
   - Stuart, T., et al. (2019). "Comprehensive Integration of Single-Cell Data." Cell, 177(7), 1888-1902. 