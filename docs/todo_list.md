# TODO List for Project Improvements

## Completed Tasks

### ✅ Implement Higher-Order Co-expression Patterns (Meta-Rules)

**Implementation Steps:**

1. **Research Phase:**
   - [ ] Research DBSCAN and HDBSCAN algorithms and their implementations in Python
   - [ ] Determine which algorithm is more suitable for our spatial data
   - [ ] Identify key parameters (eps, min_samples) for density-based clustering

2. **Development Phase:**
   - [ ] Create a new module `spatial_clustering.py` with density-based clustering functions
   - [ ] Implement methods to extract spatial points for each potential cell type
   - [ ] Apply clustering to identify cell boundaries beyond simple circles
   - [ ] Calculate confidence scores based on cluster characteristics

3. **Integration Phase:**
   - [ ] Modify `identify_astrocytes.py` to incorporate density-based clustering results
   - [ ] Create a pipeline that first uses co-expression rules, then refines with clustering
   - [ ] Develop visualization tools to show irregular cell boundaries

4. **Validation Phase:**
   - [ ] Compare results with current center-based approach
   - [ ] Validate against known cell distributions if reference data is available
   - [ ] Optimize parameters based on biological plausibility

**Expected Benefits:**
- More accurate cell boundary detection, especially for non-spherical cells
- Better handling of varying expression density within cells
- Reduced false positives from overlapping cells
- More realistic morphological representation of identified cells

**Resources Needed:**
- scikit-learn (for DBSCAN) or hdbscan library
- Reference datasets with known cell boundaries for validation
- Additional compute resources for more intensive clustering calculations

## High Priority

### 1. Enhance Rule Conflict Resolution

**Description:** Develop a more sophisticated approach for resolving conflicts when multiple cell types are identified in the same region.

**Implementation Steps:**
- [ ] Implement a scoring system that considers confidence, marker specificity, and negative markers
- [ ] Develop a cell type hierarchy based on biological knowledge
- [ ] Add visualization tools to highlight regions with type conflicts

### 2. Optimize Parameters by Cell Type

**Description:** Different cell types may require different parameter settings for optimal detection.

**Implementation Steps:**
- [ ] Create a configuration system for cell-type specific parameters
- [ ] Experiment with different parameter sets for major cell types
- [ ] Document optimal parameters for each cell type

## Lower Priority

### 4. Performance Optimization

**Description:** Optimize the pipeline for processing larger datasets more efficiently.

**Implementation Steps:**
- [ ] Profile code to identify bottlenecks
- [ ] Implement parallelization for rule generation and application
- [ ] Optimize memory usage for large spatial datasets

### 5. Interactive Visualization Tools

**Description:** Develop interactive tools for exploring identified cells.

**Implementation Steps:**
- [ ] Create web-based visualization for spatial cell distributions
- [ ] Add filtering capabilities by cell type, confidence, etc.
- [ ] Implement zoom and exploration features

## Future Directions

### 6. Implement Density-Based Clustering for Cell Detection

**Description:** Develop a system to analyze relationships between existing co-expression patterns and clique rules to create more powerful composite rules for cell identification.

**Implementation Details:**
- Created `generate_meta_rules.py` that analyzes rule co-occurrence in single-cell data
- Implemented statistical association detection with lift and Fisher's exact test
- Generated composite SPARQL rules combining multiple co-expression patterns
- Added visualization of rule relationships as a network graph
- Integrated into the pipeline with prioritization during conflict resolution

**Benefits Achieved:**
- Increased specificity through multi-pattern validation
- Better distinction between closely related cell types
- Capturing of complex biological relationships between genes
- Higher confidence scoring for cells with multiple evidence types

**Description:** Enhance cell detection by implementing density-based clustering (DBSCAN or HDBSCAN) to identify clusters of genes with irregular shapes, improving upon the current center-based approach.

### 7. Machine Learning for Rule Generation

**Description:** Explore using machine learning to generate and optimize co-expression rules.

**Implementation Steps:**
- [ ] Research appropriate ML approaches for this domain
- [ ] Collect training data from validated cell identifications
- [ ] Develop a prototype ML-based rule generator