# Astrocyte Ontology Research Test Suite

This directory contains tests for the Astrocyte Ontology Research project. The tests are organized into three categories:

- **Unit tests**: Test individual components in isolation
- **Integration tests**: Test interactions between components
- **Regression tests**: Test that previously fixed bugs remain fixed

## Test Structure

```
tests/
├── __init__.py
├── run_tests.py
├── unit/
│   ├── __init__.py
│   ├── test_astrocyte_filter.py
│   └── test_astrocyte_validation.py
├── integration/
│   ├── __init__.py
│   └── test_ontology_connections.py
└── regression/
    ├── __init__.py
    └── test_astrocyte_validation_regression.py
```

## Running Tests

You can run the tests using the `run_tests.py` script. The script supports running all tests or specific types of tests.

### Running All Tests

```bash
python tests/run_tests.py
```

### Running Specific Test Types

```bash
# Run unit tests
python tests/run_tests.py --type unit

# Run integration tests
python tests/run_tests.py --type integration

# Run regression tests
python tests/run_tests.py --type regression
```

### Verbose Output

You can add the `--verbose` or `-v` flag to get more detailed output:

```bash
python tests/run_tests.py --verbose
```

## Test Descriptions

### Unit Tests

- **test_astrocyte_filter.py**: Tests the astrocyte filtering functionality in the Single_Cell module.
- **test_astrocyte_validation.py**: Tests the astrocyte validation functionality in the Spatial_Data module.

### Integration Tests

- **test_ontology_connections.py**: Tests the connections between spatial and single-cell ontologies.

### Regression Tests

- **test_astrocyte_validation_regression.py**: Tests that previously fixed bugs in the astrocyte validation functionality remain fixed.

## Adding New Tests

When adding new tests, follow these guidelines:

1. Place unit tests in the `unit/` directory
2. Place integration tests in the `integration/` directory
3. Place regression tests in the `regression/` directory
4. Name test files with the prefix `test_`
5. Name test methods with the prefix `test_`
6. Use the `unittest` framework for all tests
7. Include docstrings for all test methods

## Test Dependencies

The tests require the following dependencies:

- Python 3.6+
- rdflib
- pandas
- numpy
- scanpy
- anndata
- scipy
- matplotlib

You can install these dependencies using pip:

```bash
pip install rdflib pandas numpy scanpy anndata scipy matplotlib
``` 