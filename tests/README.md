# Testing Framework for Astrocytes Ontology Research

This directory contains tests for the Astrocytes Ontology Research Project. The testing framework is designed to validate each component independently and also test the integration between components.

## Test Structure

- **Unit Tests**: Located in `tests/unit/`. These test individual functions and classes in isolation.
- **Integration Tests**: Located in `tests/integration/`. These test how different components work together.
- **Test Fixtures**: Common test data and fixtures are defined in `tests/conftest.py`.

## Running Tests

### Prerequisites

Make sure you have installed the testing dependencies:

```bash
pip install -r requirements-test.txt
```

### Using the Test Runner

The easiest way to run the tests is using the test runner script:

```bash
# Run all tests
./run_tests.py

# Run only unit tests
./run_tests.py --unit-only

# Run only integration tests
./run_tests.py --integration-only

# Run with verbose output
./run_tests.py -v

# Skip coverage report
./run_tests.py --no-coverage
```

### Using pytest Directly

You can also run the tests directly using pytest:

```bash
# Run all tests
pytest

# Run unit tests only
pytest tests/unit/

# Run a specific test file
pytest tests/unit/test_coexpression_rules.py

# Run a specific test class
pytest tests/unit/test_conflict_resolution.py::TestCellClass

# Run a specific test function
pytest tests/unit/test_conflict_resolution.py::TestCellClass::test_distance_calculation
```

## Coverage Reports

When running tests with the default settings, a coverage report is generated. You can view:

- **Terminal Report**: Shown in the console after the tests run
- **HTML Report**: Generated in the `coverage_report/` directory, view by opening `coverage_report/index.html` in a browser

## Writing New Tests

### Unit Tests

1. Create a new file in `tests/unit/` named `test_<component_name>.py`
2. Import the module you want to test
3. Write test classes/functions using pytest conventions
4. Use the fixtures from `conftest.py` where appropriate

### Integration Tests

1. Create a new file in `tests/integration/` named `test_<feature_name>.py`
2. Import all required components
3. Write tests that validate how components interact
4. Use mock data and fixtures to simulate the real environment

## Troubleshooting Tests

If you encounter issues with tests:

1. **Import Errors**: Make sure your test is properly importing the modules from the `scripts/` directory
2. **Module Not Found**: Check that all `__init__.py` files exist in the test directories
3. **Missing Dependencies**: Verify that all requirements in `requirements-test.txt` are installed
4. **Test Discovery Issues**: Make sure your test file and functions start with `test_`

## Code Coverage Goals

We aim for at least 80% code coverage across the codebase. Critical components like:

- Cell identification logic
- Conflict resolution
- Rule generation

Should have closer to 100% coverage due to their importance in the pipeline. 