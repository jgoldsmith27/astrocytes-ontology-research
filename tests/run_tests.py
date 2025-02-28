#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import unittest
import argparse
import time
import importlib.util
from pathlib import Path

def import_module_from_path(module_path):
    """Import a module from a file path."""
    module_name = os.path.basename(module_path).replace('.py', '')
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def discover_and_run_tests(test_type=None, verbose=False):
    """
    Discover and run tests of the specified type.
    
    Parameters:
    -----------
    test_type : str or None
        Type of tests to run: 'unit', 'integration', 'regression', or None for all tests
    verbose : bool
        Whether to show verbose output
    
    Returns:
    --------
    bool
        True if all tests passed, False otherwise
    """
    # Get the project root directory
    project_root = Path(__file__).parent.parent
    tests_dir = project_root / 'tests'
    
    # Set up test loader
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Discover tests based on test_type
    if test_type is None or test_type == 'all':
        # Discover all tests
        for test_dir in ['unit', 'integration', 'regression']:
            dir_path = tests_dir / test_dir
            if dir_path.exists():
                suite.addTests(loader.discover(str(dir_path), pattern='test_*.py'))
    else:
        # Discover tests of the specified type
        dir_path = tests_dir / test_type
        if dir_path.exists():
            suite.addTests(loader.discover(str(dir_path), pattern='test_*.py'))
    
    # Run tests
    verbosity = 2 if verbose else 1
    runner = unittest.TextTestRunner(verbosity=verbosity)
    start_time = time.time()
    result = runner.run(suite)
    duration = time.time() - start_time
    
    # Print summary
    print("\n" + "=" * 60)
    print(f"📊 TEST SUMMARY (Duration: {duration:.2f}s):")
    print(f"✅ Tests run: {result.testsRun}")
    print(f"❌ Failures: {len(result.failures)}")
    print(f"❌ Errors: {len(result.errors)}")
    print(f"⚠️ Skipped: {len(result.skipped)}")
    
    if not result.failures and not result.errors:
        print("\n🎉🎉🎉 ALL TESTS PASSED! 🎉🎉🎉")
        return True
    else:
        print("\n💣 SOME TESTS FAILED! CHECK THE LOGS!")
        return False

def main():
    """Main function to parse arguments and run tests."""
    parser = argparse.ArgumentParser(description='Run tests for the Astrocyte Ontology Research project.')
    parser.add_argument('--type', choices=['unit', 'integration', 'regression', 'all'], 
                        default='all', help='Type of tests to run')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show verbose output')
    
    args = parser.parse_args()
    
    print("\n🚀🚀🚀 RUNNING ASTROCYTE ONTOLOGY TESTS 🚀🚀🚀")
    print("=" * 60)
    print(f"Test type: {args.type}")
    print(f"Verbose: {args.verbose}")
    print("=" * 60)
    
    success = discover_and_run_tests(args.type, args.verbose)
    
    # Return exit code based on test results
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main() 