#!/usr/bin/env python3
"""
Test runner for the Astrocytes Ontology Research Project.
Runs all unit and integration tests and generates a coverage report.
"""

import os
import subprocess
import sys
import argparse
from pathlib import Path

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run tests for Astrocytes Ontology Research")
    parser.add_argument("--unit-only", action="store_true", help="Run only unit tests")
    parser.add_argument("--integration-only", action="store_true", help="Run only integration tests")
    parser.add_argument("--no-coverage", action="store_true", help="Skip coverage report")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    return parser.parse_args()

def run_tests(args):
    """Run the tests based on provided arguments."""
    # Set up base command
    cmd = ["pytest"]
    
    # Add verbosity if requested
    if args.verbose:
        cmd.append("-v")
    
    # Add coverage options unless disabled
    if not args.no_coverage:
        cmd.extend(["--cov=scripts", "--cov-report=term", "--cov-report=html:coverage_report"])
    
    # Determine which tests to run
    if args.unit_only and not args.integration_only:
        cmd.append("tests/unit/")
    elif args.integration_only and not args.unit_only:
        cmd.append("tests/integration/")
    else:
        cmd.append("tests/")
    
    # Print the command
    cmd_str = " ".join(cmd)
    print(f"Running: {cmd_str}")
    
    # Run the tests
    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"Tests failed with exit code {e.returncode}")
        return e.returncode

def main():
    """Main function to run tests."""
    args = parse_args()
    
    # Ensure we're in the project root directory
    project_root = Path(__file__).parent
    os.chdir(project_root)
    
    # Create __init__.py files if they don't exist to ensure imports work
    for directory in ["tests", "tests/unit", "tests/integration"]:
        init_file = Path(directory) / "__init__.py"
        if not init_file.exists():
            init_file.touch()
    
    # Run the tests
    return run_tests(args)

if __name__ == "__main__":
    sys.exit(main()) 