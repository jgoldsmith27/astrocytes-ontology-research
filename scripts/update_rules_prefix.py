#!/usr/bin/env python3
"""
Script to update existing SPARQL rules to use the new cell ontology prefix.

This script reads all .rq and .sparql files in the data/processed directory,
updates the prefix from astro: to cell: and updates the ontology URL from
http://example.org/ontology/astrocyte/ to http://example.org/ontology/cell/
"""

import os
import re
import glob
import argparse
from pathlib import Path

def update_file(file_path):
    """Update a single SPARQL rule file."""
    print(f"Updating file: {file_path}")
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Update the prefix declaration
    content = re.sub(
        r'PREFIX astro: <http://example.org/ontology/astrocyte/>',
        'PREFIX cell: <http://example.org/ontology/cell/>',
        content
    )
    
    # Update the astro:X references to cell:X
    content = re.sub(r'astro:', 'cell:', content)
    
    # Update astrocyte URL in URIs (if any)
    content = re.sub(
        r'<http://example.org/ontology/astrocyte/([^>]+)>',
        r'<http://example.org/ontology/cell/\1>',
        content
    )
    
    with open(file_path, 'w') as f:
        f.write(content)
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Update SPARQL rules to use the new cell ontology prefix')
    parser.add_argument('--directory', type=str, default='data/processed', 
                      help='Directory containing rule files (default: data/processed)')
    parser.add_argument('--recursive', action='store_true',
                      help='Search recursively in subdirectories')
    
    args = parser.parse_args()
    
    # Find all SPARQL rule files
    if args.recursive:
        sparql_files = glob.glob(os.path.join(args.directory, '**/*.sparql'), recursive=True)
        sparql_files += glob.glob(os.path.join(args.directory, '**/*.rq'), recursive=True)
    else:
        sparql_files = glob.glob(os.path.join(args.directory, '*.sparql'))
        sparql_files += glob.glob(os.path.join(args.directory, '*.rq'))
    
    print(f"Found {len(sparql_files)} SPARQL rule files to update")
    
    # Update each file
    updated_count = 0
    for file_path in sparql_files:
        if update_file(file_path):
            updated_count += 1
    
    print(f"Updated {updated_count} files")

if __name__ == "__main__":
    main() 