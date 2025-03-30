#!/usr/bin/env python
"""
Unified Ligand-Receptor Pair Analysis Workflow

This script provides a comprehensive workflow for analyzing ligand-receptor pairs
from the CellChat database in relation to single-cell expression data. The analysis
includes:

1. Extracting ligand-receptor pairs from the CellChat database
2. Analyzing expression data directly from the original dataset
3. Generating expression summaries for cell types and genes
4. Analyzing autocrine signaling patterns (same cell as sender and receiver)
5. Analyzing all communication patterns between cell types
6. Generating comprehensive reports and visualizations

Usage:
  python unified_lr_analysis.py [OPTIONS]

Options:
  --steps STEPS        Analysis steps to run (extract, summaries, autocrine, communication, reports, all)
  --clean              Clean output directory before running
  --threshold FLOAT    Expression threshold (default: 0.1)

Examples:
  # Run all analysis steps with default settings
  python unified_lr_analysis.py
  
  # Run only autocrine and communication analysis with a custom threshold
  python unified_lr_analysis.py --steps autocrine communication --threshold 0.2
  
  # Clean output directory and regenerate everything
  python unified_lr_analysis.py --clean

Output:
  All results are saved to the output directory. Key output files include:
  - Cell type expression summaries
  - Gene presence summaries
  - Autocrine signaling analysis
  - Cell-cell communication analysis
  - Pathway statistics
  - Summary reports in markdown format
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import argparse
from pathlib import Path

# Constants
OUTPUT_DIR = "output"
DATA_DIR = "data"

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Data loading functions
def load_cellchat_lr_pairs():
    """Load ligand-receptor pairs from CellChat database"""
    lr_pairs_file = os.path.join(OUTPUT_DIR, "cellchat_lr_pairs.csv")
    if not os.path.exists(lr_pairs_file):
        raise FileNotFoundError(f"CellChat LR pairs file not found: {lr_pairs_file}")
    return pd.read_csv(lr_pairs_file)

def load_expression_data():
    """Load expression data for ligands and receptors from original dataset"""
    # Files are in the root directory, not in a data subdirectory
    ligand_file = "ligand_expression_by_cell_type.csv"
    receptor_file = "receptor_expression_by_cell_type.csv"
    
    if not os.path.exists(ligand_file) or not os.path.exists(receptor_file):
        raise FileNotFoundError(f"Expression data files not found: {ligand_file} or {receptor_file}")
    
    # Read data - rows are cell types, columns are genes
    ligand_df_original = pd.read_csv(ligand_file, index_col=0)
    receptor_df_original = pd.read_csv(receptor_file, index_col=0)
    
    # Transpose to have genes as rows and cell types as columns
    ligand_df = ligand_df_original.transpose()
    receptor_df = receptor_df_original.transpose()
    
    return ligand_df, receptor_df

# Analysis components
def extract_cellchat_data():
    """Extract LR pairs from CellChat database"""
    print("Extracting ligand-receptor pairs from CellChat database...")
    
    # Check if output file already exists
    output_file = os.path.join(OUTPUT_DIR, "cellchat_lr_pairs.csv")
    if os.path.exists(output_file):
        print(f"CellChat data already extracted to {output_file}")
        return
    
    # Run the R script to extract data
    r_script_path = "scripts/extract_cellchat_pairs.R"
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"R script not found: {r_script_path}")
    
    import subprocess
    result = subprocess.run(["Rscript", r_script_path], capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error extracting CellChat data: {result.stderr}")
        raise RuntimeError("Failed to extract CellChat data")
    
    print(result.stdout)
    print(f"CellChat data extracted to {output_file}")

def generate_expression_summaries(expression_threshold=0.1):
    """Generate cell type and gene expression summaries"""
    print("Generating expression summaries...")
    
    # Output file paths
    cell_type_summary_file = os.path.join(OUTPUT_DIR, "cell_type_expression_summary.csv")
    gene_presence_file = os.path.join(OUTPUT_DIR, "gene_presence_summary.csv")
    gene_counts_file = os.path.join(OUTPUT_DIR, "gene_presence_counts_summary.csv")
    
    # Load expression data and CellChat pairs
    ligand_df, receptor_df = load_expression_data()
    lr_pairs = load_cellchat_lr_pairs()
    
    # Get all cell types
    cell_types = ligand_df.columns.tolist()
    
    # Get unique ligands and receptors from CellChat that exist in the expression data
    unique_ligands = sorted([gene for gene in lr_pairs['ligand'].unique() if gene in ligand_df.index])
    unique_receptors = sorted([gene for gene in lr_pairs['receptor'].unique() if gene in receptor_df.index])
    
    missing_ligands = sorted([gene for gene in lr_pairs['ligand'].unique() if gene not in ligand_df.index])
    missing_receptors = sorted([gene for gene in lr_pairs['receptor'].unique() if gene not in receptor_df.index])
    
    print(f"Found {len(unique_ligands)} ligands and {len(unique_receptors)} receptors from CellChat in expression data")
    print(f"Missing {len(missing_ligands)} ligands and {len(missing_receptors)} receptors from expression data")
    
    # Create dictionary to store cell type expression data
    cell_type_data = []
    
    # Analyze each cell type
    for cell_type in cell_types:
        # Find expressed ligands and receptors based on threshold
        expressed_ligands = ligand_df.index[ligand_df[cell_type] >= expression_threshold].tolist()
        expressed_receptors = receptor_df.index[receptor_df[cell_type] >= expression_threshold].tolist()
        
        # Filter for ligands/receptors that are in CellChat DB
        cellchat_ligands = [gene for gene in expressed_ligands if gene in unique_ligands]
        cellchat_receptors = [gene for gene in expressed_receptors if gene in unique_receptors]
        
        # Add to data
        cell_type_data.append({
            'cell_type': cell_type,
            'expressed_ligands': ','.join(expressed_ligands),
            'expressed_receptors': ','.join(expressed_receptors),
            'cellchat_ligands': ','.join(cellchat_ligands),
            'cellchat_receptors': ','.join(cellchat_receptors),
            'ligand_count': len(expressed_ligands),
            'receptor_count': len(expressed_receptors),
            'cellchat_ligand_count': len(cellchat_ligands),
            'cellchat_receptor_count': len(cellchat_receptors)
        })
    
    # Create cell type summary dataframe
    cell_type_df = pd.DataFrame(cell_type_data)
    
    # Create gene presence summary
    gene_presence_data = []
    
    # Process ligands in CellChat that exist in expression data
    for gene in unique_ligands:
        expressed_in = [cell_type for cell_type in cell_types 
                        if ligand_df.at[gene, cell_type] >= expression_threshold]
        
        gene_presence_data.append({
            'gene': gene,
            'type': 'ligand',
            'cell_types': ','.join(expressed_in),
            'cell_type_count': len(expressed_in)
        })
    
    # Process receptors in CellChat that exist in expression data
    for gene in unique_receptors:
        expressed_in = [cell_type for cell_type in cell_types 
                        if receptor_df.at[gene, cell_type] >= expression_threshold]
        
        gene_presence_data.append({
            'gene': gene,
            'type': 'receptor',
            'cell_types': ','.join(expressed_in),
            'cell_type_count': len(expressed_in)
        })
    
    # Create gene presence dataframe
    gene_presence_df = pd.DataFrame(gene_presence_data)
    
    # Create summary counts
    gene_counts = {
        'total_ligands_in_cellchat': len(lr_pairs['ligand'].unique()),
        'total_receptors_in_cellchat': len(lr_pairs['receptor'].unique()),
        'ligands_in_expression_data': len(unique_ligands),
        'receptors_in_expression_data': len(unique_receptors),
        'missing_ligands': len(missing_ligands),
        'missing_receptors': len(missing_receptors),
        'expressed_ligands': len(gene_presence_df[(gene_presence_df['type'] == 'ligand') & 
                                                 (gene_presence_df['cell_type_count'] > 0)]),
        'expressed_receptors': len(gene_presence_df[(gene_presence_df['type'] == 'receptor') & 
                                                   (gene_presence_df['cell_type_count'] > 0)])
    }
    
    # Convert to dataframe
    gene_counts_df = pd.DataFrame([gene_counts])
    
    # Save to CSV
    cell_type_df.to_csv(cell_type_summary_file, index=False)
    gene_presence_df.to_csv(gene_presence_file, index=False)
    gene_counts_df.to_csv(gene_counts_file, index=False)
    
    print(f"Cell type expression summary saved to {cell_type_summary_file}")
    print(f"Gene presence summary saved to {gene_presence_file}")
    print(f"Gene presence counts summary saved to {gene_counts_file}")
    
    return unique_ligands, unique_receptors

def analyze_autocrine_signaling(expression_threshold=0.1):
    """Analyze autocrine signaling (same cell type as sender and receiver)"""
    print("Analyzing autocrine signaling patterns...")
    
    # Output file paths
    frequency_file = os.path.join(OUTPUT_DIR, "autocrine_lr_pair_frequency_ranking.csv")
    specificity_file = os.path.join(OUTPUT_DIR, "autocrine_lr_pair_cell_type_specificity.csv")
    usage_file = os.path.join(OUTPUT_DIR, "autocrine_cell_type_lr_pair_usage.csv")
    pathway_file = os.path.join(OUTPUT_DIR, "autocrine_pathway_statistics.csv")
    summary_file = os.path.join(OUTPUT_DIR, "autocrine_ligand_receptor_summary.md")
    
    # Load data
    ligand_df, receptor_df = load_expression_data()
    lr_pairs = load_cellchat_lr_pairs()
    
    # Get all cell types
    cell_types = ligand_df.columns.tolist()
    
    # Analysis results
    autocrine_results = []
    cell_type_usage = []
    
    # Analyze each LR pair
    for _, row in lr_pairs.iterrows():
        ligand = row['ligand']
        receptor = row['receptor']
        pathway = row['pathway']
        
        # Skip if ligand or receptor not in expression data
        if ligand not in ligand_df.index or receptor not in receptor_df.index:
            continue
        
        # Check each cell type for autocrine signaling
        autocrine_cell_types = []
        
        for cell_type in cell_types:
            # Check if both ligand and receptor are expressed in this cell type
            if (ligand_df.at[ligand, cell_type] >= expression_threshold and 
                receptor_df.at[receptor, cell_type] >= expression_threshold):
                
                autocrine_cell_types.append(cell_type)
                
                # Add to cell type usage data
                cell_type_usage.append({
                    'cell_type': cell_type,
                    'ligand': ligand,
                    'receptor': receptor,
                    'pathway': pathway,
                    'ligand_expression': ligand_df.at[ligand, cell_type],
                    'receptor_expression': receptor_df.at[receptor, cell_type]
                })
        
        # Add to results if autocrine signaling is present
        if autocrine_cell_types:
            autocrine_results.append({
                'ligand': ligand,
                'receptor': receptor,
                'pathway': pathway,
                'autocrine_cell_types': ','.join(autocrine_cell_types),
                'cell_type_count': len(autocrine_cell_types)
            })
    
    # Convert to dataframes
    if autocrine_results:
        autocrine_df = pd.DataFrame(autocrine_results)
        
        # Sort by frequency (cell type count)
        autocrine_df = autocrine_df.sort_values('cell_type_count', ascending=False)
        
        # Calculate cell type specificity
        specificity_data = []
        
        for _, row in autocrine_df.iterrows():
            cell_types_list = row['autocrine_cell_types'].split(',') if row['autocrine_cell_types'] else []
            
            for cell_type in cell_types_list:
                # Count how many LR pairs use this cell type
                specificity_data.append({
                    'cell_type': cell_type,
                    'ligand': row['ligand'],
                    'receptor': row['receptor'],
                    'pathway': row['pathway']
                })
        
        specificity_df = pd.DataFrame(specificity_data)
        
        # Calculate pathway statistics
        pathway_data = []
        
        for pathway in autocrine_df['pathway'].unique():
            pathway_pairs = autocrine_df[autocrine_df['pathway'] == pathway]
            pathway_data.append({
                'pathway': pathway,
                'pair_count': len(pathway_pairs),
                'avg_cell_type_count': pathway_pairs['cell_type_count'].mean(),
                'max_cell_type_count': pathway_pairs['cell_type_count'].max()
            })
        
        pathway_df = pd.DataFrame(pathway_data)
        pathway_df = pathway_df.sort_values('pair_count', ascending=False)
        
        # Create cell type usage dataframe
        usage_df = pd.DataFrame(cell_type_usage)
        
        # Save results
        autocrine_df.to_csv(frequency_file, index=False)
        specificity_df.to_csv(specificity_file, index=False)
        usage_df.to_csv(usage_file, index=False)
        pathway_df.to_csv(pathway_file, index=False)
        
        # Generate summary report
        with open(summary_file, 'w') as f:
            f.write("# Autocrine Signaling Analysis Summary\n\n")
            f.write(f"Analysis performed on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Overview\n\n")
            f.write(f"- Total ligand-receptor pairs analyzed: {len(lr_pairs)}\n")
            f.write(f"- Pairs involved in autocrine signaling: {len(autocrine_df)}\n")
            f.write(f"- Cell types analyzed: {len(cell_types)}\n")
            f.write(f"- Expression threshold: {expression_threshold}\n\n")
            
            f.write("## Top Autocrine Signaling Pairs\n\n")
            f.write("| Ligand | Receptor | Pathway | Cell Type Count | Cell Types |\n")
            f.write("|--------|----------|---------|----------------|------------|\n")
            
            for _, row in autocrine_df.head(20).iterrows():
                cell_types_truncated = row['autocrine_cell_types']
                if len(cell_types_truncated) > 50:
                    cell_types_truncated = cell_types_truncated[:47] + "..."
                
                f.write(f"| {row['ligand']} | {row['receptor']} | {row['pathway']} | {row['cell_type_count']} | {cell_types_truncated} |\n")
            
            f.write("\n## Cell Type Usage of Autocrine Signaling\n\n")
            
            # Count pairs per cell type
            cell_type_counts = specificity_df['cell_type'].value_counts().reset_index()
            cell_type_counts.columns = ['cell_type', 'pair_count']
            cell_type_counts = cell_type_counts.sort_values('pair_count', ascending=False)
            
            f.write("| Cell Type | Autocrine Pair Count |\n")
            f.write("|-----------|----------------------|\n")
            
            for _, row in cell_type_counts.iterrows():
                f.write(f"| {row['cell_type']} | {row['pair_count']} |\n")
            
            f.write("\n## Pathway Statistics\n\n")
            f.write("| Pathway | Pair Count | Avg Cell Types | Max Cell Types |\n")
            f.write("|---------|------------|----------------|----------------|\n")
            
            for _, row in pathway_df.head(20).iterrows():
                f.write(f"| {row['pathway']} | {row['pair_count']} | {row['avg_cell_type_count']:.2f} | {row['max_cell_type_count']} |\n")
        
        print(f"Autocrine signaling analysis complete. {len(autocrine_df)} pairs found in autocrine signaling.")
        print(f"Results saved to:")
        print(f"  - {frequency_file}")
        print(f"  - {specificity_file}")
        print(f"  - {usage_file}")
        print(f"  - {pathway_file}")
        print(f"  - {summary_file}")
    else:
        print("No autocrine signaling pairs found with the current threshold.")

def analyze_all_communication(expression_threshold=0.1):
    """Analyze all possible communication patterns between cell types"""
    print("Analyzing all communication patterns between cell types...")
    
    # Output file paths
    comm_analysis_file = os.path.join(OUTPUT_DIR, "all_lr_communication_analysis.csv")
    sender_receiver_file = os.path.join(OUTPUT_DIR, "cell_sender_receiver_patterns.csv")
    network_file = os.path.join(OUTPUT_DIR, "cell_cell_communication_network.csv")
    pathway_stats_file = os.path.join(OUTPUT_DIR, "pathway_communication_statistics.csv")
    summary_file = os.path.join(OUTPUT_DIR, "cell_cell_communication_summary.md")
    
    # Load data
    ligand_df, receptor_df = load_expression_data()
    lr_pairs = load_cellchat_lr_pairs()
    
    # Get all cell types
    cell_types = ligand_df.columns.tolist()
    
    # Analysis results
    communication_results = []
    sender_receiver_patterns = []
    network_edges = []
    
    # Track total possible communication paths
    total_possible_paths = len(lr_pairs) * len(cell_types) * len(cell_types)
    active_paths = 0
    
    # Analyze each LR pair for all possible cell-cell communications
    for _, row in lr_pairs.iterrows():
        ligand = row['ligand']
        receptor = row['receptor']
        pathway = row['pathway']
        
        # Skip if ligand or receptor not in expression data
        if ligand not in ligand_df.index or receptor not in receptor_df.index:
            continue
        
        # Track cell types where this pair enables communication
        communication_cell_pairs = []
        
        # Check all possible sender-receiver combinations
        for sender in cell_types:
            # Check if sender expresses the ligand
            if ligand_df.at[ligand, sender] >= expression_threshold:
                
                for receiver in cell_types:
                    # Check if receiver expresses the receptor
                    if receptor_df.at[receptor, receiver] >= expression_threshold:
                        
                        # Record this communication path
                        communication_cell_pairs.append(f"{sender}->{receiver}")
                        
                        # Add to sender-receiver patterns
                        sender_receiver_patterns.append({
                            'ligand': ligand,
                            'receptor': receptor,
                            'pathway': pathway,
                            'sender': sender,
                            'receiver': receiver,
                            'ligand_expression': ligand_df.at[ligand, sender],
                            'receptor_expression': receptor_df.at[receptor, receiver],
                            'is_autocrine': sender == receiver
                        })
                        
                        # Add to network edges
                        network_edges.append({
                            'sender': sender,
                            'receiver': receiver,
                            'pair': f"{ligand}-{receptor}",
                            'pathway': pathway,
                            'is_autocrine': sender == receiver
                        })
                        
                        active_paths += 1
        
        # Add to results if any communication is present
        if communication_cell_pairs:
            communication_results.append({
                'ligand': ligand,
                'receptor': receptor,
                'pathway': pathway,
                'communication_pairs': ','.join(communication_cell_pairs),
                'pair_count': len(communication_cell_pairs),
                'autocrine_count': len([p for p in communication_cell_pairs if p.split('->')[0] == p.split('->')[1]]),
                'paracrine_count': len([p for p in communication_cell_pairs if p.split('->')[0] != p.split('->')[1]])
            })
    
    # Convert to dataframes
    if communication_results:
        comm_df = pd.DataFrame(communication_results)
        
        # Sort by frequency (pair count)
        comm_df = comm_df.sort_values('pair_count', ascending=False)
        
        # Create sender-receiver dataframe
        sr_df = pd.DataFrame(sender_receiver_patterns)
        
        # Create network edges dataframe
        network_df = pd.DataFrame(network_edges)
        
        # Calculate pathway statistics
        pathway_data = []
        
        for pathway in comm_df['pathway'].unique():
            pathway_pairs = comm_df[comm_df['pathway'] == pathway]
            pathway_data.append({
                'pathway': pathway,
                'lr_pair_count': len(pathway_pairs),
                'total_connections': pathway_pairs['pair_count'].sum(),
                'avg_connections_per_pair': pathway_pairs['pair_count'].mean(),
                'autocrine_connections': pathway_pairs['autocrine_count'].sum(),
                'paracrine_connections': pathway_pairs['paracrine_count'].sum()
            })
        
        pathway_df = pd.DataFrame(pathway_data)
        pathway_df = pathway_df.sort_values('total_connections', ascending=False)
        
        # Calculate cell type sender/receiver statistics
        cell_stats = []
        
        for cell_type in cell_types:
            # Count as sender
            sender_count = len(sr_df[sr_df['sender'] == cell_type])
            
            # Count as receiver
            receiver_count = len(sr_df[sr_df['receiver'] == cell_type])
            
            # Count autocrine
            autocrine_count = len(sr_df[(sr_df['sender'] == cell_type) & (sr_df['receiver'] == cell_type)])
            
            cell_stats.append({
                'cell_type': cell_type,
                'sender_count': sender_count,
                'receiver_count': receiver_count,
                'autocrine_count': autocrine_count,
                'sender_receiver_ratio': sender_count / max(receiver_count, 1)
            })
        
        cell_stats_df = pd.DataFrame(cell_stats)
        cell_stats_df = cell_stats_df.sort_values('sender_receiver_ratio', ascending=False)
        
        # Save results
        comm_df.to_csv(comm_analysis_file, index=False)
        sr_df.to_csv(sender_receiver_file, index=False)
        network_df.to_csv(network_file, index=False)
        pathway_df.to_csv(pathway_stats_file, index=False)
        
        # Generate summary report
        with open(summary_file, 'w') as f:
            f.write("# Cell-Cell Communication Analysis Summary\n\n")
            f.write(f"Analysis performed on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Overview\n\n")
            f.write(f"- Total ligand-receptor pairs analyzed: {len(lr_pairs)}\n")
            f.write(f"- Pairs enabling communication: {len(comm_df)}\n")
            f.write(f"- Cell types analyzed: {len(cell_types)}\n")
            f.write(f"- Possible communication paths: {total_possible_paths}\n")
            f.write(f"- Active communication paths: {active_paths}\n")
            f.write(f"- Expression threshold: {expression_threshold}\n\n")
            
            f.write("## Top Communication Ligand-Receptor Pairs\n\n")
            f.write("| Ligand | Receptor | Pathway | Total Connections | Autocrine | Paracrine |\n")
            f.write("|--------|----------|---------|-------------------|-----------|------------|\n")
            
            for _, row in comm_df.head(20).iterrows():
                f.write(f"| {row['ligand']} | {row['receptor']} | {row['pathway']} | {row['pair_count']} | {row['autocrine_count']} | {row['paracrine_count']} |\n")
            
            f.write("\n## Cell Type Communication Patterns\n\n")
            f.write("| Cell Type | As Sender | As Receiver | Autocrine | Sender/Receiver Ratio |\n")
            f.write("|-----------|-----------|-------------|-----------|----------------------|\n")
            
            for _, row in cell_stats_df.iterrows():
                f.write(f"| {row['cell_type']} | {row['sender_count']} | {row['receiver_count']} | {row['autocrine_count']} | {row['sender_receiver_ratio']:.2f} |\n")
            
            f.write("\n## Pathway Communication Statistics\n\n")
            f.write("| Pathway | LR Pairs | Total Connections | Avg Connections/Pair | Autocrine | Paracrine |\n")
            f.write("|---------|----------|-------------------|----------------------|-----------|------------|\n")
            
            for _, row in pathway_df.head(20).iterrows():
                f.write(f"| {row['pathway']} | {row['lr_pair_count']} | {row['total_connections']} | {row['avg_connections_per_pair']:.2f} | {row['autocrine_connections']} | {row['paracrine_connections']} |\n")
        
        print(f"Communication analysis complete. {len(comm_df)} pairs found enabling communication.")
        print(f"Results saved to:")
        print(f"  - {comm_analysis_file}")
        print(f"  - {sender_receiver_file}")
        print(f"  - {network_file}")
        print(f"  - {pathway_stats_file}")
        print(f"  - {summary_file}")
    else:
        print("No communication pairs found with the current threshold.")

def generate_reports(expression_threshold=0.1):
    """Generate comprehensive markdown and visualization reports"""
    print("Generating comprehensive reports...")
    
    # Output file paths
    report_file = os.path.join(OUTPUT_DIR, "comprehensive_lr_analysis_report.md")
    
    # Check if required data files exist
    required_files = [
        os.path.join(OUTPUT_DIR, "cell_type_expression_summary.csv"),
        os.path.join(OUTPUT_DIR, "gene_presence_summary.csv"),
        os.path.join(OUTPUT_DIR, "gene_presence_counts_summary.csv"),
        os.path.join(OUTPUT_DIR, "all_lr_communication_analysis.csv"),
        os.path.join(OUTPUT_DIR, "pathway_communication_statistics.csv")
    ]
    
    missing_files = [f for f in required_files if not os.path.exists(f)]
    if missing_files:
        print(f"Cannot generate reports. Missing files: {missing_files}")
        return
    
    # Load summary data
    gene_counts = pd.read_csv(os.path.join(OUTPUT_DIR, "gene_presence_counts_summary.csv"))
    cell_type_summary = pd.read_csv(os.path.join(OUTPUT_DIR, "cell_type_expression_summary.csv"))
    comm_analysis = pd.read_csv(os.path.join(OUTPUT_DIR, "all_lr_communication_analysis.csv"))
    pathway_stats = pd.read_csv(os.path.join(OUTPUT_DIR, "pathway_communication_statistics.csv"))
    
    # Load sender-receiver patterns if available
    sr_file = os.path.join(OUTPUT_DIR, "cell_sender_receiver_patterns.csv")
    if os.path.exists(sr_file):
        sr_patterns = pd.read_csv(sr_file)
    else:
        sr_patterns = None
    
    # Generate comprehensive report
    with open(report_file, 'w') as f:
        f.write("# Comprehensive Ligand-Receptor Analysis Report\n\n")
        f.write(f"Analysis performed on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Overview\n\n")
        f.write("This report summarizes the analysis of ligand-receptor pairs from the CellChat database in relation to single-cell expression data.\n\n")
        
        f.write("### Gene Coverage\n\n")
        f.write(f"- Total ligands in CellChat: {gene_counts['total_ligands_in_cellchat'].values[0]}\n")
        f.write(f"- Total receptors in CellChat: {gene_counts['total_receptors_in_cellchat'].values[0]}\n")
        f.write(f"- Ligands found in expression data: {gene_counts['ligands_in_expression_data'].values[0]}\n")
        f.write(f"- Receptors found in expression data: {gene_counts['receptors_in_expression_data'].values[0]}\n")
        f.write(f"- Missing ligands: {gene_counts['missing_ligands'].values[0]}\n")
        f.write(f"- Missing receptors: {gene_counts['missing_receptors'].values[0]}\n\n")
        
        f.write("### Expression Summary\n\n")
        f.write(f"- Total cell types analyzed: {len(cell_type_summary)}\n")
        f.write(f"- Expression threshold: {expression_threshold}\n")
        f.write(f"- Ligand-receptor pairs enabling communication: {len(comm_analysis)}\n\n")
        
        f.write("## Cell Type Expression Profiles\n\n")
        f.write("| Cell Type | All Ligands | All Receptors | CellChat Ligands | CellChat Receptors |\n")
        f.write("|-----------|------------|--------------|------------------|--------------------|\n")
        
        for _, row in cell_type_summary.iterrows():
            f.write(f"| {row['cell_type']} | {row['ligand_count']} | {row['receptor_count']} | {row['cellchat_ligand_count']} | {row['cellchat_receptor_count']} |\n")
        
        f.write("\n## Communication Analysis\n\n")
        
        f.write("### Top Communication Pathways\n\n")
        f.write("| Pathway | LR Pairs | Total Connections | Avg Connections/Pair | Autocrine | Paracrine |\n")
        f.write("|---------|----------|-------------------|----------------------|-----------|------------|\n")
        
        for _, row in pathway_stats.head(20).iterrows():
            f.write(f"| {row['pathway']} | {row['lr_pair_count']} | {row['total_connections']} | {row['avg_connections_per_pair']:.2f} | {row['autocrine_connections']} | {row['paracrine_connections']} |\n")
        
        f.write("\n### Top Communication Ligand-Receptor Pairs\n\n")
        f.write("| Ligand | Receptor | Pathway | Total Connections | Autocrine | Paracrine |\n")
        f.write("|--------|----------|---------|-------------------|-----------|------------|\n")
        
        for _, row in comm_analysis.head(20).iterrows():
            f.write(f"| {row['ligand']} | {row['receptor']} | {row['pathway']} | {row['pair_count']} | {row['autocrine_count']} | {row['paracrine_count']} |\n")
        
        f.write("\n## Key Findings\n\n")
        
        # Calculate top autocrine and paracrine pairs
        if len(comm_analysis) > 0:
            top_autocrine = comm_analysis.sort_values('autocrine_count', ascending=False).head(3)
            top_paracrine = comm_analysis.sort_values('paracrine_count', ascending=False).head(3)
            
            f.write("### Top Autocrine Signaling Pairs\n\n")
            for _, row in top_autocrine.iterrows():
                f.write(f"- **{row['ligand']}-{row['receptor']}** ({row['pathway']}): Used in {row['autocrine_count']} cell types\n")
            
            f.write("\n### Top Paracrine Signaling Pairs\n\n")
            for _, row in top_paracrine.iterrows():
                f.write(f"- **{row['ligand']}-{row['receptor']}** ({row['pathway']}): Enables {row['paracrine_count']} cell-cell connections\n")
        
        f.write("\n## Recommendations for Further Analysis\n\n")
        f.write("1. **Validation**: Validate top ligand-receptor pairs with experimental data\n")
        f.write("2. **Functional Analysis**: Investigate biological functions of top pathways\n")
        f.write("3. **Network Analysis**: Perform advanced network analysis of the communication patterns\n")
        f.write("4. **Spatial Context**: Consider spatial relationships between communicating cell types\n")
    
    # Generate visualizations
    try:
        # Create visualization directory
        viz_dir = os.path.join(OUTPUT_DIR, "visualizations")
        os.makedirs(viz_dir, exist_ok=True)
        
        # Plot top pathways
        plt.figure(figsize=(12, 8))
        top_pathways = pathway_stats.sort_values('total_connections', ascending=False).head(10)
        
        sns.barplot(data=top_pathways, y='pathway', x='total_connections', palette='viridis')
        plt.title('Top Communication Pathways by Connection Count')
        plt.xlabel('Number of Connections')
        plt.ylabel('Pathway')
        plt.tight_layout()
        plt.savefig(os.path.join(viz_dir, 'top_pathways.png'), dpi=300)
        plt.close()
        
        # Plot autocrine vs paracrine
        if sr_patterns is not None:
            plt.figure(figsize=(10, 6))
            autocrine_count = len(sr_patterns[sr_patterns['is_autocrine'] == True])
            paracrine_count = len(sr_patterns[sr_patterns['is_autocrine'] == False])
            
            plt.bar(['Autocrine', 'Paracrine'], [autocrine_count, paracrine_count], color=['#3498db', '#e74c3c'])
            plt.title('Autocrine vs Paracrine Signaling Paths')
            plt.ylabel('Number of Communication Paths')
            plt.tight_layout()
            plt.savefig(os.path.join(viz_dir, 'autocrine_vs_paracrine.png'), dpi=300)
            plt.close()
            
            # Cell type communication patterns
            plt.figure(figsize=(12, 8))
            
            # Count occurrences of each cell type as sender and receiver
            sender_counts = sr_patterns['sender'].value_counts().reset_index()
            sender_counts.columns = ['cell_type', 'as_sender']
            
            receiver_counts = sr_patterns['receiver'].value_counts().reset_index()
            receiver_counts.columns = ['cell_type', 'as_receiver']
            
            # Merge the dataframes
            cell_counts = pd.merge(sender_counts, receiver_counts, on='cell_type', how='outer').fillna(0)
            cell_counts = cell_counts.sort_values('as_sender', ascending=False)
            
            # Plot
            x = np.arange(len(cell_counts))
            width = 0.35
            
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.bar(x - width/2, cell_counts['as_sender'], width, label='As Sender')
            ax.bar(x + width/2, cell_counts['as_receiver'], width, label='As Receiver')
            
            ax.set_title('Cell Type Communication Patterns')
            ax.set_ylabel('Communication Path Count')
            ax.set_xticks(x)
            ax.set_xticklabels(cell_counts['cell_type'], rotation=45, ha='right')
            ax.legend()
            
            plt.tight_layout()
            plt.savefig(os.path.join(viz_dir, 'cell_type_patterns.png'), dpi=300)
            plt.close()
        
        print(f"Comprehensive report saved to {report_file}")
        print(f"Visualizations saved to {viz_dir}")
    except Exception as e:
        print(f"Error generating visualizations: {e}")
        print(f"Comprehensive report saved to {report_file}")

def clean_output_files():
    """Clean output files before running analysis"""
    output_files = [
        # List specific output files to clean
        "autocrine_lr_pair_frequency_ranking.csv",
        "autocrine_lr_pair_cell_type_specificity.csv",
        "autocrine_cell_type_lr_pair_usage.csv",
        "autocrine_pathway_statistics.csv",
        "autocrine_ligand_receptor_summary.md",
        "all_lr_communication_analysis.csv",
        "cell_sender_receiver_patterns.csv",
        "cell_cell_communication_network.csv",
        "pathway_communication_statistics.csv",
        "cell_cell_communication_summary.md",
        "cell_type_expression_summary.csv",
        "gene_presence_summary.csv",
        "gene_presence_counts_summary.csv"
    ]
    
    for file in output_files:
        file_path = os.path.join(OUTPUT_DIR, file)
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"Removed: {file_path}")

def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(description="Unified ligand-receptor pair analysis")
    parser.add_argument("--steps", nargs="+", choices=["all", "extract", "summaries", "autocrine", "communication", "reports"], 
                        default=["all"], help="Analysis steps to run")
    parser.add_argument("--clean", action="store_true", help="Clean output directory before running")
    parser.add_argument("--threshold", type=float, default=0.1, help="Expression threshold")
    args = parser.parse_args()
    
    # Get threshold from arguments
    expression_threshold = args.threshold
    
    # Clean output if requested
    if args.clean:
        clean_output_files()
    
    # Run requested steps
    steps = args.steps
    run_all = "all" in steps
    
    if run_all or "extract" in steps:
        extract_cellchat_data()
    
    if run_all or "summaries" in steps:
        generate_expression_summaries(expression_threshold)
    
    if run_all or "autocrine" in steps:
        analyze_autocrine_signaling(expression_threshold)
    
    if run_all or "communication" in steps:
        analyze_all_communication(expression_threshold)
    
    if run_all or "reports" in steps:
        generate_reports(expression_threshold)
    
    print("Analysis complete. Results saved to output directory.")

if __name__ == "__main__":
    main() 