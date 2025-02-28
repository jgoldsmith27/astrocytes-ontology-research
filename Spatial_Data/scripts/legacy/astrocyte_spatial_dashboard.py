#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, XSD, OWL
import networkx as nx
from collections import defaultdict
import argparse
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.colors as mcolors
from matplotlib.patches import Circle
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.cm as cm

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")

# Define astrocyte marker genes for different types
# Based on literature and single-cell data
ASTROCYTE_MARKERS = {
    'Protoplasmic': ['GFAP', 'S100B', 'AQP4', 'GJA1', 'SLC1A2', 'SLC1A3'],
    'Fibrous': ['GFAP', 'VIM', 'ALDH1L1', 'CD44', 'CRYAB', 'HOPX'],
    'Reactive': ['GFAP', 'VIM', 'SERPINA3', 'C3', 'CXCL10', 'STAT3', 'LCN2']
}

class AstrocyteSpatialDashboard:
    """
    A dashboard for visualizing and analyzing different types of astrocytes
    in spatial transcriptomics data.
    """
    
    def __init__(self, ttl_file=None, output_dir="../output/astrocyte_analysis"):
        """
        Initialize the dashboard with the path to the TURTLE file and output directory.
        
        Parameters:
        -----------
        ttl_file : str
            Path to the TURTLE file containing the spatial ontology
        output_dir : str
            Path to the directory where analysis results will be saved
        """
        self.ttl_file = ttl_file
        self.output_dir = output_dir
        self.graph = None
        self.gene_coordinates = {}
        self.astrocyte_data = {}
        self.current_view = "spatial"
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize UI components
        self.root = None
        self.canvas = None
        self.fig = None
        self.ax = None
        
        if ttl_file:
            self.load_data(ttl_file)
    
    def load_data(self, ttl_file):
        """
        Load the ontology data from a TURTLE file.
        
        Parameters:
        -----------
        ttl_file : str
            Path to the TURTLE file containing the spatial ontology
        """
        print(f"Loading TURTLE file: {ttl_file}")
        self.graph = Graph()
        self.graph.parse(ttl_file, format="turtle")
        print(f"Graph loaded with {len(self.graph)} triples")
        
        # Extract gene coordinates
        self.gene_coordinates = self.get_gene_coordinates()
        
        # Identify astrocyte types
        self.identify_astrocyte_types()
    
    def get_gene_coordinates(self):
        """
        Extract coordinates for all genes from the ontology.
        
        Returns:
        --------
        dict
            Dictionary mapping gene IDs to lists of (x, y, expression_level) tuples
        """
        if not self.graph:
            return {}
        
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?geneID ?x ?y ?expressionLevel
        WHERE {
            ?gene rdf:type :Gene ;
                  :expressedAt ?point ;
                  :hasExpressionLevel ?expressionLevel .
            
            ?point :hasXCoordinate ?x ;
                   :hasYCoordinate ?y .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        """
        
        results = self.graph.query(query)
        
        gene_coords = {}
        for row in results:
            gene_id = str(row.geneID)
            x = float(row.x)
            y = float(row.y)
            expr = float(row.expressionLevel)
            
            if gene_id not in gene_coords:
                gene_coords[gene_id] = []
            
            gene_coords[gene_id].append((x, y, expr))
        
        return gene_coords
    
    def identify_astrocyte_types(self):
        """
        Identify different types of astrocytes based on marker gene expression.
        """
        if not self.gene_coordinates:
            return
        
        # Create a spatial grid to aggregate gene expression
        grid_size = 50  # Size of grid cells
        grid_data = defaultdict(lambda: defaultdict(float))
        
        # Collect all spatial points
        all_points = []
        
        # Aggregate gene expression in grid cells
        for gene, coords_list in self.gene_coordinates.items():
            for x, y, expr in coords_list:
                # Discretize coordinates to grid cells
                grid_x = int(x // grid_size)
                grid_y = int(y // grid_size)
                
                # Store expression value in grid
                grid_data[(grid_x, grid_y)][gene] = max(grid_data[(grid_x, grid_y)][gene], expr)
                
                # Collect point for later clustering
                all_points.append((x, y))
        
        # Convert to DataFrame for easier manipulation
        grid_df_data = []
        for (grid_x, grid_y), gene_expr in grid_data.items():
            # Calculate center of grid cell
            center_x = grid_x * grid_size + grid_size / 2
            center_y = grid_y * grid_size + grid_size / 2
            
            # Add row for each grid cell
            row_data = {
                'grid_x': grid_x,
                'grid_y': grid_y,
                'center_x': center_x,
                'center_y': center_y
            }
            
            # Add expression values for each gene
            for gene, expr in gene_expr.items():
                row_data[gene] = expr
            
            grid_df_data.append(row_data)
        
        grid_df = pd.DataFrame(grid_df_data)
        
        # Calculate scores for each astrocyte type
        for cell_type, markers in ASTROCYTE_MARKERS.items():
            # Initialize score column
            grid_df[f'{cell_type}_score'] = 0
            
            # Sum expression of marker genes (if present in data)
            for marker in markers:
                if marker in grid_df.columns:
                    grid_df[f'{cell_type}_score'] += grid_df[marker].fillna(0)
            
            # Normalize by number of markers found
            markers_found = sum(1 for marker in markers if marker in grid_df.columns)
            if markers_found > 0:
                grid_df[f'{cell_type}_score'] /= markers_found
        
        # Determine dominant astrocyte type for each grid cell
        type_columns = [f'{cell_type}_score' for cell_type in ASTROCYTE_MARKERS.keys()]
        grid_df['dominant_type'] = grid_df[type_columns].idxmax(axis=1)
        grid_df['dominant_type'] = grid_df['dominant_type'].str.replace('_score', '')
        
        # Store results
        self.astrocyte_data = grid_df
        
        print(f"Identified {len(grid_df)} grid cells with astrocyte expression")
        for cell_type in ASTROCYTE_MARKERS.keys():
            count = sum(grid_df['dominant_type'] == cell_type)
            print(f"  - {cell_type}: {count} cells")
    
    def visualize_astrocyte_types(self):
        """
        Visualize the spatial distribution of different astrocyte types.
        
        Returns:
        --------
        matplotlib.figure.Figure
            The figure containing the visualization
        """
        if self.astrocyte_data.empty:
            return None
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Define colors for each astrocyte type
        colors = {
            'Protoplasmic': 'blue',
            'Fibrous': 'red',
            'Reactive': 'green'
        }
        
        # Plot each astrocyte type
        for cell_type, color in colors.items():
            # Filter data for this cell type
            type_data = self.astrocyte_data[self.astrocyte_data['dominant_type'] == cell_type]
            
            if not type_data.empty:
                # Plot points
                ax.scatter(
                    type_data['center_x'], 
                    type_data['center_y'],
                    c=color,
                    s=type_data[f'{cell_type}_score'] * 20,  # Size based on score
                    alpha=0.7,
                    label=f'{cell_type} Astrocytes'
                )
        
        # Add legend and labels
        ax.legend(fontsize=12)
        ax.set_xlabel('X Coordinate', fontsize=14)
        ax.set_ylabel('Y Coordinate', fontsize=14)
        ax.set_title('Spatial Distribution of Astrocyte Types', fontsize=16)
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Tight layout
        plt.tight_layout()
        
        return fig
    
    def visualize_marker_expression(self, marker_gene):
        """
        Visualize the expression of a specific marker gene.
        
        Parameters:
        -----------
        marker_gene : str
            The marker gene to visualize
            
        Returns:
        --------
        matplotlib.figure.Figure
            The figure containing the visualization
        """
        if marker_gene not in self.gene_coordinates:
            return None
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Get coordinates for the marker gene
        coords = self.gene_coordinates[marker_gene]
        
        # Extract x, y, and expression values
        x_vals = [x for x, y, expr in coords]
        y_vals = [y for x, y, expr in coords]
        expr_vals = [expr for x, y, expr in coords]
        
        # Create scatter plot with expression levels as colors
        scatter = ax.scatter(
            x_vals, 
            y_vals, 
            c=expr_vals, 
            cmap='viridis', 
            s=50, 
            alpha=0.8
        )
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(f'{marker_gene} Expression Level', fontsize=12)
        
        # Add labels and title
        ax.set_xlabel('X Coordinate', fontsize=14)
        ax.set_ylabel('Y Coordinate', fontsize=14)
        ax.set_title(f'Spatial Expression of {marker_gene}', fontsize=16)
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Tight layout
        plt.tight_layout()
        
        return fig
    
    def setup_ui(self):
        """
        Set up the user interface for the dashboard.
        """
        # Create root window
        self.root = tk.Tk()
        self.root.title("Astrocyte Spatial Dashboard")
        self.root.geometry("1200x800")
        
        # Create main frame
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create left panel for controls
        left_panel = ttk.Frame(main_frame, width=300)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # Create right panel for visualization
        right_panel = ttk.Frame(main_frame)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # Add controls to left panel
        ttk.Label(left_panel, text="Astrocyte Spatial Dashboard", font=("Arial", 14, "bold")).pack(pady=(0, 20))
        
        # Add file loading section
        file_frame = ttk.LabelFrame(left_panel, text="Data Loading")
        file_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Button(file_frame, text="Load Ontology File", command=self.load_ontology_file).pack(fill=tk.X, padx=5, pady=5)
        
        # Add visualization controls
        viz_frame = ttk.LabelFrame(left_panel, text="Visualization")
        viz_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Button(viz_frame, text="Show Astrocyte Types", command=self.show_astrocyte_types).pack(fill=tk.X, padx=5, pady=5)
        
        # Add marker gene selection
        marker_frame = ttk.LabelFrame(left_panel, text="Marker Gene Visualization")
        marker_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Flatten the marker genes list
        all_markers = sorted(set(gene for genes in ASTROCYTE_MARKERS.values() for gene in genes))
        
        self.marker_var = tk.StringVar()
        if all_markers:
            self.marker_var.set(all_markers[0])
        
        ttk.Label(marker_frame, text="Select Marker Gene:").pack(anchor=tk.W, padx=5, pady=(5, 0))
        marker_combo = ttk.Combobox(marker_frame, textvariable=self.marker_var, values=all_markers)
        marker_combo.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Button(marker_frame, text="Show Marker Expression", command=self.show_marker_expression).pack(fill=tk.X, padx=5, pady=5)
        
        # Add statistics section
        stats_frame = ttk.LabelFrame(left_panel, text="Statistics")
        stats_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.stats_text = tk.Text(stats_frame, height=10, width=30, wrap=tk.WORD)
        self.stats_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Add export section
        export_frame = ttk.LabelFrame(left_panel, text="Export")
        export_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Button(export_frame, text="Export Current View", command=self.export_current_view).pack(fill=tk.X, padx=5, pady=5)
        ttk.Button(export_frame, text="Export All Data", command=self.export_all_data).pack(fill=tk.X, padx=5, pady=5)
        
        # Add figure canvas to right panel
        self.fig = Figure(figsize=(8, 6))
        self.ax = self.fig.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=right_panel)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Update statistics
        self.update_statistics()
    
    def load_ontology_file(self):
        """
        Open a file dialog to select and load an ontology file.
        """
        file_path = filedialog.askopenfilename(
            title="Select Ontology File",
            filetypes=[("TURTLE Files", "*.ttl"), ("All Files", "*.*")]
        )
        
        if file_path:
            try:
                self.load_data(file_path)
                self.update_statistics()
                messagebox.showinfo("Success", "Ontology file loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load ontology file: {str(e)}")
    
    def show_astrocyte_types(self):
        """
        Display the spatial distribution of astrocyte types.
        """
        if not self.astrocyte_data.empty:
            # Clear current axes
            self.ax.clear()
            
            # Define colors for each astrocyte type
            colors = {
                'Protoplasmic': 'blue',
                'Fibrous': 'red',
                'Reactive': 'green'
            }
            
            # Plot each astrocyte type
            for cell_type, color in colors.items():
                # Filter data for this cell type
                type_data = self.astrocyte_data[self.astrocyte_data['dominant_type'] == cell_type]
                
                if not type_data.empty:
                    # Plot points
                    self.ax.scatter(
                        type_data['center_x'], 
                        type_data['center_y'],
                        c=color,
                        s=type_data[f'{cell_type}_score'] * 20,  # Size based on score
                        alpha=0.7,
                        label=f'{cell_type} Astrocytes'
                    )
            
            # Add legend and labels
            self.ax.legend(fontsize=12)
            self.ax.set_xlabel('X Coordinate', fontsize=14)
            self.ax.set_ylabel('Y Coordinate', fontsize=14)
            self.ax.set_title('Spatial Distribution of Astrocyte Types', fontsize=16)
            
            # Add grid
            self.ax.grid(True, linestyle='--', alpha=0.7)
            
            # Update canvas
            self.fig.tight_layout()
            self.canvas.draw()
            
            # Update current view
            self.current_view = "astrocyte_types"
        else:
            messagebox.showinfo("Info", "No astrocyte data available. Please load data first.")
    
    def show_marker_expression(self):
        """
        Display the expression of the selected marker gene.
        """
        marker_gene = self.marker_var.get()
        
        if marker_gene and marker_gene in self.gene_coordinates:
            # Clear current axes
            self.ax.clear()
            
            # Get coordinates for the marker gene
            coords = self.gene_coordinates[marker_gene]
            
            # Extract x, y, and expression values
            x_vals = [x for x, y, expr in coords]
            y_vals = [y for x, y, expr in coords]
            expr_vals = [expr for x, y, expr in coords]
            
            # Create scatter plot with expression levels as colors
            scatter = self.ax.scatter(
                x_vals, 
                y_vals, 
                c=expr_vals, 
                cmap='viridis', 
                s=50, 
                alpha=0.8
            )
            
            # Add colorbar
            cbar = self.fig.colorbar(scatter, ax=self.ax)
            cbar.set_label(f'{marker_gene} Expression Level', fontsize=12)
            
            # Add labels and title
            self.ax.set_xlabel('X Coordinate', fontsize=14)
            self.ax.set_ylabel('Y Coordinate', fontsize=14)
            self.ax.set_title(f'Spatial Expression of {marker_gene}', fontsize=16)
            
            # Add grid
            self.ax.grid(True, linestyle='--', alpha=0.7)
            
            # Update canvas
            self.fig.tight_layout()
            self.canvas.draw()
            
            # Update current view
            self.current_view = f"marker_{marker_gene}"
        else:
            messagebox.showinfo("Info", f"No data available for marker gene: {marker_gene}")
    
    def update_statistics(self):
        """
        Update the statistics text widget with current data statistics.
        """
        self.stats_text.delete(1.0, tk.END)
        
        if self.graph:
            self.stats_text.insert(tk.END, f"Ontology loaded: {len(self.graph)} triples\n\n")
        else:
            self.stats_text.insert(tk.END, "No ontology loaded\n\n")
        
        if not self.astrocyte_data.empty:
            self.stats_text.insert(tk.END, "Astrocyte Statistics:\n")
            for cell_type in ASTROCYTE_MARKERS.keys():
                count = sum(self.astrocyte_data['dominant_type'] == cell_type)
                self.stats_text.insert(tk.END, f"  - {cell_type}: {count} cells\n")
            
            # Add average expression of key markers
            self.stats_text.insert(tk.END, "\nAverage Expression:\n")
            for marker in sorted(set(gene for genes in ASTROCYTE_MARKERS.values() for gene in genes)):
                if marker in self.astrocyte_data.columns:
                    avg_expr = self.astrocyte_data[marker].mean()
                    self.stats_text.insert(tk.END, f"  - {marker}: {avg_expr:.2f}\n")
        else:
            self.stats_text.insert(tk.END, "No astrocyte data available")
    
    def export_current_view(self):
        """
        Export the current visualization to a file.
        """
        if self.current_view:
            file_path = filedialog.asksaveasfilename(
                title="Save Visualization",
                defaultextension=".png",
                filetypes=[("PNG Image", "*.png"), ("JPEG Image", "*.jpg"), ("PDF Document", "*.pdf")]
            )
            
            if file_path:
                self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Visualization saved to {file_path}")
        else:
            messagebox.showinfo("Info", "No visualization to export")
    
    def export_all_data(self):
        """
        Export all data to CSV files.
        """
        if not self.astrocyte_data.empty:
            directory = filedialog.askdirectory(title="Select Export Directory")
            
            if directory:
                # Export astrocyte data
                astrocyte_file = os.path.join(directory, "astrocyte_spatial_data.csv")
                self.astrocyte_data.to_csv(astrocyte_file, index=False)
                
                # Export marker gene statistics
                marker_stats = []
                for marker in sorted(set(gene for genes in ASTROCYTE_MARKERS.values() for gene in genes)):
                    if marker in self.astrocyte_data.columns:
                        for cell_type in ASTROCYTE_MARKERS.keys():
                            type_data = self.astrocyte_data[self.astrocyte_data['dominant_type'] == cell_type]
                            if not type_data.empty:
                                marker_stats.append({
                                    'marker_gene': marker,
                                    'astrocyte_type': cell_type,
                                    'mean_expression': type_data[marker].mean(),
                                    'median_expression': type_data[marker].median(),
                                    'max_expression': type_data[marker].max(),
                                    'cell_count': len(type_data)
                                })
                
                marker_stats_df = pd.DataFrame(marker_stats)
                marker_stats_file = os.path.join(directory, "marker_gene_statistics.csv")
                marker_stats_df.to_csv(marker_stats_file, index=False)
                
                messagebox.showinfo("Success", f"Data exported to {directory}")
        else:
            messagebox.showinfo("Info", "No data to export")
    
    def run(self):
        """
        Run the dashboard application.
        """
        self.setup_ui()
        self.root.mainloop()

def main():
    """
    Main function to run the dashboard.
    """
    parser = argparse.ArgumentParser(description="Astrocyte Spatial Dashboard")
    parser.add_argument("--ttl", type=str, help="Path to the TURTLE file containing the spatial ontology")
    parser.add_argument("--output", type=str, default="../output/astrocyte_analysis", help="Path to the output directory")
    
    args = parser.parse_args()
    
    dashboard = AstrocyteSpatialDashboard(ttl_file=args.ttl, output_dir=args.output)
    dashboard.run()

if __name__ == "__main__":
    main() 