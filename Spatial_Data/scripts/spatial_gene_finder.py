import rdflib
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, RDFS, XSD, OWL
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns
import networkx as nx
from collections import Counter, defaultdict
import argparse
import math
from matplotlib.patches import Circle
import matplotlib.colors as mcolors
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Define namespaces
ST = Namespace("http://example.org/spatial-transcriptomics#")

class SpatialGeneFinderApp:
    def __init__(self, ttl_file):
        """Initialize the Spatial Gene Finder application"""
        self.ttl_file = ttl_file
        self.graph = self.load_turtle_file(ttl_file)
        self.gene_list = self.get_gene_list()
        self.gene_coordinates = self.get_gene_coordinates()
        self.setup_ui()
    
    def load_turtle_file(self, ttl_file):
        """Load a TURTLE file into an RDF graph"""
        print(f"Loading TURTLE file: {ttl_file}")
        g = Graph()
        g.parse(ttl_file, format="turtle")
        print(f"Graph loaded with {len(g)} triples")
        return g
    
    def get_gene_list(self):
        """Get a list of all genes in the ontology"""
        query = """
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT DISTINCT ?geneID
        WHERE {
            ?gene rdf:type :Gene .
            
            # Extract gene ID from URI
            BIND(REPLACE(STR(?gene), "^.*Gene_([^_]+)_.*$", "$1") AS ?geneID)
        }
        ORDER BY ?geneID
        """
        
        results = self.graph.query(query)
        return [str(row.geneID) for row in results]
    
    def get_gene_coordinates(self):
        """Get coordinates for all genes"""
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
            x = int(row.x)
            y = int(row.y)
            expr = int(row.expressionLevel)
            
            if gene_id not in gene_coords:
                gene_coords[gene_id] = []
            
            gene_coords[gene_id].append((x, y, expr))
        
        return gene_coords
    
    def find_genes_within_distance(self, target_gene, distance_threshold):
        """Find all genes within a specific distance of the target gene"""
        if target_gene not in self.gene_coordinates:
            return pd.DataFrame()
        
        # Get target gene coordinates
        target_coords = self.gene_coordinates[target_gene]
        
        # Prepare results
        results = []
        
        # Check distance to all other genes
        for gene_id, coords_list in self.gene_coordinates.items():
            if gene_id == target_gene:
                continue
            
            for tx, ty, texpr in target_coords:
                for ox, oy, oexpr in coords_list:
                    # Calculate Euclidean distance
                    distance = math.sqrt((tx - ox)**2 + (ty - oy)**2)
                    
                    if distance <= distance_threshold:
                        results.append({
                            'target_gene': target_gene,
                            'nearby_gene': gene_id,
                            'target_x': tx,
                            'target_y': ty,
                            'nearby_x': ox,
                            'nearby_y': oy,
                            'distance': distance,
                            'target_expression': texpr,
                            'nearby_expression': oexpr
                        })
        
        # Convert to DataFrame
        df = pd.DataFrame(results)
        
        if not df.empty:
            # Sort by distance
            df = df.sort_values('distance')
            
            # Get unique gene pairs with minimum distance
            df = df.loc[df.groupby(['target_gene', 'nearby_gene'])['distance'].idxmin()]
        
        return df
    
    def find_co_expressed_genes(self, target_gene):
        """Find genes that are co-expressed with the target gene"""
        query = f"""
        PREFIX : <http://example.org/spatial-transcriptomics#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        
        SELECT ?gene2ID (COUNT(*) AS ?coExpressionCount)
        WHERE {{
            ?gene1 rdf:type :Gene ;
                   :expressedAt ?point1 .
            
            ?gene2 rdf:type :Gene ;
                   :expressedAt ?point2 .
            
            ?point1 :locatedNear ?point2 .
            
            # Extract gene IDs from URIs
            BIND(REPLACE(STR(?gene1), "^.*Gene_([^_]+)_.*$", "$1") AS ?gene1ID)
            BIND(REPLACE(STR(?gene2), "^.*Gene_([^_]+)_.*$", "$1") AS ?gene2ID)
            
            # Filter for target gene
            FILTER(?gene1ID = "{target_gene}" && ?gene1 != ?gene2)
        }}
        GROUP BY ?gene2ID
        ORDER BY DESC(?coExpressionCount)
        """
        
        results = self.graph.query(query)
        
        # Convert to DataFrame
        data = []
        for row in results:
            data.append({
                'gene': str(row.gene2ID),
                'co_expression_count': int(row.coExpressionCount)
            })
        
        return pd.DataFrame(data)
    
    def visualize_gene_spatial_relationships(self, target_gene, nearby_genes_df):
        """Visualize the spatial relationships between genes"""
        if target_gene not in self.gene_coordinates or nearby_genes_df.empty:
            return None
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot target gene
        for x, y, expr in self.gene_coordinates[target_gene]:
            ax.scatter(x, y, color='red', s=100, label=f'{target_gene} (Target)')
            ax.add_patch(Circle((x, y), radius=200, fill=False, color='red', linestyle='--', alpha=0.5))
        
        # Plot nearby genes
        plotted_genes = set()
        colors = list(mcolors.TABLEAU_COLORS)
        color_idx = 0
        
        for _, row in nearby_genes_df.iterrows():
            gene = row['nearby_gene']
            
            if gene not in plotted_genes:
                for x, y, expr in self.gene_coordinates[gene]:
                    ax.scatter(x, y, color=colors[color_idx % len(colors)], s=50, label=gene)
                
                plotted_genes.add(gene)
                color_idx += 1
        
        # Draw connections
        for _, row in nearby_genes_df.iterrows():
            ax.plot([row['target_x'], row['nearby_x']], 
                    [row['target_y'], row['nearby_y']], 
                    'k-', alpha=0.3)
            
            # Add distance label
            mid_x = (row['target_x'] + row['nearby_x']) / 2
            mid_y = (row['target_y'] + row['nearby_y']) / 2
            ax.text(mid_x, mid_y, f"{row['distance']:.1f}", fontsize=8)
        
        # Set labels and title
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title(f'Spatial Relationships of {target_gene} with Nearby Genes')
        
        # Add legend
        ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1))
        
        plt.tight_layout()
        return fig
    
    def visualize_co_expression_network(self, target_gene, co_expression_df):
        """Visualize the co-expression network for a target gene"""
        if co_expression_df.empty:
            return None
        
        # Create graph
        G = nx.Graph()
        
        # Add target gene
        G.add_node(target_gene, size=500, color='red')
        
        # Add co-expressed genes (limit to top 15)
        top_genes = co_expression_df.head(15)
        
        for _, row in top_genes.iterrows():
            gene = row['gene']
            count = row['co_expression_count']
            
            # Add node
            G.add_node(gene, size=300, color='lightblue')
            
            # Add edge
            G.add_edge(target_gene, gene, weight=count)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Set up layout
        pos = nx.spring_layout(G, seed=42)
        
        # Get node attributes
        node_sizes = [G.nodes[node].get('size', 300) for node in G.nodes()]
        node_colors = [G.nodes[node].get('color', 'lightblue') for node in G.nodes()]
        
        # Get edge weights for line thickness
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        
        # Normalize weights for visualization
        max_weight = max(edge_weights) if edge_weights else 1
        normalized_weights = [1 + 5 * (w / max_weight) for w in edge_weights]
        
        # Draw the network
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=0.8, ax=ax)
        nx.draw_networkx_edges(G, pos, width=normalized_weights, alpha=0.5, edge_color='gray', ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif', ax=ax)
        
        # Add edge labels (co-expression counts)
        edge_labels = {(u, v): G[u][v]['weight'] for u, v in G.edges()}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8, ax=ax)
        
        ax.set_title(f'Co-expression Network for {target_gene}')
        ax.axis('off')
        
        plt.tight_layout()
        return fig
    
    def setup_ui(self):
        """Set up the user interface"""
        self.root = tk.Tk()
        self.root.title("Spatial Gene Finder")
        self.root.geometry("1200x800")
        
        # Create main frame
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create control panel
        control_frame = ttk.LabelFrame(main_frame, text="Controls", padding=10)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        
        # Gene selection
        ttk.Label(control_frame, text="Select Gene:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.gene_var = tk.StringVar()
        self.gene_combo = ttk.Combobox(control_frame, textvariable=self.gene_var, values=self.gene_list)
        self.gene_combo.grid(row=0, column=1, sticky=tk.W, pady=5)
        self.gene_combo.set(self.gene_list[0] if self.gene_list else "")
        
        # Distance threshold
        ttk.Label(control_frame, text="Distance Threshold:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.distance_var = tk.StringVar(value="200")
        distance_entry = ttk.Entry(control_frame, textvariable=self.distance_var, width=10)
        distance_entry.grid(row=1, column=1, sticky=tk.W, pady=5)
        
        # Search button
        search_button = ttk.Button(control_frame, text="Find Nearby Genes", command=self.on_search)
        search_button.grid(row=2, column=0, columnspan=2, pady=10)
        
        # Co-expression button
        coexpr_button = ttk.Button(control_frame, text="Find Co-expressed Genes", command=self.on_coexpr_search)
        coexpr_button.grid(row=3, column=0, columnspan=2, pady=10)
        
        # Results frame
        results_frame = ttk.Notebook(main_frame)
        results_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Spatial tab
        self.spatial_frame = ttk.Frame(results_frame)
        results_frame.add(self.spatial_frame, text="Spatial Relationships")
        
        # Co-expression tab
        self.coexpr_frame = ttk.Frame(results_frame)
        results_frame.add(self.coexpr_frame, text="Co-expression Network")
        
        # Results table tab
        self.table_frame = ttk.Frame(results_frame)
        results_frame.add(self.table_frame, text="Results Table")
        
        # Initialize empty canvas for plots
        self.spatial_canvas = None
        self.coexpr_canvas = None
        
        # Initialize empty table
        self.results_table = ttk.Treeview(self.table_frame)
        self.results_table.pack(fill=tk.BOTH, expand=True)
    
    def on_search(self):
        """Handle search button click"""
        target_gene = self.gene_var.get()
        try:
            distance_threshold = float(self.distance_var.get())
        except ValueError:
            distance_threshold = 200
        
        # Find nearby genes
        nearby_genes_df = self.find_genes_within_distance(target_gene, distance_threshold)
        
        # Update results table
        self.update_results_table(nearby_genes_df)
        
        # Update spatial visualization
        self.update_spatial_visualization(target_gene, nearby_genes_df)
    
    def on_coexpr_search(self):
        """Handle co-expression search button click"""
        target_gene = self.gene_var.get()
        
        # Find co-expressed genes
        co_expression_df = self.find_co_expressed_genes(target_gene)
        
        # Update co-expression visualization
        self.update_coexpr_visualization(target_gene, co_expression_df)
    
    def update_results_table(self, df):
        """Update the results table with new data"""
        # Clear existing table
        for item in self.results_table.get_children():
            self.results_table.delete(item)
        
        # Configure columns
        if not df.empty:
            columns = ['target_gene', 'nearby_gene', 'distance', 'target_expression', 'nearby_expression']
            self.results_table['columns'] = columns
            
            # Format columns
            self.results_table.column('#0', width=0, stretch=tk.NO)
            for col in columns:
                self.results_table.column(col, anchor=tk.CENTER, width=100)
                self.results_table.heading(col, text=col.replace('_', ' ').title())
            
            # Add data
            for idx, row in df.iterrows():
                values = [row[col] for col in columns]
                self.results_table.insert('', tk.END, values=values)
    
    def update_spatial_visualization(self, target_gene, df):
        """Update the spatial visualization"""
        # Clear existing canvas
        if self.spatial_canvas:
            self.spatial_canvas.get_tk_widget().destroy()
        
        # Create new visualization
        fig = self.visualize_gene_spatial_relationships(target_gene, df)
        
        if fig:
            # Create canvas
            self.spatial_canvas = FigureCanvasTkAgg(fig, master=self.spatial_frame)
            self.spatial_canvas.draw()
            self.spatial_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def update_coexpr_visualization(self, target_gene, df):
        """Update the co-expression visualization"""
        # Clear existing canvas
        if self.coexpr_canvas:
            self.coexpr_canvas.get_tk_widget().destroy()
        
        # Create new visualization
        fig = self.visualize_co_expression_network(target_gene, df)
        
        if fig:
            # Create canvas
            self.coexpr_canvas = FigureCanvasTkAgg(fig, master=self.coexpr_frame)
            self.coexpr_canvas.draw()
            self.coexpr_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def run(self):
        """Run the application"""
        self.root.mainloop()

def main():
    parser = argparse.ArgumentParser(description='Spatial Gene Finder')
    parser.add_argument('--ttl', default="../data/spatial_transcriptomics_advanced.ttl",
                        help='Path to TURTLE file')
    args = parser.parse_args()
    
    app = SpatialGeneFinderApp(args.ttl)
    app.run()

if __name__ == "__main__":
    main() 