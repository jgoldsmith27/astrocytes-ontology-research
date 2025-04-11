import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from pathlib import Path

# Set up paths
PROJECT_ROOT = Path(".")
OUTPUT_DIR = PROJECT_ROOT / "results"
OUTPUT_DIR.mkdir(exist_ok=True)

print("Loading data...")
# Load data
ligand_data = pd.read_csv("data/spatial/processed/ligand_spots.csv")
receptor_data = pd.read_csv("data/spatial/processed/receptor_spots.csv")
interactions = pd.read_csv("data/cellchat/processed/interaction.csv")

print("Filtering and normalizing expression data...")
# 1. Set expression threshold (filter out noise)
min_expression = 0.001  # 0.1% of spots
ligand_data_filtered = ligand_data[ligand_data['percent_spots'] > min_expression]
receptor_data_filtered = receptor_data[receptor_data['percent_spots'] > min_expression]

print(f"Filtered ligands: {len(ligand_data_filtered)}/{len(ligand_data)} ({len(ligand_data_filtered)/len(ligand_data):.1%})")
print(f"Filtered receptors: {len(receptor_data_filtered)}/{len(receptor_data)} ({len(receptor_data_filtered)/len(receptor_data):.1%})")

# 2. Normalize expression values (0-1 scale)
l_scaler = MinMaxScaler()
r_scaler = MinMaxScaler()
ligand_data_filtered['normalized_expression'] = l_scaler.fit_transform(
    ligand_data_filtered[['percent_spots']])
receptor_data_filtered['normalized_expression'] = r_scaler.fit_transform(
    receptor_data_filtered[['percent_spots']])

# Create lookups with original values and normalized values
ligand_expr = {row['gene']: {
    'percent': row['percent_spots'],
    'normalized': row['normalized_expression']
} for _, row in ligand_data_filtered.iterrows()}

receptor_expr = {row['gene']: {
    'percent': row['percent_spots'], 
    'normalized': row['normalized_expression']
} for _, row in receptor_data_filtered.iterrows()}

print("Analyzing interaction pairs...")
# 3. Analyze interaction pairs
valid_pairs = []

for _, row in interactions.iterrows():
    ligand = row['ligand']
    receptor_complex = row['receptor']
    pathway = row['pathway_name']
    
    # Check if ligand exists in our filtered data
    if ligand not in ligand_expr:
        continue
        
    # Handle receptor complexes
    receptors = receptor_complex.split('_')
    receptor_scores = {'percent': [], 'normalized': []}
    valid_receptors = []
    
    for receptor in receptors:
        if receptor in receptor_expr:
            receptor_scores['percent'].append(receptor_expr[receptor]['percent'])
            receptor_scores['normalized'].append(receptor_expr[receptor]['normalized'])
            valid_receptors.append(receptor)
    
    # Skip if no receptors found
    if not valid_receptors:
        continue
    
    # Calculate average receptor expression (both raw and normalized)
    avg_receptor_percent = sum(receptor_scores['percent']) / len(receptor_scores['percent'])
    avg_receptor_normalized = sum(receptor_scores['normalized']) / len(receptor_scores['normalized'])
    
    # Calculate different types of scores
    raw_score = ligand_expr[ligand]['percent'] * avg_receptor_percent
    normalized_score = ligand_expr[ligand]['normalized'] * avg_receptor_normalized
    
    valid_pairs.append({
        'ligand': ligand,
        'receptor': receptor_complex,
        'valid_receptors': '_'.join(valid_receptors),
        'pathway': pathway,
        'ligand_percent': ligand_expr[ligand]['percent'],
        'receptor_percent': avg_receptor_percent,
        'ligand_normalized': ligand_expr[ligand]['normalized'],
        'receptor_normalized': avg_receptor_normalized,
        'raw_score': raw_score,
        'normalized_score': normalized_score
    })

print(f"Found {len(valid_pairs)} valid ligand-receptor pairs")

# 4. Convert to DataFrame and sort by normalized score
pairs_df = pd.DataFrame(valid_pairs)
pairs_df = pairs_df.sort_values('normalized_score', ascending=False)

# Save full results
pairs_df.to_csv(OUTPUT_DIR / "ligand_receptor_pairs.csv", index=False)

print("Analyzing pathways...")
# 5. Group by pathway and calculate summary statistics
pathway_summary = pairs_df.groupby('pathway').agg({
    'normalized_score': ['mean', 'sum', 'count'],
    'ligand': 'nunique',
    'valid_receptors': 'nunique'
})

# Flatten column names
pathway_summary.columns = ['avg_score', 'total_score', 'pair_count', 'unique_ligands', 'unique_receptors']

# Calculate a weighted score that considers both score and diversity
pathway_summary['weighted_score'] = pathway_summary['total_score'] * np.sqrt(pathway_summary['pair_count'])
pathway_summary = pathway_summary.sort_values('weighted_score', ascending=False)

# Save pathway summary
pathway_summary.to_csv(OUTPUT_DIR / "pathway_summary.csv")

# 6. Print top results
print("\nTop 10 Ligand-Receptor Pairs (Normalized Score):")
print(pairs_df[['ligand', 'receptor', 'pathway', 'normalized_score']].head(10))

print("\nTop 5 Signaling Pathways:")
print(pathway_summary[['weighted_score', 'pair_count', 'unique_ligands', 'unique_receptors']].head(5))

print("\nCreating visualizations...")
# 7. Create visualizations

# Top pairs bar chart
plt.figure(figsize=(12, 8))
top_pairs = pairs_df.head(15)
sns.barplot(x='normalized_score', y=top_pairs['ligand'] + ' - ' + top_pairs['receptor'], data=top_pairs)
plt.title('Top 15 Ligand-Receptor Pairs by Normalized Expression Score')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "top_pairs.png", dpi=300)

# Top pathways bar chart
plt.figure(figsize=(12, 8))
top_pathways = pathway_summary.head(10)
sns.barplot(x='weighted_score', y=top_pathways.index, data=top_pathways)
plt.title('Top 10 Pathways by Weighted Score')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "top_pathways.png", dpi=300)

# Heatmap of top pathways
plt.figure(figsize=(12, 8))
metrics = ['avg_score', 'total_score', 'pair_count', 'unique_ligands', 'unique_receptors']
sns.heatmap(pathway_summary[metrics].head(15), annot=True, cmap='viridis', fmt='.1f')
plt.title('Pathway Metrics Heatmap (Top 15)')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "pathway_metrics_heatmap.png", dpi=300)

# Pathway composition: For the top pathway, show the contributions of different pairs
top_pathway = pathway_summary.index[0]
top_pathway_pairs = pairs_df[pairs_df['pathway'] == top_pathway].sort_values('normalized_score', ascending=False)

plt.figure(figsize=(12, 8))
sns.barplot(x='normalized_score', y=top_pathway_pairs['ligand'] + ' - ' + top_pathway_pairs['receptor'], 
            data=top_pathway_pairs.head(10))
plt.title(f'Top Pairs in {top_pathway} Pathway')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "top_pathway_composition.png", dpi=300)

print(f"Results saved to {OUTPUT_DIR}") 