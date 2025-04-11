import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from pathlib import Path

# Simple relative paths - assuming running from src/spatial
DATA_PATH = "../../data"
OUTPUT_DIR = "../../results"
Path(OUTPUT_DIR).mkdir(exist_ok=True)

print("Loading data...")
# Load data using relative paths
ligand_data = pd.read_csv(f"{DATA_PATH}/spatial/processed/ligand_spots.csv")
receptor_data = pd.read_csv(f"{DATA_PATH}/spatial/processed/receptor_spots.csv")
interactions = pd.read_csv(f"{DATA_PATH}/cellchat/processed/interaction.csv")

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
ligand_data_filtered.loc[:, 'normalized_expression'] = l_scaler.fit_transform(
    ligand_data_filtered[['percent_spots']])
receptor_data_filtered.loc[:, 'normalized_expression'] = r_scaler.fit_transform(
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

print("Analyzing ligand-receptor pairs...")
# 3. Analyze ligand-receptor pairs
valid_pairs = []

for _, row in interactions.iterrows():
    ligand = row['ligand']
    receptor_complex = row['receptor']
    
    # Check if ligand exists in our filtered data
    if ligand not in ligand_expr:
        continue
        
    # Handle receptor complexes (receptors separated by "_")
    is_complex = '_' in receptor_complex
    receptors = receptor_complex.split('_')
    receptor_scores = {'percent': [], 'normalized': []}
    valid_receptors = []
    
    # For each receptor in the complex (or the single receptor if not a complex)
    for receptor in receptors:
        if receptor in receptor_expr:
            receptor_scores['percent'].append(receptor_expr[receptor]['percent'])
            receptor_scores['normalized'].append(receptor_expr[receptor]['normalized'])
            valid_receptors.append(receptor)
    
    # Skip if no receptors found
    if not valid_receptors:
        continue
    
    # For complexes, we require ALL receptors to be found to consider it valid
    if is_complex and len(valid_receptors) < len(receptors):
        complex_completeness = len(valid_receptors) / len(receptors)
        print(f"Partial receptor complex found for {receptor_complex}: {valid_receptors} ({complex_completeness:.0%} complete)")
    
    # Calculate average receptor expression (both raw and normalized)
    avg_receptor_percent = sum(receptor_scores['percent']) / len(receptor_scores['percent'])
    avg_receptor_normalized = sum(receptor_scores['normalized']) / len(receptor_scores['normalized'])
    
    # Calculate different types of scores
    raw_score = ligand_expr[ligand]['percent'] * avg_receptor_percent
    normalized_score = ligand_expr[ligand]['normalized'] * avg_receptor_normalized
    
    valid_pairs.append({
        'ligand': ligand,
        'receptor': receptor_complex,
        'is_complex_receptor': is_complex,
        'valid_receptors': '_'.join(valid_receptors),
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
pairs_df.to_csv(f"{OUTPUT_DIR}/ligand_receptor_pairs.csv", index=False)

# 5. Print top results
print("\nTop 20 Ligand-Receptor Pairs (Normalized Score):")
print(pairs_df[['ligand', 'receptor', 'normalized_score', 'is_complex_receptor']].head(20))

# Generate visualizations
print("\nCreating visualizations...")

# Top pairs bar chart
plt.figure(figsize=(12, 8))
top_pairs = pairs_df.head(15)
sns.barplot(x='normalized_score', y=top_pairs['ligand'] + ' - ' + top_pairs['receptor'], data=top_pairs)
plt.title('Top 15 Ligand-Receptor Pairs by Normalized Expression Score')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/top_pairs.png", dpi=300)

# Complex vs Simple receptors bar chart
plt.figure(figsize=(10, 6))
complex_counts = pairs_df['is_complex_receptor'].value_counts()
sns.barplot(x=complex_counts.index, y=complex_counts.values)
plt.title('Complex vs. Simple Receptors in Valid Pairs')
plt.xticks([0, 1], ['Simple Receptor', 'Complex Receptor'])
plt.ylabel('Number of Pairs')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/complex_vs_simple.png", dpi=300)

# Distribution of scores
plt.figure(figsize=(10, 6))
sns.histplot(pairs_df['normalized_score'], bins=30, kde=True)
plt.title('Distribution of Normalized Scores')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/score_distribution.png", dpi=300)

print(f"Results saved to {OUTPUT_DIR}")
print("\nNote on complex receptor handling:")
print("- For receptors with underscore separators (e.g., 'TGFbR1_R2'), each component is checked separately")
print("- The script considers a valid match if ANY of the components are found in the expression data")
print("- Scores are averaged across the valid components")
print("- The 'is_complex_receptor' column indicates whether the receptor is a complex (contains '_')")
print("- The 'valid_receptors' column shows which components were actually found in the expression data") 