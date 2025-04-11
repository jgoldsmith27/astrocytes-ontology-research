import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Load data
print("Loading data...")
module_data = pd.read_csv("../../data/cell_modules/human_cerebral_cortex_condition_control_integrated_zpower_matrix.csv", 
                         header=None, names=['gene', 'ensembl', 'score1', 'score2', 'module'])
pairs_data = pd.read_csv("../../results/ligand_receptor_pairs.csv")

# Get top 15 pairs
top_pairs = pairs_df = pairs_data.sort_values('normalized_score', ascending=False).head(15)

# Create module context analysis
module_context = []
for _, pair in top_pairs.iterrows():
    ligand = pair['ligand']
    receptor = pair['receptor'].split('_')  # Handle complex receptors
    
    # Find ligand module
    ligand_info = module_data[module_data['gene'] == ligand]
    ligand_module = ligand_info['module'].values[0] if not ligand_info.empty else None
    
    # Find receptor modules (might be multiple for complexes)
    receptor_modules = []
    for r in receptor:
        r_info = module_data[module_data['gene'] == r]
        if not r_info.empty:
            receptor_modules.append(r_info['module'].values[0])
    
    module_context.append({
        'ligand': ligand,
        'receptor': pair['receptor'],
        'normalized_score': pair['normalized_score'],
        'ligand_module': ligand_module,
        'receptor_modules': receptor_modules,
        'is_same_module': all(m == ligand_module for m in receptor_modules) if receptor_modules else None,
        'interaction_type': 'Intra-module' if (receptor_modules and all(m == ligand_module for m in receptor_modules)) 
                          else 'Inter-module' if receptor_modules else 'Unknown'
    })

# Convert to DataFrame
context_df = pd.DataFrame(module_context)

# Save results to CSV
context_df.to_csv('../../results/module_context_analysis.csv', index=False)
print("Detailed results saved to module_context_analysis.csv")

# Print results
print("\nModule Context Analysis of Top 15 Ligand-Receptor Pairs:")
print(context_df[['ligand', 'receptor', 'ligand_module', 'receptor_modules', 'interaction_type', 'normalized_score']])

# Visualizations
print("\nCreating visualizations...")

# 1. Module interaction heatmap
plt.figure(figsize=(10, 8))
interactions = []
for _, row in context_df.iterrows():
    if row['receptor_modules']:
        for r_module in row['receptor_modules']:
            interactions.append((row['ligand_module'], r_module, row['normalized_score']))

if interactions:
    interaction_df = pd.DataFrame(interactions, columns=['ligand_module', 'receptor_module', 'score'])
    pivot_df = interaction_df.pivot_table(values='score', 
                                        index='ligand_module', 
                                        columns='receptor_module', 
                                        aggfunc='sum')
    
    # Save interaction matrix to CSV
    pivot_df.to_csv('../../results/module_interaction_matrix.csv')
    print("Module interaction matrix saved to module_interaction_matrix.csv")
    
    sns.heatmap(pivot_df, annot=True, cmap='viridis', fmt='.2f')
    plt.title('Module Interaction Strength')
    plt.xlabel('Receptor Module')
    plt.ylabel('Ligand Module')
else:
    plt.text(0.5, 0.5, 'No module interactions found', ha='center', va='center')

plt.tight_layout()
plt.savefig('../../results/module_interactions.png', dpi=300)

# 2. Interaction type distribution
plt.figure(figsize=(8, 6))
interaction_counts = context_df['interaction_type'].value_counts()
sns.barplot(x=interaction_counts.index, y=interaction_counts.values)
plt.title('Distribution of Interaction Types')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('../../results/interaction_types.png', dpi=300)

# Summary statistics
print("\nSummary Statistics:")
print(f"Total pairs analyzed: {len(context_df)}")
print("\nInteraction type distribution:")
print(context_df['interaction_type'].value_counts())
print("\nMost common module pairs:")
module_pairs = [(row['ligand_module'], m) for _, row in context_df.iterrows() 
                for m in row['receptor_modules'] if row['receptor_modules']]
if module_pairs:
    pair_counts = pd.Series(module_pairs).value_counts().head()
    print(pair_counts)

# Save summary statistics to CSV
summary_stats = {
    'Metric': ['Total Pairs', 'Intra-module Pairs', 'Inter-module Pairs', 'Unknown'],
    'Count': [
        len(context_df),
        sum(context_df['interaction_type'] == 'Intra-module'),
        sum(context_df['interaction_type'] == 'Inter-module'),
        sum(context_df['interaction_type'] == 'Unknown')
    ]
}
pd.DataFrame(summary_stats).to_csv('../../results/module_analysis_summary.csv', index=False)
print("\nSummary statistics saved to module_analysis_summary.csv")

print("\nAll results saved to ../../results/") 