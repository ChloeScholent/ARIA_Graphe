import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

summary_df = pd.read_csv("csv/drug_set_stats_summary.csv")


summary_df = summary_df.sort_values("file").reset_index(drop=True)

x = np.arange(len(summary_df)) 
width = 0.35  

fig, ax = plt.subplots(2, 1, figsize=(16, 10))

# Avg shared proteins 
ax[0].bar(x - width/2, summary_df['avg_shared_proteins_side_effect'], width, label='Side Effect', color='red')
ax[0].bar(x + width/2, summary_df['avg_shared_proteins_no_side_effect'], width, label='No Side Effect', color='gray')
ax[0].set_ylabel('Average Shared Proteins')
ax[0].set_xticks(x)
ax[0].set_xticklabels(summary_df['file'], rotation=45, ha='right')
ax[0].set_title('Average Shared Proteins per Drug Set')
ax[0].legend()

# Avg PPI interactions
ax[1].bar(x - width/2, summary_df['avg_ppi_side_effect'], width, label='Side Effect', color='red')
ax[1].bar(x + width/2, summary_df['avg_ppi_no_side_effect'], width, label='No Side Effect', color='gray')
ax[1].set_ylabel('Average PPI Interactions')
ax[1].set_xticks(x)
ax[1].set_xticklabels(summary_df['file'], rotation=45, ha='right')
ax[1].set_title('Average Protein–Protein Interactions per Drug Set')
ax[1].legend()

plt.tight_layout()
plt.show()

plt.figure(figsize=(12,6))
plt.scatter(summary_df['avg_shared_proteins_side_effect'], summary_df['avg_ppi_side_effect'],
            color='red', label='Side Effect', s=80)
plt.scatter(summary_df['avg_shared_proteins_no_side_effect'], summary_df['avg_ppi_no_side_effect'],
            color='gray', label='No Side Effect', s=80)
plt.xlabel('Average Shared Proteins')
plt.ylabel('Average PPI Interactions')
plt.title('Shared Proteins vs PPI Interactions per Drug Set')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


metrics = ['avg_shared_proteins_side_effect', 'avg_ppi_side_effect',
           'avg_shared_proteins_no_side_effect', 'avg_ppi_no_side_effect']

plt.figure(figsize=(16,6))
sns.heatmap(summary_df[metrics].T, annot=True, fmt=".2f", cmap="coolwarm", cbar=True,
            yticklabels=metrics, xticklabels=summary_df['file'])
plt.title("Heatmap of Drug Set Metrics")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
