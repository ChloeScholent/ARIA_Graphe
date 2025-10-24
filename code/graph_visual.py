import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# --- Load CSVs ---
drug_protein_df = pd.read_csv("csv/bio-decagon-targets-all.csv")
protein_protein_df = pd.read_csv("csv/short_protein.csv")

# --- Find sets of proteins ---
proteins_in_drug_links = set(drug_protein_df["Gene"])
proteins_in_protein_links = set(protein_protein_df["Gene 1"]) | set(protein_protein_df["Gene 2"])

# Keep only proteins that appear in BOTH sets
proteins_to_keep = proteins_in_drug_links & proteins_in_protein_links

print(f"Proteins in drug links: {len(proteins_in_drug_links)}")
print(f"Proteins in protein links: {len(proteins_in_protein_links)}")
print(f"Proteins kept (drug + protein): {len(proteins_to_keep)}")

# --- Filter data ---
filtered_drug_protein_df = drug_protein_df[drug_protein_df["Gene"].isin(proteins_to_keep)]
filtered_protein_protein_df = protein_protein_df[
    protein_protein_df["Gene 1"].isin(proteins_to_keep) &
    protein_protein_df["Gene 2"].isin(proteins_to_keep)
]

# --- Build graph ---
G = nx.Graph()

# Add Drug–Protein edges
for _, row in filtered_drug_protein_df.iterrows():
    drug, protein = row["STITCH"], row["Gene"]
    G.add_node(drug, mode="drug")
    G.add_node(protein, mode="protein")
    G.add_edge(drug, protein, relation="drug-protein")

# Add Protein–Protein edges
for _, row in filtered_protein_protein_df.iterrows():
    p1, p2 = row["Gene 1"], row["Gene 2"]
    G.add_node(p1, mode="protein")
    G.add_node(p2, mode="protein")
    G.add_edge(p1, p2, relation="protein-protein")

# --- Separate node sets ---
drug_nodes = [n for n, d in G.nodes(data=True) if d["mode"] == "drug"]
protein_nodes = [n for n, d in G.nodes(data=True) if d["mode"] == "protein"]

# --- Define cloud layout ---
def cloud_positions(center_x, center_y, n_nodes, radius=2.0):
    angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    positions = {}
    for i, angle in enumerate(angles):
        x = center_x + np.cos(angle) * radius * np.random.uniform(0.7, 1.2)
        y = center_y + np.sin(angle) * radius * np.random.uniform(0.7, 1.2)
        positions[i] = (x, y)
    return positions

# Generate cloud coordinates
drug_pos = cloud_positions(0, 2, len(drug_nodes), radius=2)
protein_pos = cloud_positions(0, -2, len(protein_nodes), radius=2.5)

# Merge into a single position dictionary
pos = {drug: drug_pos[i] for i, drug in enumerate(drug_nodes)}
pos.update({protein: protein_pos[i] for i, protein in enumerate(protein_nodes)})

# --- Draw the network ---
plt.figure(figsize=(10, 8))

nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes, node_color="skyblue", node_shape="s", label="Drugs", node_size=600)
nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color="lightgreen", node_shape="o", label="Proteins", node_size=400)

# Edges
drug_protein_edges = [(u, v) for u, v, d in G.edges(data=True) if d["relation"] == "drug-protein"]
protein_protein_edges = [(u, v) for u, v, d in G.edges(data=True) if d["relation"] == "protein-protein"]

nx.draw_networkx_edges(G, pos, edgelist=drug_protein_edges, edge_color="gray", alpha=0.7)
nx.draw_networkx_edges(G, pos, edgelist=protein_protein_edges, edge_color="orange", style="dashed", alpha=0.8)

nx.draw_networkx_labels(G, pos, font_size=8)

plt.title("Filtered Drug–Protein Network (Proteins linked to Drugs AND Proteins)", fontsize=12)
plt.legend(scatterpoints=1)
plt.axis("off")
plt.tight_layout()
plt.show()
