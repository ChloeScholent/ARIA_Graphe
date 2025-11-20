import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


drug_protein_df = pd.read_csv("csv/short_bio-decagon-targets-all.csv")
protein_protein_df = pd.read_csv("csv/bio-decagon-ppi.csv")
drug_combo_df = pd.read_csv("csv/bio-decagon-combo.csv")  

proteins_in_drug_links = set(drug_protein_df["STITCH"])
proteins_in_protein_links = set(protein_protein_df["Gene 1"]) | set(protein_protein_df["Gene 2"])
proteins_to_keep = proteins_in_drug_links & proteins_in_protein_links

filtered_drug_protein_df = drug_protein_df[drug_protein_df["Gene"].isin(proteins_to_keep)]
filtered_protein_protein_df = protein_protein_df[
    protein_protein_df["Gene 1"].isin(proteins_to_keep) &
    protein_protein_df["Gene 2"].isin(proteins_to_keep)
]


G = nx.Graph()

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

drug_nodes = {n for n, d in G.nodes(data=True) if d["mode"] == "drug"}

for _, row in drug_combo_df.iterrows():
    d1, d2 = row["STITCH 1"], row["STITCH 2"]
    # Add edge only if both drugs are already in the graph
    if d1 in drug_nodes and d2 in drug_nodes:
        G.add_edge(d1, d2, relation="drug-drug")


drug_nodes = [n for n, d in G.nodes(data=True) if d["mode"] == "drug"]
protein_nodes = [n for n, d in G.nodes(data=True) if d["mode"] == "protein"]


def cloud_positions(center_x, center_y, n_nodes, radius=2.0):
    angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    positions = {}
    for i, angle in enumerate(angles):
        x = center_x + np.cos(angle) * radius * np.random.uniform(0.7, 1.2)
        y = center_y + np.sin(angle) * radius * np.random.uniform(0.7, 1.2)
        positions[i] = (x, y)
    return positions

drug_pos = cloud_positions(0, 2, len(drug_nodes), radius=2)
protein_pos = cloud_positions(0, -2, len(protein_nodes), radius=2.5)

pos = {drug: drug_pos[i] for i, drug in enumerate(drug_nodes)}
pos.update({protein: protein_pos[i] for i, protein in enumerate(protein_nodes)})


plt.figure(figsize=(10, 8))

# Nodes
nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes, node_color="skyblue", node_shape="s", label="Drugs", node_size=600)
nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color="lightgreen", node_shape="o", label="Proteins", node_size=400)

# Edges
drug_protein_edges = [(u, v) for u, v, d in G.edges(data=True) if d["relation"] == "drug-protein"]
protein_protein_edges = [(u, v) for u, v, d in G.edges(data=True) if d["relation"] == "protein-protein"]
drug_drug_edges = [(u, v) for u, v, d in G.edges(data=True) if d["relation"] == "drug-drug"]

nx.draw_networkx_edges(G, pos, edgelist=drug_protein_edges, edge_color="gray", alpha=0.7)
nx.draw_networkx_edges(G, pos, edgelist=protein_protein_edges, edge_color="orange", style="dashed", alpha=0.8)
nx.draw_networkx_edges(G, pos, edgelist=drug_drug_edges, edge_color="red", width=2.0, style="solid")

# Labels
nx.draw_networkx_labels(G, pos, font_size=8)

plt.title("Drug–Protein–Protein–Combination Network", fontsize=12)
plt.axis("off")
plt.legend(scatterpoints=1)
plt.tight_layout()
plt.show()
