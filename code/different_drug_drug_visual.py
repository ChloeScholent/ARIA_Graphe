import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from drug_drug_effect import diff_drug_pair

# Load drug–protein interactions
dpi = pd.read_csv("csv/different_drug_protein.csv", header=None, names=["drug", "protein"])

# Drug–Drug side effects as list of sets
drug_side_effects = diff_drug_pair

# Initialize graph
G = nx.Graph()

# --- Add nodes ---

# Add drug nodes
drugs = set(dpi['drug']).union({d for pair in drug_side_effects for d in pair})
for d in drugs:
    G.add_node(d, type='drug')

# Add protein nodes
proteins = set(dpi['protein'])
for p in proteins:
    G.add_node(p, type='protein')

# --- Add edges ---

# Add drug–protein edges
for _, row in dpi.iterrows():
    G.add_edge(row['drug'], row['protein'], type='drug-protein')

# Add drug–drug edges (side effects)
for pair in drug_side_effects:
    if len(pair) == 2:  # ensure it's a pair
        d1, d2 = tuple(pair)
        G.add_edge(d1, d2, type='side-effect', effect='cancer')
    else:
        print(f"Skipping invalid pair: {pair}")

# --- Visualization ---

# Separate node types
drug_nodes = [n for n, attr in G.nodes(data=True) if attr['type'] == 'drug']
protein_nodes = [n for n, attr in G.nodes(data=True) if attr['type'] == 'protein']

# Layout
pos = nx.spring_layout(G, seed=42)

# Draw nodes
nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes, node_color='lightcoral', node_shape='o', label='Drugs')
nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color='lightblue', node_shape='s', label='Proteins')

# Draw edges by type
side_effect_edges = [(u, v) for u, v, e in G.edges(data=True) if e['type'] == 'side-effect']
dpi_edges = [(u, v) for u, v, e in G.edges(data=True) if e['type'] == 'drug-protein']

nx.draw_networkx_edges(G, pos, edgelist=side_effect_edges, edge_color='red', style='solid', width=2, label='Presence of side effect')
nx.draw_networkx_edges(G, pos, edgelist=dpi_edges, edge_color='gray', style='dashed', label='Drug–Protein interaction')

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=8)

# Final touches
plt.legend()
plt.title("Drug–Drug Side Effects and Drug–Protein Interactions")
plt.axis("off")
plt.show()
