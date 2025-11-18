import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from extract_drug_drug_effect import drug_pair
from pprint import pprint as print

# ==============================
# 1️⃣ Load Data
# ==============================

# Load the drug–protein interactions CSV (no header)
dpi = pd.read_csv(
    "csv/similar_drug_protein.csv",
    header=None,
    names=["drug", "protein"]  # assign column names
)

# Example: list of sets for drug–drug side effects (side effect: "cancer")
drug_side_effects = drug_pair

# ==============================
# 2️⃣ Filter Drug–Drug Pairs
# ==============================

# Keep only pairs where both drugs exist in the drug–protein file
valid_drugs = set(dpi["drug"])
filtered_pairs = [pair for pair in drug_side_effects if pair.issubset(valid_drugs)]

# Collect drugs involved in at least one valid side-effect pair
drugs_with_side_effects = {d for pair in filtered_pairs for d in pair}

# Filter the drug–protein data to only include these drugs
dpi_filtered = dpi[dpi["drug"].isin(drugs_with_side_effects)]

print(f"Original drug–drug pairs: {len(drug_side_effects)}")
print(f"Valid pairs (both drugs in CSV): {len(filtered_pairs)}")
print(f"Drugs with at least one side effect: {len(drugs_with_side_effects)}")

# ==============================
# 3️⃣ Build Graph
# ==============================

G = nx.Graph()

# Add drug nodes (only those with side effects)
for d in drugs_with_side_effects:
    G.add_node(d, type='drug')

# Add protein nodes connected to those drugs
proteins = set(dpi_filtered["protein"])
for p in proteins:
    G.add_node(p, type='protein')

# Add drug–protein edges
for _, row in dpi_filtered.iterrows():
    G.add_edge(row["drug"], row["protein"], type='drug-protein')

# Add drug–drug side-effect edges
for pair in filtered_pairs:
    d1, d2 = tuple(pair)
    G.add_edge(d1, d2, type='side-effect', effect='cancer')

# ==============================
# 4️⃣ Visualization
# ==============================

# Separate node types
drug_nodes = [n for n, attr in G.nodes(data=True) if attr["type"] == "drug"]
protein_nodes = [n for n, attr in G.nodes(data=True) if attr["type"] == "protein"]

# Layout
pos = nx.spring_layout(G, seed=42)

# Draw nodes
nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes, node_color="lightcoral", node_shape="o", label="Drugs")
nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color="lightblue", node_shape="s", label="Proteins")

# Draw edges
side_effect_edges = [(u, v) for u, v, e in G.edges(data=True) if e["type"] == "side-effect"]
dpi_edges = [(u, v) for u, v, e in G.edges(data=True) if e["type"] == "drug-protein"]

nx.draw_networkx_edges(G, pos, edgelist=side_effect_edges, edge_color="red", width=2, label="Side effect: cancer")
nx.draw_networkx_edges(G, pos, edgelist=dpi_edges, edge_color="gray", style="dashed", label="Drug–Protein interaction")

# Labels
nx.draw_networkx_labels(G, pos, font_size=8)

# Final touches
plt.legend()
plt.title("Drug–Drug Side Effects (Cancer) and Drug–Protein Interactions")
plt.axis("off")
plt.show()
