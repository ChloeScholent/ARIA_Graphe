import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from extract_drug_drug_effect import diff_drug_pair 

# ==============================
# 1️⃣ Load Data
# ==============================

# Load drug–protein interactions (no header)
dpi = pd.read_csv(
    "csv/different_drug_protein.csv",
    header=None,
    names=["drug", "protein"]
)

# Load protein–protein interactions (with header: "Gene 1", "Gene 2")
ppi = pd.read_csv("csv/bio-decagon-ppi.csv")  # headers are "Gene 1" and "Gene 2"

# List of sets for drug–drug side effects (side effect: "cancer")
drug_side_effects = diff_drug_pair

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

# ==============================
# 3️⃣ Remove proteins linked to only one drug
# ==============================

# Count how many unique drugs interact with each protein
protein_counts = dpi_filtered.groupby("protein")["drug"].nunique()

# Keep only proteins that have interactions with ≥ 2 drugs
multi_drug_proteins = protein_counts[protein_counts >= 2].index

# Filter the dataframe
dpi_filtered = dpi_filtered[dpi_filtered["protein"].isin(multi_drug_proteins)]

# ==============================
# 4️⃣ Filter Protein–Protein Interactions
# ==============================

# Keep only protein–protein edges where both proteins are already in our selected set
selected_proteins = set(dpi_filtered["protein"])

ppi_filtered = ppi[
    ppi["Gene 1"].isin(selected_proteins) & ppi["Gene 2"].isin(selected_proteins)
]

print(f"Original drug–drug pairs: {len(drug_side_effects)}")
print(f"Valid pairs (both drugs in CSV): {len(filtered_pairs)}")
print(f"Drugs with at least one side effect: {len(drugs_with_side_effects)}")
print(f"Proteins with ≥ 2 drug interactions: {len(multi_drug_proteins)}")
print(f"Filtered protein–protein interactions: {len(ppi_filtered)}")

# ==============================
# 5️⃣ Build Graph
# ==============================

G = nx.Graph()

# Add drug nodes
for d in set(dpi_filtered["drug"]):
    G.add_node(d, type="drug")

# Add protein nodes
for p in multi_drug_proteins:
    G.add_node(p, type="protein")

# Add drug–protein edges
for _, row in dpi_filtered.iterrows():
    G.add_edge(row["drug"], row["protein"], type="drug-protein")

# Add drug–drug side-effect edges
for pair in filtered_pairs:
    d1, d2 = tuple(pair)
    if d1 in G.nodes and d2 in G.nodes:
        G.add_edge(d1, d2, type="side-effect", effect="cancer")

# Add protein–protein edges
for _, row in ppi_filtered.iterrows():
    G.add_edge(row["Gene 1"], row["Gene 2"], type="protein-protein")

# ==============================
# 6️⃣ Visualization (Natural Layout)
# ==============================

# Use spring layout for natural structure
pos = nx.spring_layout(G, seed=42)

drug_nodes = [n for n, attr in G.nodes(data=True) if attr["type"] == "drug"]
protein_nodes = [n for n, attr in G.nodes(data=True) if attr["type"] == "protein"]

# Draw nodes
plt.figure(figsize=(12, 7))
nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes, node_color="lightcoral", node_shape="o", label="Drugs")
nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color="lightblue", node_shape="s", label="Proteins")

# Draw edges by type
side_effect_edges = [(u, v) for u, v, e in G.edges(data=True) if e["type"] == "side-effect"]
dpi_edges = [(u, v) for u, v, e in G.edges(data=True) if e["type"] == "drug-protein"]
ppi_edges = [(u, v) for u, v, e in G.edges(data=True) if e["type"] == "protein-protein"]

nx.draw_networkx_edges(G, pos, edgelist=side_effect_edges, edge_color="red", width=2, label="Side effect: cancer")
nx.draw_networkx_edges(G, pos, edgelist=dpi_edges, edge_color="gray", style="dashed", label="Drug–Protein interaction")
nx.draw_networkx_edges(G, pos, edgelist=ppi_edges, edge_color="blue", style="dotted", label="Protein–Protein interaction")

# Labels
nx.draw_networkx_labels(G, pos, font_size=8)

# Legend and formatting
plt.legend()
plt.title("Drug–Drug Side Effects (Cancer) with Drug–Protein and Protein–Protein Interactions")
plt.axis("off")
plt.tight_layout()
plt.show()
