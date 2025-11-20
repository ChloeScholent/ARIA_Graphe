import csv
import networkx as nx
import matplotlib.pyplot as plt


combo_path = "csv/bio-decagon-combo.csv"
drug_protein_path = "csv/fever_hyperT_drug_protein.csv"  
ppi_path = "csv/bio-decagon-ppi.csv"


FEVER_CODE = "C0018621"
HYPERTENSION_CODE = "C0020542"

# Drug-drug pair both side effetcs
pair_to_effects = {}
with open(combo_path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        pair = tuple(sorted([row["STITCH 1"], row["STITCH 2"]]))
        effect = row["Polypharmacy Side Effect"]
        pair_to_effects.setdefault(pair, set()).add(effect)

# Filter for fever + hypertension
selected_pairs = [
    pair for pair, effects in pair_to_effects.items()
    if FEVER_CODE in effects and HYPERTENSION_CODE in effects
]

# Extract all drugs from those pairs
drug_list = sorted(set([d for pair in selected_pairs for d in pair]))

print(f"✅ {len(selected_pairs)} drug pairs with both fever and hypertension.")
print(f"✅ {len(drug_list)} unique drugs involved.")

# rug–protein interactions
drug_proteins = []
with open(drug_protein_path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        drug = row["Drug"] if "Drug" in row else row["col1"]
        protein = str(row["Protein"] if "Protein" in row else row["col2"])
        if drug in drug_list:
            drug_proteins.append((drug, protein))

# protein–protein interactions
ppi_edges = []
with open(ppi_path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        p1 = str(row["Protein 1"]) if "Protein 1" in row else str(list(row.values())[0])
        p2 = str(row["Protein 2"]) if "Protein 2" in row else str(list(row.values())[1])
        ppi_edges.append((p1, p2))

# Build graph
G = nx.Graph()

# Add drug nodes
for d in drug_list:
    G.add_node(d, type="drug")

# Add protein nodes
protein_nodes = {p for _, p in drug_proteins} | {p for pair in ppi_edges for p in pair}
for p in protein_nodes:
    G.add_node(p, type="protein")

# Add edges
for d1, d2 in selected_pairs:
    G.add_edge(d1, d2, type="drug_drug", label="Fever+Hypertension")

for d, p in drug_proteins:
    G.add_edge(d, p, type="drug_protein")

for p1, p2 in ppi_edges:
    G.add_edge(p1, p2, type="ppi")

print(f"✅ Graph built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")



plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, k=0.3, seed=42)

drug_nodes = [n for n, attr in G.nodes(data=True) if attr["type"] == "drug"]
protein_nodes = [n for n, attr in G.nodes(data=True) if attr["type"] == "protein"]

nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes, node_color="skyblue", label="Drugs", node_size=300)
nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color="lightgreen", label="Proteins", node_size=200)


edge_colors = []
for u, v, data in G.edges(data=True):
    if data["type"] == "drug_drug":
        edge_colors.append("red")
    elif data["type"] == "drug_protein":
        edge_colors.append("gray")
    else:
        edge_colors.append("lightgray")

nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=1.2, alpha=0.7)
nx.draw_networkx_labels(G, pos, labels={d: d for d in drug_nodes}, font_size=8)

plt.legend(scatterpoints=1)
plt.title("Drug–Protein Interaction Network (Fever & Hypertension)")
plt.axis("off")
plt.tight_layout()
plt.show()
