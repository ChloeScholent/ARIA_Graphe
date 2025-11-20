import networkx as nx
import random
import itertools

num_drugs = 1774
num_proteins = 19081

# Random edge counts
num_drug_drug_edges = random.randint(14256, 14257)
num_drug_protein_edges = random.randint(400, 100000)
num_protein_protein_edges = random.randint(2000, 100000)

print("Edge counts chosen:")
print(f"  Drug–Drug edges: {num_drug_drug_edges}")
print(f"  Drug–Protein edges: {num_drug_protein_edges}")
print(f"  Protein–Protein edges: {num_protein_protein_edges}\n")

# Create graph
G = nx.Graph()

# Add nodes
for i in range(num_drugs):
    G.add_node(f"Drug_{i+1}", type="drug")
for i in range(num_proteins):
    G.add_node(f"Protein_{i+1}", type="protein")

drugs = [n for n, d in G.nodes(data=True) if d["type"] == "drug"]
proteins = [n for n, d in G.nodes(data=True) if d["type"] == "protein"]

def add_random_edges(G, nodes_a, nodes_b, edge_type, num_edges, allow_same_set=False):
    possible_edges = (
        list(itertools.combinations(nodes_a, 2)) if allow_same_set
        else list(itertools.product(nodes_a, nodes_b))
    )
    random.shuffle(possible_edges)
    added = 0
    for edge in possible_edges:
        if added >= num_edges:
            break
        if G.has_edge(*edge):
            continue
        G.add_edge(edge[0], edge[1], type=edge_type)
        added += 1

add_random_edges(G, drugs, drugs, "drug_drug", num_drug_drug_edges, allow_same_set=True)
add_random_edges(G, drugs, proteins, "drug_protein", num_drug_protein_edges)
add_random_edges(G, proteins, proteins, "protein_protein", num_protein_protein_edges, allow_same_set=True)



# Map drug → connected proteins
drug_to_proteins = {
    d: {nbr for nbr in G.neighbors(d) if G.nodes[nbr]["type"] == "protein"}
    for d in drugs
}

linked_drug_pairs = [
    (d1, d2) for d1, d2 in itertools.combinations(drugs, 2)
    if G.has_edge(d1, d2) and G.edges[d1, d2]["type"] == "drug_drug"
]

# Average shared proteins between linked drugs 
shared_counts = []
for d1, d2 in linked_drug_pairs:
    shared = len(drug_to_proteins[d1].intersection(drug_to_proteins[d2]))
    shared_counts.append(shared)
avg_shared_linked = sum(shared_counts) / len(shared_counts) if shared_counts else 0

# Metric 2: Average PPI between linked drugs
ppi_counts = []
for d1, d2 in linked_drug_pairs:
    ppi_count = 0
    for p1 in drug_to_proteins[d1]:
        for p2 in drug_to_proteins[d2]:
            if G.has_edge(p1, p2) and G.edges[p1, p2]["type"] == "protein_protein":
                ppi_count += 1
    ppi_counts.append(ppi_count)
avg_ppi_linked = sum(ppi_counts) / len(ppi_counts) if ppi_counts else 0

# Network-level stats
density = nx.density(G)
degrees = dict(G.degree())
avg_degree = sum(degrees.values()) / len(degrees)
max_degree = max(degrees.values())
min_degree = min(degrees.values())


print("=== Graph Statistics ===")
print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")
print(f"Density: {density:.4f}")
print(f"Average degree: {avg_degree:.2f}")
print(f"Maximum degree: {max_degree}")
print(f"Minimum degree: {min_degree}")

print("\n=== Drug Interaction Metrics ===")
print(f"Average shared proteins between linked drugs: {avg_shared_linked:.2f}")
print(f"Average protein–protein interactions between linked drugs: {avg_ppi_linked:.2f}")
