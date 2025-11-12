import random
import networkx as nx
import matplotlib.pyplot as plt

# Exemple : garder 5000 arêtes au hasard
sample_edges = random.sample(list(G.edges()), 5000)

H = nx.Graph()
H.add_edges_from(sample_edges)

plt.figure(figsize=(10,8))
pos = nx.spring_layout(H, seed=42)
nx.draw(H, pos, node_size=10, node_color='lightblue', edge_color='gray', with_labels=False)
plt.show()