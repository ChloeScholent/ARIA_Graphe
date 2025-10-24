import networkx as nx
import matplotlib.pyplot as plt
import scipy as sp

n = 7236  #nodes
m = 22270 #edges

G = nx.gnm_random_graph(n, m)

# some properties
print("node degree clustering")
for v in nx.nodes(G):
    print(f"{v} {nx.degree(G, v)} {nx.clustering(G, v)}")

print("the adjacency list")
with open('erdos-renyi_edges.txt', "w") as f:
    for line in nx.edges(G):
        print(line)
        f.write(f'{line}\n')

# pos = nx.spring_layout(G)
# nx.draw(G, pos=pos)
# plt.show()

