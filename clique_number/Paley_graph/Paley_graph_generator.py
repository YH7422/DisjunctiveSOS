import numpy as np
import networkx as nx

qvec = [17, 29, 37, 41, 61]

for q in qvec:
    G = nx.paley_graph(q)
    A = nx.to_numpy_array(G)
    filename = 'Paley_q' + str(q) + '.csv'
    np.savetxt(filename, A)
