import sys
import os
import pickle
import networkx as nx
from parse_rich_newick2 import parse_rich_newick

network_string = sys.argv[1]
path = os.path.expanduser(sys.argv[2])

network = parse_rich_newick(network_string)

print(f'Dumping network {network_string} to path {path}...')

pickle.dump(network, open(path, 'wb'))

viz = nx.nx_pydot.to_pydot(network)
print(viz)
