import sys
import os
import pickle
import networkx as nx
from parse_rich_newick2 import parse_rich_newick
from run_PhyloNet import extract_network_string

in_path = os.path.expanduser(sys.argv[1])
out_path = os.path.expanduser(sys.argv[2])

network_string = extract_network_string(in_path)
network = parse_rich_newick(network_string)

print(f'Dumping network {network_string} to output path {out_path}...')

pickle.dump(network, open(out_path, 'wb'))

viz = nx.nx_pydot.to_pydot(network)
print(viz)
