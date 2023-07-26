import sys
import pickle
import networkx as nx
from compare_PhyloNet import compare_networks
from write_rich_newick import write_rich_newick

network_path = sys.argv[1]

network = pickle.load(open(network_path, "rb"))

viz = nx.nx_pydot.to_pydot(network)

print(viz)
print(write_rich_newick(network))

viz.write_pdf("viz2.pdf")