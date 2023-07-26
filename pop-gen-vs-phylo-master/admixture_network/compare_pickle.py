import sys
import pickle
import networkx as nx
from compare_PhyloNet import compare_networks
from write_rich_newick import write_rich_newick

network1_path = sys.argv[1]
network2_path = sys.argv[2]

network1 = pickle.load(open(network1_path, "rb"))
network2 = pickle.load(open(network2_path, "rb"))

# for attr in network2.edges.values():
#     if "weight" in attr:
#         attr["length"] = attr["weight"]
#         attr.pop("weight")

# viz1 = nx.nx_pydot.to_pydot(network1)
# viz2 = nx.nx_pydot.to_pydot(network2)

# print(viz1)
# print(viz2)

# viz1.write_pdf("easier1.pdf")
# viz2.write_pdf("harder1.pdf")

print(compare_networks(network1, network2))

# print(write_rich_newick(network1))
# print(write_rich_newick(network2))
