import re
import os
import time
from write_rich_newick import write_rich_newick

def strip_edge_data(network_string):
    """
    Strips out the edge data from a Rich Newick string, leaving just the network topology remaining.
    """
    return re.sub(r':[:.0-9Ee-]+', '', network_string)

# print(strip_edge_data("(1:0.054091306690912076,(((2:0.03644309641657073,(5:0.03079916221791329,(6:0.003999897752669873)I14#H14:0.026799264465243416::0.2555694689261563)I7:0.005643934198657441)I5:0.006323279507830765,(4:0.031166933924481208,3:0.031166933924481208)I11:0.011599441999920287)I4:9.03E-3,I14#H14:0.04778727413298191::0.7444305310738437)I3:0.0023041348052602953)I1;"))

def compare_networks(network1, network2, method="luay", phylonet_prefix="java -jar PhyloNet_3.8.2.jar"):
    """
    Computes the distance between two phylogenetic networks, based on their topologies, using PhyloNet.
    """
    # Validate method string
    if method not in ["tree", "tri", "cluster", "luay"]:
        return None

    # Convert NetworkX DiGraphs to Rich Newick strings
    network1_string = write_rich_newick(network1)
    network2_string = write_rich_newick(network2)

    # Strip edge data from Rich Newick strings
    network1_string = strip_edge_data(network1_string)
    network2_string = strip_edge_data(network2_string)

    # network1_string = "((a,(b,(c)x#1)),((d,x#1),e));"
    # network2_string = "((((a, (c)x#1), d), (b, x)), e);"

    # Build PhyloNet input
    input_str = (
        f'#NEXUS\n'
        f'BEGIN NETWORKS;\n'
        f'Network net1 = {network1_string}\n'
        f'Network net2 = {network2_string}\n'
        f'END;\n'
        f'BEGIN PHYLONET;\n'
        f'Cmpnets net1 net2 -m {method};\n'
        f'END;\n'
    )

    # Write PhyloNet input to file
    with open("compare.nex", "w") as f:
        f.write(input_str)

    # Call PhyloNet, and extract the result
    stream = os.popen(f'{phylonet_prefix} compare.nex')
    for line in stream:
        line = line.strip()
        if line.startswith("The"):
            result = tuple(float(x) for x in line.split(":")[1].split())
    
    # Clean up input file
    os.remove("compare.nex")
    
    return result if len(result) > 1 else result[0]
