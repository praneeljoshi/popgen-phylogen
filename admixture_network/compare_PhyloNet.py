import os
import time
from write_rich_newick import write_rich_newick

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
