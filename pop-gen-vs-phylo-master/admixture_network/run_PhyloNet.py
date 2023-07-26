import os
from parse_rich_newick2 import parse_rich_newick

def extract_network_string(file):
    """
    Given a PhyloNet output file, extracts the Rich Newick string corresponding to the best network it contains.
    """
    with open(file, "r") as f:
        flag = False # Used to get the line after a match is found
        for line in f:
            if flag:
                return line.strip()
            elif line.startswith("Inferred Network #1:") or line.startswith("Overall MAP"):
                flag = True # We want the next line
    with open(file, "r") as f:
        best_network = ""
        flag = False
        for line in f:
            if flag:
                best_network = line.strip().split("]")[-1]
                flag = False
            elif line.startswith("Likelihood : Topology : Full network string"):
                flag = True # We want the next line
        return best_network # Ends up being the final best network, which is what we want
    return None

def run_PhyloNet(file, output=None, phylonet_prefix="PhyloNet_3.8.2.jar", max_heap=None):
    """
    Runs PhyloNet on an input file, and returns the network that it produces.
    Note that some of the command prefixes may have to be altered depending on your OS.
    """
    if output is None:
        output = os.path.splitext(file)[0] + ".out"
    os.system(f'java{" -Xmx" + max_heap + " " if max_heap is not None else ""} -jar {phylonet_prefix} {file} > {output}')
    network_string = extract_network_string(output)
    return parse_rich_newick(network_string)

if __name__ == "__main__":
    # test = run_PhyloNet("test.nex")

    InferNetwork_ML = parse_rich_newick(extract_network_string('Example_Output/InferNetwork_ML.txt'))
    MCMC_GT = parse_rich_newick(extract_network_string('Example_Output/MCMC_GT.txt'))
    MCMC_SEQ = parse_rich_newick(extract_network_string('Example_Output/MCMC_SEQ.txt'))
    MLE_BiMarkers = parse_rich_newick(extract_network_string('Example_Output/MLE_BiMarkers.txt'))
    InferNetwork_MP = parse_rich_newick(extract_network_string('Example_Output/InferNetwork_MP.txt'))
    InferNetwork_MPL = parse_rich_newick(extract_network_string('Example_Output/InferNetwork_MPL.txt'))
    MCMC_BiMarkers = parse_rich_newick(extract_network_string('Example_Output/MCMC_BiMarkers.txt'))
    