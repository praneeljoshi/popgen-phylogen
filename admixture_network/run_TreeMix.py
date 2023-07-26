import os
import gzip
import networkx as nx

def run_TreeMix(file, admixture_count, outgroup=None, snp_group_size=None, output="out_stem", treemix_prefix="treemix"):
    """
    Runs TreeMix on an input file, and returns the network that it produces.
    Note that some of the command prefixes may have to be altered depending on your OS.
    """
    # Run TreeMix with given parameters
    cmd = f'{treemix_prefix} -i {file} -o {output} -m {admixture_count} '
    if outgroup is not None:
        cmd += f'-root {outgroup} '
    if snp_group_size is not None:
        cmd += f'-k {snp_group_size} '
    os.system(cmd)
    
    # Parse output graph
    with gzip.open(f'{output}.vertices.gz', 'rt') as vertices_file, gzip.open(f'{output}.edges.gz', 'rt') as edges_file:
        # Initalize network
        network = nx.DiGraph()

        # Add nodes
        for node_line in [x.split() for x in vertices_file.readlines()]:
            node_id = int(node_line[0])
            network.add_node(node_id)
            if node_line[2] == "ROOT":
                network.nodes[node_id]["type"] = "root"
            elif node_line[4] == "TIP":
                network.nodes[node_id]["type"] = "leaf"
                network.nodes[node_id]["population"] = int(node_line[1])
            else:
                network.nodes[node_id]["type"] = "internal"
        
        # Add edges
        for edge_line in [x.split() for x in edges_file.readlines()]:
            u_id = int(edge_line[0])
            v_id = int(edge_line[1])
            network.add_edge(u_id, v_id)
            network.edges[u_id, v_id]["weight"] = float(edge_line[2])
            if edge_line[4] == "MIG":
                network.edges[u_id, v_id]["admixture"] = float(edge_line[3]) # This seems to be right for the proportion, but I'm not totally sure

        # Alter admixture events so that there is an admixture node for each one
        for mix_edge in [edge for edge in network.edges if "admixture" in network.edges[edge]]: # Iterate over admixture edges
            mix_parent, mix_target = mix_edge
            mix_prop = network.edges[mix_edge]["admixture"]
            mix_node = max(network.nodes) + 1 # Create an admixture node
            network.add_node(mix_node)
            network.nodes[mix_node]["type"] = "admixture"
            for parent_node in list(network.predecessors(mix_target)): # Iterate over edges incoming to target node
                network.add_edge(parent_node, mix_node) # Add edge to admixture node
                for key, value in network.edges[parent_node, mix_target].items(): # Copy attributes from previous edge
                    network.edges[parent_node, mix_node][key] = value
                network.remove_edge(parent_node, mix_target) # Delete previous edge
            network.add_edge(mix_node, mix_target) # Connect admixture node and target node (left unweighted)
            network.nodes[mix_node]["mix_parent"] = mix_parent # Record node from which the admixture edge originated
            network.nodes[mix_node]["proportion"] = mix_prop # Record admixture "proportion": note that the TreeMix paper says that this value is *correlated with* the proportion, and is not perfectly scaled
            network.edges[mix_parent, mix_node].pop("admixture", None) # Clean up admixture attribute

    return network

if __name__ == "__main__":
    inferred_network = run_TreeMix("TreeMix_input.gz", 1, 5)

    test_viz = nx.nx_pydot.to_pydot(inferred_network)
    test_viz.write_pdf("inferred.pdf")
    print(test_viz)
