import newick
import networkx as nx
import re

# Notes: some processing might need to be done to the population labels after they are read in (int() call, etc.)

def replace_match(match):
    """
    Handles regular expression matches.
    """
    match = match.group()
    contents = match.split(":")
    return f'{"@".join(contents[:1] + contents[2:])}:{contents[1]}'

def modify_newick_str(string):
    """
    Modifies a Rich Newick string in order to parse it correctly with the newick module.
    """
    return re.sub(r'[a-zA-Z0-9]*(:[.0-9Ee-]*){2,3}', replace_match, string)

def parse_with_newick(string):
    """
    Modifies a string, and then parses it with the newick module.
    """
    return newick.loads(modify_newick_str(string), strip_comments=True)

def newick_to_nx_helper(node, parent_id, network):
    """
    Recursive helper function for newick_to_nx.
    """
    # Store and increment node ID counter
    node_id = network.graph["counter"]
    network.graph["counter"] += 1

    # Add current node to graph
    network.add_node(node_id)

    # Read attributes out of node name
    if node.name is None:
        name = None
        support = None
        proportion = None
    else:
        name_data = node.name.split("@")
        name = name_data[0]
        support = float(name_data[1]) if len(name_data) > 1 and name_data[1] != "" else None
        proportion = float(name_data[2]) if len(name_data) > 2 and name_data[2] != "" else None

    # Add attributes to current node
    if name is None or name.startswith("I") and "#" not in name:
        if node_id == 0: # Root node
             network.nodes[node_id]["type"] = "root"
        else: # Internal node
            network.nodes[node_id]["type"] = "internal"
    elif "#" in name: # Proto-admixture node
        network.nodes[node_id]["type"] = "proto_admixture"
        network.nodes[node_id]["proto_name"] = name
        network.nodes[node_id]["proportion"] = proportion
    else: # Leaf node
        network.nodes[node_id]["type"] = "leaf"
        network.nodes[node_id]["population"] = name # An int() call might be necessary here, eventually

    # Connect node to its parent
    if parent_id is not None:
        network.add_edge(parent_id, node_id)
        network.edges[parent_id, node_id]["weight"] = float(node.length)
        network.edges[parent_id, node_id]["support"] = support
    
    # Call recursively for node's children
    for child in node.descendants:
        newick_to_nx_helper(child, node_id, network)
    
    return

def newick_to_nx(tree):
    """
    Converts a tree generated by the newick module into a NetworkX DiGraph.
    """
    # Extract node from tree list
    node = tree[0] # Is it possible for this list to have more than one element?

    # Initalize network object
    network = nx.DiGraph()
    network.graph["counter"] = 0

    # Call recursive helper function
    newick_to_nx_helper(node, None, network)

    # Remove leftover counter graph attribute
    network.graph.pop("counter")

    return network

def merge_proto_admixture(network):
    """
    Merge together matching proto-admixture nodes into proper admixture nodes. Mutates the input network.
    """
    # Generate list of all proto-admixture nodes
    proto_admixture_nodes = [node for node, attrs in network.nodes.items() if attrs["type"] == "proto_admixture"]
    
    # Sort list in order to group nodes with same proto_name and always choose the higher proportion as node1
    proto_admixture_nodes.sort(key=lambda node: (network.nodes[node]["proto_name"], network.nodes[node]["proportion"]), reverse=True)
    
    # Iterate over pairs of nodes
    for i in range(0, len(proto_admixture_nodes), 2):
        node1, node2 = proto_admixture_nodes[i:i+2]

        # Record mix_parent
        network.nodes[node1]["mix_parent"] = next(network.predecessors(node1))

        # Merge edges incoming to node2
        for u in network.predecessors(node2):
            network.add_edge(u, node1)
            for key, value in network.edges[u, node2].items():
                network.edges[u, node1][key] = value

        # Merge edges outgoing from node2
        for v in network.neighbors(node2):
            network.add_edge(node1, v)
            for key, value in network.edges[node2, v].items():
                network.edges[node1, v][key] = value

        # Remove node2
        network.remove_node(node2)

        # Remove proto_name attribute and change node type
        network.nodes[node1].pop("proto_name")
        network.nodes[node1]["type"] = "admixture"

    return

def parse_rich_newick(string):
    """
    Parses a network represented as a Rich Newick string, and returns it as a NetworkX DiGraph.
    """
    network = newick_to_nx(parse_with_newick(string))
    merge_proto_admixture(network)
    return network

if __name__ == "__main__":
    # test = "(((B:0.0)#H1:0.05::0.8,(C:0.002,#H1:0.002::0.2):0.048):0.01,A:0.06);"
    # test = "(((3, 4)Z#H1, 1), (Z#H1, 2));"
    # test = "[0.016870977597323037](((R:0.037090964181698674,(Q:0.0032885622841274066,((L:0.0012618693588044786,A:0.0012618693588044786)I0:7.205400554551246E-4)I4#H1:0.0013061528698678033::0.8806236173358228)I2:0.03380240189757127)I5:0.051455374837895,C:0.08854633901959368)I3:0.014942950793962975,I4#H1:0.10150688039929705::0.11937638266417716)I1;"
    test = "((((L:0.014588518345685399,(Q:0.0030585579153499486,(R:0.0029652591560123936)I2#H1:9.329875933755498E-5::0.8662633897731737)I4:0.01152996043033545)I0:0.009561663688947442,A:0.02415018203463284)I5:0.056399078448244724,C:0.08054926048287757)I1:8.753025849028901,I2#H1:8.830609850355767::0.13373661022682626)I3;"

    print(modify_newick_str(test))

    test_network = parse_rich_newick(test)

    test_viz = nx.nx_pydot.to_pydot(test_network)
    print(test_viz)
    test_viz.write_pdf("test.pdf")
