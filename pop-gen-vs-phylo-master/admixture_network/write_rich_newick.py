import networkx as nx
import copy

def format_edge_data(data):
    """
    Helper function to format edge data for Rich Newick strings.
    """
    # Remove all trailing None values
    while len(data) > 0 and data[-1] is None:
        data.pop()

    # Convert elements to strings, with None values becoming empty strings
    data = [str(x) if x is not None else "" for x in data]

    # Assemble data strings
    data_string = ":".join(data)

    # Only add the begenning colon if the string is not empty
    if len(data_string) > 0:
        data_string = ":" + data_string

    return data_string

def write_rich_newick_helper(network, node, parent):
    """
    Recursive helper function for write_rich_newick.
    """
    # Call recursively to get child Rich Newick strings, unless we are revisiting an admixture node
    child_strings = []
    if "skip_children" not in network.nodes[node]:
        for child in network.neighbors(node):
            child_strings.append(write_rich_newick_helper(network, child, node))

    # Form this node's string out of the child strings
    if len(child_strings) > 0:
        output = '(' + ",".join(child_strings) + ')'
    else:
        output = ''

    # Add node label to string
    node_type = network.nodes[node]["type"]
    if node_type in ["leaf", "outgroup"]:
        output += str(network.nodes[node]["population"])
    elif node_type in ["internal", "root"]:
        output += "I" + str(node)
    elif node_type == "admixture":
        output += "I" + str(node) + "#H" + str(node)
        network.nodes[node]["skip_children"] = True

    # Add edge data to string (specific keys might need to be tweaked here)
    if node_type != "root":
        length = network.edges[parent, node].get("length")
        support = network.edges[parent, node].get("support")
        if node_type == "admixture":
            if network.nodes[node]["mix_parent"] == parent: # This will produce an error for the GTmix output unless a mix_parent is chosen somehow
                probability = network.nodes[node]["proportion"]
            else:
                probability = 1 - network.nodes[node]["proportion"]
        else:
            probability = None
        
        output += format_edge_data([length, support, probability])

    return output

def write_rich_newick(network):
    """
    Given a properly formatted NetworkX DiGraph as an input, returns the Rich Newick string correpsonding to it.
    """
    root = next(node for node, attrs in network.nodes.items() if attrs["type"] == "root")

    output = write_rich_newick_helper(network, root, None) + ";"

    for attrs in network.nodes.values():
        if "skip_children" in attrs:
            attrs.pop("skip_children")

    return output

if __name__ == "__main__":
    from admixture_network import generate_admixture_networks

    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500000, True, 3)
    test_network, _ = test_results[0]

    # test_network = nx.DiGraph()
    # test_network.add_nodes_from(range(7))
    # test_network.add_edges_from([(0,1), (0,2), (1,3), (1,4), (2,5), (2,6)])

    # test_network.nodes[0]["type"] = "root"
    # test_network.nodes[1]["type"] = "internal"
    # test_network.nodes[2]["type"] = "internal"
    # test_network.nodes[3]["type"] = "leaf"
    # test_network.nodes[4]["type"] = "leaf"
    # test_network.nodes[5]["type"] = "leaf"
    # test_network.nodes[6]["type"] = "leaf"

    # test_network.nodes[3]["population"] = "A"
    # test_network.nodes[4]["population"] = "B"
    # test_network.nodes[5]["population"] = "C"
    # test_network.nodes[6]["population"] = "D"

    test_viz = nx.nx_pydot.to_pydot(test_network)
    test_viz.write_pdf("test.pdf")
    print(test_viz)

    test_string = write_rich_newick(test_network)
    print(test_string)