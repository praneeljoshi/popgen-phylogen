import networkx as nx

# TODO
# Quoted strings "" for node names (not super important for now)

def unpack_node_data(node_string):
    """
    Unpacks the data associated with a node string, and returns it as a list.
    """
    data = node_string.split(":")
    data = [None if x == '' else x for x in data] + [None] * (4 - len(data))
    for i in range(len(data)):
        if data[i] is not None:
            if i == 0:
                data[i] = int(data[i]) if data[i].isnumeric() else str(data[i])
            else:
                data[i] = float(data[i])
    return data

def merge_nodes(node1, node2, network):
    """
    Merges node1 into node2. Overwrites conflicting information attached to node2 with that of node1.
    """
    # Merge attributes
    for key, value in network.nodes[node1].items():
        network.nodes[node2][key] = value

    # Merge incoming edges
    for u in network.predecessors(node1):
        network.add_edge(u, node2)
        for key, value in network.edges[u, node1].items():
            network.edges[u, node2][key] = value

    # Merge outgoing edges
    for v in network.neighbors(node1):
        network.add_edge(node2, v)
        for key, value in network.edges[node1, v].items():
            network.edges[node2, v][key] = value

    # Delete node1
    network.remove_node(node1)

    return

def merge_admixture_using_tag(node, mix_tag, network):
    """
    Merges node with another node that has the hybrid tag mix_tag, if such an other node exists. Otherwise, records that node has the hybrid tag mix_tag.
    """
    mix_tags = network.graph["mix_tags"]
    if mix_tag in mix_tags: # Other node exists
        other_node = mix_tags[mix_tag]
        merge_nodes(node, other_node, network)
    else: # Other node does not exist
        mix_tags[mix_tag] = node
    return

def parse_rich_newick_helper(string, parent, network):
    """
    Recursive helper function for parse_rich_newick.
    """
    # Create node and connect it to parent
    node = network.graph["counter"]
    network.graph["counter"] += 1
    network.add_node(node)
    if parent is not None:
        network.add_edge(parent, node)

    # Parse string
    if len(string) > 0 and string[0] == "(": # Recursive case: list of child nodes
        # Initalize variables
        chunk = ""
        paren_depth = 0
        internal_data = None

        # Iterate over characters in string
        for i, char in enumerate(string):
            if char == "," and paren_depth == 1: # A full sublist has been captured in the chunk
                parse_rich_newick_helper(chunk, node, network) # Call recursively for current chunk
                chunk = ""
            elif char == "(":
                paren_depth += 1
                if paren_depth != 1: # Skip outer list open parenthesis
                    chunk += char
            elif char == ")":
                if paren_depth != 1: # Skip outer list close parenthesis
                    chunk += char
                paren_depth -= 1
            elif paren_depth == 0: # End of the outer list has been reached
                internal_data = string[i:]
                break
            elif char != " ": # Skip space characters
                chunk += char
        
        # Call recursively for final chunk
        parse_rich_newick_helper(chunk, node, network)

        # Record node types and fill placeholder length and support values
        if parent is None:
            network.nodes[node]["type"] = "root"
        else:
            network.nodes[node]["type"] = "internal"
            network.edges[parent, node]["length"] = None
            network.edges[parent, node]["support"] = None

        # Record internal node data if it is present
        if internal_data is not None:
            data = unpack_node_data(internal_data)
            if parent is not None:
                network.edges[parent, node]["length"] = data[1]
                network.edges[parent, node]["support"] = data[2]
                if isinstance(data[0], str) and "#" in data[0]: # Internal admixture node
                    network.nodes[node]["type"] = "admixture"
                    network.nodes[node]["mix_parent"] = parent
                    network.nodes[node]["proportion"] = data[3]
                    mix_tag = data[0].split("#")[1]
                    merge_admixture_using_tag(node, mix_tag, network)
                    
    else: # Base case: leaf node or leaf admixture node
        data = unpack_node_data(string)
        if parent is not None:
            network.edges[parent, node]["length"] = data[1]
            network.edges[parent, node]["support"] = data[2]
        if isinstance(data[0], str) and "#" in data[0]: # Leaf admixture node
            network.nodes[node]["type"] = "admixture"
            network.nodes[node]["mix_parent"] = parent
            network.nodes[node]["proportion"] = data[3]
            mix_tag = data[0].split("#")[1]
            merge_admixture_using_tag(node, mix_tag, network)
        else: # Leaf node
            network.nodes[node]["type"] = "leaf"
            network.nodes[node]["population"] = data[0]
    return

def strip_bracket_comments(string):
    """
    Strips comments, enclosed in square brackets, from string.
    """
    new_string = ""
    comment = False
    for i, char in enumerate(string):
        if char == "[":
            comment = True
        elif char == "]":
            comment = False
        elif comment:
            continue
        else:
            new_string += char
    return new_string

def parse_rich_newick(string, strip_comments=True):
    """
    Parses a network in the Rich Newick format into a NetworkX DiGraph.
    """
    # Strip comments from string, if requested
    if strip_comments:
        string = strip_bracket_comments(string)

    # Remove trailing semicolon if it is present
    if string[-1] == ";":
        string = string[:-1]

    # Initalize DiGraph object
    network = nx.DiGraph()
    network.graph["counter"] = 0
    network.graph["mix_tags"] = {}

    # Call recursive helper function
    parse_rich_newick_helper(string, None, network)

    # Clean up counter attribute
    network.graph.pop("counter")
    network.graph.pop("mix_tags")

    return network

if __name__ == "__main__":
    # test = "(B,(A,C,E),D);"
    # test = "(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"
    # test = "(B:6.0,(A:5.0,C:3.0,E:4.0)Ancestor1:5.0,D:11.0);"
    # test = "((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);"
    # test = "(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;"
    # test = "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);"
    # test = "A;"
    # test = "((A,B),(C,D));"
    # test = "(Alpha,Beta,Gamma,Delta,,Epsilon,,,);"

    # test = "(((B:0.0)#H1:0.05::0.8,(C:0.002,#H1:0.002::0.2):0.048):0.01,A:0.06);"
    # test = "(((3, 4)Z#H1, 1), (Z#H1, 2));"
    test = "[0.016870977597323037](((R:0.037090964181698674,(Q:0.0032885622841274066,((L:0.0012618693588044786,A:0.0012618693588044786)I0:7.205400554551246E-4)I4#H1:0.0013061528698678033::0.8806236173358228)I2:0.03380240189757127)I5:0.051455374837895,C:0.08854633901959368)I3:0.014942950793962975,I4#H1:0.10150688039929705::0.11937638266417716)I1;"
    # test = "((((L:0.014588518345685399,(Q:0.0030585579153499486,(R:0.0029652591560123936)I2#H1:9.329875933755498E-5::0.8662633897731737)I4:0.01152996043033545)I0:0.009561663688947442,A:0.02415018203463284)I5:0.056399078448244724,C:0.08054926048287757)I1:8.753025849028901,I2#H1:8.830609850355767::0.13373661022682626)I3;"

    test_network = parse_rich_newick(test)

    test_viz = nx.nx_pydot.to_pydot(test_network)
    print(test_viz)
    test_viz.write_pdf("test.pdf")
