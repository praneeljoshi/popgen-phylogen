# https://www.math.utah.edu/mathcircle/tessler-catalan-notes.pdf (slide 14)

from dyck_word import generate_dyck_word
import networkx as nx

def _generate_random_fbt_helper(w, G, parent):
    """
    Recursive helper function for generate_random_fbt. Given a Dyck word w and a target graph G, constructs a FBT.
    """
    node_id = G.graph["counter"] # Store current counter value
    G.graph["counter"] += 1 # Increment counter

    G.add_node(node_id)
    if parent is not None:
        G.add_edge(parent, node_id)

    if len(w) == 0:
         G.nodes[node_id]["type"] = "leaf"
         return
    else:
        # Initialize variables that persist across loop iterations
        x = ""
        y = ""
        zigzag_height = 0
        for i, char in enumerate(w):
            # Keep track of zigzag line height in order to determine if substring up to i is balanced
            if char == "1":
                zigzag_height += 1
            elif char == "2":
                zigzag_height -= 1
        
            if zigzag_height == 0: # Finds the first balanced substring, form is 1x2y
                x = w[1:i] # extract substring between first matching 1 and 2
                y = w[i+1:] # extract the rest of the string, after that matching 2
                break
        
        G.nodes[node_id]["type"] = "internal"
        _generate_random_fbt_helper(x, G, node_id) # Recursive call for left subtree
        _generate_random_fbt_helper(y, G, node_id) # Recursive call for right subtree
        return

def generate_random_fbt(l):
    """
    Randomly generates a full binary tree with l leaf nodes, following a uniform distribution.

    Total nodes n = 2 * l - 1
    """
    G = nx.DiGraph()
    G.graph["counter"] = 0
    w = generate_dyck_word(l - 1)
    _generate_random_fbt_helper(w, G, None)
    G.graph.pop("counter")
    return G

if __name__ == "__main__":
    # import matplotlib.pyplot as plt

    # test = generate_random_fbt(5)
    # nx.draw(test, with_labels=True, font_weight='bold')
    # plt.show()

    test = generate_random_fbt(5)
    test_viz = nx.nx_pydot.to_pydot(test)
    
    print(test_viz)
    test_viz.write_pdf("test_fbt.pdf")