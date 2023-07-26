import random
import copy
import networkx as nx
from random_fbt import generate_random_fbt

# Define degree conditions

# Root node --> indegree = 0, outdegree = 2
# Internal nodes --> indegree = 1, outdegree = 2
# Admixture nodes --> indegree = 2, outdegree = 1
# Leaf nodes --> indegree = 1, outdegree = 0

max_outdegree = {
    "root": 2,
    "internal": 2,
    "admixture": 1,
    "leaf": 0
}
max_indegree = {
    "root": 0,
    "internal": 1,
    "admixture": 2,
    "leaf": 1
}

def generate_topology_only_extant(pop_count, admixture_count):
    """
    Generates the topology of an admixture network where only extant populations are admixed.
    """
    network = generate_random_fbt(pop_count + admixture_count)
    leaves = {x for x in network.nodes if network.nodes[x]["type"] == "leaf"}

    # Label root node
    network.nodes[0]["type"] = "root"

    # Select pairs of leaves that do not share a parent node
    while True:
        combined_leaves = random.sample(leaves, 2 * admixture_count)
        no_siblings = True
        for i in range(0, len(combined_leaves), 2):
            if next(network.predecessors(combined_leaves[i])) == next(network.predecessors(combined_leaves[i+1])):
                no_siblings = False
                break
        if no_siblings:
            break

    # Form admixture nodes
    for i in range(0, len(combined_leaves), 2):
        # Merge combined_leaves[i] and combined_leaves[i+1] into combined_leaves[i], and make it an admixture node
        network.add_edge(next(network.predecessors(combined_leaves[i+1])), combined_leaves[i])
        network.nodes[combined_leaves[i]]["type"] = "admixture"
        leaves.remove(combined_leaves[i])

        # Make combined_leaves[i+1] the new leaf node
        network.remove_edge(next(network.predecessors(combined_leaves[i+1])), combined_leaves[i+1])
        network.add_edge(combined_leaves[i], combined_leaves[i+1])
    
    return network

# https://www.liebertpub.com/doi/pdf/10.1089/cmb.2015.0228
# http://phylnet.univ-mlv.fr/tools/randomNtkGenerator.php
# https://en.wikipedia.org/wiki/Handshaking_lemma
# https://stackoverflow.com/questions/20246417/how-to-detect-if-adding-an-edge-to-a-directed-graph-results-in-a-cycle

def generate_topologies_any(pop_count, admixture_count, n):
    """
    Generates the topologies of n admixture networks where any populations can be admixed.
    """
    internal_count = admixture_count + pop_count - 2 # Root node is not an internal node for counting purposes (-1 -1)
    total_count = admixture_count + pop_count + internal_count + 1

    base_network = nx.DiGraph()
    base_network.add_nodes_from(range(total_count))

    # Initalize sets that keep track of which nodes are eligible to be the tails and heads of new edges, respectively
    base_eligible_u = set(base_network.nodes)
    base_eligible_v = set(base_network.nodes)

    # Tag nodes and update starting sets appropriately
    for node, attributes in base_network.nodes.items():
        if node == 0: # Root node
            base_network.nodes[node]["type"] = "root"
            base_eligible_v.remove(node)
        elif 0 < node <= internal_count: # Internal node
            base_network.nodes[node]["type"] = "internal"
        elif internal_count < node <= admixture_count + internal_count: # Admixture node
            base_network.nodes[node]["type"] = "admixture"
        else: # Leaf node
            base_network.nodes[node]["type"] = "leaf"
            base_eligible_u.remove(node)

    edge_count = (2 + 3 * (internal_count + admixture_count) + pop_count) // 2 # Handshaking lemma

    networks = []

    while len(networks) < n:
        # Copy the initial network setup for this iteration
        network = copy.deepcopy(base_network)
        eligible_u = copy.deepcopy(base_eligible_u)
        eligible_v = copy.deepcopy(base_eligible_v)

        dead_end = False

        while len(network.edges) < edge_count: # Keep adding edges until we have enough
            while not dead_end:
                v = random.choice(tuple(eligible_v))

                reachable_from_v = [v] + [k for h, k in nx.bfs_edges(network, v)]                

                u_choices = tuple(eligible_u.difference(reachable_from_v).difference(network.predecessors(v))) # Don't choose u's that are reachable from v or already point to v

                if len(u_choices) == 0: # If we are in a situation where no valid edge can be added that ends at v
                    eligible_v.remove(v) # v is no longer eligible
                    if len(eligible_v) == 0 and len(network.edges) != edge_count: # If this empties out our set of eligible v's before we are finished, we must start over
                        dead_end = True
                else:
                    u = random.choice(u_choices)
                    break

            if dead_end: # This allows the inner loop to break all the way out
                break

            # Add edge from u to v
            network.add_edge(u, v)

            # Update sets
            if (network.nodes[u]["type"], network.out_degree(u)) in max_outdegree.items():
                eligible_u.remove(u)
            if (network.nodes[v]["type"], network.in_degree(v)) in max_indegree.items():
                eligible_v.remove(v)
                if len(eligible_v) == 0 and len(network.edges) != edge_count: # If this empties out our set of eligible v's before we are finished, we must start over
                        dead_end = True

        if not dead_end:
            networks.append(network) # Only save the previous network if it is actually complete

    return networks

def add_outgroup_time_mix_tags(network, time_interval, outgroup_time_bonus, admixture_prop):
    """
    Takes in the topology of an admixture network, and adds an outgroup as well as time and admixture proportion tags.
    Modifies the network in place as well as returning it.
    """
    # Scan for leaves and old root, and also add in admixture proportions
    leaves = set()
    old_root = 0

    for node, attributes in network.nodes.items():
        if attributes["type"] == "leaf":
            leaves.add(node)
        elif attributes["type"] == "root":
            old_root = node
        elif attributes["type"] == "admixture":
            attributes["proportion"] = admixture_prop
            attributes["mix_parent"] = min(network.predecessors(node)) # generate_ms_command follows this convention

    # Add the outgroup
    outgroup = max(network.nodes) + 1
    new_root = outgroup + 1

    network.add_node(outgroup)
    network.nodes[outgroup]["type"] = "outgroup"

    network.add_node(new_root)
    network.nodes[new_root]["type"] = "root"

    network.add_edge(new_root, outgroup)
    network.add_edge(new_root, old_root)

    network.nodes[old_root]["type"] = "internal"

    leaves.add(outgroup)

    # Assign times to each internal/admixture node
    network_topo_sort = list(nx.topological_sort(network))
    network_topo_sort_reversed = reversed(network_topo_sort)
    time_order = [x for x in network_topo_sort_reversed if x not in leaves]

    for i, leaf in enumerate(leaves):
        network.nodes[leaf]["time"] = 0
        network.nodes[leaf]["population"] = i + 1

    time = time_interval
    for node in time_order:
        network.nodes[node]["time"] = time
        time += time_interval # Floating point errors?

    network.nodes[new_root]["time"] += outgroup_time_bonus

    return network

def generate_ms_command(network, alleles_per_pop, loci_count, mutation, recombination, locus_length, ms_prefix="ms"):
    """
    Generates a command that uses the ms software package to generate random haplotypes based on a given admixture network
    and other parameters.
    """
    # Scan for leaves
    leaves = set()

    for node, attributes in network.nodes.items():
        if attributes["type"] in ["leaf", "outgroup"]:
            leaves.add(node)

    # Build a list of nodes in time order, excluding leaves
    time_order_with_leaves = sorted(network.nodes, key=lambda x: network.nodes[x]["time"])
    time_order = [x for x in time_order_with_leaves if x not in leaves]

    # Write the beginning of the command (could differ slightly depending on OS and name of ms executable)
    pop_count = len(leaves) # Includes the outgroup
    cmd = f"{ms_prefix} {pop_count * alleles_per_pop} {loci_count} -t {mutation} -r {recombination} {locus_length} -I {pop_count} {(str(alleles_per_pop) + ' ') * pop_count}"

    # Initalize the edges that are connected to leaf nodes with population ids
    for edge in network.edges.data():
        if edge[1] in leaves:
            edge[2]["population"] = network.nodes[edge[1]]["population"]

    # Build up the command by moving up the tree
    new_pop_counter = 1
    for focus_node in time_order:
        focus_type = network.nodes[focus_node]["type"]
        lower_edges = [x for x in network.edges.data() if x[0] == focus_node]
        higher_edges = sorted([x for x in network.edges.data() if x[1] == focus_node], key=lambda y: y[0]) # Define an edge order for admixture
        if focus_type in ["internal", "root"]: # Merge the two lower edge populations into one higher edge population
            cmd += f"-ej {network.nodes[focus_node]['time']} {lower_edges[0][2]['population']} {lower_edges[1][2]['population']} "
            if focus_type != "root": # Don't propagate upwards from root
                higher_edges[0][2]['population'] = lower_edges[1][2]['population'] # must be this population, not the other one
        elif focus_type == "admixture": # Split the one lower edge population into two higher edge populations
            cmd += f"-es {network.nodes[focus_node]['time']} {lower_edges[0][2]['population']} {network.nodes[focus_node]['proportion']} "
            higher_edges[0][2]['population'] = lower_edges[0][2]['population'] # proportion --> parent with lower id
            higher_edges[1][2]['population'] = pop_count + new_pop_counter # (1 - proportion) --> parent with higher id
            new_pop_counter += 1

    # Remove leftover "population" edge attributes
    for edge in network.edges.data():
        edge[2].pop("population", None)

    return cmd

def generate_admixture_networks(pop_count, admixture_count, time_interval, outgroup_time_bonus, admixture_prop, alleles_per_pop, loci_count, mutation, recombination, locus_length, only_extant_admixture, n, ms_prefix="ms"):
    """
    Generates n random admixture networks according to the provided parameters. 
    Returns a list of tuples, where the first element in each tuple is an admixture network, and the second element is a command for ms.
    """
    if only_extant_admixture:
        networks = [generate_topology_only_extant(pop_count, admixture_count) for _ in range(n)]
    else:
        networks = generate_topologies_any(pop_count, admixture_count, n)
    networks_modified = [add_outgroup_time_mix_tags(network, time_interval, outgroup_time_bonus, admixture_prop) for network in networks]
    return list([(network, generate_ms_command(network, alleles_per_pop, loci_count, mutation, recombination, locus_length, ms_prefix)) for network in networks_modified])

if __name__ == "__main__":
    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500000, True, 3)
    print(test_results)
    test_network, test_cmd = test_results[0]

    test_viz = nx.nx_pydot.to_pydot(test_network)
    test_viz.write_pdf("test.pdf")
    print(test_viz)
    print(test_cmd)
