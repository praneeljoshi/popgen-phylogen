import random
import copy
import networkx as nx
import xml.etree.ElementTree as ET
import os
from random_fbt import generate_random_fbt
from parse_rich_newick2 import parse_rich_newick, modify_BirthDeath_str

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

def BirthHybrid(taxonset_count, taxa_per_set, iterations, origin=0.1, birth_rate=20, hybrid_rate=10, simtag=0, work_path='', beast_prefix='beast'):
    '''
    Runs the BirthHybrid model through Beast2. Formats XML file using given inputs.
    Command prefix will need to be modified if Beast isn't in your PATH.

    Returns an array of networks in Newick format.
    '''
    if work_path != '' and work_path[-1] != '/':
        work_path += '/'

    base_sim_xml = """
    <beast namespace="beast.core:beast.evolution.alignment:beast.app" version="2.4">
    <run id="networkSimulator" outputFileName="" iterations="" spec="speciesnetwork.simulator.BirthHybridSimulator">
        <speciesNetwork id="network:species" spec="speciesnetwork.Network">
            <taxonset id="taxonSuperset" spec="TaxonSet">
            </taxonset>
        </speciesNetwork>
        <parameter id="time.origin" name="origin"></parameter>
        <parameter id="rate.birth" name="birthRate"></parameter>
        <parameter id="rate.hybrid" name="hybridRate"></parameter>
        </run>
    </beast>"""

    network_root = ET.fromstring(base_sim_xml)

    # Do the job of BEAUti by configuring relevant XML files
    network_run = network_root.find('run')
    network_run.set('iterations', str(iterations)) # How many trees
    network_run.set('outputFileName', f'{work_path}sim{simtag}.trees') # For calling function in a loop and storing sims or similar uses
    network_run.find('parameter[@name="origin"]').text = str(origin)
    network_run.find('parameter[@name="birthRate"]').text = str(birth_rate)
    network_run.find('parameter[@name="hybridRate"]').text = str(hybrid_rate)

    taxonsets = network_run.find('speciesNetwork').find('taxonset')
    for set_id in range(1, taxonset_count + 1):
        tset = ET.SubElement(taxonsets, 'taxon', {'id': str(set_id), 'spec': 'TaxonSet'})
        for taxon in range(taxa_per_set):
            ET.SubElement(tset, 'taxon', {'id': f'{set_id}.{taxon}', 'spec': 'Taxon'})

    pretty_xml(network_root)

    with open(f'{work_path}simulation.xml', 'w') as f:
        f.write(ET.tostring(network_root, encoding="unicode"))

    bubble_pop_xml = """
    <beast namespace="beast.core:beast.app" version="2.4">
        <run id="network.summary" spec="speciesnetwork.tools.SummarizePosterior" inputFileName="netsim/sim0.trees" outputFileName="netsim/popped0.trees" burnin="0" />
    </beast>
    """

    bubble_pop_root = ET.fromstring(bubble_pop_xml)
    bubble_pop_run = bubble_pop_root.find('run')
    bubble_pop_run.set('inputFileName', f'{work_path}sim{simtag}.trees')
    bubble_pop_run.set('outputFileName', f'{work_path}popped{simtag}.trees')

    pretty_xml(bubble_pop_root)

    with open(f'{work_path}bubblepop.xml', 'w') as f:
        f.write(ET.tostring(bubble_pop_root, encoding="unicode"))

    os.system(f'{beast_prefix} {work_path}simulation.xml')
    os.system(f'{beast_prefix} {work_path}bubblepop.xml')

    with open(f'{work_path}popped{simtag}.trees') as f:
        output = [line[:-2] for line in f.readlines()]

    return output

def pretty_xml(element):
    '''
    fixes ET's lazy spacing, indenting of xmls for readability.
    takes in a start node and formats the clade that it roots.
    '''
    indent = '    '
    queue = [(element,
              0)]  # node, depth
    while queue:
        node, depth = queue.pop(0)
        kids = [(child, depth + 1) for child in node]
        if kids:
            node.text = '\n' + indent * (depth + 1)
        if queue:
            node.tail = '\n' + indent * queue[0][1]
        else:
            node.tail = '\n' + indent * (depth - 1)
        queue[0:0] = kids  # kids before siblings

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
            attributes["mix_parent"] = min(network.predecessors(node)) # Choose arbitrarily

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

def add_BirthDeath_time_tags(network):
    """
    Takes in an admixture network generated by BEAST, and adds time tags based on branch lengths.
    Modifies the network in place as well as returning it.
    """
    # Iterate over nodes in reverse topological order
    for node in reversed(list(nx.topological_sort(network))):
        if network.nodes[node]['type'] == 'leaf': # Assign leaf nodes a time of 0
            network.nodes[node]['time'] = 0
        else: # Assign other nodes a time equal to their child's time plus the length of the branch between them and their child
            neighbor = next(network.neighbors(node))
            network.nodes[node]['time'] = network.nodes[neighbor]['time'] + network.edges[node, neighbor]['length']
    
    return network

def strip_extra_root(network):
    """
    Removes the extra root node generated by BEAST.
    """
    old_root = next(node for node, attr in network.nodes.items() if attr['type'] == 'root')
    new_root = next(network.neighbors(old_root))
    network.remove_node(old_root)
    network.nodes[new_root]['type'] = 'root'
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
    # time_order_with_leaves = reversed(list(nx.topological_sort(network)))
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
        higher_edges = [x for x in network.edges.data() if x[1] == focus_node]
        if focus_type in ["internal", "root"]: # Merge the two lower edge populations into one higher edge population
            cmd += f"-ej {network.nodes[focus_node]['time']} {lower_edges[0][2]['population']} {lower_edges[1][2]['population']} "
            if focus_type != "root": # Don't propagate upwards from root
                higher_edges[0][2]['population'] = lower_edges[1][2]['population'] # must be this population, not the other one
        elif focus_type == "admixture": # Split the one lower edge population into two higher edge populations
            cmd += f"-es {network.nodes[focus_node]['time']} {lower_edges[0][2]['population']} {network.nodes[focus_node]['proportion']} "
            if network.nodes[focus_node]['mix_parent'] == higher_edges[1][0]: # If necessary, swap higher edges so that mix_parent is always first
                higher_edges[0], higher_edges[1] = higher_edges[1], higher_edges[0]
            higher_edges[0][2]['population'] = lower_edges[0][2]['population'] # proportion --> mix_parent
            higher_edges[1][2]['population'] = pop_count + new_pop_counter # (1 - proportion) --> non mix_parent
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

# def generate_BirthHybrid_networks(n, pop_count, alleles_per_pop, loci_count, mutation, recombination, locus_length, sim_path, bubble_pop_path, netcount=1, simtag='0', beast_dir='', ms_prefix="ms"):
def generate_BirthHybrid_networks(n, pop_count, alleles_per_pop, loci_count, mutation, recombination, locus_length, origin=0.1, birth_rate=20, hybrid_rate=10, work_path="", simtag=0, beast_prefix="beast", ms_prefix="ms"):
    """
    Generates n random admixture networks according to the provided parameters. 
    Returns a list of tuples, where the first element in each tuple is an admixture network, and the second element is a command for ms.
    """
    # network_strings = BirthHybrid(pop_count, alleles_per_pop, n, sim_path, bubble_pop_path, netcount=netcount, simtag=simtag, beast_dir=beast_dir)
    network_strings = BirthHybrid(pop_count, alleles_per_pop, n, origin=origin, birth_rate=birth_rate, hybrid_rate=hybrid_rate, simtag=simtag, work_path=work_path, beast_prefix=beast_prefix)
    networks = []
    for network_str in network_strings:
        network_str = modify_BirthDeath_str(network_str)
        network = parse_rich_newick(network_str)
        if not all(attr["length"] >= 0 for attr in network.edges.values()): # Skip networks that have negative-length branches (very unlikely), should probably be handled more robustly in the future
            continue
        add_BirthDeath_time_tags(network)
        strip_extra_root(network)
        cmd = generate_ms_command(network, alleles_per_pop, loci_count, mutation, recombination, locus_length, ms_prefix=ms_prefix)
        networks.append((network, cmd))
    return networks

def generate_BirthHybrid_networks_admixture_target(n, admixture_count, pop_count, alleles_per_pop, loci_count, mutation, recombination, locus_length, origin=0.1, birth_rate=20, hybrid_rate=10, work_path="", simtag=0, beast_prefix="beast", ms_prefix="ms"):
    """
    Generates n random admixture networks according to the provided parameters. All networks have exactly admixture_count admixture nodes.
    Returns a list of tuples, where the first element in each tuple is an admixture network, and the second element is a command for ms.

    This function is fairly slow, as it just uses rejection sampling.
    """
    networks = []
    while len(networks) < n:
        new_networks = generate_BirthHybrid_networks(n, pop_count, alleles_per_pop, loci_count, mutation, recombination, locus_length, origin=origin, birth_rate=birth_rate, hybrid_rate=hybrid_rate, work_path=work_path, simtag=simtag, beast_prefix=beast_prefix, ms_prefix=ms_prefix)
        for network_pair in new_networks:
            network_admixture_count = len([node for node, attrs in network_pair[0].nodes.items() if attrs['type'] == 'admixture']) # Count admixture nodes
            if network_admixture_count == admixture_count:
                networks.append(network_pair)
        print(f'*** FINISHED GENERATING {len(networks)} / {n} NETWORKS SO FAR ***')
        simtag += 1
    return networks[:n] # Returns exactly n networks (there could be some left over because we're working in batches)
        
if __name__ == "__main__":
    # test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500000, True, 3)
    # print(test_results)
    # test_network, test_cmd = test_results[0]
    
    # test_viz = nx.nx_pydot.to_pydot(test_network)
    # test_viz.write_pdf("test.pdf")
    # print(test_viz)
    # print(test_cmd)

    # trees = BirthHybrid(3, 1, 100, beast_prefix="/home/ehs3/pop-gen-vs-phylo/bin/beast/bin/beast", work_path="/home/ehs3/pop-gen-vs-phylo/repo/admixture_network/generation_work/")
    # trees = BirthHybrid(3, 1, 100)
    # for t in trees:
    #     print(t)

    # test_tree = "(((0[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.009349672795262914)#H2[&gamma=0.786270348702741,gamma_95%HPD={0.786270348702741,0.786270348702741},gamma_mean=0.786270348702741,gamma_median=0.786270348702741,gamma_range={0.786270348702741,0.786270348702741},height_95%HPD={0.009349672795262914,0.009349672795262914},height_mean=0.009349672795262914,height_median=0.009349672795262914,height_range={0.009349672795262914,0.009349672795262914}]:0.025087984589743942,(#H2[&height_95%HPD={0.009349672795262914,0.009349672795262914},height_mean=0.009349672795262914,height_median=0.009349672795262914,height_range={0.009349672795262914,0.009349672795262914}]:0.019226199833411728,((1[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.006027444538420607,(2[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:4.8142157626229753E-4)#H1[&gamma=0.35082257770178016,gamma_95%HPD={0.35082257770178016,0.35082257770178016},gamma_mean=0.35082257770178016,gamma_median=0.35082257770178016,gamma_range={0.35082257770178016,0.35082257770178016},height_95%HPD={4.8142157626229753E-4,4.8142157626229753E-4},height_mean=4.8142157626229753E-4,height_median=4.8142157626229753E-4,height_range={4.8142157626229753E-4,4.8142157626229753E-4}]:0.005546022962158309)S4[&height_95%HPD={0.006027444538420607,0.006027444538420607},height_mean=0.006027444538420607,height_median=0.006027444538420607,height_range={0.006027444538420607,0.006027444538420607}]:0.012261954585493728,#H1[&height_95%HPD={4.8142157626229753E-4,4.8142157626229753E-4},height_mean=4.8142157626229753E-4,height_median=4.8142157626229753E-4,height_range={4.8142157626229753E-4,4.8142157626229753E-4}]:0.017807977547652037)S3[&height_95%HPD={0.018289399123914335,0.018289399123914335},height_mean=0.018289399123914335,height_median=0.018289399123914335,height_range={0.018289399123914335,0.018289399123914335}]:0.010286473504760307)S2[&height_95%HPD={0.028575872628674642,0.028575872628674642},height_mean=0.028575872628674642,height_median=0.028575872628674642,height_range={0.028575872628674642,0.028575872628674642}]:0.005861784756332214)S1[&height_95%HPD={0.034437657385006856,0.034437657385006856},height_mean=0.034437657385006856,height_median=0.034437657385006856,height_range={0.034437657385006856,0.034437657385006856}]:0.06556234261499315)[&height_95%HPD={0.1,0.1},height_mean=0.1,height_median=0.1,height_range={0.1,0.1},topologySupport=0.01];"
    # test_tree = "((((((2[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.019615663453484963)#H1[&gamma=0.2994302774252735,gamma_95%HPD={0.2994302774252735,0.2994302774252735},gamma_mean=0.2994302774252735,gamma_median=0.2994302774252735,gamma_range={0.2994302774252735,0.2994302774252735},height_95%HPD={0.019615663453484963,0.019615663453484963},height_mean=0.019615663453484963,height_median=0.019615663453484963,height_range={0.019615663453484963,0.019615663453484963}]:0.006387872793067476,1[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02600353624655244)S4[&height_95%HPD={0.02600353624655244,0.02600353624655244},height_mean=0.02600353624655244,height_median=0.02600353624655244,height_range={0.02600353624655244,0.02600353624655244}]:0.013530595820480301,(#H1[&height_95%HPD={0.019615663453484963,0.019615663453484963},height_mean=0.019615663453484963,height_median=0.019615663453484963,height_range={0.019615663453484963,0.019615663453484963}]:0.00837804404913188)#H2[&gamma=0.20703895116445048,gamma_95%HPD={0.20703895116445048,0.20703895116445048},gamma_mean=0.20703895116445048,gamma_median=0.20703895116445048,gamma_range={0.20703895116445048,0.20703895116445048},height_95%HPD={0.027993707502616844,0.027993707502616844},height_mean=0.027993707502616844,height_median=0.027993707502616844,height_range={0.027993707502616844,0.027993707502616844}]:0.011540424564415896)S3[&height_95%HPD={0.03953413206703274,0.03953413206703274},height_mean=0.03953413206703274,height_median=0.03953413206703274,height_range={0.03953413206703274,0.03953413206703274}]:0.005084950676059741,((0[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.03134086959764923,#H2[&height_95%HPD={0.027993707502616844,0.027993707502616844},height_mean=0.027993707502616844,height_median=0.027993707502616844,height_range={0.027993707502616844,0.027993707502616844}]:0.0033471620950323855)S5[&height_95%HPD={0.03134086959764923,0.03134086959764923},height_mean=0.03134086959764923,height_median=0.03134086959764923,height_range={0.03134086959764923,0.03134086959764923}]:0.0045555258200908055)#H3[&gamma=0.8844087766171624,gamma_95%HPD={0.8844087766171624,0.8844087766171624},gamma_mean=0.8844087766171624,gamma_median=0.8844087766171624,gamma_range={0.8844087766171624,0.8844087766171624},height_95%HPD={0.035896395417740035,0.035896395417740035},height_mean=0.035896395417740035,height_median=0.035896395417740035,height_range={0.035896395417740035,0.035896395417740035}]:0.008722687325352446)S2[&height_95%HPD={0.04461908274309248,0.04461908274309248},height_mean=0.04461908274309248,height_median=0.04461908274309248,height_range={0.04461908274309248,0.04461908274309248}]:0.04070431210196185,#H3[&height_95%HPD={0.035896395417740035,0.035896395417740035},height_mean=0.035896395417740035,height_median=0.035896395417740035,height_range={0.035896395417740035,0.035896395417740035}]:0.049426999427314294)S1[&height_95%HPD={0.08532339484505433,0.08532339484505433},height_mean=0.08532339484505433,height_median=0.08532339484505433,height_range={0.08532339484505433,0.08532339484505433}]:0.014676605154945677)[&height_95%HPD={0.1,0.1},height_mean=0.1,height_median=0.1,height_range={0.1,0.1},topologySupport=0.01];"
    # test_tree = "((((((1[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0030930701801671068,2[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0030930701801671068)S6[&height_95%HPD={0.0030930701801671068,0.0030930701801671068},height_mean=0.0030930701801671068,height_median=0.0030930701801671068,height_range={0.0030930701801671068,0.0030930701801671068}]:0.02749025893309541)#H2[&gamma=0.46626070679030673,gamma_95%HPD={0.46626070679030673,0.46626070679030673},gamma_mean=0.46626070679030673,gamma_median=0.46626070679030673,gamma_range={0.46626070679030673,0.46626070679030673},height_95%HPD={0.030583329113262517,0.030583329113262517},height_mean=0.030583329113262517,height_median=0.030583329113262517,height_range={0.030583329113262517,0.030583329113262517}]:0.0012576821181869563)#H3[&gamma=0.46536867634862245,gamma_95%HPD={0.46536867634862245,0.46536867634862245},gamma_mean=0.46536867634862245,gamma_median=0.46536867634862245,gamma_range={0.46536867634862245,0.46536867634862245},height_95%HPD={0.03184101123144947,0.03184101123144947},height_mean=0.03184101123144947,height_median=0.03184101123144947,height_range={0.03184101123144947,0.03184101123144947}]:0.019412901704320487,((0[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.02680972598865039)#H1[&gamma=0.9857035230934854,gamma_95%HPD={0.9857035230934854,0.9857035230934854},gamma_mean=0.9857035230934854,gamma_median=0.9857035230934854,gamma_range={0.9857035230934854,0.9857035230934854},height_95%HPD={0.02680972598865039,0.02680972598865039},height_mean=0.02680972598865039,height_median=0.02680972598865039,height_range={0.02680972598865039,0.02680972598865039}]:0.02226357894377356,((#H1[&height_95%HPD={0.02680972598865039,0.02680972598865039},height_mean=0.02680972598865039,height_median=0.02680972598865039,height_range={0.02680972598865039,0.02680972598865039}]:0.005983566708209542)#H4[&gamma=0.25927936779745997,gamma_95%HPD={0.25927936779745997,0.25927936779745997},gamma_mean=0.25927936779745997,gamma_median=0.25927936779745997,gamma_range={0.25927936779745997,0.25927936779745997},height_95%HPD={0.03279329269685993,0.03279329269685993},height_mean=0.03279329269685993,height_median=0.03279329269685993,height_range={0.03279329269685993,0.03279329269685993}]:0.005866800414770876,(#H3[&height_95%HPD={0.03184101123144947,0.03184101123144947},height_mean=0.03184101123144947,height_median=0.03184101123144947,height_range={0.03184101123144947,0.03184101123144947}]:0.0015088325152052967,#H4[&height_95%HPD={0.03279329269685993,0.03279329269685993},height_mean=0.03279329269685993,height_median=0.03279329269685993,height_range={0.03279329269685993,0.03279329269685993}]:5.565510497948373E-4)S5[&height_95%HPD={0.03334984374665477,0.03334984374665477},height_mean=0.03334984374665477,height_median=0.03334984374665477,height_range={0.03334984374665477,0.03334984374665477}]:0.0053102493649760385)S4[&height_95%HPD={0.03866009311163081,0.03866009311163081},height_mean=0.03866009311163081,height_median=0.03866009311163081,height_range={0.03866009311163081,0.03866009311163081}]:0.01041321182079314)S3[&height_95%HPD={0.04907330493242395,0.04907330493242395},height_mean=0.04907330493242395,height_median=0.04907330493242395,height_range={0.04907330493242395,0.04907330493242395}]:0.0021806080033460115)S2[&height_95%HPD={0.05125391293576996,0.05125391293576996},height_mean=0.05125391293576996,height_median=0.05125391293576996,height_range={0.05125391293576996,0.05125391293576996}]:0.03611494504029575,#H2[&height_95%HPD={0.030583329113262517,0.030583329113262517},height_mean=0.030583329113262517,height_median=0.030583329113262517,height_range={0.030583329113262517,0.030583329113262517}]:0.05678552886280319)S1[&height_95%HPD={0.08736885797606571,0.08736885797606571},height_mean=0.08736885797606571,height_median=0.08736885797606571,height_range={0.08736885797606571,0.08736885797606571}]:0.012631142023934297)[&height_95%HPD={0.1,0.1},height_mean=0.1,height_median=0.1,height_range={0.1,0.1},topologySupport=0.01];"
    # test_tree = "(((0[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0247188415732171)#H1[&gamma=0.6220467377408075,gamma_95%HPD={0.04550077336740799,0.9908656749707229},gamma_mean=0.6220467377408075,gamma_median=0.6628491132514791,gamma_range={0.04550077336740799,0.9908656749707229},height_95%HPD={0.007054305844583511,0.04847366730947777},height_mean=0.0247188415732171,height_median=0.021120238695547452,height_range={0.007054305844583511,0.04847366730947777}]:0.034975525070592414,(#H1[&height_95%HPD={0.007054305844583511,0.04847366730947777},height_mean=0.0247188415732171,height_median=0.021120238695547452,height_range={0.007054305844583511,0.04847366730947777}]:0.009061583480214078,(2[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.004827462020965706,1[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.004827462020965706)S3[&height_95%HPD={1.315765215474246E-4,0.015828858294814743},height_mean=0.004827462020965706,height_median=0.002857200413255538,height_range={1.315765215474246E-4,0.015828858294814743}]:0.028952963032465473)S2[&height_95%HPD={0.01227276414513978,0.05100759338750334},height_mean=0.03378042505343118,height_median=0.03565383382155617,height_range={0.01227276414513978,0.05100759338750334}]:0.02591394159037834)S1[&height_95%HPD={0.035951787949442024,0.08229907315787746},height_mean=0.05969436664380952,height_median=0.058019182871923494,height_range={0.035951787949442024,0.08229907315787746}]:0.04030563335619049)[&height_95%HPD={0.1,0.1},height_mean=0.1,height_median=0.1,height_range={0.1,0.1},topologySupport=0.06];"

    # test_tree = "((((2[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.04183317364224939,(3[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.012679861136794798)#H1[&gamma=0.8531915506540236,gamma_95%HPD={0.8531915506540236,0.8531915506540236},gamma_mean=0.8531915506540236,gamma_median=0.8531915506540236,gamma_range={0.8531915506540236,0.8531915506540236},height_95%HPD={0.012679861136794798,0.012679861136794798},height_mean=0.012679861136794798,height_median=0.012679861136794798,height_range={0.012679861136794798,0.012679861136794798}]:0.02915331250545459)S4[&height_95%HPD={0.04183317364224939,0.04183317364224939},height_mean=0.04183317364224939,height_median=0.04183317364224939,height_range={0.04183317364224939,0.04183317364224939}]:0.018680420989285296,((1[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.04396102038812405)#H2[&gamma=0.9927926132498487,gamma_95%HPD={0.9927926132498487,0.9927926132498487},gamma_mean=0.9927926132498487,gamma_median=0.9927926132498487,gamma_range={0.9927926132498487,0.9927926132498487},height_95%HPD={0.04396102038812405,0.04396102038812405},height_mean=0.04396102038812405,height_median=0.04396102038812405,height_range={0.04396102038812405,0.04396102038812405}]:0.009369711597019957,#H1[&height_95%HPD={0.012679861136794798,0.012679861136794798},height_mean=0.012679861136794798,height_median=0.012679861136794798,height_range={0.012679861136794798,0.012679861136794798}]:0.04065087084834921)S3[&height_95%HPD={0.05333073198514401,0.05333073198514401},height_mean=0.05333073198514401,height_median=0.05333073198514401,height_range={0.05333073198514401,0.05333073198514401}]:0.007182862646390675)S2[&height_95%HPD={0.06051359463153468,0.06051359463153468},height_mean=0.06051359463153468,height_median=0.06051359463153468,height_range={0.06051359463153468,0.06051359463153468}]:0.010013396176792082,#H2[&height_95%HPD={0.04396102038812405,0.04396102038812405},height_mean=0.04396102038812405,height_median=0.04396102038812405,height_range={0.04396102038812405,0.04396102038812405}]:0.026565970420202714)S1[&height_95%HPD={0.07052699080832676,0.07052699080832676},height_mean=0.07052699080832676,height_median=0.07052699080832676,height_range={0.07052699080832676,0.07052699080832676}]:0.02947300919167324)[&height_95%HPD={0.1,0.1},height_mean=0.1,height_median=0.1,height_range={0.1,0.1},topologySupport=0.5];"

    # test_tree = "((((((3[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:4.568612263172317E-4)#H1[&gamma=0.9080881870624686,gamma_95%HPD={0.9080881870624686,0.9080881870624686},gamma_mean=0.9080881870624686,gamma_median=0.9080881870624686,gamma_range={0.9080881870624686,0.9080881870624686},height_95%HPD={4.568612263172317E-4,4.568612263172317E-4},height_mean=4.568612263172317E-4,height_median=4.568612263172317E-4,height_range={4.568612263172317E-4,4.568612263172317E-4}]:0.013370237330424917)#H2[&gamma=0.8079722236814135,gamma_95%HPD={0.8079722236814135,0.8079722236814135},gamma_mean=0.8079722236814135,gamma_median=0.8079722236814135,gamma_range={0.8079722236814135,0.8079722236814135},height_95%HPD={0.013827098556742148,0.013827098556742148},height_mean=0.013827098556742148,height_median=0.013827098556742148,height_range={0.013827098556742148,0.013827098556742148}]:6.632734260142437E-4)#H3[&gamma=0.012442598670674654,gamma_95%HPD={0.012442598670674654,0.012442598670674654},gamma_mean=0.012442598670674654,gamma_median=0.012442598670674654,gamma_range={0.012442598670674654,0.012442598670674654},height_95%HPD={0.014490371982756392,0.014490371982756392},height_mean=0.014490371982756392,height_median=0.014490371982756392,height_range={0.014490371982756392,0.014490371982756392}]:0.0326150065405114,(((#H1[&height_95%HPD={4.568612263172317E-4,4.568612263172317E-4},height_mean=4.568612263172317E-4,height_median=4.568612263172317E-4,height_range={4.568612263172317E-4,4.568612263172317E-4}]:0.0023833173954552617,2[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.0028401786217724934)S6[&height_95%HPD={0.0028401786217724934,0.0028401786217724934},height_mean=0.0028401786217724934,height_median=0.0028401786217724934,height_range={0.0028401786217724934,0.0028401786217724934}]:0.018820165517727597)#H4[&gamma=0.6419741250638837,gamma_95%HPD={0.6419741250638837,0.6419741250638837},gamma_mean=0.6419741250638837,gamma_median=0.6419741250638837,gamma_range={0.6419741250638837,0.6419741250638837},height_95%HPD={0.02166034413950009,0.02166034413950009},height_mean=0.02166034413950009,height_median=0.02166034413950009,height_range={0.02166034413950009,0.02166034413950009}]:0.013204351096084047,(#H4[&height_95%HPD={0.02166034413950009,0.02166034413950009},height_mean=0.02166034413950009,height_median=0.02166034413950009,height_range={0.02166034413950009,0.02166034413950009}]:0.0041704256169054765,1[&height_95%HPD={0.0,0.0},height_mean=0.0,height_median=0.0,height_range={0.0,0.0}]:0.025830769756405567)S5[&height_95%HPD={0.025830769756405567,0.025830769756405567},height_mean=0.025830769756405567,height_median=0.025830769756405567,height_range={0.025830769756405567,0.025830769756405567}]:0.00903392547917857)S4[&height_95%HPD={0.03486469523558414,0.03486469523558414},height_mean=0.03486469523558414,height_median=0.03486469523558414,height_range={0.03486469523558414,0.03486469523558414}]:0.012240683287683657)S2[&height_95%HPD={0.047105378523267794,0.047105378523267794},height_mean=0.047105378523267794,height_median=0.047105378523267794,height_range={0.047105378523267794,0.047105378523267794}]:0.011843023302296779,(#H3[&height_95%HPD={0.014490371982756392,0.014490371982756392},height_mean=0.014490371982756392,height_median=0.014490371982756392,height_range={0.014490371982756392,0.014490371982756392}]:0.008301750933149013,#H2[&height_95%HPD={0.013827098556742148,0.013827098556742148},height_mean=0.013827098556742148,height_median=0.013827098556742148,height_range={0.013827098556742148,0.013827098556742148}]:0.008965024359163257)S3[&height_95%HPD={0.022792122915905405,0.022792122915905405},height_mean=0.022792122915905405,height_median=0.022792122915905405,height_range={0.022792122915905405,0.022792122915905405}]:0.03615627890965917)S1[&height_95%HPD={0.05894840182556457,0.05894840182556457},height_mean=0.05894840182556457,height_median=0.05894840182556457,height_range={0.05894840182556457,0.05894840182556457}]:0.04105159817443545)[&height_95%HPD={0.10000000000000002,0.10000000000000002},height_mean=0.10000000000000002,height_median=0.10000000000000002,height_range={0.10000000000000002,0.10000000000000002},topologySupport=0.5];"

    # test_tree = modify_BirthDeath_str(test_tree)
    # print(test_tree)

    # test_network = parse_rich_newick(test_tree)

    # add_BirthDeath_time_tags(test_network)
    # strip_extra_root(test_network)

    # test_viz = nx.nx_pydot.to_pydot(test_network)
    # test_viz.write_pdf("test.pdf")
    # print(test_viz)

    # test_cmd = generate_ms_command(test_network, 4, 2, 50, 50, 500000)

    # print(test_cmd)

    # test_networks = generate_BirthHybrid_networks(2, 3, 1, 2, 50, 50, 500000, work_path="/home/ehs3/pop-gen-vs-phylo/repo/admixture_network/generation_work/", beast_prefix="/home/ehs3/pop-gen-vs-phylo/bin/beast/bin/beast")
    test_networks = generate_BirthHybrid_networks(100, 3, 1, 2, 50, 50, 500000, work_path="netsim/")

    # test_networks = generate_BirthHybrid_networks_admixture_target(50, 7, 3, 1, 2, 50, 50, 500000, work_path="/home/ehs3/pop-gen-vs-phylo/repo/admixture_network/generation_work/", beast_prefix="/home/ehs3/pop-gen-vs-phylo/bin/beast/bin/beast")
    # test_networks = generate_BirthHybrid_networks_admixture_target(10, 1, 3, 1, 2, 50, 50, 500000, hybrid_rate=3, work_path="/home/ehs3/pop-gen-vs-phylo/repo/admixture_network/generation_work/", beast_prefix="/home/ehs3/pop-gen-vs-phylo/bin/beast/bin/beast")

    print(test_networks)