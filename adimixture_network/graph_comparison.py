#TO DO:
#Implement compute_inference_error
#thouroughly test perc_correct_admixed and especially compute_RF_distance


import networkx as nx
import dendropy
from itertools import permutations, product
from copy import deepcopy

test_graph_1 = nx.DiGraph();
# test_graph_2 = nx.DiGraph();
# test_graph_3 = nx.DiGraph();
#
# #              G1                       G2                     G3
# #
# #               0                       0                       0
# #             /   \                   /   \                   /  \
# #            1     2                 |     1                1     4
# #          /   \  /  \               |    /  \
# #         |     3    |               4   5    6
# #         |     |    |
# #         4     5    6
#
test_graph_1.add_nodes_from(range(7))
test_graph_1.add_edges_from([(0,1),(0,2),(1,4),(1,3),(2,3),(2,6),(3,5)])
#
# test_graph_2.add_nodes_from([0,1,4,5,6])
# test_graph_2.add_edges_from([(0,1),(0,4),(1,5),(1,6)])
#
# test_graph_3.add_nodes_from([0,1,4])
# test_graph_3.add_edges_from([(0,1),(0,4)])





def admixture_error(pred_graph, actual_graph):
    """
    computes the percentage of correctly predicted admixed nodes
    note this measure is asymetric, and has to be entered in the order above
    (this implementation only works for extant admixtures, I'm not sure what it would mean otherwise)

    for comparison of admixture proportion, assumes proportion is recorderd in the node's attributes

    Returns percent admixed, list of proportion differences
    """
    total_correct = 0
    total_admixed = 0
    prop_diffs = []
    #iterate through leaves whose parents are admixed
    for node in actual_graph.nodes:
        if not any(actual_graph.neighbors(node)):
            #iterators are annoying and I'm too lazy to write good code
            #checks that node's parent has two parents
            if sum(1 for temp in actual_graph.predecessors(next(actual_graph.predecessors(node)))) == 2:
                total_admixed += 1;
                #check if the same node is admixed in pred_graph
                if sum(1 for temp in pred_graph.predecessors(next(pred_graph.predecessors(node)))) == 2:
                    total_correct += 1;
                    prop_diffs.append(pred_graph.nodes[node]['proportion'] - actual_graph.nodes[node]['proportion'])


    if total_admixed == 0:  #relatively reasonable behavior to avoid div by 0
        print("WARNING: No admixed populations found! Percentage correct set to 0 by default.")
        total_admixed = 1

    return total_correct/total_admixed,







def compute_newick_tuple(dtree, root):
    """
    Recursively translates from nx DiGraph to newick tree
    """
    #base case root is a leaf
    if sum(1 for temp in dtree.neighbors(root)) == 0: #because iterators are annoying
        return root;
    else:
        arr = [] #otherwise break into a tree for each child
        for child in dtree.neighbors(root):
            arr.append(compute_newick_tuple(dtree, child));
        return tuple(arr)





def compute_RF_distance(dtree1, dtree2):
    #It seems that debdropy may use a different notation that gives double the distance as compared to some references I've seen.
    #   something to keep in mind for later
    """
    Converts from an nx DiGraph representation of the two trees to a dendropy Tree representation,
    then computes Robinson Foulds distance

    Could I have made this prettier with a loop around it? Yes. Would that have been harder to read? Also yes.
    """

    if not nx.is_tree(dtree1) or not nx.is_tree(dtree2):
        print("WARNING: dtree1 or dtree2 is not a tree.")

    #find root of each tree
    for root1 in dtree1:
        if sum(1 for temp in dtree1.predecessors(root1)) == 0:
            break;
    for root2 in dtree2:
        if sum(1 for temp in dtree2.predecessors(root2)) == 0:
            break;

    #this is needed for some unknown reason
    tns = dendropy.TaxonNamespace();

    #reformat nx to dendropy string
    s = str(compute_newick_tuple(dtree1, root1)) +";"
    dendro_tree1 = dendropy.Tree.get(data=s, schema="newick", taxon_namespace=tns) #create dendro tree
    print(dendro_tree1.as_string(schema="newick"));

    #reformat nx to dendropy string
    s = str(compute_newick_tuple(dtree2, root2)) +";"
    dendro_tree2 = dendropy.Tree.get(data=s, schema="newick", taxon_namespace=tns) #create dendro tree
    print(dendro_tree2.as_string(schema="newick"));

    return dendropy.calculate.treecompare.symmetric_difference(dendro_tree1,dendro_tree2) #the whole reason for converting



def network_trees(network):
    '''
    takes in a network and returns a list of 2^a trees contained in it, where a is # of admixture nodes
    '''
    trees = []
    a_parents = [] # array of lists of potential incoming edges of each admixture node
    frame = deepcopy(network) # a graph that does not specify lineage of ad nodes upon which trees will be built

    for node in network.nodes:
        if len(list(network.predecessors(node))) == 2: #ID ad node
            in_edges = [(u,node) for u in list(network.predecessors(node))]
            a_parents.append(in_edges)
            frame.remove_edges_from(in_edges)

    # print(a_parents)
    n = len(a_parents) # number of ad nodes
    # makes bit arrays. each bit is a choice of one of two potential incoming edges for each ad node
    possible_choices = set(product({0,1}, repeat=n))
    for choice in possible_choices:
        # print(choice)
        # print(type(deepcopy(frame)))
        t = deepcopy(frame)
        t.add_edges_from([a_parents[i][choice[i]] for i in range(n)])
        trees.append(t)

    # a failsafe i put together in case product() doensnt work how I think it does...

    # edge_choice = [0 for i in range(n)]
    # while (1):
    #     trees.append(deepcopy(frame).add_edges_from({edges[i][edge_choice[i]] for i in range(n)}))
    #     #iterate through every possible bit array of size n (all possible combinations of a_parents)
    #     next = 0
    #     while (next < n and edge_choice[next] == 1):
    #         next += 1
    #     if next == n:
    #         break #we will have accounted for all trees
    #     edge_choice[next] += 1
    #     for i in range(next):
    #         edge_choice[i] = 0
    #     '''
    #     for n = 3, algorithm above will iterate over bit arrays as follows:
    #     0 0 0
    #     1 0 0   next = 0
    #     0 1 0   next = 1
    #     1 1 0   next = 0
    #     0 0 1   next = 2
    #     1 0 1   next = 0
    #     0 1 1   next = 1
    #     1 1 1   next = 0
    #     '''
    return trees


def compute_inference_error(predicted, actual):
    """
    Creates T, T'
    computes RF for all 2^a! possible matching functions m: T -> T'
    "normalizes" RF for each
    returns smallest error over all possible matchings. (best match inference error)
    """
    RF_averages = {}
    T_pred  = network_trees(predicted)
    T_actual = network_trees(actual)
    t_count = len(T_pred)

    # f: T -> [0,t_count-1]
    possible_matchings = set(permutations(range(t_count)))
    for matching in possible_matchings:
        # (T_pred[i], T_actual[matching[i]]) for i in range(t_count) is g: [0,t_count-1] -> T'
        RFs = {compute_RF_distance(T_pred[i], T_actual[matching[i]]) for i in range(t_count)}
        RF_averages.add(sum(RFs)/t_count)

    return min(RF_averages) #BMIE

n = 0
T = network_trees(test_graph_1)
for t in T:
    n += 1
    test_viz = (nx.nx_pydot.to_pydot(t))
    print(test_viz)
    test_viz.write_pdf("test" + str(n) + ".pdf")