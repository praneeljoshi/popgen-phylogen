import os
from skbio import DistanceMatrix
from skbio.tree import nj
from write_Structure import write_Structure
import dendropy


def run_Structure(params, dir=''):
    os.system(f'structure {params} -m {dir}Structure_input/mainparams -e {dir}Structure_input/extraparams -i {dir}Structure_input/structure_input.txt '
              f'-o {dir}Structure_input/structure_output')

    pop_count = int(params.split()[3])
    dists = []
    ids = []

    with open('Structure_input/structure_output_f') as outfile:
        output = outfile.readlines()[37 + pop_count: 37 + 2 * pop_count]
        for x in output:
            csv = x.replace('-', '0').split()
            ids.append(csv.pop(0))
            dists.append([float(y) for y in csv])

    dm = DistanceMatrix(dists, ids)
    newick = nj(dm, result_constructor=str)
    print(nj(dm).ascii_art())

def newick_to_nx(newick):
    '''
    takes in a weighted newick tuple and returns a
    '''
if __name__ == '__main__':
    import networkx as nx
    from call_ms import call_ms
    from admixture_network import generate_admixture_networks

    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500, True, 1)
    test_network, test_cmd = test_results[0]

    test_viz = nx.nx_pydot.to_pydot(test_network)
    print(test_viz)

    test_data = call_ms(test_cmd)
    param = write_Structure(test_data, True)
    run_Structure(param)

    # pop_count = 5
    # dists = []
    # ids = []
    # with open('Structure_input/structure_output_f') as outfile:
    #     output = outfile.readlines()[37 + pop_count: 37 + 2 * pop_count]
    #     for x in output:
    #         csv = x.replace('-', '0').split()
    #         ids.append(csv.pop(0))
    #         dists.append([float(y) for y in csv])
    # dm = DistanceMatrix(dists, ids)
    # newick = nj(dm, result_constructor=str)
    #
    # print(ids)
    # print(dists)
    # print(newick)
    # print()