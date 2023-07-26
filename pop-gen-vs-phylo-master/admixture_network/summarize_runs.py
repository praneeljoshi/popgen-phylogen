import sys
import os
import pickle
import statistics
from compare_PhyloNet import compare_networks

def summarize_results(base_dir, compare_method="luay", phylonet_prefix="java -jar PhyloNet_3.8.2.jar"):
    """
    Given the path to a directory of results, outputs a dictionary where the keys are methods and the values are lists of distances.
    """
    methods = ["gtmix", "phylonet_mcmc_bimarkers", "phylonet_mle_bimarkers", "treemix"]

    results = {method:[] for method in methods}

    for method, method_data in results.items():
        output_network_files = [f for f in os.listdir(f'{base_dir}/output/{method}') if os.path.isfile(f'{base_dir}/output/{method}/{f}')]

        for output_network_file in output_network_files:
            output_network = pickle.load(open(f'{base_dir}/output/{method}/{output_network_file}', "rb"))
            input_network = pickle.load(open(f'{base_dir}/networks/{output_network_file}', "rb"))

            distance = compare_networks(output_network, input_network, method=compare_method, phylonet_prefix=phylonet_prefix)
            method_data.append(distance)

    return results

if __name__ == "__main__":
    data_dir = os.path.expanduser(sys.argv[1])

    pop_counts = [4, 6, 8, 10]

    all_results = {}
    # all_results = {pop_count:summarize_results(f'{data_dir}/{pop_count}_pop/') for pop_count in pop_counts}

    # for pop_count in pop_counts:
    #     all_results[pop_count] = summarize_results(f'{data_dir}/{pop_count}_pop/')

    # print(all_results)

    all_results = {4: {'gtmix': [6.0, 6.0, 3.0, 2.0, 6.0, 6.0, 6.0, 6.0, 6.0, 3.0], 'phylonet_mcmc_bimarkers': [3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.0, 3.0, 4.0, 3.0], 'phylonet_mle_bimarkers': [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.0, 3.0, 3.0, 3.0], 'treemix': [6.0, 6.0, 6.0, 6.0, 4.0, 3.0, 6.0, 6.0, 6.0, 6.0]}, 6: {'gtmix': [6.0, 8.0, 8.0, 8.0, 7.0, 8.0, 8.0, 8.0, 7.0, 8.0], 'phylonet_mcmc_bimarkers': [6.0, 6.0, 6.0, 6.0, 4.0, 7.0, 6.0, 6.0, 4.0, 6.0], 'phylonet_mle_bimarkers': [8.0, 6.0, 8.0, 8.0, 7.0, 8.0, 8.0, 8.0, 3.0, 8.0], 'treemix': [8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 5.0, 8.0, 8.0]}, 8: {'gtmix': [7.0, 8.0, 9.0, 9.0, 5.0, 10.0, 9.0, 10.0, 9.0, 9.0], 'phylonet_mcmc_bimarkers': [7.0, 5.0, 7.0, 10.0, 9.0, 8.0, 8.0, 8.0, 6.0, 6.0], 'phylonet_mle_bimarkers': [9.0, 8.0, 9.0, 10.0, 7.0, 10.0, 9.0, 10.0, 10.0, 6.0], 'treemix': [9.0, 10.0, 8.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0, 10.0]}, 10: {'gtmix': [10.0, 11.0, 10.0, 11.0, 12.0, 12.0, 12.0, 10.0, 12.0, 9.0], 'phylonet_mcmc_bimarkers': [10.0, 10.0, 10.0, 11.0, 10.0], 'phylonet_mle_bimarkers': [10.0, 8.0, 11.0, 12.0, 12.0, 12.0, 11.0, 8.0, 9.0], 'treemix': [12.0, 9.0, 10.0, 9.0, 12.0, 12.0, 12.0, 11.0, 12.0, 9.0]}}

    print("method: average distance (standard deviation)")

    for pop_count, results in all_results.items():
        print(f'\nPopulation count: {pop_count}\n')
        for method, distances in results.items():
            print(f'"{method}": {statistics.mean(distances)} ({statistics.stdev(distances)})')
