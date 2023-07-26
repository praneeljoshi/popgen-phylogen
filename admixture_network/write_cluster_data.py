import pickle
from admixture_network import generate_admixture_networks
from call_ms import call_ms
from write_GTmix_input import write_GTmix_input
from write_TreeMix_input import write_TreeMix_input
from write_PhyloNet_input import write_input, build_bimarker_nexus, build_MCMC_BiMarkers_input, build_MLE_BiMarkers_input

base = "C:\\Users\\Eliot Solomon\\Documents\\Research\\pop-gen-vs-phylo\\data"

count = 2

networks = [list(x) for x in generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500000, True, count)]

networks = [x + [call_ms(x[1])] for x in networks] # Will have to be changed for gene trees

for i, network_data in enumerate(networks):
    # Save network using pickle
    pickle.dump(network_data[0], open(f'{base}/networks/{i}.p', 'wb'))

    # Write GTmix input
    write_GTmix_input(network_data[2], f'{base}/input/gtmix/{i}/')

    # Write TreeMix input
    write_TreeMix_input(network_data[2], f'{base}/input/treemix/{i}.gz')

    # Write PhyloNet input
    bimarker_nexus, bimarker_taxa, bimarker_taxon_map = build_bimarker_nexus(network_data[2])
    MCMC_BiMarkers_input = build_MCMC_BiMarkers_input(bimarker_nexus, bimarker_taxa, bimarker_taxon_map)
    MLE_BiMarkers_input = build_MLE_BiMarkers_input(bimarker_nexus, bimarker_taxon_map)
    write_input(MCMC_BiMarkers_input, f'{base}/input/phylonet_mcmc_bimarkers/{i}.nex')
    write_input(MLE_BiMarkers_input, f'{base}/input/phylonet_mle_bimarkers/{i}.nex')