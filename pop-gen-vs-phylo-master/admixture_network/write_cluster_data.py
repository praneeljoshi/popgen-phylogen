import pickle
import os
import sys
from admixture_network import generate_admixture_networks, generate_BirthHybrid_networks, generate_BirthHybrid_networks_admixture_target
from call_ms import call_ms
from write_GTmix_input import write_GTmix_input
from write_TreeMix_input import write_TreeMix_input
from write_PhyloNet_input import write_input, build_bimarker_nexus, build_MCMC_BiMarkers_input, build_MLE_BiMarkers_input

def create_dir(path):
    """
    Creates a directory if it does not already exist.
    """
    if not os.path.exists(path):
        os.makedirs(path)
    return

# base = "C:\\Users\\Eliot Solomon\\Documents\\Research\\pop-gen-vs-phylo\\data"
# base = "/home/ehs3/pop-gen-vs-phylo/data"
base = sys.argv[1]

# Build directory structure
top_level = ["generation_work", "input", "networks", "output"]
methods = ["gtmix", "treemix", "phylonet_mcmc_bimarkers", "phylonet_mle_bimarkers"]
for name in top_level:
    create_dir(base + "/" + name)
for name in methods:
    create_dir(base + "/input/" + name)
    create_dir(base + "/output/" + name)

count = 2

# networks = [list(x) for x in generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500000, True, count)]
# networks = [list(x) for x in generate_BirthHybrid_networks_admixture_target(50, 7, 3, 1, 2, 50, 50, 500000, work_path=f'{base}/generation_work/', beast_prefix="/home/ehs3/pop-gen-vs-phylo/bin/beast/bin/beast")]
networks = [list(x) for x in generate_BirthHybrid_networks_admixture_target(10, int(sys.argv[2]), 6, 1, 50, 50, 50, 500000, hybrid_rate=3, work_path=f'{base}/generation_work/', beast_prefix="/home/ehs3/pop-gen-vs-phylo/bin/beast/bin/beast", ms_prefix="/home/ehs3/pop-gen-vs-phylo/bin/ms/ms")]

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
    MCMC_BiMarkers_input = build_MCMC_BiMarkers_input(bimarker_nexus, bimarker_taxa, bimarker_taxon_map, max_reticulation=int(sys.argv[2]))
    MLE_BiMarkers_input = build_MLE_BiMarkers_input(bimarker_nexus, bimarker_taxon_map, max_reticulation=int(sys.argv[2]))
    write_input(MCMC_BiMarkers_input, f'{base}/input/phylonet_mcmc_bimarkers/{i}.nex')
    write_input(MLE_BiMarkers_input, f'{base}/input/phylonet_mle_bimarkers/{i}.nex')