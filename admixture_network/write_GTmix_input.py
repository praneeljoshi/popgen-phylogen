import os
from admixture_network import generate_admixture_networks
from call_ms import call_ms

def write_GTmix_input(data, path):
    """
    Given a data tuple returned by call_ms, writes out the GTmix input that corresponds to that data.
    path refers to the directory that will be filled with locus folders: it will be created if necessary.
    """
    if not os.path.exists(path):
        os.makedirs(path)

    # Write out the population information file
    with open(f'{path}/listPopInfo-all.txt', 'w') as info_file:
        haplotype_counter = 1

        # Iterate through the data for the first locus (assumes that all loci are similarly organized)
        for population_id, population_haplotypes in enumerate(data[0][1], 1):
            haplotype_str = " ".join([str(x) for x in range(haplotype_counter, haplotype_counter + len(population_haplotypes))])
            haplotype_counter += len(population_haplotypes)

            info_file.write(f'{population_id} {len(population_haplotypes)} {haplotype_str}\n')

    # Write out the haplotype files for each locus
    for locus_id, locus in enumerate(data):
        locus_path = f'{path}/{locus_id}'

        if not os.path.exists(locus_path):
            os.makedirs(locus_path)
        
        # Extract position information
        positions = locus[0]

        with open(f'{locus_path}/locus-{locus_id}.hap', 'w') as haplotype_file:
            # Write position information
            positions_str = " ".join([str(x) for x in positions])
            haplotype_file.write(f'{positions_str}\n')

            # Write haplotypes in order
            for population_haplotypes in locus[1]:
                for haplotype in population_haplotypes:
                    haplotype_file.write(f'{haplotype}\n')

if __name__ == "__main__":
    import networkx as nx
    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500000, True, 1)

    test_network, test_cmd = test_results[0]

    test_viz = nx.nx_pydot.to_pydot(test_network)
    test_viz.write_pdf("test.pdf")

    test_data = call_ms(test_cmd)

    write_GTmix_input(test_data, "GTmix_input")
