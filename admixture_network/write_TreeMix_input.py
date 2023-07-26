import os
import gzip
from admixture_network import generate_admixture_networks
from call_ms import call_ms

def write_TreeMix_input(data, path):
    """
    Given a data tuple returned by call_ms, writes out the TreeMix input that corresponds to that data.
    """
    # with open(path, "w") as f:
    with gzip.open(path, "wt") as f: # File must be gzipped!
        # Count populations, and write the first line of the file based on this count
        pop_count = len(data[0][1])
        f.write(" ".join(str(x) for x in range(1, pop_count + 1)) + "\n")

        # Count alleles per population
        alleles_per_pop = len(data[0][1][0])

        for locus in data: # Iterate through loci
            haplotype_length = len(locus[1][0][0]) # Sample access: [1][population][member]
            for snp_idx in range(haplotype_length): # Iterate through SNP positions
                line = ""
                for pop_idx in range(pop_count): # Iterate through populations, building up a line for the SNP
                    # Count 1's at a given SNP position for this population
                    count = sum(int(locus[1][pop_idx][allele_idx][snp_idx]) for allele_idx in range(alleles_per_pop))
                    line += f"{count},{alleles_per_pop - count} " # Subtract to find the complement
                f.write(line + "\n") # Write SNP line

if __name__ == "__main__":
    import networkx as nx
    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500, True, 1)
    # test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 100, 500, 50, 50, 500000, True, 1)

    test_network, test_cmd = test_results[0]

    test_viz = nx.nx_pydot.to_pydot(test_network)
    print(test_viz)
    test_viz.write_pdf("test.pdf")

    test_data = call_ms(test_cmd)

    # write_TreeMix_input(test_data, "TreeMix_input.txt")
    write_TreeMix_input(test_data, "TreeMix_input.gz")