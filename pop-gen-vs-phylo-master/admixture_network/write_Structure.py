import os

def write_Structure(data, use_pop_data, dir='Structure_input'):
    """
    sets up the command for running structure and also prepares the necessary data files
    takes in data, and boolean use_pop_data
    returns run command
    """
    #a dictionary (no values at the moment). may be useful if expanding functionality to construct mainparams ground-up
    #param = {"MAXPOPS" , "BURNIN", "NUMREPS", "INFILE", "OUTFILE", "NUMINDS", "NUMLOCI", "PLOIDY", "MISSING", "ONEROWPERIND", "LABEL", "POPDATA", "POPFLAG", "LOCDATA", "PHENOTYPE", "EXTRACOLS", "MARKERNAMES", "RECESSIVEALLELES", "MAPDISTANCES", "PHASED", "MARKOVPHASE", "NOTAMBIGUOUS"}

    # if not os.path.exists(dir):
    #     os.makedirs(dir)

    #get relevant information for specific run
    loci_count = len(data)
    snps = sum([len(loc[1][0][0]) for loc in data])
    pop_count = len(data[0][1])
    ind_count = len(data[0][1][0])
    #initialize with variable porameters command

    with open(f'{dir}/structure_input.txt', 'w') as input:

        for pop in range(pop_count):
            for indiv in range (ind_count):
                input.write(f'{pop}.{indiv}') #initialize data line
                if use_pop_data: #option to not specify population may be useful (see Structure documentation)
                    input.write(f' {pop} {1}')

                for locus in data:
                    codes = code_alleles(locus[1][pop][indiv])
                    for code in codes:
                        input.write(f' {code}')
                input.write('\n')

    var_inst = f'-L {snps} -K {pop_count} -N {ind_count*pop_count}'

    return var_inst

def code_alleles(hap):
    """
    takes in set of alleles and a start value from which to code
    returns ordered codings of the same alleles
    - at the moment, coding is biallelic and simply parses haplotypes into SNPs
        (another approach could be coding entire haplotypes. assing a code to each allele in the set of unique
        alleles at each loci)
    """
    return [snp for snp in hap]

if __name__ == '__main__':
    import networkx as nx
    from call_ms import call_ms
    from admixture_network import generate_admixture_networks

    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500, True, 1)
    # test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 100, 500, 50, 50, 500000, True, 1)

    test_network, test_cmd = test_results[0]

    test_viz = nx.nx_pydot.to_pydot(test_network)
    # print(test_viz)

    test_data = call_ms(test_cmd)
    # write_TreeMix_input(test_data, "TreeMix_input.txt")
    print(write_Structure(test_data, True))
