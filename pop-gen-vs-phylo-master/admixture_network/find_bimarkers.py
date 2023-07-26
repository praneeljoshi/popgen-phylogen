import itertools

# data[locus][population][individual] --> bimarkers
# data[locus][locations/haplotypes][population][individual] --> sequences

def find_bimarkers(data):
    """
    Given a data tuple returned from call_seq_gen, returns a data list corresponding to the biallelic markers in that data. The format of the
    returned list is identical to that of the tuple returned by call_ms.
    """
    # Initalize output list
    output = []

    # Iterate over loci
    for i, locus in enumerate(data):
        # Initalize new locus structure
        new_locus = ["", [[""] * len(pop) for pop in locus]]

        # Count the length of the sequence at this locus
        n = len(locus[0][0])

        # Iterate over all of the base positions, along with a tuple of every base that exists at that position within the locus
        for j, bases in enumerate(zip(*itertools.chain(*locus))):
            # Count the number of unique bases
            unique_bases = len(set(bases))

            if unique_bases == 2: # Position corresponds to a SNP
                # Count the sizes of each population, which are used to reassemble the new locus
                pop_sizes = [len(x) for x in locus]

                # Print out some information; if this is removed, the outermost enumerate() can be as well
                # print(f'SNP found at locus {i}, position {j} = {j / n}')

                # Record the position of the SNP
                new_locus[0] += f'{j / n} '

                # Initalize a counter to keep track of the base index
                k = 0

                # Iterate over populations
                for pop, pop_size in enumerate(pop_sizes):
                    # Iterate over individuals within this population
                    for indiv in range(pop_size):
                        # Identify the individual's base at this location
                        base = bases[k + indiv]

                        # Convert the base to a "0" or "1"
                        new_locus[1][pop][indiv] += "0" if base == bases[0] else "1"
                    
                    # Update the counter
                    k += pop_size

        # Record the new locus
        new_locus[0] = new_locus[0].strip()
        output.append(new_locus)
    
    return output

if __name__ == "__main__":
    import networkx as nx
    from admixture_network import generate_admixture_networks
    from call_ms import call_ms
    from call_seq_gen import call_seq_gen
    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 70, True, 1)

    test_network, test_ms_cmd = test_results[0]

    test_viz = nx.nx_pydot.to_pydot(test_network)
    test_viz.write_pdf("test1.pdf")

    _, test_tree_lines = call_ms(test_ms_cmd, True)

    test_seqs = call_seq_gen(test_tree_lines, 70, 4, 4, seq_gen_prefix=".\seq-gen.exe")

    # test_seqs = ((('CAAACCGCCCAGACTCTACCTTCGGTTTGAGTTCAATGTCGGGCCAAACCATCCGCGACCCTAGTAAAGT', 'CGACCCTCTGAGACTCAACGAGCGGTCTGGGGTCGAAGTAGAACGAAACTATCGGCGACTCTGATACAGT', 'CGCATACCTGATATTCTGCTTGGTGTCTGAGGTCACTGTAGGGCGAGACTAACCCACACGCTAGTAGGGT', 'CAAACCGCCCAGACTTTACCTTCGGTTTGAGTTCAATGTTGGGCCAAACCATCCGCGACCCTAGTAAAGT'), ('CTAGCCTCTGGTATTCACCCTGCGGTCTTGGTTTACTGCAGGGCAACACTATCGGCTACCCAGATAGAGT', 'CGTACCGTTGGTACTCAACCTGGGGTTTGAGTTCACTACAGGGCGAAGCTAACGGCTACCCAGATAGAGT', 'CCAACCGCTGGAACCCAACCTTCGGTCTGAGGTTACTGTATGGCGAAGCTATCCGCGTCCCAGCTAGCGT', 'CGCACACGTGATATTCTGCCTGCGGTCTGAGATGACTGCAGGGCAAAGCTATCCGCGTGTCAGGTAGAGT'), ('CGACCATCTAAGACTCAACCTGGTGCCTGAGGTCGCTGTAAGGCGGGACGTACCTGCACGCTAATTTAGT', 'CATACCGCTGGAATTCATCGTGGGGTCGTAGCTTTCTGCAGGGCAAAACTAACTCACACGAACAGAAGGT', 'CATACCGCTGGAATTCAACGTGGGGTCGTAGTTTTCCGCAGGGCAAAACTAACTCACACGAACATAGAGC', 'CGACCCTGGGGAACTGTTCGAGCGGTCGTAGTTTACGGCAGCGCTAAACTACCAGCGTCCCTAATTTAGT'), ('CCATCCGCTCAGACTCTAACTGCGGTGTGAGGTCAGTGTATGGCAAAATTATCGGCTACCCTAATAGGGT', 'CATACCGCTTAGACTCATCGAGCGGTCCTAGTTTACTGCAGCGCTAAACTATCCCACACGCATGCAAAGT', 'AATACCCCCAGAACTCACCCTGCGGTCTGAGGTCACTGTATGGCAAAATTATCGGCTACCCTAACAGGGT', 'CGCATAACTGAAACCCTACCTGCGGAGTGAGTTCACTCTAGGGCGAGACTAACCCACACGCAAGTAAAGT'), ('ACCTCACACAAACTCCACGGATAGTGCGGTATCCGCATTGGAACAGTTCACATTGACACCTTAACGGGGC', 'ACCTCACACCAACTCGACTGCAATCGCGACATTGGCCTAGACAGGGCCCTAAATTACCGCTTAATAAGTT', 'CCATGACGCCTGCTCGACTAATGTGGCGACATTGGCATAGTCAGGAGCCTAACGTAGACCATTGTCAGTC', 'CCATAACGCGAGATCCATCCCAGTCGCGACTCAGGCATGGTCAGGGGCCTAATTTTCCCCTTAATAAGTT')), (('TGCAAGATGATATATTATCATATTCTATTTACATTCGCTGGACTTCGCGTATCGCGTTACGTCGCCTGTG', 'AGTACCTCGTCCACCGATCATGCAGTCTTTAAGTTCGTGGTTTATAGCTCGTATAGATATGTAACCTGTG', 'AGTACCTCGGTCACCGTTCATGCACTGTTTAAATACAATCTGCACAGCGTATAGCGATACGGCGTCGGGG', 'TGCATGTCAATCACCGATCATGTTGTAGACAATTTCTCTGGTCATAACGGCTCGGGTTAAGCCGTCTGGG'), ('TGCAAGTATGTTACCGAAGCCAATTTTTTTAAATTCGCTGTGTATAGCGCTTATTGATATCCCACCTGGG', 'TGTAACTACGTCCATAGTCTTCACCTGTAGAAGATCCCTGGGTCTGTACCTTCTGGAGACGGCGTCAGGC', 'TGCATGTTTGTCTTCTTTGACATCGTTATCAAATTCCCCGGGTCTGTTCACTCCGGTGTCCGAGCCAGGG', 'GGCATTTCATTCACGGCTCCTCCACTATTTAAGTTCAAGGTGCACAACGTATATTGATATGGCACCTGTG'), ('TGAAAGTATGTCATCGAGTATATTCTCTTTCAATTCGTTGTGCTTCGCGTGTCTGATGTCAACGCCCGGG', 'GGCACCTATATCATCGATCCTGCAGTCTTGAGGTTCGTTGTTTATAGCCCTTATTGATATGCCGCCTTGC', 'TGCAAGTATGTTACGGACTTTATTGTAGTCAAATTCGTTGTGCTTCGCGTGTCTGTTGTCAACGTCTTTG', 'GGCAACCAGATCCATGATCATGCGGTAAAGAAGATGCCTGGGTCTGTACCTACTGGACACGGCGTCGCGG'), ('TGTAACTTTGTCACCGATCATGAGGTAAAGAAGATTCTTGGGTCTGTACTTTCTGGACATGGCACCTGGC', 'TGCAAGTTGGTCAATAATTATATACTCTTTACTTCCGCTGTGCTCTGCGCTTATGGACACGGCGTCGCGG', 'TGCAAGTTTGTCCATAATTTTCAACTGTTGCAATTCGCTGTGTATTAAGAATCTGATGTCGGCGTCGTGC', 'TGTATGTCTTTGACGGAGGCCAAGCTGCTTAAATTAGCTGCATACAACGCCTCTGATGTTGGCGTCTGTG'), ('CGCACCGAAAGTATACATTCTTTTCTGAATGCATCCGTTTAGGACTTTCAGTCTCGTAAAGGGCCCAGGC', 'TGCACCGAAGCGATACATTGTTATCTCAGATTGTCCTTGTAGTATGTTAAGTCTCGTTAAGGGGACAAGC', 'TGCACCGAAGCGATACATTGTTATGTTAAAGGATCCTCTGGGCACTGTAAGCCAGGCTAAGGGTGCAGGC', 'TAAGCGCAAGCTTAAAAGCCGTCTATGCTGAACTCCGTTTAACATTGTAAGCCACGGTGAGGGCGCGGGA')))

    # print(test_seqs)
    
    test_bimarkers = find_bimarkers(test_seqs)

    print(test_bimarkers)

