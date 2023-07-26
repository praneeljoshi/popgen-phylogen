import collections
import os
import re

# Eventually you should seperate out the MrBayes functions into another file

def build_bimarker_nexus(data):
    """
    Given a data tuple returned from call_ms, returns a multiline string corrsponding to the NEXUS file format representation of that data.
    Also returns the taxa and taxon map corresponding to the data.
    """
    # data[locus][locations/haplotypes][population][individual]

    # Count number of individuals, and total length of each individual's bitstrings
    ntax = sum(len(x) for x in data[0][1])
    nchar = sum(len(x[1][0][0]) for x in data)

    # Build the beginning of the string
    nexus = (
        f'#NEXUS\n'
        f'Begin data;\n'
        f'Dimensions ntax={ntax} nchar={nchar};\n'
        f'Format datatype=dna symbols="012" missing=? gap=-;\n'
        f'Matrix\n'
    )

    taxa = []
    taxon_map = collections.defaultdict(list)

    # Build data lines by iterating through loci, holding individuals constant; also records information in taxa and taxon_map
    for pop_idx in range(len(data[0][1])):
        for indiv_idx in range(len(data[0][1][0])):
            taxon_str = f'I{pop_idx + 1}-{indiv_idx + 1} '
            taxa.append(taxon_str.strip())
            taxon_map[pop_idx + 1].append(taxon_str.strip())
            for locus in data:
                taxon_str += locus[1][pop_idx][indiv_idx]
            nexus += taxon_str + "\n"

    # Add the end of the NEXUS block
    nexus += ";End;\n"

    return nexus, taxa, dict(taxon_map)

def PhyloNet_taxa_loci_string(taxa_or_loci):
    """
    Returns the PhyloNet taxa/loci string corresponding to a list of taxa/loci.
    """
    return f'({",".join(taxa_or_loci)})'

def PhyloNet_taxon_map_string(taxon_map):
    """
    Returns the PhyloNet taxon map string string corresponding to a taxon map.
    """
    return f'<{";".join(str(pop) + ":" + ",".join(lst) for pop, lst in taxon_map.items())}>'

def build_MCMC_BiMarkers_input(nexus, taxa, taxon_map, max_reticulation=1, chain_length=500000, burn_in_length="200000", sample_freq=500, seed=12345678, threads=None):
    """
    Given the results of build_bimarker_nexus, returns a multiline string representing the MCMC_BiMarkers input corresponding to that data.
    """
    # Use the NEXUS file as a starting point
    input_str = nexus

    # Convert taxa and taxon_map to strings
    taxa_string = PhyloNet_taxa_loci_string(taxa)
    taxon_map_string = PhyloNet_taxon_map_string(taxon_map)

    # Build the rest of the input based on the given parameters
    input_str += (
        f'BEGIN PHYLONET;\n'
        f'MCMC_BiMarkers -cl {chain_length} -bl {burn_in_length} -sf {sample_freq} -mr {max_reticulation}{ "-pl "+ str(threads) if threads is not None else ""}\n'
        f'-sd {seed}\n'
        f'\n'
        f'-taxa {taxa_string}\n'
        f'-tm {taxon_map_string};\n'
        f'\n'
        f'END;\n'
    )

    return input_str

def build_MLE_BiMarkers_input(nexus, taxon_map, max_reticulation=1, max_runs=100, max_examinations=50000, num_optimums=10, max_failures=50, pseudo=True, seed=12345678, threads=None):
    """
    Given the results of build_bimarker_nexus, returns a multiline string representing the MLE_BiMarkers input corresponding to that data.
    """
    # Use the NEXUS file as a starting point
    input_str = nexus

    # Convert taxon_map to string
    taxon_map_string = PhyloNet_taxon_map_string(taxon_map)

    # Build the rest of the input based on the given parameters
    input_str += (
        f'BEGIN PHYLONET;\n'
        f'MLE_BiMarkers -mnr {max_runs} -mec {max_examinations} -mno {num_optimums} -mf {max_failures} {"-pseudo" if pseudo else ""} -mr {max_reticulation}{" -pl "+ str(threads) if threads is not None else ""} -sd {seed} -tm {taxon_map_string};\n'
        f'END;\n'
    )

    return input_str

def build_alignment_nexus(data):
    """
    Given a data tuple returned from call_seq_gen, returns a multiline string corrsponding to the NEXUS file format representation of that data.
    Also returns the taxon map and loci list corresponding to the data.
    """
    # data[locus][population][individual]

    # Count number of individuals, and total length of each individual's bitstrings
    ntax = sum(len(x) for x in data[0])
    nchar = sum(len(x[0][0]) for x in data)

    # Build the beginning of the string
    nexus = (
        f'#NEXUS\n'
        f'Begin data;\n'
        f'Dimensions ntax={ntax} nchar={nchar};\n'
        f'Format datatype=dna symbols="ACTG" missing=? gap=- interleave=yes;\n' # interleave=yes might cause problems with PhyloNet?
        f'Matrix\n'
    )

    loci = []
    taxon_map = collections.defaultdict(list)

    # Build data lines by iterating through loci, populations, and individuals; also records information in loci and taxon_map
    for locus_idx, locus in enumerate(data):
        locus_str = f'L{locus_idx}'
        loci.append(locus_str)
        locus_len = len(locus[0][0])
        nexus += f'[{locus_str}, {locus_len}]\n'
        for pop_idx, population in enumerate(locus):
            for indiv_idx, individual in enumerate(population):
                taxon_str = f'I{pop_idx + 1}-{indiv_idx + 1}'
                nexus += f'{taxon_str} {individual}\n'
                if taxon_str not in taxon_map[pop_idx + 1]:
                    taxon_map[pop_idx + 1].append(taxon_str)

    # Add the end of the NEXUS block
    nexus += ";End;\n"

    return nexus, loci, dict(taxon_map)

def build_MCMC_SEQ_input(nexus, loci, taxon_map, max_reticulation=4, chain_length=10000000, burn_in_length=2000000, sample_freq=5000, seed=12345678, out_directory=None, threads=None):
    """
    Given the results of build_alignment_nexus, returns a multiline string representing the MCMC_SEQ input corresponding to that data.
    """
    # Use the NEXUS file as a starting point
    input_str = nexus

    # Remove the " interleave=yes" string from the prebuilt NEXUS file
    input_str = input_str.replace(" interleave=yes", "")

    # Convert loci and taxon_map to strings
    loci_string = PhyloNet_taxa_loci_string(loci)
    taxon_map_string = PhyloNet_taxon_map_string(taxon_map)

    # Build the rest of the input based on the given parameters
    input_str += (
        f'BEGIN PHYLONET;\n'
        f'MCMC_SEQ -cl {chain_length} -bl {burn_in_length} -sf {sample_freq} -mr {max_reticulation}{" -pl "+ str(threads) if threads is not None else ""} -sd {seed} -tm {taxon_map_string}{" -pl "+ str(threads) if threads is not None else ""}{" -dir "+ out_directory if out_directory is not None else ""};\n'
        f'END;\n'
    )

    return input_str

def build_InferNetwork_helper(command, nexus, taxon_map, loci_spec_string, max_reticulation, runs, threads):
    """
    Given the results of build_alignment_nexus (taxon_map) and run_MrBayes (nexus, loci_spec_string), returns a multiline string representing the InferNetwork input corresponding 
    to that data and the specified command.
    """
    # Use the NEXUS file as a starting point
    input_str = nexus

    # Convert taxon_map to string
    taxon_map_string = PhyloNet_taxon_map_string(taxon_map)

    # Build the rest of the input based on the given parameters
    input_str += (
        f'BEGIN PHYLONET;\n'
        f'{command} {loci_spec_string} {max_reticulation} -a {taxon_map_string} -x {runs} -pl {threads};\n'
        f'END;\n'
    )

    return input_str

def build_InferNetwork_MP_input(nexus, taxon_map, loci_spec_string, max_reticulation=1, runs=5, threads=1):
    """
    Given the results of build_alignment_nexus (taxon_map) and run_MrBayes (nexus, loci_spec_string), returns a multiline string representing the InferNetwork_MP input corresponding 
    to that data.
    """
    return build_InferNetwork_helper("InferNetwork_MP", nexus, taxon_map, loci_spec_string, max_reticulation, runs, threads)

def build_InferNetwork_ML_input(nexus, taxon_map, loci_spec_string, max_reticulation=1, runs=5, threads=1):
    """
    Given the results of build_alignment_nexus (taxon_map) and run_MrBayes (nexus, loci_spec_string), returns a multiline string representing the InferNetwork_ML input corresponding 
    to that data.
    """
    return build_InferNetwork_helper("InferNetwork_ML", nexus, taxon_map, loci_spec_string, max_reticulation, runs, threads)

def build_InferNetwork_MPL_input(nexus, taxon_map, loci_spec_string, max_reticulation=1, runs=5, threads=1):
    """
    Given the results of build_alignment_nexus (taxon_map) and run_MrBayes (nexus, loci_spec_string), returns a multiline string representing the InferNetwork_MPL input corresponding 
    to that data.
    """
    return build_InferNetwork_helper("InferNetwork_MPL", nexus, taxon_map, loci_spec_string, max_reticulation, runs, threads)

def build_MCMC_GT_input(nexus, taxon_map, loci_spec_string, max_reticulation=None, chain_length=1100000, burn_in_length=100000, sample_freq=1000, seed=12345678, pseudo=False, threads=1):
    """
    Given the results of build_alignment_nexus (taxon_map) and run_MrBayes (nexus, loci_spec_string), returns a multiline string representing the MCMC_GT input corresponding 
    to that data.
    """
    # Use the NEXUS file as a starting point
    input_str = nexus

    # Convert taxon_map to string
    taxon_map_string = PhyloNet_taxon_map_string(taxon_map)

    # Build the rest of the input based on the given parameters
    input_str += (
        f'BEGIN PHYLONET;\n'
        f'MCMC_GT {loci_spec_string}{" -mr " + str(max_reticulation) if max_reticulation is not None else ""} -cl {chain_length} -bl {burn_in_length} -sf {sample_freq} -sd {seed} -tm {taxon_map_string} -pl {threads}{" -pseudo " if pseudo else ""};\n'
        f'END;\n'
    )

    return input_str

def write_input(input_str, file):
    """
    Writes out an input NEXUS string to a file.
    """
    with open(file, "w") as f:
        f.write(input_str)
    return

if __name__ == "__main__":
    import networkx as nx
    from admixture_network import generate_admixture_networks
    from call_ms import call_ms
    from call_seq_gen import call_seq_gen
    from call_MrBayes import build_MrBayes_input, run_MrBayes

    # Generate test admixture network
    # test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500000, True, 1)
    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 70, True, 1)

    test_network, test_ms_cmd = test_results[0]

    test_viz = nx.nx_pydot.to_pydot(test_network)
    test_viz.write_pdf("test.pdf")

    # Call ms to generate haplotypes and trees
    test_data, test_tree_lines = call_ms(test_ms_cmd, True)

    # Call Seq-Gen to generate sequence alignments
    test_seqs = call_seq_gen(test_tree_lines, 70, 4, 4, seq_gen_prefix=".\seq-gen.exe")

    # Generate base bimarker nexus, taxa list, and taxon map
    bimarker_nexus, bimarker_taxa, bimarker_taxon_map = build_bimarker_nexus(test_data)

    # Generate base alignment nexus, loci list, and taxa map
    alignment_nexus, alignment_loci, alignment_taxon_map = build_alignment_nexus(test_seqs)

    # Generate bimarker PhyloNet inputs
    MCMC_BiMarkers_input = build_MCMC_BiMarkers_input(bimarker_nexus, bimarker_taxa, bimarker_taxon_map)
    MLE_BiMarkers_input = build_MLE_BiMarkers_input(bimarker_nexus, bimarker_taxon_map)

    # Generate direct sequence inference PhyloNet input
    MCMC_SEQ_input = build_MCMC_SEQ_input(alignment_nexus, alignment_loci, alignment_taxon_map)

    # Generate MrBayes input for inferring gene trees from sequences
    MrBayes_input = build_MrBayes_input(alignment_nexus, alignment_loci, alignment_taxon_map)

    # Write MrBayes input to file
    with open("MrBayes\\test.nex", "w") as f:
        f.write(MrBayes_input)

    # Run MrBayes on input
    trees_nexus, trees_loci_spec = run_MrBayes("MrBayes\\test.nex", ".\mb.exe", 50)

    # Generate gene tree PhyloNet inputs
    InferNetwork_MP_input = build_InferNetwork_MP_input(trees_nexus, alignment_taxon_map, trees_loci_spec)
    InferNetwork_ML_input = build_InferNetwork_ML_input(trees_nexus, alignment_taxon_map, trees_loci_spec)
    InferNetwork_MPL_input = build_InferNetwork_MPL_input(trees_nexus, alignment_taxon_map, trees_loci_spec)
    MCMC_GT_input = build_MCMC_GT_input(trees_nexus, alignment_taxon_map, trees_loci_spec)

    # Write out test PhyloNet inputs
    write_input(MCMC_BiMarkers_input, "test_inputs/MCMC_BiMarkers_input.nex")
    write_input(MLE_BiMarkers_input, "test_inputs/MLE_BiMarkers_input.nex")
    write_input(MCMC_SEQ_input, "test_inputs/MCMC_SEQ_input.nex")
    write_input(InferNetwork_MP_input, "test_inputs/InferNetwork_MP_input.nex")
    write_input(InferNetwork_ML_input, "test_inputs/InferNetwork_ML_input.nex")
    write_input(InferNetwork_MPL_input, "test_inputs/InferNetwork_MPL_input.nex")
    write_input(MCMC_GT_input, "test_inputs/MCMC_GT_input.nex")