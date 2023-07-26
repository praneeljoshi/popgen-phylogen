import os
import re

def build_MrBayes_input(nexus, loci, taxon_map, generations=1000000, chains=1):
    """
    Given the results of build_alignment_nexus, returns a multiline string representing the MrBayes input corresponding to that data.
    """
    # See manual p76
    # Known bug (doesn't affect us, because we only care about gene trees, not the species tree): https://www.gitmemory.com/issue/NBISweden/MrBayes/116/514587164

    # Use the NEXUS file as a starting point
    input_str = nexus

    # Add the beginning of the MrBayes block
    input_str += (
        f'Begin mrbayes;\n'
        f'set autoclose=yes nowarn=yes;\n'
    )

    # Scan for loci delimeters, and build MrBayes CHARSETs from them; also record locus_names as we go
    counter = 0
    locus_names = []
    for line in nexus.splitlines():
        if line.startswith("[") and line.endswith("]"):
            locus_name, locus_length = line[1:-1].split(", ")
            locus_length = int(locus_length)
            input_str += f'CHARSET {locus_name} = {counter + 1}-{counter + locus_length};\n'
            counter += locus_length
            locus_names.append(locus_name)

    # Specify partition information
    input_str += (
        f'partition genes = {len(locus_names)}: {",".join(locus_names)};\n'
        f'set partition = genes;\n'
    )

    # Specify taxon list (map)
    input_str += "speciespartition species = " + ", ".join("P" + str(pop) + ":" + " ".join(lst) for pop, lst in taxon_map.items()) + ";\n"
    input_str += "set speciespartition = species;\n"

    # Unlink partition's topology parameter, invoke the multispecies coalescent, set the model, and run MrBayes
    input_str += (
        f'unlink topology = (all);\n'
        f'prset topologypr = speciestree;\n'
        f'prset brlenspr = clock:speciestree;\n'
        f'lset nst=2 rates=gamma;\n'
        f'mcmc ngen={generations} nchains={chains};\n'
        f'sumt;\n'
        f'End;\n'
    )

    return input_str

def substitute_individual_names(string, map):
    """
    Substitutes individual names into a tree string.
    """
    for key, value in map.items():
        # Only replaces numbers that are in between parentheses or commas
        string = re.sub(f'[(),]{key}[(),]', lambda match: match.group()[0] + value + match.group()[-1], string)
    return string

def run_MrBayes(file, mrbayes_prefix="mb", max_per_locus=float('inf')):
    """
    Runs MrBayes on an input file, and returns a string corresponding to the NEXUS file representation of the trees that it produces, as well as a string
    corresponding to the loci that the different trees represent.
    """
    # Call MrBayes
    os.system(f'{mrbayes_prefix} {file}')

    tree_lines_string = (
        '#NEXUS\n'
        'BEGIN TREES;\n'
    )

    loci_spec_list = []

    # Find parent directory path
    directory = os.path.dirname(file)

    # Find tree file names, and make sure they are in order
    treefiles = [filename for filename in os.listdir(directory) if filename.endswith(".trprobs")] # This would have to be altered to not return the final file if the aformentioned bug in MrBayes is ever fixed
    treefiles.sort(key=lambda filename: int("".join(char for char in filename if char.isdigit())))

    # Iterate through files, reading in trees
    overall_counter = 0 # Counts the number of trees read overall
    for treefile_name in treefiles:
        path = f'{directory}/{treefile_name}'
        with open(path, "r") as treefile:
            substitution_map = {}
            locus_counter = 0 # Counts the number of trees read for this locus
            loci_start = overall_counter
            for line in treefile:
                line = line.strip()
                if len(line) == 0:
                    continue
                elif line[0].isdigit():
                    key, value = line.strip(",;").split()
                    substitution_map[key] = value
                elif line.startswith("tree"):
                    tree_data = line.split()[-3:]
                    tree_data[2] = substitute_individual_names(tree_data[2], substitution_map)
                    tree_string = " ".join(tree_data)
                    tree_lines_string += f'Tree gt{overall_counter} = {tree_string}\n'
                    overall_counter += 1
                    locus_counter += 1
                if locus_counter >= max_per_locus:
                    break
            loci_stop = overall_counter - 1
            loci_spec_list.append(f'{{gt{loci_start}-gt{loci_stop}}}')
    
    tree_lines_string += 'END;\n'

    loci_spec_string = f'({", ".join(loci_spec_list)})'

    return tree_lines_string, loci_spec_string