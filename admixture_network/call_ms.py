import os
import networkx as nx
from admixture_network import generate_admixture_networks

def call_ms(cmd, trees=False):
    """
    Calls ms, given a command and two parameters that the command contains.
    """
    # Extract relevant parameters from command
    pop_count = int(cmd.split()[9]) # Includs the outgroup
    alleles_per_pop = int(cmd.split()[1]) // pop_count

    # Run ms (with an additional argument requesting gene trees added if necessary), and store its result (excluding the first two lines)
    stream = os.popen(cmd + (" -T" if trees else ""))
    output_lines = list([x.strip() for x in stream.readlines()])[3:]

    # Seperate out gene trees if necessary
    if trees:
        tree_lines = []
        new_output_lines = []
        for line in output_lines:
            if line.startswith("[") and line.endswith(";"):
                tree_lines.append(line)
            else:
                new_output_lines.append(line)
        output_lines = new_output_lines

    # Determine step for breaking output into locus chunks
    step = pop_count * alleles_per_pop + 4

    output = []

    for i in range(0, len(output_lines), step):
        # Filter out empty lines
        chunk = [x for x in output_lines[i+2:i+step] if x != ""]

        # Take the first line, split it, and store it as the locus's positions
        positions = tuple([float(x) for x in chunk.pop(0).split()[1:]])

        alleles = []

        # Break the remainder of the chunk into groups of alleles_per_pop lines
        for j in range(0, len(chunk), alleles_per_pop):
            alleles.append(tuple(chunk[j:j+alleles_per_pop]))

        output.append((positions, tuple(alleles)))

    # Return appropriate results
    if trees:
        return tuple(output), tuple(tree_lines)
    else:
        return tuple(output)

if __name__ == "__main__":
    test_results = generate_admixture_networks(4, 1, 0.02, 0.3, 0.5, 4, 2, 50, 50, 500000, True, 1)

    test_network, test_cmd = test_results[0]

    # test_viz = nx.nx_pydot.to_pydot(test_network)
    # test_viz.write_pdf("test_ms.pdf")

    test_haplotypes, test_tree_lines = call_ms(test_cmd, True)

    print(test_haplotypes)
    # print(test_tree_lines)
    print(test_cmd)
