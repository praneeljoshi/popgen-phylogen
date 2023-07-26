import os
import networkx as nx

def run_GTmix(directory, trees_per_locus, admixture_count, outgroup=None, max_trees=500, rent_prefix="java -jar RentPlus.jar", treepicker_prefix="./treepicker-linux64", gtmix_prefix="./gtmix-linux64"):
    """
    Runs GTmix on an input directory, and returns the network that it produces.
    Note that some of the command prefixes may have to be altered depending on your OS.
    """
    # Initalize variables
    subdir_count = len(next(os.walk(directory))[1])
    all_chosen_trees = ""

    # Loop through locus directories, computing and selecting gene trees for each one
    for i in range(subdir_count):
        print(i)
        os.system(f'{rent_prefix} {directory}/{i}/locus-{i}.hap')
        os.system(f'{treepicker_prefix} {directory}/{i}/locus-{i}.hap.trees {directory}/{i}/locus-{i}.hap {trees_per_locus} > {directory}/{i}/locus-{i}.hap.trees.chosen')
        with open(f'{directory}/{i}/locus-{i}.hap.trees.chosen', 'r') as f:
            data = f.read()
            if not data.startswith("Can not open gene tree"):
                all_chosen_trees += data

    # Write out all of the chosen trees to an input file
    with open(f'{directory}/locus-all.trees.chosen', 'w') as f:
        f.write(all_chosen_trees)

    # Call GTmix appropriately
    if outgroup is None:
        os.system(f'{gtmix_prefix} -n {admixture_count} -T {max_trees} -P {directory}/listPopInfo-all.txt {directory}/locus-all.trees.chosen')
    else:
        os.system(f'{gtmix_prefix} -n {admixture_count} -T {max_trees} -P {directory}/listPopInfo-all.txt -r {outgroup} {directory}/locus-all.trees.chosen')

    # Read in output network
    inferred_network = nx.read_gml("optimal-network.gml", label="id")

    # Modify network so that it is consistent with the input network
    for node, attributes in inferred_network.nodes.items():
        if inferred_network.in_degree(node) == 0: # Root node
            attributes["type"] = "root"
        elif attributes["label"] == " ": # Internal node
            attributes["type"] = "internal"
        elif attributes["label"].startswith("Mix"): # Admixture node
            attributes["type"] = "admixture"
            attributes["proportion"] = float(attributes["label"][5:])
        else: # Leaf node
            attributes["type"] = "leaf"
            attributes["population"] = int(attributes["label"])
        attributes.pop("label", None)
        attributes.pop("defaultAtrribute", None)
        attributes.pop("vgj", None)

    for edge, attributes in inferred_network.edges.items():
        attributes["weight"] = float(attributes["label"])
        attributes.pop("label", None)

    return inferred_network

if __name__ == "__main__":
    inferred_network = run_GTmix("GTmix_input", 20, 1, treepicker_prefix="wsl ./treepicker-linux64", gtmix_prefix="wsl ./gtmix-linux64")

    test_viz = nx.nx_pydot.to_pydot(inferred_network)
    test_viz.write_pdf("inferred.pdf")
    print(test_viz)
