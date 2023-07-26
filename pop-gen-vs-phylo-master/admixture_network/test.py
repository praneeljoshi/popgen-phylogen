from admixture_network import generate_admixture_networks
from call_ms import call_ms
from write_GTmix_input import write_GTmix_input
from run GTmix import run_GTmix
from graph_comparison import admixture_error, compute_inference_error

def test(algorithm, pop_count, ad_count, time_interval, outgroup_tbonus, ad_prop, alleles_per_pop, loci_count, mutation,
         recombination, locus_len, extant_admixture, net_count, trees_per_locus):
    """
    takes in algorithm and relevant parameters, and spits out measures of performance
    """
    simulations = generate_admixture_networks(pop_count, ad_count, time_interval, outgroup_tbonus, ad_prop,
                                               alleles_per_pop, loci_count, mutation, recombination, locus_len,
                                               extant_admixture, net_count)
    sim_haps = []
    sim_nets = []
    predictions = []
    bmie = []
    ad_errors = []

    n = 0
    for sim in simulations:
        sim_nets.append = sim[0]
        hap = call_ms(sim[1])
        sim_haps.append(hap)
        if algorithm == "gt":
            write_GTmix_input(hap, "gt" + str(n))
        n += 1

    for run in range(n):
        #run predictions
        if algorithm == "gt":
            predict = run_GTmix("gt" + str(n), trees_per_locus, ad_count)
        predictions.append(predict)

        #error calculations
        bmie.append(compute_inference_error(predict, sim_nets[n]), pop_count)
        ad_errors.append(admixture_error(predict, sim_nets))

    return bmie, ad_errors