"""
Uses the graph_comparison methods to benchmark the (currently two) algorithms
for generating networks

This does mostly the same thing as test.py, so if this gets deleted I won't be offended

The inputs I've written work for my filepaths;
talk to me (Ian) if you have issues changing the paths for your system

"""

from admixture_network import generate_admixture_networks
from call_ms import call_ms
from write_GTmix_input import write_GTmix_input
from write_TreeMix_input import write_TreeMix_input
from run_GTmix import run_GTmix
from run_TreeMix import run_TreeMix
from graph_comparison import *



#range of values for testing                  #currently just fixed values
pop_count = 4
admixture_count = 1
time_interval = 0.02
outgroup_time_bonus = 0.3
admixture_prop = 0.5
alleles_per_pop = 4
loci_count = 2
mutation = 50
recombination = 50
locus_length = 500000
only_extant_admixture = True
n = 10                                      #If you're testing set this to 1 or 2. This takes a while to run
trees_per_locus = 20
max_trees = 500

rent_prefix = "java -jar ~/Desktop/GTmix-master/RentPlus.jar"       #  <--  Change these for yourself
treepicker_prefix = "~/Desktop/GTmix-master/treepicker-mac"         #  <--
gtmix_prefix = "~/Desktop/GTmix-master/gtmix-mac"                   #  <--
treemix_prefix= "treemix"                                           #  <--

tm_output_base= "./tmfiles/output"
base_gtpath = "./gtfiles"
base_tmpath = "./tmfiles"





#create test network
ntwrk_cmd_tpls = generate_admixture_networks(pop_count,
                                            admixture_count,
                                            time_interval,
                                            outgroup_time_bonus,
                                            admixture_prop,
                                            alleles_per_pop,
                                            loci_count,
                                            mutation,
                                            recombination,
                                            locus_length,
                                            only_extant_admixture,
                                            n)

#call ms
ms_outs = [];
for ntwrk, cmd in ntwrk_cmd_tpls:
    ms_outs.append(call_ms(cmd))

#print(ms_outs)


gt_avg_adm_err = 0
tm_avg_adm_err = 0

gt_avg_inf_err = 0
tm_avg_inf_err = 0
#actual benchmarking
for idx, data in enumerate(ms_outs):
    gtpath = base_gtpath + "/" + str(idx)
    tmpath = base_tmpath + "/" + str(idx) + ".gz"
    tm_output = tm_output_base + "/" + str(idx)

    ntwrk = ntwrk_cmd_tpls[idx][0]

    #run TreeMix and GTmix
    write_GTmix_input(data, gtpath)
    gt_ntwrk = run_GTmix(gtpath, trees_per_locus, admixture_count, None, max_trees, rent_prefix, treepicker_prefix, gtmix_prefix)

    write_TreeMix_input(data, tmpath)
    tm_ntwrk = run_TreeMix(tmpath, admixture_count, None, None, tm_output, treemix_prefix)


    if only_extant_admixture: #admixture error
        gt_adm_err = admixture_error(gt_ntwrk, ntwrk)
        tm_adm_err = admixture_error(tm_ntwrk, ntwrk)

        gt_avg_adm_err += gt_adm_err[0]/n
        tm_avg_adm_err += tm_adm_err[0]/n

    #rf error
    gt_inf_err = compute_inference_error(gt_ntwrk, ntwrk, pop_count)
    tm_inf_err = compute_inference_error(tm_ntwrk, ntwrk, pop_count)

    gt_avg_inf_err += gt_inf_err/n
    tm_avg_inf_err += tm_inf_err/n

if only_extant_admixture:
    print(f"GTmix Error: admixture_err = {gt_avg_adm_err}, inference_err = {gt_avg_inf_err}")
    print(f"TreeMix Error: admixture_err = {tm_avg_adm_err}, inference_err = {tm_avg_inf_err}")

else:
    print(f"GTmix Error: inference_err = {gt_avg_inf_err}")
    print(f"TreeMix Error: inference_err = {tm_avg_inf_err}")











#
