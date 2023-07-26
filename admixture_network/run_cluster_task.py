import os
import sys
import pickle
from run_GTmix import run_GTmix
from run_PhyloNet import run_PhyloNet
from run_TreeMix import run_TreeMix

base = "/home/ehs3/pop-gen-vs-phylo/data"

# $env:SLURM_PROCID = '0'
rank = int(os.environ.get("SLURM_PROCID"))

task = sys.argv[1]

print(f'Running with rank {rank} and task {task}...')

if task == "gtmix":
    inferred_network = run_GTmix(f'{base}/input/gtmix/{rank}/', 20, 1, treepicker_prefix="/home/ehs3/pop-gen-vs-phylo/bin/gtmix/treepicker", gtmix_prefix="/home/ehs3/pop-gen-vs-phylo/bin/gtmix/gtmix", rent_prefix="java -jar /home/ehs3/pop-gen-vs-phylo/bin/gtmix/RentPlus.jar", output=f'{base}/input/gtmix/{rank}/optimal-network.gml') # Be careful with the GTmix output file location
    pickle.dump(inferred_network, open(f'{base}/output/gtmix/{rank}.p', 'wb'))
elif task == "treemix":
    inferred_network = run_TreeMix(f'{base}/input/treemix/{rank}.gz', 1, treemix_prefix="/home/ehs3/pop-gen-vs-phylo/bin/treemix/treemix", output=f'/home/ehs3/pop-gen-vs-phylo/data/input/treemix/out_stem{rank}')
    pickle.dump(inferred_network, open(f'{base}/output/treemix/{rank}.p', 'wb'))
elif task == "phylonet_mcmc_bimarkers":
    inferred_network = run_PhyloNet(f'{base}/input/phylonet_mcmc_bimarkers/{rank}.nex', phylonet_prefix="/home/ehs3/pop-gen-vs-phylo/bin/phylonet/PhyloNet_3.8.2.jar")
    pickle.dump(inferred_network, open(f'{base}/output/phylonet_mcmc_bimarkers/{rank}.p', 'wb'))
elif task == "phylonet_mle_bimarkers":
    inferred_network = run_PhyloNet(f'{base}/input/phylonet_mle_bimarkers/{rank}.nex', phylonet_prefix="/home/ehs3/pop-gen-vs-phylo/bin/phylonet/PhyloNet_3.8.2.jar")
    pickle.dump(inferred_network, open(f'{base}/output/phylonet_mle_bimarkers/{rank}.p', 'wb'))