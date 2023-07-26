import os
import sys
import pickle
from run_GTmix import run_GTmix
from run_PhyloNet import run_PhyloNet
from run_TreeMix import run_TreeMix

# $env:SLURM_PROCID = '0'
rank = int(os.environ.get("SLURM_PROCID"))

task = sys.argv[1]

base = sys.argv[2] if len(sys.argv) > 2 else "/home/ehs3/pop-gen-vs-phylo/data"

mix_count = int(sys.argv[3]) if len(sys.argv) > 3 else 1

n = int(sys.argv[4]) if len(sys.argv) > 4 else 1

for i in range(rank * n, (rank + 1) * n):

    print(f'Running with rank {rank} and task {task} on input {i} if input exists...')

    if task == "gtmix" and os.path.exists(f'{base}/input/gtmix/{i}/'):
        inferred_network = run_GTmix(f'{base}/input/gtmix/{i}/', 10, mix_count, treepicker_prefix="/home/ehs3/pop-gen-vs-phylo/bin/gtmix/treepicker", gtmix_prefix="/home/ehs3/pop-gen-vs-phylo/bin/gtmix/gtmix", rent_prefix="java -jar /home/ehs3/pop-gen-vs-phylo/bin/gtmix/RentPlus.jar", output=f'{base}/input/gtmix/{rank}/optimal-network.gml') # Be careful with the GTmix output file location
        pickle.dump(inferred_network, open(f'{base}/output/gtmix/{i}.p', 'wb'))
    elif task == "treemix" and os.path.exists(f'{base}/input/treemix/{i}.gz'):
        inferred_network = run_TreeMix(f'{base}/input/treemix/{i}.gz', mix_count, treemix_prefix="/home/ehs3/pop-gen-vs-phylo/bin/treemix/treemix", output=f'{base}/input/treemix/out_stem{rank}')
        pickle.dump(inferred_network, open(f'{base}/output/treemix/{i}.p', 'wb'))
    elif task == "phylonet_mcmc_bimarkers" and os.path.exists(f'{base}/input/phylonet_mcmc_bimarkers/{i}.nex'):
        inferred_network = run_PhyloNet(f'{base}/input/phylonet_mcmc_bimarkers/{i}.nex', phylonet_prefix="/home/ehs3/pop-gen-vs-phylo/bin/phylonet/PhyloNet_3.8.2.jar")
        pickle.dump(inferred_network, open(f'{base}/output/phylonet_mcmc_bimarkers/{i}.p', 'wb'))
    elif task == "phylonet_mle_bimarkers" and os.path.exists(f'{base}/input/phylonet_mle_bimarkers/{i}.nex'):
        inferred_network = run_PhyloNet(f'{base}/input/phylonet_mle_bimarkers/{i}.nex', phylonet_prefix="/home/ehs3/pop-gen-vs-phylo/bin/phylonet/PhyloNet_3.8.2.jar")
        pickle.dump(inferred_network, open(f'{base}/output/phylonet_mle_bimarkers/{i}.p', 'wb'))