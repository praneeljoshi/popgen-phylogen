# pop-gen-vs-phylo
Code for the pop gen vs. phylo research project. 

We've built out a pipeline to automatically generate random admixture networks; generate multiple sequence aligment and biallelic marker data based on these networks; run GTmix, TreeMix, Structure, and PhyloNet on this data; and finally evaluate the results of each method by computing relevant statistics. 

## Python Dependencies
* `networkx`
* `dendropy`
* `newick`
* `pydot` (optional, needed to output PDFs)

These can all be installed easily with pip (or pip3, depending on your Python installation): `pip install networkx dendropy newick pydot`.

## Other Dependencies
* [`ms`](http://home.uchicago.edu/rhudson1/source/mksamples.html)
* [`Seq-Gen`](http://tree.bio.ed.ac.uk/software/seqgen/)
* [`GTmix`](https://github.com/yufengwudcs/GTmix)
* [`PhyloNet`](https://bioinfocs.rice.edu/phylonet)
* [`TreeMix`](https://bitbucket.org/nygcresearch/treemix/wiki/Home)
* [`Structure`](https://web.stanford.edu/group/pritchardlab/structure.html)
* [`Graphviz`](https://graphviz.org/) (optional, needed to output PDFs)
