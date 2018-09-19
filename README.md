# mutation-paper

Scripts for performing and analyzing evolutionary simulations in manuscript titled "Exonic splicing enhancers are highly mutable and maintained by selection on the protein code" by Stephen Rong\*, Christy L. Rhine\*, Jing Wang, Kamil J. Cygan, Luke Buerer, and William G. Fairbrother (\*contributed equally to this work). Simulations are based on either random nucleotide sequences or human protein-coding exonic sequences evolved under mutation rate profiles with varying levels of mutational bias and purifying selection scenarios with varying levels of protein-coding constraint on the translated sequence. 

Custom scripts written by Stephen Rong (Fairbrother Lab, Brown University). Have questions? Contact stephen[underscore]rong[at]brown[dot]edu or post a git issue.

### Contents:

data/ folder contains all required data files for performing simulations and supporting files for plotting results, with the exception of the file data/hg19-unipAliSwissprot-introns.txt (974 Mb), which can be downloaded from the UCSC Table Browser (GRCh37/hg19 assembly, Genes and Gene predictions, UniProt, SwissProt Aln., introns only, accessed on 2018/07/05).

results/ folder contains placeholder directories for simulations and figures.

scripts/ folder contains the following scripts:



### Dependencies:
Python (>=2.7.14), NumPy (>=1.14.2), pandas (>=0.22.0), Biopython (>=1.68)

R (>=3.4.4), tidyverse (>=1.2.1), ggthemr (>=1.1.0), patchwork (>=0.0.1), biomaRt (>=2.34.2), seqinr (>=3.4-5), SDMTools (>=1.1-221)

Simulations were performed on a HPC cluster with Slurm. Figures were generated on a MacBook Pro (Early 2015) running High Sierra (10.13.1).
