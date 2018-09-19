# mutation-paper

Scripts for performing and analyzing evolutionary simulations in manuscript titled "Exonic splicing enhancers are highly mutable and maintained by selection on the protein code" by Stephen Rong\*, Christy L. Rhine\*, Jing Wang, Kamil J. Cygan, Luke Buerer, and William G. Fairbrother (\*contributed equally to this work). Simulations are based on either random nucleotide sequences or human protein-coding exonic sequences evolved under mutation rate profiles with varying levels of mutational bias and purifying selection scenarios with varying levels of protein-coding constraint on the translated sequence. Details of simulations are described in manuscript.

Custom scripts written by Stephen Rong (Fairbrother Lab, Brown University). Have questions? Contact stephen[underscore]rong[at]brown[dot]edu or post a git issue.

### Contents:

data/ folder contains all required data files for performing simulations and supporting files for plotting results, with the exception of the file data/hg19-unipAliSwissprot-introns.txt (too large at 974 Mb), which can be downloaded from the UCSC Table Browser (GRCh37/hg19 assembly, Genes and Gene Predictions group, UniProt track, SwissProt Aln. table, introns only, accessed on 2018/07/05).

results/ folder contains placeholder directories for simulations and figures.

scripts/ folder contains the following scripts:

- mu_sims_matrix.R is used to preprocess ERM rates from Carlson et al. Nat. Comm. 2018 into mutation rate profiles with varying levels of mutational bias. Output files are already included in data/.

- mu_sims_module.py is a module containing functions for loading requisite data of ERM rates, EI scores, and Grantham scores, generating random sequences, performing simulations, and tracking summary statistics.

- mu_sims_random.py, mu_sims_exon.py, mu_sims_intron.py specify the simulations based on random sequences, human exonic sequences, and human intronic sequences (only used to compute intronic mean ERM rate and mean EI score).

- mu_sims_exon_sh.py, mu_sims_intron_sh.py, mu_sims_random_sh.py are used to generate Slurm job files that partition simulation iterations into many jobs, and also generate mu_sims_exon.sh, mu_sims_intron.sh, mu_sims_random.sh files for submitting all jobs to Slurm. Run \*\_sh.py files and then \*.sh files to output simulation statistics to results/simulations/.

- mu_sims_random_scaled_0.py, mu_sims_random_scaled_50.py, mu_sims_random_scaled_100.py, mu_sims_random_scaled_200.py, mu_sims_random_scaled_0_sh.py, mu_sims_random_scaled_50_sh.py, mu_sims_random_scaled_100_sh.py, mu_sims_random_scaled_200_sh.py are variations on the above mu_sims_random.py and mu_sims_random_sh.py scripts used to generate results for varying levels of mutational bias. Run \*\_sh.py files and then \*.sh files to output simulation statistics to results/simulations/.

- mu_sims_figures.R is used to generate all simulation related figures in manuscript. mu_sims_intermediate.R is used to generate the figure for intermediate constraint simulations. Run these files to output simulation related figures to results/figures/.

### Dependencies:
Python (>=2.7.14), with NumPy (>=1.14.2), pandas (>=0.22.0), and Biopython (>=1.68)

R (>=3.4.4), with tidyverse (>=1.2.1), ggthemr (>=1.1.0), patchwork (>=0.0.1), biomaRt (>=2.34.2), seqinr (>=3.4-5), and SDMTools (>=1.1-221)

Simulations were performed on a HPC cluster running Slurm. Figures were generated on a MacBook Pro (Early 2015) running High Sierra (10.13.1).

### License:

TBD

