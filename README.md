# mutation-paper

Scripts for performing and analyzing evolutionary simulations in paper titled "Mutational bias and the co-evolution of protein and splicing code" by Stephen Rong, Luke Buerer, Christy L. Rhine, Jing Wang, Kamil J. Cygan, and William G. Fairbrother. Simulations are initialized with either random sequences or human exonic sequences evolved under mutation rate profiles with varying levels of mutational bias and purifying selection. Details of simulations are described in paper.

Scripts written by Stephen Rong (Fairbrother Lab, Brown University). Have questions? Contact stephen[underscore]rong[at]brown[dot]edu or post a git issue.

Last updated: Oct 16, 2019

### Contents:

data/ contains all required data files for performing simulations.

results/ contains placeholder directories for simulation output.

scripts/ contains the following:

- mu_sims_matrix.R is used to preprocess ERM rates from Carlson et al. Nat. Comm. (2018) into mutation rate profiles with varying levels of mutational bias. Outputs are already included in data/.

- mu_sims_module.py is a module containing functions for loading  ERM rates, EI scores, and Grantham scores; generating random sequences; performing simulations; and tracking summary statistics.

- mu_sims_random.py, mu_sims_exon.py, mu_sims_intron.py specify the simulations initialized with random sequences, human exonic sequences, and human intronic sequences (only used to compute genome-wide means for introns), respectively.

- mu_sims_exon_sh.py, mu_sims_intron_sh.py, mu_sims_random_sh.py are used to generate Slurm job files in scripts/temp_scripts/ that split runs into many jobs, and also to generate mu_sims_exon.sh, mu_sims_intron.sh, mu_sims_random.sh for submitting jobs. Run \*\_sh.py and then \*.sh to output simulations to results/simulations/.

- mu_sims_random_scaled_0.py, mu_sims_random_scaled_50.py, mu_sims_random_scaled_100.py, mu_sims_random_scaled_200.py, mu_sims_random_scaled_0_sh.py, mu_sims_random_scaled_50_sh.py, mu_sims_random_scaled_100_sh.py, mu_sims_random_scaled_200_sh.py are variations on the above mu_sims_random.py and mu_sims_random_sh.py used to generate results for varying levels of mutational bias. Run \*\_sh.py and then \*.sh to output simulations to results/simulations/.

- scripts/temp_scripts/test.sh ... 

- results/simulations/\*_0_2_\* files ... 

### Dependencies:
Python (>=2.7.14), with NumPy (>=1.14.2), pandas (>=0.22.0), and Biopython (>=1.68).

R (>=3.4.4), with tidyverse (>=1.2.1), ggthemr (>=1.1.0).

Simulations were performed on a HPC cluster running Slurm.
