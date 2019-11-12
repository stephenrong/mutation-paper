# Mutational bias and the protein-code shape the evolution of splicing enhancers

Stephen Rong, Luke Buerer, Christy L. Rhine, Jing Wang, Kamil J. Cygan, and William G. Fairbrother

Scripts for running the simulations as described in the paper. Simulations are initialized on a 'genome' of random sequences or human exon sequences. Substitutions are drawn in proportion to estimated relative mutation (ERM) rates from Carlson et al. (Nature Communiations 2018) based on 7-mer sequence context. Simulations are run under models with (analogous to exons) or without purifying selection against changes to the protein coding sequence (analogous to introns).

Scripts written by Stephen Rong (PhD Candidate, Fairbrother Lab, Brown University). Have questions? Contact stephen[underscore]rong[at]brown[dot]edu or post a git issue.

Last updated: Oct 28, 2019

### Contents:

data/ contains data files required for performing simulations with the exception of the file data/hg19-unipAliSwissprot-introns.txt.gz (292 Mb), which can instead be downloaded and gzipped from the UCSC Table Browser (GRCh37/hg19 assembly, Genes and Gene Predictions group, UniProt track, SwissProt Aln. table, introns only, accessed on 2018/07/05).

results/ contains simulation output:

- figures/ contains plots of rescaled ERM rates generated by scripts/mu_sims_matrix.R

- simulations/ is a placeholder for simulation output, and includes example simulation output from running scripts/temp_scripts/test.sh.

scripts/ contains the following files:

- mu_sims_matrix.R is used to rescale ERM rates from Carlson et al. (Nat. Comm. 2018) in order to vary overall levels of mutational bias as described in the paper. Output already included in data/.

- mu_sims_module.py is a module containing functions for loading ERM rates, EI scores, Rosenberg scores, and Grantham scores, initializing random sequences, initializing from user-supplied sequences, running the simulations, and tracking summary statistics.

- mu_sims_random.py, mu_sims_exon.py, mu_sims_intron.py specify the simulations initialized with random sequences, human exonic sequences, and human intronic sequences, respectively (mu_sims_intron.py is only used to compute genome-wide means for introns).

- mu_sims_exon_sh.py, mu_sims_intron_sh.py, mu_sims_random_sh.py are used to generate slurm job files in scripts/temp_scripts/, splitting runs into multiple jobs, and also to generate mu_sims_exon.sh, mu_sims_intron.sh, mu_sims_random.sh for submitting all jobs. Run \*\_sh.py and then \*.sh to output simulations to results/simulations/.

- mu_sims_random_scaled_0.py, mu_sims_random_scaled_50.py, mu_sims_random_scaled_100.py, mu_sims_random_scaled_200.py, mu_sims_random_scaled_0_sh.py, mu_sims_random_scaled_50_sh.py, mu_sims_random_scaled_100_sh.py, mu_sims_random_scaled_200_sh.py are variations of the mu_sims_random.py and mu_sims_random_sh.py files for vaying levels of mutational bias. Run \*\_sh.py and then \*.sh to output simulations to results/simulations/.

- temp_scripts/test.sh produces example simulation output. Output already included in results/simulations/.

- analytical_model_double.R solves for the mutation-selection balance of exon vs intron motif ratio for randomly sampled mutation rate matrices with different mutational bias, and for different levels of purifying selection, and plots Supplementary Fig. 7.

### Dependencies:
Python (>=2.7.14), with NumPy (>=1.14.2), pandas (>=0.22.0), and Biopython (>=1.68)

R (>=3.4.4), with tidyverse (>=1.2.1), ggthemr (>=1.1.0)

Simulations were performed on the Fairbrother Lab server running slurm.
