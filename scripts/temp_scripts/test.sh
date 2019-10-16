#!/bin/bash

#SBATCH -n 1
time python mu_sims_exon.py 0 2
time python mu_sims_intron.py 0 2
time python mu_sims_random.py 0 2
time python mu_sims_random_scaled_0.py 0 2
time python mu_sims_random_scaled_50.py 0 2
time python mu_sims_random_scaled_100.py 0 2
time python mu_sims_random_scaled_200.py 0 2
