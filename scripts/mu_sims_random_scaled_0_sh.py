#!/usr/bin/env python

with open("mu_sims_random_scaled_0.sh", "w+") as g:
	g.write("#!/bin/bash\n")
	for i in range(1, 26):
		with open("./temp_scripts/mu_sims_random_scaled_0_" + str(i*200) + ".sh", "w+") as f:
			f.write("#!/bin/bash\n")
			f.write("\n")
			f.write("#SBATCH -n 1\n")
			f.write("time python mu_sims_random_scaled_0.py " + str((i-1)*200) + " " + str(i*200))
		g.write("sbatch ./temp_scripts/mu_sims_random_scaled_0_" + str(i*200) + ".sh\n")
