#!/usr/bin/env python

with open("mu_sims_exon.sh", "w+") as g:
	g.write("#!/bin/bash\n")
	for i in range(1, 107):
		with open("./temp_scripts/mu_sims_exon_" + str(i*200) + ".sh", "w+") as f:
			f.write("#!/bin/bash\n")
			f.write("\n")
			f.write("#SBATCH -n 1\n")
			f.write("time python mu_sims_exon.py " + str((i-1)*200) + " " + str(i*200))
		g.write("sbatch ./temp_scripts/mu_sims_exon_" + str(i*200) + ".sh\n")

	with open("./temp_scripts/mu_sims_exon_" + str(21328) + ".sh", "w+") as f:
		f.write("#!/bin/bash\n")
		f.write("\n")
		f.write("#SBATCH -n 1\n")
		f.write("time python mu_sims_exon.py " + str(i*200) + " " + str(21328))
	g.write("sbatch ./temp_scripts/mu_sims_exon_" + str(21328) + ".sh\n")
