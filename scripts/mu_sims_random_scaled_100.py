#!/usr/bin/env python

from __future__ import division
import sys
from mu_sims_module import *
from numpy.random import randint

# # # random start
if __name__ == '__main__':
	iter_start = int(sys.argv[1])
	iter_end = int(sys.argv[2])
	np.random.seed(iter_start+1)
	
	# import mut_matrix
	mut_file = "../data/mu-matrix-7mers-scaled_100.txt"
	mut_matrix = get_mut_matrix(mut_file)

	# output files
	run_file = "../results/simulations/random_scaled_100_"+str(iter_start)+"_"+str(iter_end)
	with open(run_file+"_track_seq.txt", "w+") as track_seq, open(run_file+"_track_mut.txt", "w+") as track_mut, open(run_file+"_track_verbose.txt", "w+") as track_verbose:

		# add header
		prefix = "init_cond"+"\t"+"constr_cond"+"\t"+"run_number"
		track_seq.write(prefix+"\t"+"seq_length"+"\t"+"mutation"+"\t"+"mut_scaled"+"\t"+"mutation_tot"+"\t"+"mut_scaled_tot"+"\t"+"sequence"+"\n")
		track_mut.write(prefix+"\t"+"seq_length"+"\t"+"mutation"+"\t"+"mut_scaled"+"\t"+"mutation_tot"+"\t"+"mut_scaled_tot"+"\t"+"mut_mean"+"\t"+"esr_mean"+"\n")
		track_verbose.write(prefix+"\t"+"seq_length"+"mutation"+"\t"+"mut_scaled"+"\t"+"mutation_tot"+"\t"+"mut_scaled_tot"+"\t"+"codonBefore"+"\t" "codonAfter"+"\t"+"wtMotif"+"\t"+"mtMotif"+"\t"+"mutRate"+"\t"+"aminoBefore"+"\t"+"aminoAfter"+"\t"+"aminoCheck"+"\n")

		for i in range(iter_start, iter_end):
			init_cond = "random"
			constr_cond = "neutral"
			run_number = str(i)
			print init_cond+"_"+constr_cond+"_"+run_number
			seq_sequence = initiate_seq_sequence(999)
			muSimsConstraintGranthamSimulation(mut_matrix, track_seq, track_mut, track_verbose, init_cond, constr_cond, run_number, seq_sequence, 9990, 50, 300, False)

		for i in range(iter_start, iter_end):
			init_cond = "random"
			constr_cond = "identity"
			run_number = str(i)
			print init_cond+"_"+constr_cond+"_"+run_number
			seq_sequence = initiate_seq_sequence(999)
			muSimsConstraintGranthamSimulation(mut_matrix, track_seq, track_mut, track_verbose, init_cond, constr_cond, run_number, seq_sequence, 9990, 50, 0, False)

		# for i in range(iter_start, iter_end):
		# 	init_cond = "random"
		# 	constr_cond = "grantham"
		# 	run_number = str(i)
		# 	print init_cond+"_"+constr_cond+"_"+run_number
		# 	seq_sequence = initiate_seq_sequence(999)
		# 	muSimsConstraintGranthamSimulation(mut_matrix, track_seq, track_mut, track_verbose, init_cond, constr_cond, run_number, seq_sequence, 9990, 50, 30, False)

	# track kmer freqs
	for kmer_size in range(6, 7):
		out_seq_kmers(run_file+"_track_seq.txt", run_file+"_track_kmers"+str(kmer_size)+".txt", kmer_size, True)

	# track mut, esr, and rosenberg scores
	out_seq_mut(run_file+"_track_seq.txt", run_file+"_track_mut.txt")
