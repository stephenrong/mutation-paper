#!/usr/bin/env python

from __future__ import division
import sys
import gzip
from mu_sims_module import *
from numpy.random import randint

# # # real sequences
if __name__ == '__main__':
	iter_start = int(sys.argv[1])
	iter_end = int(sys.argv[2])
	np.random.seed(iter_start+2)

	# import mut_matrix
	mut_file = "../data/mu-matrix-7mers.txt"
	mut_matrix = get_mut_matrix(mut_file)

	# import sequences
	with gzip.open("../data/hg19-unipAliSwissprot-cds_genes.txt.gz", "rt") as handle:
		seq_exon = [record.seq for record in SeqIO.parse(handle, "fasta")]

	# output files
	run_file = "../results/simulations/exon_"+str(iter_start)+"_"+str(iter_end)
	with open(run_file+"_track_seq.txt", "w+") as track_seq, open(run_file+"_track_mut.txt", "w+") as track_mut, open(run_file+"_track_verbose.txt", "w+") as track_verbose:

		# add header
		prefix = "init_cond"+"\t"+"constr_cond"+"\t"+"run_number"
		track_seq.write(prefix+"\t"+"seq_length"+"\t"+"mutation"+"\t"+"mut_scaled"+"\t"+"mutation_tot"+"\t"+"mut_scaled_tot"+"\t"+"sequence"+"\n")
		track_mut.write(prefix+"\t"+"seq_length"+"\t"+"mutation"+"\t"+"mut_scaled"+"\t"+"mutation_tot"+"\t"+"mut_scaled_tot"+"\t"+"mut_mean"+"\t"+"esr_mean"+"\n")
		track_verbose.write(prefix+"\t"+"seq_length"+"mutation"+"\t"+"mut_scaled"+"\t"+"mutation_tot"+"\t"+"mut_scaled_tot"+"\t"+"codonBefore"+"\t" "codonAfter"+"\t"+"wtMotif"+"\t"+"mtMotif"+"\t"+"mutRate"+"\t"+"aminoBefore"+"\t"+"aminoAfter"+"\t"+"aminoCheck"+"\n")

		for i in range(iter_start, iter_end):
			init_cond = "exon"
			constr_cond = "neutral"
			run_number = str(i)
			print init_cond+"_"+constr_cond+"_"+run_number
			seq_sequence = seq_exon[i]
			mut_step = 50
			mut_final = 10*len(seq_sequence)
			muSimsConstraintGranthamSimulation(mut_matrix, track_seq, track_mut, track_verbose, init_cond, constr_cond, run_number, seq_sequence, mut_final, mut_step, 300, False)

		for i in range(iter_start, iter_end):
			init_cond = "exon"
			constr_cond = "identity"
			run_number = str(i)
			print init_cond+"_"+constr_cond+"_"+run_number
			seq_sequence = seq_exon[i]
			mut_step = 50
			mut_final = 10*len(seq_sequence)
			muSimsConstraintGranthamSimulation(mut_matrix, track_seq, track_mut, track_verbose, init_cond, constr_cond, run_number, seq_sequence, mut_final, mut_step, 0, False)

		# for i in range(iter_start, iter_end):
		# 	init_cond = "exon"
		# 	constr_cond = "grantham"
		# 	run_number = str(i)
		# 	print init_cond+"_"+constr_cond+"_"+run_number
		# 	seq_sequence = seq_exon[i]
		# 	mut_step = 50
		# 	mut_final = 10*len(seq_sequence)
		# 	muSimsConstraintGranthamSimulation(mut_matrix, track_seq, track_mut, track_verbose, init_cond, constr_cond, run_number, seq_sequence, mut_final, mut_step, 30, False)

	# track kmer freqs
	for kmer_size in range(6, 7):
		out_seq_kmers(run_file+"_track_seq.txt", run_file+"_track_kmers"+str(kmer_size)+".txt", kmer_size, True)

	# track mut, esr, and rosenberg scores
	out_seq_mut(run_file+"_track_seq.txt", run_file+"_track_mut.txt")
