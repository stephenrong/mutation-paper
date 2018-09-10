#!/usr/bin/env python

from __future__ import division
import copy as cp
import numpy as np
import pandas as pd
import random as rd
import bisect as bisect
from math import floor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
from itertools import product

# Note 2018/06/22: 
# 1. All three constraints now based on Grantham, thus, all three constraints now disallow premature stop codons.
# 2. Human protein coding sequences now based on UniProt SwissProt alignment to hg19, manually curated non-redundant.
# 3. Now performing simulations on all human gene sequences, reformatted for running on lab linux cluster.

# # # helper functions
def get_mut_matrix(file_name):
	mu_temp = pd.read_table(file_name)
	mut_matrix = dict()
	for i in range(0, len(mu_temp)):
		mut_matrix[(mu_temp["wtMotif"][i], mu_temp["mtMotif"][i])] = mu_temp["ERV_rel_rate"][i]
	return mut_matrix

def get_esr_matrix(file_name):
	esr_temp = pd.read_table(file_name, header=None, names=["wtMotif", "ESR_chasin_score", "ESR_type"])
	esr_matrix = dict()
	for i in range(0, len(esr_temp)):
		esr_matrix[esr_temp["wtMotif"][i]] = esr_temp["ESR_chasin_score"][i]
	return esr_matrix

def get_middle_mutate(sevenmer):
	sevenmer_A = sevenmer[0:3]+"A"+sevenmer[4:7]
	sevenmer_C = sevenmer[0:3]+"C"+sevenmer[4:7]
	sevenmer_G = sevenmer[0:3]+"G"+sevenmer[4:7]
	sevenmer_T = sevenmer[0:3]+"T"+sevenmer[4:7]
	sevenmer_list = [sevenmer_A, sevenmer_C, sevenmer_G, sevenmer_T]
	return sevenmer_list

def get_cyclic_substring(seq_sequence, start, end):
	len_seq = len(seq_sequence)
	start_mod = start % len_seq
	end_mod = end % len_seq
	if start_mod < end_mod:
		return seq_sequence[start_mod:end_mod]
	else:
		return seq_sequence[start_mod:len_seq]+seq_sequence[0:end_mod]

def initiate_seq_sequence(seq_length=1000):
	seq_sequence = "".join([np.random.choice(["A", "C", "G", "T"]) for i in range(0, seq_length)])
	seq_sequence = Seq(seq_sequence, IUPAC.unambiguous_dna)
	return seq_sequence

def initiate_seq_mutrates(mut_matrix, seq_sequence):
	len_seq = len(seq_sequence)
	seq_mutrates = np.zeros([4, len_seq])
	for i in range(0, len_seq):
		wtMotif = get_cyclic_substring(seq_sequence, i-3, i+4)
		mtMotif_list = get_middle_mutate(wtMotif)
		for j in range(0, len(mtMotif_list)):
			mtMotif = mtMotif_list[j]
			if wtMotif != mtMotif:
				seq_mutrates[j, i] = mut_matrix[(wtMotif, mtMotif)]
	return seq_mutrates

def initate_seq_sumrates(seq_mutrates):
	seq_sumrates = np.sum(seq_mutrates, axis=0)
	return seq_sumrates

def weighted_choice(weights):
	# https://stackoverflow.com/questions/24140114/
	# fast-way-to-obtain-a-random-index-from-an-array-of-weights-in-python
	# way faster than using np.random.choice, but slowed by np.cumsum
    cs_weights = np.cumsum(weights)
    return bisect.bisect(cs_weights, np.random.random()*cs_weights[-1])

def get_mutation_event(seq_sequence, seq_mutrates, seq_sumrates):
	j = weighted_choice(seq_sumrates)
	i = weighted_choice(seq_mutrates[0:4,j])
	return (i, j)

def update_constraint_none(mut_matrix, seq_sequence, seq_mutrates, seq_sumrates, temp_index):
	len_seq = len(seq_sequence)
	update_middle = ["A", "C", "G", "T"]
	# update sequence by-copy
	seq_sequence_update = seq_sequence[0:temp_index[1]]+\
		update_middle[temp_index[0]]+seq_sequence[temp_index[1]+1:len_seq]
	# update mutrates in-place
	seq_mutrates_update = seq_mutrates
	start_mod = (temp_index[1]-3) % len_seq
	end_mod = (temp_index[1]+4) % len_seq
	if start_mod < end_mod:
		replace_indices = range(start_mod, end_mod)
	else:
		replace_indices = range(start_mod, len_seq)+range(0, end_mod)
	for i in replace_indices:
		wtMotif = get_cyclic_substring(seq_sequence_update, i-3, i+4)
		mtMotif_list = get_middle_mutate(wtMotif)
		for j in range(0, len(mtMotif_list)):
			mtMotif = mtMotif_list[j]
			if wtMotif != mtMotif:
				seq_mutrates_update[j, i] = mut_matrix[(wtMotif, mtMotif)]
			else:
				seq_mutrates_update[j, i] = 0
	# update sumrates in-place
	seq_sumrates_update = seq_sumrates
	if start_mod < end_mod:
		seq_sumrates_update[start_mod:end_mod] = np.sum(seq_mutrates_update[0:4, start_mod:end_mod], axis=0)
	else:
		seq_sumrates_update[start_mod:len_seq] = np.sum(seq_mutrates_update[0:4, start_mod:len_seq], axis=0)
		seq_sumrates_update[0:end_mod] = np.sum(seq_mutrates_update[0:4, 0:end_mod], axis=0)
	return seq_sequence_update, seq_mutrates_update, seq_sumrates_update

def get_kmer_profile(seq_sequence, kmer_size):
	kmer_profile = [seq_sequence[i:i+kmer_size] for i in range(0, len(seq_sequence)-kmer_size+1)]
	return kmer_profile

def get_grantham_matrix():
	# (Grantham 1974) DOI:10.1126/science.185.4154.862
	# https://github.com/ashutoshkpandey/Annotation/blob/master/Grantham_score_calculator.py
	grantham_matrix = {
		'S':{'R':110, 'L':145, 'P':74, 'T':58, 'A':99, 'V':124, 'G':56, 'I':142, 'F':155, 'Y':144, 'C':112, 'H':89, 'Q':68, 'N':46, 'K':121, 'D':65, 'E':80, 'M':135, 'W':177},
		'R':{'R':0, 'L':102, 'P':103, 'T':71, 'A':112, 'V':96, 'G':125, 'I':97, 'F':97, 'Y':77, 'C':180, 'H':29, 'Q':43, 'N':86, 'K':26, 'D':96, 'E':54, 'M':91, 'W':101, 'S':110},
		'L':{'R':102, 'L':0, 'P':98, 'T':92, 'A':96, 'V':32, 'G':138, 'I':5, 'F':22, 'Y':36, 'C':198, 'H':99, 'Q':113, 'N':153, 'K':107, 'D':172, 'E':138, 'M':15, 'W':61, 'S':145},
		'P':{'R':103, 'L':98, 'P':0, 'T':38, 'A':27, 'V':68, 'G':42, 'I':95, 'F':114, 'Y':110, 'C':169, 'H':77, 'Q':76, 'N':91, 'K':103, 'D':108, 'E':93, 'M':87, 'W':147, 'S':74},
		'T':{'R':71, 'L':92, 'P':38, 'T':0, 'A':58, 'V':69, 'G':59, 'I':89, 'F':103, 'Y':92, 'C':149, 'H':47, 'Q':42, 'N':65, 'K':78, 'D':85, 'E':65, 'M':81, 'W':128, 'S':58},
		'A':{'R':112, 'L':96, 'P':27, 'T':58, 'A':0, 'V':64, 'G':60, 'I':94, 'F':113, 'Y':112, 'C':195, 'H':86, 'Q':91, 'N':111, 'K':106, 'D':126, 'E':107, 'M':84, 'W':148, 'S':99},
		'V':{'R':96, 'L':32, 'P':68, 'T':69, 'A':64, 'V':0, 'G':109, 'I':29, 'F':50, 'Y':55, 'C':192, 'H':84, 'Q':96, 'N':133, 'K':97, 'D':152, 'E':121, 'M':21, 'W':88, 'S':124},
		'G':{'R':125, 'L':138, 'P':42, 'T':59, 'A':60, 'V':109, 'G':0, 'I':135, 'F':153, 'Y':147, 'C':159, 'H':98, 'Q':87, 'N':80, 'K':127, 'D':94, 'E':98, 'M':127, 'W':184, 'S':56},
		'I':{'R':97, 'L':5, 'P':95, 'T':89, 'A':94, 'V':29, 'G':135, 'I':0, 'F':21, 'Y':33, 'C':198, 'H':94, 'Q':109, 'N':149, 'K':102, 'D':168, 'E':134, 'M':10, 'W':61, 'S':142},
		'F':{'R':97, 'L':22, 'P':114, 'T':103, 'A':113, 'V':50, 'G':153, 'I':21, 'F':0, 'Y':22, 'C':205, 'H':100, 'Q':116, 'N':158, 'K':102, 'D':177, 'E':140, 'M':28, 'W':40, 'S':155},
		'Y':{'R':77, 'L':36, 'P':110, 'T':92, 'A':112, 'V':55, 'G':147, 'I':33, 'F':22, 'Y':0, 'C':194, 'H':83, 'Q':99, 'N':143, 'K':85, 'D':160, 'E':122, 'M':36, 'W':37, 'S':144},
		'C':{'R':180, 'L':198, 'P':169, 'T':149, 'A':195, 'V':192, 'G':159, 'I':198, 'F':205, 'Y':194, 'C':0, 'H':174, 'Q':154, 'N':139, 'K':202, 'D':154, 'E':170, 'M':196, 'W':215, 'S':112},
		'H':{'R':29, 'L':99, 'P':77, 'T':47, 'A':86, 'V':84, 'G':98, 'I':94, 'F':100, 'Y':83, 'C':174, 'H':0, 'Q':24, 'N':68, 'K':32, 'D':81, 'E':40, 'M':87, 'W':115, 'S':89},
		'Q':{'R':43, 'L':113, 'P':76, 'T':42, 'A':91, 'V':96, 'G':87, 'I':109, 'F':116, 'Y':99, 'C':154, 'H':24, 'Q':0, 'N':46, 'K':53, 'D':61, 'E':29, 'M':101, 'W':130, 'S':68},
		'N':{'R':86, 'L':153, 'P':91, 'T':65, 'A':111, 'V':133, 'G':80, 'I':149, 'F':158, 'Y':143, 'C':139, 'H':68, 'Q':46, 'N':0, 'K':94, 'D':23, 'E':42, 'M':142, 'W':174, 'S':46},
		'K':{'R':26, 'L':107, 'P':103, 'T':78, 'A':106, 'V':97, 'G':127, 'I':102, 'F':102, 'Y':85, 'C':202, 'H':32, 'Q':53, 'N':94, 'K':0, 'D':101, 'E':56, 'M':95, 'W':110, 'S':121},
		'D':{'R':96, 'L':172, 'P':108, 'T':85, 'A':126, 'V':152, 'G':94, 'I':168, 'F':177, 'Y':160, 'C':154, 'H':81, 'Q':61, 'N':23, 'K':101, 'D':0, 'E':45, 'M':160, 'W':181, 'S':65},
		'E':{'R':54, 'L':138, 'P':93, 'T':65, 'A':107, 'V':121, 'G':98, 'I':134, 'F':140, 'Y':122, 'C':170, 'H':40, 'Q':29, 'N':42, 'K':56, 'D':45, 'E':0, 'M':126, 'W':152, 'S':80},
		'M':{'R':91, 'L':15, 'P':87, 'T':81, 'A':84, 'V':21, 'G':127, 'I':10, 'F':28, 'Y':36, 'C':196, 'H':87, 'Q':101, 'N':142, 'K':95, 'D':160, 'E':126, 'M':0, 'W':67, 'S':135}, 
		'W':{'R':101, 'L':61,  'P':147,  'T':128,  'A':148,  'V':88,  'G':184,  'I':61,  'F':40,  'Y':37,  'C':215,  'H':115,  'Q':130,  'N':174,  'K':110,  'D':181,  'E':152,  'M':67,  'W':0 , 'S':177}
	}
	return grantham_matrix

def check_constraint_identity(seq_sequence, temp_index):
	# get codon change
	update_change = ["A", "C", "G", "T"]
	codon_start = int(floor(temp_index[1]/3)*3)
	codon_end = codon_start+3
	codon_change = temp_index[1] - codon_start
	codon_before = seq_sequence[codon_start:codon_end]
	codon_after = codon_before[0:codon_change]+update_change[temp_index[0]]+codon_before[codon_change+1:3]
	# get amino acid change
	amino_before = codon_before.translate()
	amino_after = codon_after.translate()
	# check amino acid change
	amino_check = amino_before == amino_after
	return codon_before, codon_after, amino_check

def check_constraint_grantham(seq_sequence_init, seq_sequence, temp_index, threshold):
	# get codon change
	update_change = ["A", "C", "G", "T"]
	codon_start = int(floor(temp_index[1]/3)*3)
	codon_end = codon_start+3
	codon_change = temp_index[1] - codon_start
	codon_before = seq_sequence_init[codon_start:codon_end]
	codon_temp = seq_sequence[codon_start:codon_end]
	codon_after = codon_temp[0:codon_change]+update_change[temp_index[0]]+codon_temp[codon_change+1:3]
	# get amino acid change
	amino_before = codon_before.translate()
	amino_after = codon_after.translate()
	# check amino acid change
	# print (amino_before, amino_after)
	if amino_before == amino_after:
		amino_check = True
	elif amino_before == "*" or amino_after == "*":
		amino_check = False
	else:
		amino_grantham = grantham_matrix[amino_before][amino_after]
		amino_check = amino_grantham <= threshold
	return codon_before, codon_after, amino_check

def get_esrmean(seq_sequence):
	seq_kmers = [str(seq_sequence[i:i+6]) for i in range(0, len(seq_sequence)-5)]
	esr_mean = np.mean([esr_matrix[i] for i in seq_kmers])
	return esr_mean

def get_kmerfreq(seq_sequence, kmer_size):
	seq_kmers = [str(seq_sequence[i:i+kmer_size]) for i in range(0, len(seq_sequence)-kmer_size+1)]
	count_kmers = Counter(seq_kmers)
	return count_kmers

def out_seq_kmers(in_file, out_file, kmer_size, verbose):
	seq_temp = pd.read_table(in_file)
	max_value = np.max(seq_temp["mut_scaled"])
	with open(out_file, "w+") as f:
		count_keys = sorted(["".join(k) for k in list(product(["A", "C", "T", "G"], repeat=kmer_size))])
		prefix = "init_cond"+"\t"+"constr_cond"+"\t"+"run_number"
		f.write(prefix+"\t"+"seq_length"+"\t"+"mutation"+"\t"+"mut_scaled"+"\t"+"mutation_tot"+"\t"+"mut_scaled_tot"+"\t"+"\t".join(count_keys)+"\n")
		for i in range(0, len(seq_temp)):
			prefix = str(seq_temp["init_cond"][i])+"\t"+str(seq_temp["constr_cond"][i])+"\t"+str(seq_temp["run_number"][i])
			if verbose:
				count_kmers = get_kmerfreq(seq_temp["sequence"][i], kmer_size)
				temp_kmers = [count_kmers.get(j, 0) for j in count_keys]
				f.write(prefix+"\t"+str(seq_temp["seq_length"][i])+"\t"+str(seq_temp["mutation"][i])+"\t"+str(seq_temp["mut_scaled"][i])+"\t"+str(seq_temp["mutation_tot"][i])+"\t"+str(seq_temp["mut_scaled_tot"][i])+"\t"+"\t".join([str(j) for j in temp_kmers])+"\n")
			else:
				if seq_temp["mut_scaled"][i] == 0 or seq_temp["mut_scaled"][i] == max_value:
					count_kmers = get_kmerfreq(seq_temp["sequence"][i], kmer_size)
					temp_kmers = [count_kmers.get(j, 0) for j in count_keys]
					f.write(prefix+"\t"+str(seq_temp["seq_length"][i])+"\t"+str(seq_temp["mutation"][i])+"\t"+str(seq_temp["mut_scaled"][i])+"\t"+str(seq_temp["mutation_tot"][i])+"\t"+str(seq_temp["mut_scaled_tot"][i])+"\t"+"\t".join([str(j) for j in temp_kmers])+"\n")

# # # global variables
# mut_file = "../data/mu-matrix-7mers.txt"
# mut_matrix = get_mut_matrix(mut_file)
esr_file = "../data/chasin_ESS_ESE_numbers.txt"
esr_matrix = get_esr_matrix(esr_file)
grantham_matrix = get_grantham_matrix()
update_middle = ["A", "C", "G", "T"]

# # # simulation functions
def muSimsConstraintGranthamSimulation(mut_matrix, track_seq, track_mut, track_verbose, init_cond, constr_cond, run_number, seq_sequence, mut_final, mut_step, threshold, verbose):
	# initialize sequence
	# seq_sequence = initiate_seq_sequence(seq_length)
	seq_mutrates = initiate_seq_mutrates(mut_matrix, seq_sequence)
	seq_sumrates = initate_seq_sumrates(seq_mutrates)
	seq_sequence_init = cp.deepcopy(seq_sequence)
	# save information
	prefix = init_cond+"\t"+constr_cond+"\t"+run_number
	iter_print = [np.int(j) for j in np.floor(np.arange(0, mut_step+1)*mut_final/mut_step)]
	i = 0
	j = 0
	track_seq.write(prefix+"\t"+str(len(seq_sequence))+"\t"+str(i)+"\t"+str(i/len(seq_sequence))+"\t"+str(j)+"\t"+str(j/len(seq_sequence))+"\t"+str(seq_sequence)+"\n")
	track_mut.write(prefix+"\t"+str(len(seq_sequence))+"\t"+str(i)+"\t"+str(i/len(seq_sequence))+"\t"+str(j)+"\t"+str(j/len(seq_sequence))+"\t"+str(np.mean(seq_sumrates)/3)+"\t"+str(get_esrmean(seq_sequence)) +"\n")
	for i in range(1, mut_final+1):
		# get mutation
		amino_check = False
		while not amino_check:
			j += 1
			temp_index = get_mutation_event(seq_sequence, seq_mutrates, seq_sumrates)
			codon_before, codon_after, amino_check = check_constraint_grantham(seq_sequence_init, seq_sequence, temp_index, threshold)
			if verbose and str(codon_after.translate())=="*":
				wtMotif = get_cyclic_substring(seq_sequence, temp_index[1]-3, temp_index[1]+4)
				mtMotif = wtMotif[0:3]+update_middle[temp_index[0]]+wtMotif[4:7]
				track_verbose.write(prefix+"\t"+str(len(seq_sequence))+"\t"+str(i)+"\t"+str(i/len(seq_sequence))+"\t"+str(j)+"\t"+str(j/len(seq_sequence))+"\t"+\
					str(codon_before)+"\t"+str(codon_after)+"\t"+str(wtMotif)+"\t"+str(mtMotif)+"\t"+str(mut_matrix[(wtMotif, mtMotif)])+"\t"+\
					str(codon_before.translate())+"\t"+str(codon_after.translate())+"\t"+str(amino_check)+"\n")
		# update sequence
		seq_sequence, seq_mutrates, seq_sumrates = update_constraint_none(mut_matrix, seq_sequence, seq_mutrates, seq_sumrates, temp_index)
		# save information
		# if (i % np.ceil(mut_final/mut_step) == 0):
		if i in iter_print:
			track_seq.write(prefix+"\t"+str(len(seq_sequence))+"\t"+str(i)+"\t"+str(i/len(seq_sequence))+"\t"+str(j)+"\t"+str(j/len(seq_sequence))+"\t"+str(seq_sequence)+"\n")
			track_mut.write(prefix+"\t"+str(len(seq_sequence))+"\t"+str(i)+"\t"+str(i/len(seq_sequence))+"\t"+str(j)+"\t"+str(j/len(seq_sequence))+"\t"+str(np.mean(seq_sumrates)/3)+"\t"+str(get_esrmean(seq_sequence)) +"\n")
	return seq_sequence
