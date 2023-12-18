#!/usr/bin/env python3

import glob
from Bio import SeqIO
from itertools import product
import pandas as pd


#Read the abi file-----record = SeqIO.read("*.ab1", "abi")

#Show the annotations attribute of the record-----print(record.annotations.keys())
#Under abif_raw is another dictionary of data-----print(record.annotations["abif_raw"].keys())

#DATA9-Channel one analyzed data (G)-----record.annotations['abif_raw']['DATA9']
#DATA10-Channel two analyzed data (A)-----record.annotations['abif_raw']['DATA10']
#DATA11-Channel three analyzed data (T)-----record.annotations['abif_raw']['DATA11']
#DATA12-Channel four analyzed data (C)-----record.annotations['abif_raw']['DATA12']

#Base location-----record.annotations['abif_raw']['PLOC2']
#Basecalled sequence-----record.annotations['abif_raw']['PBAS2']


#Load the CRISPR library data for reference.
def get_lib():

	lib_path = glob.glob('*.txt')
	f = open(lib_path[0])
	list_lib = f.readlines()
	f.close()
	dict_lib = {}
	for line in list_lib:
		dict_lib.update({line.split('\t',1)[0]:line.split('\t',1)[1]})#20-bp target recognition sequence : target genes.

	return dict_lib

#Get the sequences of sgRNAs from abi file.
def get_seq():

	file = open('Results','w')
	abi_path = glob.glob('*.ab1')
	reference_lib = get_lib()

	for path in abi_path:

		sample_name = path.split('.')[0].split('_')[0]
		record = SeqIO.read(path, "abi")
		output_seq = str(record.annotations['abif_raw']['PBAS2'])
		flag = output_seq.index('GTTTTAGAG')#Locate the position of the end base. Change the sequence when using distinct scaffold.
		base_location = record.annotations['abif_raw']['PLOC2'][flag-22:flag-2]#Locate the position of 20-bp CRISPR-Cas9 sgRNA.
		DATA9 = record.annotations['abif_raw']['DATA9']#Channel one
		DATA10 = record.annotations['abif_raw']['DATA10']#Channel two
		DATA11 = record.annotations['abif_raw']['DATA11']#Channel three
		DATA12 = record.annotations['abif_raw']['DATA12']#Channel four
		Value_G = [DATA9[index] for index in base_location]
		Value_A = [DATA10[index] for index in base_location]
		Value_T = [DATA11[index] for index in base_location]
		Value_C = [DATA12[index] for index in base_location]
		df = pd.DataFrame({"G":Value_G,"A":Value_A,"T":Value_T,"C":Value_C})
		seq_len = len(base_location)
		Seq_pocket = []

		for loc_index in range(seq_len):
			base_call = ""
			value_list = df.loc[loc_index]
			value_sum = sum(value_list)
			G_ratio = value_list[0]/value_sum
			A_ratio = value_list[1]/value_sum
			T_ratio = value_list[2]/value_sum
			C_ratio = value_list[3]/value_sum

			if G_ratio > 0.10:#Assume the cut-off ratio.
				base_call = base_call + "G"
			if A_ratio > 0.10:
				base_call = base_call + "A"
			if T_ratio > 0.10:
				base_call = base_call + "T"
			if C_ratio > 0.10:
				base_call = base_call + "C"					
			Seq_pocket.append(base_call)
						
		#Deconvolve sgRNAs
		degenerate_seq1 = Seq_pocket[:10]
		degenerate_seq2 = Seq_pocket[10:]
		car_prod1 = product(*degenerate_seq1)#Cartesian product
		car_prod2 = product(*degenerate_seq2)#Cartesian product
		total_target = reference_lib.keys()
		total_target_half1 = [i[:10] for i in total_target]
		total_target_half2 = [i[10:] for i in total_target]
		half1,half2 = [],[]
		for candidate_seq1 in car_prod1:
			candidate_seq1 = "".join(candidate_seq1)
			if candidate_seq1 in total_target_half1:
				half1.append(candidate_seq1)
		for candidate_seq2 in car_prod2:
			candidate_seq2 = "".join(candidate_seq2)
			if candidate_seq2 in total_target_half2:
				half2.append(candidate_seq2)

		for seq_half1 in half1:
			for seq_half2 in half2:
				candidate_seq = seq_half1 + seq_half2
				if candidate_seq in total_target:
					info = reference_lib[candidate_seq]
					file.write(sample_name+'\t'+candidate_seq+'\t'+info)

	file.close()

	return None

if __name__ == '__main__':
	get_seq()
	