#!/usr/bin/python
import sys
f1=open(sys.argv[1],"r")

# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

# this script use the python script to identify the possible haplotype variables and stimate the frequency of them.
# Example: python 05_haplotypes_ind_raw_num.py $fileind

dico_haplotypes={}
total_reads=float(0)

from Bio import SeqIO
for seq_record in SeqIO.parse(f1, "fasta"):
	total_reads+=1
	if str(seq_record.seq) in dico_haplotypes:
		dico_haplotypes[str(seq_record.seq)]+=1
	else:
		dico_haplotypes[str(seq_record.seq)]=float(1)

Haplotype1=['haplotype1',0]
Haplotype2=['haplotype2',0]

for haplotype in dico_haplotypes:
	if  dico_haplotypes[haplotype] > 10:
		if dico_haplotypes[haplotype] > Haplotype1[1]:
			if Haplotype1[1] > Haplotype2[1]:
				Haplotype2=Haplotype1
			Haplotype1=[haplotype ,dico_haplotypes[haplotype]]
		else:
			if dico_haplotypes[haplotype] > Haplotype2[1]:
				Haplotype2=[haplotype ,dico_haplotypes[haplotype]]

genome_ind=[]

if (Haplotype2[1]/total_reads) < ((Haplotype1[1]/total_reads)+10):
	if (Haplotype2[1]/total_reads) > ((Haplotype1[1]/total_reads)-10):
		genome_ind=['Heterozygous']
		print sys.argv[1][0:(len(sys.argv[1])-22)], genome_ind[0], Haplotype1[1]/total_reads, Haplotype1[1], Haplotype2[1]/total_reads, Haplotype2[1], Haplotype1[0], Haplotype2[0]
else:
	genome_ind=['Homozygous']
	print sys.argv[1][0:(len(sys.argv[1])-22)], genome_ind[0], Haplotype1[1]/total_reads, Haplotype1[1], Haplotype2[1]/total_reads, Haplotype2[1], Haplotype1[0], Haplotype2[0]

	


