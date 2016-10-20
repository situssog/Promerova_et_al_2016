#! /bin/bash -l
#SBATCH -A b2014144
#SBATCH -p core
#SBATCH -J haplotype_ind
#SBATCH -t 48:00:00

# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

# this script use the python script to identify the possible haplotype variables and stimate the frequency of them.
# it is design to use al the fasta files in which the tags have been already removed.
# Example: sbash ../05_haplotypes_ind_raw_num.sh

ls [0-9]*_tagcleaner_*.fasta > list_files.txt

while read fileind
	do
		ind=$(echo  ${fileind%.fasta} | cut -f1 -d'_' )
		ind_code=$(echo  ${fileind%_tagcleaner_????.fasta} )
		python ../05_haplotypes_ind_raw_num.py $fileind >> haplotypes.txt		
done < list_files.txt

rm -f list_files.txt

