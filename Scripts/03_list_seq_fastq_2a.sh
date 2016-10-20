#! /bin/bash -l
#SBATCH -A b2014144
#SBATCH -p core
#SBATCH -J take_seq
#SBATCH -t 48:00:00

# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

# this script take the indviduals sequences in the fasta file, and use the names of the sequences to filter the fastq file. With this, we will have a fastq file only with the sequences for each individual.
# example: sbatch ../03_list_seq_fastq_2a.sh ../../INBOX/lib1_R1.fastq ../../INBOX/lib1_R2.fastq

lib_f="$1"
lib_r="$2"

ls [0-9]*.fas > list_files.txt
while read fileind
	do
		sbatch ../03_list_seq_fastq_2b.sh $fileind $lib_f $lib_r
done < list_files.txt

rm -f list_files.txt


