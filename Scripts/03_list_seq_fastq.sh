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
# example: sbatch ../03_list_seq_fastq.sh ../lib2.fastq

lib_f="$1"

ls [0-9]*.fas > list_files.txt
while read fileind
	do
		echo $fileind
		grep '^>' $fileind > list_seq"$fileind".txt
		sed -i -e 's/>//g' list_seq"$fileind".txt
		grep -A 3 -f list_seq"$fileind".txt $lib_f > "$fileind"tq
		rm -f list_seq"$fileind".txt
done < list_files.txt

rm -f list_files.txt

