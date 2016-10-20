#! /bin/bash -l
#SBATCH -A b2014144
#SBATCH -p core
#SBATCH -J take_seq_2
#SBATCH -t 24:00:00

# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

# this script take the indviduals sequences in the fasta file, and use the names of the sequences to filter the fastq file. With this, we will have a fastq file only with the sequences for each individual.
# example: sbatch ../03_list_seq_fastq_2b.sh $fileind ../../INBOX/lib1_R1.fastq ../../INBOX/lib1_R2.fastq

file_ind="$1"
lib_f="$2"
lib_r="$3"

echo $file_ind
grep '^>' $file_ind > list_seq"$file_ind".txt
sed -i -e 's/>//g' list_seq"$file_ind".txt
grep -A 3 -f list_seq"$file_ind".txt $lib_f > ${file_ind%.fas}'_R1.fastq'
sed -i -e 's/ 1:N:/ 2:N:/g' list_seq"$file_ind".txt
grep -A 3 -f list_seq"$file_ind".txt $lib_r > ${file_ind%.fas}'_R2.fastq'
rm -f list_seq"$file_ind".txt

