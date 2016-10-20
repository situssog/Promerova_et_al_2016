#! /bin/bash -l
#SBATCH -A b2014144
#SBATCH -p core
#SBATCH -J trim_seq
#SBATCH -t 48:00:00

# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

# this script will take a fastq file and remove the teg and primer sequences of all reads
# example: sbatch ../04_remove_tags_primers_seq_test.sh ../complete_table_editted_tab.csv

module load perl

table="$1"

ls 4801_191_1.fas > list_files.txt

while read fileind
	do
		ind=$(echo  ${fileind%.fas} | cut -f1 -d'_' )
		for_tag=$(awk -F '\t' '$1 ~ '${ind}' {print $10}' $table )
		rev_tag=$(awk -F '\t' '$1 ~ '${ind}' {print $18}' $table )'NNN'
		ind_code=$(echo  ${fileind%.fas} )
		perl /proj/b2014144/bin/tagcleaner-standalone-0.16/tagcleaner.pl -fasta $fileind -out_format 1 -64 -line_width 0 -nomatch 3 -tag5 $for_tag -tag3 $rev_tag
done < list_files.txt

rm -f list_files.txt


