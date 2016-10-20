#! /bin/bash -l
#SBATCH -A b2014144
#SBATCH -p core
#SBATCH -J split_seq
#SBATCH -t 48:00:00

# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

# this script generate individual fasta files for each individual in each library
# the script can be run with the following comand:
# 02_split_seq.sh library_file.fas number_lib table_with_samples_information.txt script_sort_by_TAGs.pl
# 02_split_seq.sh: corresponds to this script
# library_file.fas: fasta file with all the sequences, this is the result of transforming the FastQ file to FASTA
# number_lib: number for 1 to 6 depending to the lib file to use
# table_with_samples_information.txt: in this case conrresponds to the file "complete_table.txt"
# script_sort_by_TAGs.pl: script file "sort_by_TAGs.pl" with path of location
# example: sbatch 02_split_seq.sh ./lib1.fas 1 ./complete_table_editted_lib1.txt ./sort_by_TAGs.pl

module load perl

lib_f="$1"
lib_num="$2"
table="$3"
script_sort="$4"

while read line
	do
		IFS=$'\t' read -a array <<< "$line"
		if [ ${array[2]} -eq $lib_num ]
			echo $line
			echo ${array[2]}
			echo $lib_num
			then
				echo ${array[12]}
				echo ${array[17]}
				perl $script_sort -i $lib_f -t ${array[12]}
				perl $script_sort -i 'outfile_'${array[12]}'.fas' -t ${array[17]}
				rm 'outfile_'${array[12]}'.fas'
				rm outfile_NO_TAG.fas
				mv 'outfile_'${array[17]}'.fas' ${array[0]}'_'${array[1]}'.fas'
		fi
done < $table

