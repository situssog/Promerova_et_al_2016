# ==================================================
# Documents and scripts were written by: Sergio Tusso
# 2016
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# +++++++++++++++++++++++++++++++++++++++++++++++++

# +++++++++++++++++++++++++++++++++++++++++++++++++

# ================== Protocol for pape: =================
# "No evidence for MHC class II-based non-random mating at the gametic haplotype in Atlantic salmon"

# This is a step by step protocol for genotyping
# All scripts are in the folder Scripts

########

# Merging paired reads
# Software FLASH 
# http://ccb.jhu.edu/software/FLASH/

flash /raw_data/140619_M00485_0137_000000000-A9G08/Sample_Library_1/Library_1_TGACCA_L001_R1_001.fastq /raw_data/140619_M00485_0137_000000000-A9G08/Sample_Library_1/Library_1_TGACCA_L001_R2_001.fastq
flash /raw_data/140619_M00485_0137_000000000-A9G08/Sample_Library_2/Library_2_ACAGTG_L001_R1_001.fastq /raw_data/140619_M00485_0137_000000000-A9G08/Sample_Library_2/Library_2_ACAGTG_L001_R2_001.fastq
flash /raw_data/140619_M00485_0137_000000000-A9G08/Sample_Library_3/Library_3_GCCAAT_L001_R1_001.fastq /raw_data/140619_M00485_0137_000000000-A9G08/Sample_Library_3/Library_3_GCCAAT_L001_R2_001.fastq
flash /raw_data/140630_M00629_0010_000000000-AA24Y/Sample_Library_4/Library_4_CAGATC_L001_R1_001.fastq /raw_data/140630_M00629_0010_000000000-AA24Y/Sample_Library_4/Library_4_CAGATC_L001_R2_001.fastq
flash /raw_data/140630_M00629_0010_000000000-AA24Y/Sample_Library_5/Library_5_ACTTGA_L001_R1_001.fastq /raw_data/140630_M00629_0010_000000000-AA24Y/Sample_Library_5/Library_5_ACTTGA_L001_R2_001.fastq
flash /raw_data/140630_M00629_0010_000000000-AA24Y/Sample_Library_6/Library_6_GATCAG_L001_R1_001.fastq /raw_data/140630_M00629_0010_000000000-AA24Y/Sample_Library_6/Library_6_GATCAG_L001_R2_001.fastq

# the output of FLASH is a FASTQ file with whole sequences
# A Perl script was used to transform the FASTQ file to FASTA file
# the script was obtained from the link:
# http://www.uni-mainz.de/FB/Biologie/Anthropologie/487_ENG_HTML.php
# the output is the FASTA file called out.fa that I transformed to lib1.fa

perl FASTQ_to_FASTA.pl -i out.extendedFrags.fastq
mv out.fa lib1.fa

perl ./FASTQ_to_FASTA.pl -i lib2.fastq
perl ./FASTQ_to_FASTA.pl -i lib3.fastq
perl ./FASTQ_to_FASTA.pl -i lib4.fastq
perl ./FASTQ_to_FASTA.pl -i lib5.fastq
perl ./FASTQ_to_FASTA.pl -i lib6.fastq


# From the same package of scripts, sort_by_TAGs.pl was used to filter the fasta file by tags (barcodes).
# to find the sequences of each individual, both the primer and the tag must be used. For that the forward tag shown in the file "complete_table_editted_tab.csv" (in folder Samples_info) and the reverse primer needs to be transformed to reverse complement before using the script
# example: perl /home/evuser/Downloads/Software/sort_by_TAGs.pl -i lib1.fas -t gatcgcgaGTGTTTTATTGGGTTTCTTTTCTC
# the original sequence in the table was:
# tagctagtCTCTAAATTACTTCTCTCTCTTAC and it was transformed to GTAAGAGAGAGAAGTAATTTAGAGactagcta before using it

perl sort_by_TAGs.pl -i outfile_gatcgcgaGTGTTTTATTGGGTTTCTTTTCTC.fas -t GTAAGAGAGAGAAGTAATTTAGAGactagcta

# The list of sequences for the reverse primer was obtained using the following commands:
# these commands produce a list in text file that were included in the original table with all de individuals.


awk '//{print $17}' complete_table.txt > seq_ori_t.txt
grep -v w_seq_r seq_ori_t.txt > seq_ori.txt
awk 'BEGIN {
  j = n = split("A C G T a c g t", t)
  for (i = 0; ++i <= n;)
    map[t[i]] = t[j--]
  }
{
  if (/^>/) print
  else {
    for (i = length; i; i--)
    printf "%s", map[substr($0, i, 1)]
    print x
  }   
  }' seq_ori.txt > rev_seq_ori.txt

  
# The individuals were sorted by library
awk -F $' ' '$3 ~ 1' complete_table_editted_tab.csv > complete_table_editted_lib1.txt
awk -F $' ' '$3 ~ 2' complete_table_editted_tab.csv > complete_table_editted_lib2.txt
awk -F $' ' '$3 ~ 3' complete_table_editted_tab.csv > complete_table_editted_lib3.txt
awk -F $' ' '$3 ~ 4' complete_table_editted_tab.csv > complete_table_editted_lib4.txt
awk -F $' ' '$3 ~ 5' complete_table_editted_tab.csv > complete_table_editted_lib5.txt
awk -F $' ' '$3 ~ 6' complete_table_editted_tab.csv > complete_table_editted_lib6.txt

# those tables were used to run the script 02_split_seq.sh and filter all read for each individual in a single file

02_split_seq.sh lib1.fas 1 complete_table_editted_lib1.txt sort_by_TAGs.pl
02_split_seq_2.sh lib2.fas 2 complete_table_editted_lib2.txt sort_by_TAGs.pl
02_split_seq_2.sh lib3.fas 3 complete_table_editted_lib3.txt sort_by_TAGs.pl
02_split_seq_2.sh lib4.fas 4 complete_table_editted_lib4.txt sort_by_TAGs.pl
02_split_seq_2.sh lib5.fas 5 complete_table_editted_lib5.txt sort_by_TAGs.pl
02_split_seq_2.sh lib6.fas 6 complete_table_editted_lib6.txt sort_by_TAGs.pl

# tagcleaner.pl was used to remove the tag and primer sequences of each read in each individual file
# 04_remove_tags_primers_seq.sh was used for that
# This will remove the tagges and primer sequences and create a new fasta
# The script is automatised, so it will do the job for all the individuals inside of one library.

04_remove_tags_primers_seq.sh complete_table_editted_tab.csv

# Script 05_haplotypes_ind_raw_num.sh was used to identify the diferent possible haplotypes for each individual.
# This script uses the python script 05_haplotypes_ind_raw_num.py
# The outcome will be a single file for all the individuals of the library with the two more dominant haplotypes and the sequences.
# The script estimates the freqeuncy for each the the possible haplotypes and determines if the individual is homozygote or heterozygote base on the frequencies.
# The output is in the file haplotypes.txt


# merge the information of the libraries from 1 to 3
cat lib[0-3]/haplotype* > ADD_haplotypes.txt

# the same for the libraries 4 to 6
cat lib[4-6]/haplotype* > DAB_haplotypes.txt

# A list with all candidate haplotypes was created to filter it and check how many unique haplotypes there are.
# the output file shows unique haplotype

awk '//{print $9"\n"$10}' ADD_haplotypes.txt | sort -u | grep -v haplotype1 | grep -v haplotype2 > ADD_list_raw_haplotypes.txt
awk '//{print $9"\n"$10}' DAB_haplotypes.txt | sort -u | grep -v haplotype1 | grep -v haplotype2 > DAB_list_raw_haplotypes.txt

nl DAA_list_raw_haplotypes.txt | awk '//{print ">haplotype_"$1"\n"$2}' - > DAA_list_raw_haplotypes.fasta
nl DAB_list_raw_haplotypes.txt | awk '//{print ">haplotype_"$1"\n"$2}' - > DAB_list_raw_haplotypes.fasta

# offspring genotypes were compared with parental genotypes by blast
# The output is the file with the name of the parental haplotypes with the corresponding offspring haplotype name


makeblastdb -in DAA_alleles_parents.fas -dbtype 'nucl'
makeblastdb -in DAB_alleles_parents.fas -dbtype 'nucl'

blastn -query DAA_list_raw_haplotypes.fasta -db DAA_alleles_parents.fas -outfmt 6 -max_target_seqs 1 > DAA_haplotypes_parent_offspring.txt
blastn -query DAB_list_raw_haplotypes.fasta -db DAB_alleles_parents.fas -outfmt 6 -max_target_seqs 1 > DAB_haplotypes_parent_offspring.txt

# The R script "07_table_haplotypes.R" was used to convine all the tables. Transform the sequences to names (the names of the parental haplotypes) and to integrate the information of each family to produce a table with the frequencies of each haplotype in each family.
# The script uses several tables containes in the folder 05_offspring_haplotypes. Many of them were produced previously and other were created manually.
