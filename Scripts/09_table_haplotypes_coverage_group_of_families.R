# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls(all=TRUE))

setwd("/home/sergio/Dropbox/Uppsala/Salmon")

DAA_haplotypes_parent_offspring<-read.table("DAA_haplotypes_parent_offspring.txt", F)
DAB_haplotypes_parent_offspring<-read.table("DAB_haplotypes_parent_offspring.txt", F)

name_change_toParental_DAA<-c()
for(i in seq(1,length(DAA_haplotypes_parent_offspring[,1]))){
  if (DAA_haplotypes_parent_offspring[i,3]==100.00){
    name_change_toParental_DAA[length(name_change_toParental_DAA)+1]<-as.character(DAA_haplotypes_parent_offspring[i,2])
  } else {
    name_change_toParental_DAA[length(name_change_toParental_DAA)+1]<-paste(DAA_haplotypes_parent_offspring[i,1],"_",DAA_haplotypes_parent_offspring[i,2], sep="")
  }
} 


name_change_toParental_DAB<-c()
for(i in seq(1,length(DAB_haplotypes_parent_offspring[,1]))){
  if (DAB_haplotypes_parent_offspring[i,3]==100.00){
    name_change_toParental_DAB[length(name_change_toParental_DAB)+1]<-as.character(DAB_haplotypes_parent_offspring[i,2])
  } else {
    name_change_toParental_DAB[length(name_change_toParental_DAB)+1]<-paste(DAB_haplotypes_parent_offspring[i,1],"_",DAB_haplotypes_parent_offspring[i,2], sep="")
  }
} 


# before running the following lines you must run these lines in the terminal:
# awk 'NR%2{printf $0"\t";next;}1' DAA_list_raw_haplotypes.fasta | sed 's/>//g' - > DAA_list_raw_haplotypes_forR.txt
# awk 'NR%2{printf $0"\t";next;}1' DAB_list_raw_haplotypes.fasta | sed 's/>//g' - > DAB_list_raw_haplotypes_forR.txt


DAA_haplotypes_parent_offspring_seq<-read.table("DAA_list_raw_haplotypes_forR.txt", F)
DAB_haplotypes_parent_offspring_seq<-read.table("DAB_list_raw_haplotypes_forR.txt", F)
head(DAA_haplotypes_parent_offspring_seq)
head(DAA_haplotypes_parent_offspring)

table_names_toParental_DAA<-data.frame(name=c(as.character(DAA_haplotypes_parent_offspring[,1]),"No_haplotype"), change_name=c(name_change_toParental_DAA,"No_haplotype"), seq=c(as.character(DAA_haplotypes_parent_offspring_seq[,2]),"haplotype1"))
table_names_toParental_DAB<-data.frame(name=c(as.character(DAB_haplotypes_parent_offspring[,1]),"No_haplotype"), change_name=c(name_change_toParental_DAB,"No_haplotype"), seq=c(as.character(DAB_haplotypes_parent_offspring_seq[,2]),"haplotype1"))

table_names_toParental_DAA
head(table_names_toParental_DAA)
table_names_toParental_DAB

# DAA_haplotypes<-read.table("DAA_haplotypes.txt", F)
# DAB_haplotypes<-read.table("DAB_haplotypes.txt", F)

# this are the files used after filtering the haplotypes tables, deleting the samples without enough reads.

DAA_haplotypes<-read.table("DAA_haplotypes.txt", F)
DAB_haplotypes<-read.table("DAB_haplotypes.txt", F)


head(DAA_haplotypes)

names(DAA_haplotypes)<-c("ind_code","het","freq1","freq2","num1","num2","num_hap","total_reads","haplotype1_seq","haplotype2_seq")
names(DAB_haplotypes)<-c("ind_code","het","freq1","freq2","num1","num2","num_hap","total_reads","haplotype1_seq","haplotype2_seq")

table_names_toParental_DAA<-cbind(table_names_toParental_DAA[,1:2],haplotype1_seq=table_names_toParental_DAA[,3], haplotype2_seq=c(as.character(table_names_toParental_DAA[-(length(table_names_toParental_DAA[,3])),3]),"haplotype2"))
table_names_toParental_DAB<-cbind(table_names_toParental_DAB[,1:2],haplotype1_seq=table_names_toParental_DAB[,3], haplotype2_seq=c(as.character(table_names_toParental_DAB[-(length(table_names_toParental_DAB[,3])),3]),"haplotype2"))


head(table_names_toParental_DAA)

temporal<-merge(DAA_haplotypes,table_names_toParental_DAA, by="haplotype1_seq")
temporal<-temporal[,c(c(2:9),12,10)]
names(temporal)[9]<-"haplotype1"
names(temporal)[10]<-"haplotype2_seq"
temporal2<-merge(temporal,table_names_toParental_DAA, by="haplotype2_seq")
final_haplotypes_offspring_DAA<-temporal2[,c(c(2:10),12)]
names(final_haplotypes_offspring_DAA)[10]<-"haplotype2"


for (i in seq(1,length(final_haplotypes_offspring_DAA[,1]))){
  if (final_haplotypes_offspring_DAA[i,2]=="Homozygous"){
    final_haplotypes_offspring_DAA[i,10]<-final_haplotypes_offspring_DAA[i,9]
  }
}


temporal<-merge(DAB_haplotypes,table_names_toParental_DAB, by="haplotype1_seq")
temporal<-temporal[,c(c(2:9),12,10)]
names(temporal)[9]<-"haplotype1"
names(temporal)[10]<-"haplotype2_seq"
temporal2<-merge(temporal,table_names_toParental_DAB, by="haplotype2_seq")
final_haplotypes_offspring_DAB<-temporal2[,c(c(2:10),12)]
names(final_haplotypes_offspring_DAB)[10]<-"haplotype2"

for (i in seq(1,length(final_haplotypes_offspring_DAB[,1]))){
  if (final_haplotypes_offspring_DAB[i,2]=="Homozygous"){
    final_haplotypes_offspring_DAB[i,10]<-final_haplotypes_offspring_DAB[i,9]
  }
}


head(final_haplotypes_offspring_DAA)




























##### analyses for families

code_ind<-c()
family<-c()
ind_fam<-c()
for (i in seq(1,length(final_haplotypes_offspring_DAA[,1]))){
  code_ind[length(code_ind)+1]<-strsplit(as.character(final_haplotypes_offspring_DAA[i,1]),"_")[[1]][1]
  family[length(family)+1]<-strsplit(as.character(final_haplotypes_offspring_DAA[i,1]),"_")[[1]][2]
  ind_fam[length(ind_fam)+1]<-strsplit(as.character(final_haplotypes_offspring_DAA[i,1]),"_")[[1]][3]
}

head(final_haplotypes_offspring_DAA)

head(data.frame(code_ind,family,ind_fam,final_haplotypes_offspring_DAA[,c(2:10)]))

final_haplotypes_offspring_DAA<-data.frame(code_ind,family,ind_fam,final_haplotypes_offspring_DAA[,c(2:10)])


selected_families<-c(190:195)

DAA_0201_table<-final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$family %in% selected_families,]

head(DAA_0201_table)
DAA_0201_table


DAA_0201_table<-DAA_0201_table[with(DAA_0201_table,order(family)),]

freq_DAA_0201<-c()
freq_other<-c()
num_reads_DAA_0201<-c()
num_reads_other<-c()

for (i in c(1:length(DAA_0201_table[,1]))){
  if(DAA_0201_table[i,4]=="Homozygous"){
    freq_DAA_0201<-c(freq_DAA_0201,"NA")
    freq_other<-c(freq_other,"NA")
    num_reads_DAA_0201<-c(num_reads_DAA_0201, "NA")
    num_reads_other<-c(num_reads_other, "NA")
  } else {
    if (DAA_0201_table[i,11]=="DAA-0201"){
      freq_DAA_0201<-c(freq_DAA_0201,DAA_0201_table[i,5])
      freq_other<-c(freq_other,DAA_0201_table[i,6])
      num_reads_DAA_0201<-c(num_reads_DAA_0201, DAA_0201_table[i,7])
      num_reads_other<-c(num_reads_other, DAA_0201_table[i,8])
    } else if (DAA_0201_table[i,12]=="DAA-0201"){
      freq_DAA_0201<-c(freq_DAA_0201,DAA_0201_table[i,6])
      freq_other<-c(freq_other,DAA_0201_table[i,5])
      num_reads_DAA_0201<-c(num_reads_DAA_0201, DAA_0201_table[i,8])
      num_reads_other<-c(num_reads_other, DAA_0201_table[i,7])
    } else {
      freq_DAA_0201<-c(freq_DAA_0201,DAA_0201_table[i,5])
      freq_other<-c(freq_other,DAA_0201_table[i,6])
      num_reads_DAA_0201<-c(num_reads_DAA_0201, DAA_0201_table[i,7])
      num_reads_other<-c(num_reads_other, DAA_0201_table[i,8])
    }
  }
}


freq_DAA_0201
freq_other
num_reads_DAA_0201
num_reads_other

head(DAA_0201_table)

final_table_DAA_0201<-data.frame(DAA_0201_table[,c(1:4,11,12,10)],freq_DAA_0201, freq_other, 
                                 num_reads_DAA_0201, num_reads_other)

final_table_DAA_0201

# 
# write.table(final_table_DAA_0201, file = "summary_table_DAA_0201.txt", append = FALSE, quote = F, 
#             col.names = T, row.names = F)







