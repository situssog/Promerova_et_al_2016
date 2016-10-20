# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls(all=TRUE))

setwd("/home/evuser/Desktop/Salmon/05_offspring_haplotypes")
# setwd("C:/Users/SITG/Dropbox/Uppsala/Salmon/temporal_05_offspring_haplotypes")

DAA_haplotypes_parent_offspring<-read.table("DAA_haplotypes_parent_offspring_eddited.txt", T)
DAB_haplotypes_parent_offspring<-read.table("DAB_haplotypes_parent_offspring_eddited.txt", T)

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

table_names_toParental_DAA<-data.frame(name=c(as.character(DAA_haplotypes_parent_offspring[,1]),"No_haplotype"), change_name=c(name_change_toParental_DAA,"No_haplotype"), seq=c(as.character(DAA_haplotypes_parent_offspring_seq[,2]),"haplotype1"))
table_names_toParental_DAB<-data.frame(name=c(as.character(DAB_haplotypes_parent_offspring[,1]),"No_haplotype"), change_name=c(name_change_toParental_DAB,"No_haplotype"), seq=c(as.character(DAB_haplotypes_parent_offspring_seq[,2]),"haplotype1"))

table_names_toParental_DAA
head(table_names_toParental_DAA)
table_names_toParental_DAB

# DAA_haplotypes<-read.table("DAA_haplotypes.txt", F)
# DAB_haplotypes<-read.table("DAB_haplotypes.txt", F)

# this are the files used after filtering the haplotypes tables, deleting the samples without enough reads.

DAA_haplotypes<-read.table("DAA_haplotypes_sperm_samples.txt", F)
DAB_haplotypes<-read.table("DAB_haplotypes_sperm_samples.txt", F)


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


(final_haplotypes_offspring_DAB)

# 
# write.table(final_haplotypes_offspring_DAA, file = "summary_sperm_samples_DAA.txt", append = FALSE, quote = F, 
#             col.names = T, row.names = F)
# 
# write.table(final_haplotypes_offspring_DAB, file = "summary_sperm_samples_DAB.txt", append = FALSE, quote = F, 
#             col.names = T, row.names = F)
# 

