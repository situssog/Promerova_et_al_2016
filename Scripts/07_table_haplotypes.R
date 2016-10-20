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
table_names_toParental_DAB

# DAA_haplotypes<-read.table("DAA_haplotypes.txt", F)
# DAB_haplotypes<-read.table("DAB_haplotypes.txt", F)

# this are the files used after filtering the haplotypes tables, deleting the samples without enough reads.

DAA_haplotypes<-read.table("DAA_haplotypes_editted_01_10.txt", F)
DAB_haplotypes<-read.table("DAB_haplotypes_editted_01_10.txt", F)


head(DAA_haplotypes)

names(DAA_haplotypes)<-c("ind_code","het","freq1","freq2","num1","num2","num_hap","total_reads","haplotype1_seq","haplotype2_seq")
names(DAB_haplotypes)<-c("ind_code","het","freq1","freq2","num1","num2","num_hap","total_reads","haplotype1_seq","haplotype2_seq")

table_names_toParental_DAA<-cbind(table_names_toParental_DAA[,1:2],haplotype1_seq=table_names_toParental_DAA[,3], haplotype2_seq=c(as.character(table_names_toParental_DAA[-(length(table_names_toParental_DAA[,3])),3]),"haplotype2"))
table_names_toParental_DAB<-cbind(table_names_toParental_DAB[,1:2],haplotype1_seq=table_names_toParental_DAB[,3], haplotype2_seq=c(as.character(table_names_toParental_DAB[-(length(table_names_toParental_DAB[,3])),3]),"haplotype2"))




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

final_haplotypes_offspring_DAA<-data.frame(code_ind,family,ind_fam,final_haplotypes_offspring_DAA[,c(2,9,10)])

table_families<-read.table("sample_info.csv", T)
head(table_families)
names(table_families)[3]<-"family"

table_parental_males_DAA<-read.table("males_summary_DAA.csv", T)
head(table_parental_males_DAA)
names(table_parental_males_DAA)<-c("male_ID", "haplotype_parentalMale_1", "haplotype_parentalMale_2")
names(table_parental_males_DAA)

######## the following line is to check the genotypes of the males #########

# head(final_haplotypes_offspring_DAA[,2])
# head(table_families[,3])
# 
# !(final_haplotypes_offspring_DAA[,2]%in%table_families[,3])
# 
# head(final_haplotypes_offspring_DAA[(final_haplotypes_offspring_DAA[,2]%in%table_families[,3]),])
# 
# family_num_test=15
# (final_haplotypes_offspring_DAA[!(final_haplotypes_offspring_DAA[,2]%in%table_families[,3]) & final_haplotypes_offspring_DAA$family==family_num_test,])

####################################

temporal<-merge(final_haplotypes_offspring_DAA, table_families, by="family")
final_haplotypes_offspring_DAA<-merge(temporal, table_parental_males_DAA, by="male_ID")
head(final_haplotypes_offspring_DAA)

dim(final_haplotypes_offspring_DAA)

paternal_allele_DAA<-c()
maternal_allele_DAA<-c()
for (i in seq(1,length(final_haplotypes_offspring_DAA[,1]))){
  if (as.character(final_haplotypes_offspring_DAA[i,5])=="Homozygous"){
    paternal_allele_DAA[length(paternal_allele_DAA)+1]<-as.character(final_haplotypes_offspring_DAA[i,6])
    maternal_allele_DAA[length(maternal_allele_DAA)+1]<-as.character(final_haplotypes_offspring_DAA[i,6])
  } else {
    haplotype1<-as.character(final_haplotypes_offspring_DAA[i,6])
    haplotype2<-as.character(final_haplotypes_offspring_DAA[i,7])
    haplotype1_male<-as.character(final_haplotypes_offspring_DAA[i,9])
    haplotype2_male<-as.character(final_haplotypes_offspring_DAA[i,10])
      if (haplotype1==haplotype1_male){
      paternal_allele_DAA[length(paternal_allele_DAA)+1]<-haplotype1
      maternal_allele_DAA[length(maternal_allele_DAA)+1]<-haplotype2
      } else if (haplotype1==haplotype2_male){
      paternal_allele_DAA[length(paternal_allele_DAA)+1]<-haplotype1
      maternal_allele_DAA[length(maternal_allele_DAA)+1]<-haplotype2
      } else if (haplotype2==haplotype1_male){
      paternal_allele_DAA[length(paternal_allele_DAA)+1]<-haplotype2
      maternal_allele_DAA[length(maternal_allele_DAA)+1]<-haplotype1
      } else if (haplotype2==haplotype2_male){
      paternal_allele_DAA[length(paternal_allele_DAA)+1]<-haplotype2
      maternal_allele_DAA[length(maternal_allele_DAA)+1]<-haplotype1
      } else {
      paternal_allele_DAA[length(paternal_allele_DAA)+1]<-"No_match"
      maternal_allele_DAA[length(maternal_allele_DAA)+1]<-"No_match"
      }
  }
}

final_haplotypes_offspring_DAA<-data.frame(final_haplotypes_offspring_DAA,paternal_allele_DAA,maternal_allele_DAA)


code_ind<-c()
family<-c()
ind_fam<-c()
for (i in seq(1,length(final_haplotypes_offspring_DAB[,1]))){
  code_ind[length(code_ind)+1]<-strsplit(as.character(final_haplotypes_offspring_DAB[i,1]),"_")[[1]][1]
  family[length(family)+1]<-strsplit(as.character(final_haplotypes_offspring_DAB[i,1]),"_")[[1]][2]
  ind_fam[length(ind_fam)+1]<-strsplit(as.character(final_haplotypes_offspring_DAB[i,1]),"_")[[1]][3]
}

head(final_haplotypes_offspring_DAB)

final_haplotypes_offspring_DAB<-data.frame(code_ind,family,ind_fam,final_haplotypes_offspring_DAB[,c(2,9,10)])

table_families<-read.table("sample_info.csv", T)
head(table_families)
names(table_families)[3]<-"family"

table_parental_males_DAB<-read.table("males_summary_DAB.csv", T)
head(table_parental_males_DAB)
names(table_parental_males_DAB)<-c("male_ID", "haplotype_parentalMale_1", "haplotype_parentalMale_2")
names(table_parental_males_DAB)


temporal<-merge(final_haplotypes_offspring_DAB, table_families, by="family")
final_haplotypes_offspring_DAB<-merge(temporal, table_parental_males_DAB, by="male_ID")
head(final_haplotypes_offspring_DAB)

dim(final_haplotypes_offspring_DAB)

paternal_allele_DAB<-c()
maternal_allele_DAB<-c()
for (i in seq(1,length(final_haplotypes_offspring_DAB[,1]))){
  if (as.character(final_haplotypes_offspring_DAB[i,5])=="Homozygous"){
    paternal_allele_DAB[length(paternal_allele_DAB)+1]<-as.character(final_haplotypes_offspring_DAB[i,6])
    maternal_allele_DAB[length(maternal_allele_DAB)+1]<-as.character(final_haplotypes_offspring_DAB[i,6])
  } else {
    haplotype1<-as.character(final_haplotypes_offspring_DAB[i,6])
    haplotype2<-as.character(final_haplotypes_offspring_DAB[i,7])
    haplotype1_male<-as.character(final_haplotypes_offspring_DAB[i,9])
    haplotype2_male<-as.character(final_haplotypes_offspring_DAB[i,10])
    if (haplotype1==haplotype1_male){
      paternal_allele_DAB[length(paternal_allele_DAB)+1]<-haplotype1
      maternal_allele_DAB[length(maternal_allele_DAB)+1]<-haplotype2
    } else if (haplotype1==haplotype2_male){
      paternal_allele_DAB[length(paternal_allele_DAB)+1]<-haplotype1
      maternal_allele_DAB[length(maternal_allele_DAB)+1]<-haplotype2
    } else if (haplotype2==haplotype1_male){
      paternal_allele_DAB[length(paternal_allele_DAB)+1]<-haplotype2
      maternal_allele_DAB[length(maternal_allele_DAB)+1]<-haplotype1
    } else if (haplotype2==haplotype2_male){
      paternal_allele_DAB[length(paternal_allele_DAB)+1]<-haplotype2
      maternal_allele_DAB[length(maternal_allele_DAB)+1]<-haplotype1
    } else {
      paternal_allele_DAB[length(paternal_allele_DAB)+1]<-"No_match"
      maternal_allele_DAB[length(maternal_allele_DAB)+1]<-"No_match"
    }
  }
}


final_haplotypes_offspring_DAB<-data.frame(final_haplotypes_offspring_DAB,paternal_allele_DAB,maternal_allele_DAB)





final_haplotypes_offspring_DAA<-final_haplotypes_offspring_DAA[,c(1:10)]
final_haplotypes_offspring_DAB<-final_haplotypes_offspring_DAB[,c(1:10)]



head(final_haplotypes_offspring_DAA)
head(final_haplotypes_offspring_DAB)


table_parental_females_DAA<-read.table("females_summary_DAA", T)
head(table_parental_females_DAA)
names(table_parental_females_DAA)<-c("female_ID", "haplotype_maternalFemale_1", "haplotype_maternalFemale_2")
names(table_parental_females_DAA)


final_haplotypes_offspring_DAA<-merge(final_haplotypes_offspring_DAA, table_parental_females_DAA, by="female_ID")
head(final_haplotypes_offspring_DAA)



table_parental_females_DAB<-read.table("females_summary_DAB", T)
head(table_parental_females_DAB)
names(table_parental_females_DAB)<-c("female_ID", "haplotype_maternalFemale_1", "haplotype_maternalFemale_2")
names(table_parental_females_DAB)


final_haplotypes_offspring_DAB<-merge(final_haplotypes_offspring_DAB, table_parental_females_DAB, by="female_ID")
head(final_haplotypes_offspring_DAB)


paternal_allele_DAA<-c()
maternal_allele_DAA<-c()
for (i in seq(1,length(final_haplotypes_offspring_DAA[,1]))){
  haplotype1<-as.character(final_haplotypes_offspring_DAA[i,7])
  haplotype2<-as.character(final_haplotypes_offspring_DAA[i,8])
  paternal_hap_1<-as.character(final_haplotypes_offspring_DAA[i,9])
  paternal_hap_2<-as.character(final_haplotypes_offspring_DAA[i,10])
  maternal_hap_1<-as.character(final_haplotypes_offspring_DAA[i,11])
  maternal_hap_2<-as.character(final_haplotypes_offspring_DAA[i,12])
  if ((haplotype1 %in% c(paternal_hap_1,paternal_hap_2)) & (haplotype2 %in% c(maternal_hap_1,maternal_hap_2))){
    paternal_allele_DAA[length(paternal_allele_DAA)+1]<-haplotype1
    maternal_allele_DAA[length(maternal_allele_DAA)+1]<-haplotype2
  } else if ((haplotype2 %in% c(paternal_hap_1,paternal_hap_2)) & (haplotype1 %in% c(maternal_hap_1,maternal_hap_2))) {
    paternal_allele_DAA[length(paternal_allele_DAA)+1]<-haplotype2
    maternal_allele_DAA[length(maternal_allele_DAA)+1]<-haplotype1    
  } else {
    paternal_allele_DAA[length(paternal_allele_DAA)+1]<-"No_match"
    maternal_allele_DAA[length(maternal_allele_DAA)+1]<-"No_match"
  }
}
  


paternal_allele_DAB<-c()
maternal_allele_DAB<-c()
for (i in seq(1,length(final_haplotypes_offspring_DAB[,1]))){
  haplotype1<-as.character(final_haplotypes_offspring_DAB[i,7])
  haplotype2<-as.character(final_haplotypes_offspring_DAB[i,8])
  paternal_hap_1<-as.character(final_haplotypes_offspring_DAB[i,9])
  paternal_hap_2<-as.character(final_haplotypes_offspring_DAB[i,10])
  maternal_hap_1<-as.character(final_haplotypes_offspring_DAB[i,11])
  maternal_hap_2<-as.character(final_haplotypes_offspring_DAB[i,12])
  if ((haplotype1 %in% c(paternal_hap_1,paternal_hap_2)) & (haplotype2 %in% c(maternal_hap_1,maternal_hap_2))){
    paternal_allele_DAB[length(paternal_allele_DAB)+1]<-haplotype1
    maternal_allele_DAB[length(maternal_allele_DAB)+1]<-haplotype2
  } else if ((haplotype2 %in% c(paternal_hap_1,paternal_hap_2)) & (haplotype1 %in% c(maternal_hap_1,maternal_hap_2))) {
    paternal_allele_DAB[length(paternal_allele_DAB)+1]<-haplotype2
    maternal_allele_DAB[length(maternal_allele_DAB)+1]<-haplotype1    
  } else {
    paternal_allele_DAB[length(paternal_allele_DAB)+1]<-"No_match"
    maternal_allele_DAB[length(maternal_allele_DAB)+1]<-"No_match"
  }
}


final_haplotypes_offspring_DAA<-data.frame(final_haplotypes_offspring_DAA,paternal_allele_DAA,maternal_allele_DAA)
final_haplotypes_offspring_DAB<-data.frame(final_haplotypes_offspring_DAB,paternal_allele_DAB,maternal_allele_DAB)


head(final_haplotypes_offspring_DAA)
head(final_haplotypes_offspring_DAB)


final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$paternal_allele_DAA=="No_match",]
final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$paternal_allele_DAB=="No_match",]

final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$family=="11",]
final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$male_ID=="12",]
final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$female_ID=="10",]

length(final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$family=="181",1])


### exclude repited samples

head(final_haplotypes_offspring_DAA)

f=181
head(famility_table)
i=1
c(unique(famility_table$ind_fam))
unique_ind_codes_DAA

famility_table
as.vector(famility_table[famility_table$ind_fam==i,4][1])

unique_ind_codes_DAA<-c()
for(f in c(levels(final_haplotypes_offspring_DAA$family))){
  if (length(final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$family==f,1])>0){
    famility_table<-final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$family==f,]
    for (i in c(unique(famility_table$ind_fam))){
      unique_ind_codes_DAA<-c(unique_ind_codes_DAA,as.vector(famility_table[famility_table$ind_fam==i,4][1]))
    }
  }
}

unique_ind_codes_DAB<-c()
for(f in c(levels(final_haplotypes_offspring_DAB$family))){
  if (length(final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$family==f,1])>0){
    famility_table<-final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$family==f,]
    for (i in c(unique(famility_table$ind_fam))){
      unique_ind_codes_DAB<-c(unique_ind_codes_DAB,as.vector(famility_table[famility_table$ind_fam==i,4][1]))
    }
  }
}

final_haplotypes_offspring_DAA<-final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$code_ind %in% unique_ind_codes_DAA,]
final_haplotypes_offspring_DAB<-final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$code_ind %in% unique_ind_codes_DAB,]



##### Information for families

female_ID<-c()
male_ID<-c()
family<-c()
hap_pat_1<-c()
hap_pat_2<-c()
hap_mat_1<-c()
hap_mat_2<-c()
hap_pat_1_num<-c()
hap_pat_2_num<-c()
hap_mat_1_num<-c()
hap_mat_2_num<-c()
gen_1_1_num<-c()
gen_1_2_num<-c()
gen_2_1_num<-c()
gen_2_2_num<-c()

for(i in c(levels(final_haplotypes_offspring_DAA$family))){
  if (length(final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$family==i,1])>0){
    female_ID<-c(female_ID,final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$family==i,1][1])
    male_ID<-c(male_ID,final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$family==i,2][1])
    family<-c(family,i)
    table_tem<-final_haplotypes_offspring_DAA[final_haplotypes_offspring_DAA$family==i,]
    hap_pat_1<-c(hap_pat_1,as.character(table_tem[,9][1]))
    hap_pat_2<-c(hap_pat_2,as.character(table_tem[,10][1]))
    hap_mat_1<-c(hap_mat_1,as.character(table_tem[,11][1]))
    hap_mat_2<-c(hap_mat_2,as.character(table_tem[,12][1]))
    if (as.character(table_tem[,9][1])!=as.character(table_tem[,10][1])){
      if (as.character(table_tem[,11][1])!=as.character(table_tem[,12][1])){
        hap_pat_1_num<-c(hap_pat_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]),1]))
        hap_pat_2_num<-c(hap_pat_2_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,10][1]),1]))
        hap_mat_1_num<-c(hap_mat_1_num, length(table_tem[table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        hap_mat_2_num<-c(hap_mat_2_num, length(table_tem[table_tem$maternal_allele_DAA==as.character(table_tem[,12][1]),1]))
        gen_1_1_num<-c(gen_1_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        gen_1_2_num<-c(gen_1_2_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,12][1]),1]))
        gen_2_1_num<-c(gen_2_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,10][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        gen_2_2_num<-c(gen_2_2_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,10][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,12][1]),1]))
      } else {
        hap_pat_1_num<-c(hap_pat_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]),1]))
        hap_pat_2_num<-c(hap_pat_2_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,10][1]),1]))
        hap_mat_1_num<-c(hap_mat_1_num, length(table_tem[table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        hap_mat_2_num<-c(hap_mat_2_num, 0)
        gen_1_1_num<-c(gen_1_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        gen_1_2_num<-c(gen_1_2_num, 0)
        gen_2_1_num<-c(gen_2_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,10][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        gen_2_2_num<-c(gen_2_2_num, 0)
      }
    } else {
      if (as.character(table_tem[,11][1])!=as.character(table_tem[,12][1])){
        hap_pat_1_num<-c(hap_pat_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]),1]))
        hap_pat_2_num<-c(hap_pat_2_num, 0)
        hap_mat_1_num<-c(hap_mat_1_num, length(table_tem[table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        hap_mat_2_num<-c(hap_mat_2_num, length(table_tem[table_tem$maternal_allele_DAA==as.character(table_tem[,12][1]),1]))
        gen_1_1_num<-c(gen_1_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        gen_1_2_num<-c(gen_1_2_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,12][1]),1]))
        gen_2_1_num<-c(gen_2_1_num, 0)
        gen_2_2_num<-c(gen_2_2_num, 0)
      } else {
        hap_pat_1_num<-c(hap_pat_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]),1]))
        hap_pat_2_num<-c(hap_pat_2_num, 0)
        hap_mat_1_num<-c(hap_mat_1_num, length(table_tem[table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        hap_mat_2_num<-c(hap_mat_2_num, 0)
        gen_1_1_num<-c(gen_1_1_num, length(table_tem[table_tem$paternal_allele_DAA==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAA==as.character(table_tem[,11][1]),1]))
        gen_1_2_num<-c(gen_1_2_num, 0)
        gen_2_1_num<-c(gen_2_1_num, 0)
        gen_2_2_num<-c(gen_2_2_num, 0)
      }
    }    
  }
}


summary_genotypes_fam_DAA<-data.frame(female_ID, male_ID, family, 
                                      hap_pat_1, hap_pat_2, hap_mat_1, hap_mat_2, 
                                      hap_pat_1_num, freq_hap_pat_1=c(hap_pat_1_num/(hap_pat_1_num+hap_pat_2_num)),
                                      hap_pat_2_num, freq_hap_pat_2=c(hap_pat_2_num/(hap_pat_1_num+hap_pat_2_num)),
                                      hap_mat_1_num, freq_hap_mat_1=c(hap_mat_1_num/(hap_mat_1_num+hap_mat_2_num)),
                                      hap_mat_2_num, freq_hap_mat_2=c(hap_mat_2_num/(hap_mat_1_num+hap_mat_2_num)),
                                      gen_1_1_num, gen_1_2_num, gen_2_1_num, gen_2_2_num, 
                                      freq_1_1=c(gen_1_1_num/(gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)), 
                                      freq_1_2=c(gen_1_2_num/(gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)), 
                                      freq_2_1=c(gen_2_1_num/(gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)), 
                                      freq_2_2=c(gen_2_2_num/(gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)), 
                                      total_Num_ind=c((gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)))



female_ID<-c()
male_ID<-c()
family<-c()
hap_pat_1<-c()
hap_pat_2<-c()
hap_mat_1<-c()
hap_mat_2<-c()
hap_pat_1_num<-c()
hap_pat_2_num<-c()
hap_mat_1_num<-c()
hap_mat_2_num<-c()
gen_1_1_num<-c()
gen_1_2_num<-c()
gen_2_1_num<-c()
gen_2_2_num<-c()


for(i in c(levels(final_haplotypes_offspring_DAB$family))){
  if (length(final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$family==i,1])>0){
    female_ID<-c(female_ID,final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$family==i,1][1])
    male_ID<-c(male_ID,final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$family==i,2][1])
    family<-c(family,i)
    table_tem<-final_haplotypes_offspring_DAB[final_haplotypes_offspring_DAB$family==i,]
    hap_pat_1<-c(hap_pat_1,as.character(table_tem[,9][1]))
    hap_pat_2<-c(hap_pat_2,as.character(table_tem[,10][1]))
    hap_mat_1<-c(hap_mat_1,as.character(table_tem[,11][1]))
    hap_mat_2<-c(hap_mat_2,as.character(table_tem[,12][1]))
    if (as.character(table_tem[,9][1])!=as.character(table_tem[,10][1])){
      if (as.character(table_tem[,11][1])!=as.character(table_tem[,12][1])){
        hap_pat_1_num<-c(hap_pat_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]),1]))
        hap_pat_2_num<-c(hap_pat_2_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,10][1]),1]))
        hap_mat_1_num<-c(hap_mat_1_num, length(table_tem[table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        hap_mat_2_num<-c(hap_mat_2_num, length(table_tem[table_tem$maternal_allele_DAB==as.character(table_tem[,12][1]),1]))
        gen_1_1_num<-c(gen_1_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        gen_1_2_num<-c(gen_1_2_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,12][1]),1]))
        gen_2_1_num<-c(gen_2_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,10][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        gen_2_2_num<-c(gen_2_2_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,10][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,12][1]),1]))
      } else {
        hap_pat_1_num<-c(hap_pat_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]),1]))
        hap_pat_2_num<-c(hap_pat_2_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,10][1]),1]))
        hap_mat_1_num<-c(hap_mat_1_num, length(table_tem[table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        hap_mat_2_num<-c(hap_mat_2_num, 0)
        gen_1_1_num<-c(gen_1_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        gen_1_2_num<-c(gen_1_2_num, 0)
        gen_2_1_num<-c(gen_2_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,10][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        gen_2_2_num<-c(gen_2_2_num, 0)
      }
    } else {
      if (as.character(table_tem[,11][1])!=as.character(table_tem[,12][1])){
        hap_pat_1_num<-c(hap_pat_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]),1]))
        hap_pat_2_num<-c(hap_pat_2_num, 0)
        hap_mat_1_num<-c(hap_mat_1_num, length(table_tem[table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        hap_mat_2_num<-c(hap_mat_2_num, length(table_tem[table_tem$maternal_allele_DAB==as.character(table_tem[,12][1]),1]))
        gen_1_1_num<-c(gen_1_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        gen_1_2_num<-c(gen_1_2_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,12][1]),1]))
        gen_2_1_num<-c(gen_2_1_num, 0)
        gen_2_2_num<-c(gen_2_2_num, 0)
      } else {
        hap_pat_1_num<-c(hap_pat_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]),1]))
        hap_pat_2_num<-c(hap_pat_2_num, 0)
        hap_mat_1_num<-c(hap_mat_1_num, length(table_tem[table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        hap_mat_2_num<-c(hap_mat_2_num, 0)
        gen_1_1_num<-c(gen_1_1_num, length(table_tem[table_tem$paternal_allele_DAB==as.character(table_tem[,9][1]) & table_tem$maternal_allele_DAB==as.character(table_tem[,11][1]),1]))
        gen_1_2_num<-c(gen_1_2_num, 0)
        gen_2_1_num<-c(gen_2_1_num, 0)
        gen_2_2_num<-c(gen_2_2_num, 0)
      }
    }    
  }
}


summary_genotypes_fam_DAB<-data.frame(female_ID, male_ID, family, 
                                      hap_pat_1, hap_pat_2, hap_mat_1, hap_mat_2, 
                                      hap_pat_1_num, freq_hap_pat_1=c(hap_pat_1_num/(hap_pat_1_num+hap_pat_2_num)),
                                      hap_pat_2_num, freq_hap_pat_2=c(hap_pat_2_num/(hap_pat_1_num+hap_pat_2_num)),
                                      hap_mat_1_num, freq_hap_mat_1=c(hap_mat_1_num/(hap_mat_1_num+hap_mat_2_num)),
                                      hap_mat_2_num, freq_hap_mat_2=c(hap_mat_2_num/(hap_mat_1_num+hap_mat_2_num)),
                                      gen_1_1_num, gen_1_2_num, gen_2_1_num, gen_2_2_num, 
                                      freq_1_1=c(gen_1_1_num/(gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)), 
                                      freq_1_2=c(gen_1_2_num/(gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)), 
                                      freq_2_1=c(gen_2_1_num/(gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)), 
                                      freq_2_2=c(gen_2_2_num/(gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)), 
                                      total_Num_ind=c((gen_1_1_num+gen_1_2_num+gen_2_1_num+gen_2_2_num)))



summary_genotypes_fam_DAA
summary_genotypes_fam_DAB

# write.table(summary_genotypes_fam_DAA, file = "summary_genotypes_fam_DAA.txt", append = FALSE, quote = F, 
#             col.names = T, row.names = F)
# 
# write.table(summary_genotypes_fam_DAB, file = "summary_genotypes_fam_DAB.txt", append = FALSE, quote = F, 
#             col.names = T, row.names = F)



