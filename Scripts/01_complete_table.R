# ==================================================
# Sergio Tusso
# Simone Immler Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014
# +++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls(all=TRUE))
# setwd("/home/evuser/Dropbox/Uppsala/Salmon/Samples info")
setwd("/home/evuser/Desktop/Salmon/Sample_information")
# setwd("C:/Users/SITG/Dropbox/Uppsala/salmon/Samples\ info")

Ind_data<-read.table("Table_ind_code.txt", T)
Primer_data<-read.table("Table_primers.txt", T)
Tags_data<-read.table("Table_tags.txt", T)

head(Ind_data)
head(Primer_data)
head(Tags_data)
str(Ind_data)
str(Primer_data)
str(Tags_data)

# these lines change the class data in each column

Ind_data[,1]<-as.factor(Ind_data[,1])
Ind_data[,3]<-as.factor(Ind_data[,3])
Ind_data[,6]<-as.factor(Ind_data[,6])

Primer_data[,2]<-as.factor(Primer_data[,2])

# this line deletes all lines without sample
Ind_data<-Ind_data[Ind_data$Ind_code!="0",]

# this loop identifies the forward and reverse primer for each sample
forward_code<-c()
reverse_code<-c()
for (i in seq(1, length(Ind_data[,1]))){
  forward_code<-c(forward_code,as.character(Primer_data[Primer_data$Loc==as.character(Ind_data[i,4]) & Primer_data$Plate==as.character(Ind_data[i,5]) & Primer_data$Col==as.character(Ind_data[i,6]) & Primer_data$Row==as.character(Ind_data[i,7]),7]))
  reverse_code<-c(reverse_code,as.character(Primer_data[Primer_data$Loc==as.character(Ind_data[i,4]) & Primer_data$Plate==as.character(Ind_data[i,5]) & Primer_data$Col==as.character(Ind_data[i,6]) & Primer_data$Row==as.character(Ind_data[i,7]),8]))
}

Ind_data<-cbind(Ind_data, forward_code, reverse_code)




# this loop identifies the forward and reverse code for each sample and give the tag and primer sequences
sequences_F<-c("SFRD","SFRD","SFRD")
sequences_R<-c("SFRD","SFRD","SFRD")
for (i in seq(1, length(Ind_data[,1]))){
  sequences_F<-rbind(sequences_F,(Tags_data[Tags_data$Locus==as.character(Ind_data[i,4]) & Tags_data$primer==as.character(Ind_data[i,8]),seq(2,4)]))
  sequences_R<-rbind(sequences_R,(Tags_data[Tags_data$Locus==as.character(Ind_data[i,4]) & Tags_data$primer==as.character(Ind_data[i,9]),seq(2,4)]))
}

sequences_F<-sequences_F[-1,]
sequences_R<-sequences_R[-1,]

Ind_data<-cbind(Ind_data, 
                w_seq_f=sequences_F[,1], tag_f=sequences_F[,2], locus_seq_f=sequences_F[,3], w_seq_f_editted=paste(tag_f=sequences_F[,2], locus_seq_f=sequences_F[,3], sep=""),
                w_seq_r=sequences_R[,1], tag_r=sequences_R[,2], locus_seq_r=sequences_R[,3], w_seq_r_concatenated=paste(tag_r=sequences_R[,2], locus_seq_r=sequences_R[,3], sep=""))

head(Ind_data)
  
# this line save the table in a text file

write.table(Ind_data, file = "complete_table.txt", append = FALSE, quote = F, 
            col.names = T, row.names = F, sep="  ")

# the following lines should be run in batch before continue!!!
# 
# awk '//{print $17}' complete_table.txt > seq_ori_t.txt
# grep -v w_seq_r seq_ori_t.txt > seq_ori.txt
# awk 'BEGIN {
#   j = n = split("A C G T a c g t", t)
#   for (i = 0; ++i <= n;)
#     map[t[i]] = t[j--]
#   }
# {
#   if (/^>/) print
#   else {
#     for (i = length; i; i--)
#     printf "%s", map[substr($0, i, 1)]
#     print x
# 	}	  
#   }' seq_ori.txt > rev_seq_ori.txt


head(Ind_data)
# 
# REVERSE_COM_w_seq_r_concatenated<-read.table("rev_seq_ori.txt", F)
# Ind_data<-data.frame(Ind_data, REVERSE_COM_w_seq_r_concatenated)
# 
# write.table(Ind_data, file = "complete_table_editted_tab.csv", append = FALSE, quote = F, 
#             col.names = T, row.names = F)

