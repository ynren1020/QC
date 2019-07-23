#########2019-07-23##########################################
#########################QC BMI1 data########################
##separate BMI1 expressed genes(rpkm>1) into four groups#####
#####################pie chart ##############################
#############################################################

#RPKM - Reads per kilo base per million mapped reads
#Formula
#RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )

#numReads - number of reads mapped to a gene sequence
#geneLength - length of the gene sequence
#totalNumReads - total number of mapped reads of a sample

library("biomaRt")
library(dplyr)
library(tidyr)
library(stringr)

BMI1<-read.delim("BMI1_id_counts.txt",stringsAsFactors = FALSE,skip = 1,header = TRUE) #60211
##RPKM##
BMI1<-BMI1%>%
  mutate(BMI1_1_rpkm=BMI1_1_sorted.bam/(Length/1000*sum(BMI1_1_sorted.bam)/1000000),BMI1_2_rpkm=BMI1_2_sorted.bam/(Length/1000*sum(BMI1_2_sorted.bam)/1000000))%>%
  arrange(BMI1_1_rpkm,BMI1_2_rpkm)

BMI1sub<-filter(BMI1,BMI1_1_rpkm>1&BMI1_2_rpkm>1) #10600 all expressed genes (include novel lncrna)

##10570 known genes##30 G.*** novel lncrnas##
BMI1known<-BMI1sub %>%
  filter(str_detect(Geneid, "ENSG"))
BMI1known<-BMI1known%>%
  mutate(gene_name=str_split_fixed(BMI1known$Geneid,"[.]",n=2)[,1])

BMI1_expressed_knowngenes<-as.data.frame(BMI1known$gene_name)

write.table(BMI1_expressed_knowngenes,"BMI1_expressed_knowngenes_0723.txt",sep="\t",quote = FALSE,col.names = FALSE,row.names = FALSE)

##9315 proteincoding;779 known lncrnas; 110 small rnas; 10570-10204=366 others;30 novel lncRNAs.

types<-c("protein coding","known lncRNAs","small RNAs","novel lncRNAs","other")
counts<-c(9315,779,110,30,366)

##combine novel lncRNA and known lncRNAs##
types<-c("protein coding","lncRNAs","small RNAs","pseudogene")
counts<-c(9315,809,110,366)

piechart<-data.frame(genetype=types,counts=counts,lbs=paste0(types, ",",counts))

pct<-round(counts/sum(counts)*100,1)

pie(piechart$counts,labels = piechart$lbs, col=rainbow(length(piechart$counts)),
    main="Pie Chart of gene types")


##other types##
input<-"BMI1_expressed_knowngenes_0723_others_geneid.gtf"
others<-read.delim(input,header = FALSE,stringsAsFactors = FALSE)
datsub<-separate(others,V9,c(paste0("V",9:20)),sep=";")





                    