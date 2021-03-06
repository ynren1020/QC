############################2019-04-30################################
##choose top20 genes and other 7 lncrnas for test#####################
##prad,and su2c#######################################################
######################################################################

#library(tidyverse)
#library(ggpubr)
library(dplyr)
library(tidyr)

#################some genes' junction cannot match to first junction position annotation####################

##first junction annotation file##
#args<-commandArgs(TRUE)
#args[1]<-

hg38dancr<-read.delim("hg38_junction_HAVANA.ENSEMBL.EEF2.txt",header = FALSE,stringsAsFactors = FALSE)
#hg38malat1<-read.delim("hg38_junction_HAVANA.ENSEMBL.MALAT1.txt",header = FALSE,stringsAsFactors = FALSE)
##su2c_hybrid janno files only contain DANCR##
hyb_dancr<-read.delim("su2c.janno.EEF2.txt",header = FALSE,stringsAsFactors = FALSE)
#hyb_malat1<-read.delim("su2c_hybrid_MALAT1.txt",header=FALSE,stringsAsFactors = FALSE)

hyb_dancr<-separate(hyb_dancr,V1,c("sample","chrom"),":")
#hyb_malat1<-separate(hyb_malat1,V1,c("sample","chrom"),":")

##join##
hyb_dancr.join<-left_join(hg38dancr,hyb_dancr,by=c("V1"="chrom","V4"="V2","V5"="V3"))
#hyb_malat1.join<-left_join(hg38malat1,hyb_malat1,by=c("V1"="chrom","V4"="V2","V5"="V3"))

##select sample:V16#
hyb_dancr.join<-hyb_dancr.join[,11:24]
hyb_dancr.join<-unique(hyb_dancr.join)
hyb_dancr.join<-hyb_dancr.join%>%group_by(sample)%>%
  summarize(score=sum(V5.y))

##total junction##
hybtotal<-read.delim("su2c.total.janno.txt",header = FALSE,stringsAsFactors = FALSE)
##join##
hyb_dancr.join.total<-full_join(hyb_dancr.join,hybtotal,by=c("sample"="V1"))

hyb_dancr.join.total<-hyb_dancr.join.total%>%
  mutate(scoreN=score/V2*1000000)%>%
  mutate(scoreNlog=log10(scoreN))

hyb_dancr.join.total<-na.omit(hyb_dancr.join.total)
write.table(hyb_dancr.join.total,"su2c_hyb_dancr.join.total.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")


