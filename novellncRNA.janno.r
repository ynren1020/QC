###################2019-05-01#######################
##5 novel lncRNAs first exon annotation#############
##+ exon1end-exon2start#############################
##- exon.(last-1).end-exon.last.start###############
####################################################
library(dplyr)
library(tidyr)

novel<-read.delim("novel5lncrna_hg38_filter.gtf",header = FALSE,stringsAsFactors = FALSE)

novel<-novel%>%separate("V9",c("V9","V10","V11"),";") #" exon_number 3"

##count exons##
novel.sum<-novel%>%group_by(V9,V10)%>%
  summarize(exon=n_distinct(V11))

novel<-left_join(novel,novel.sum,by=c("V9","V10"))

##2 exons##
novel2<-filter(novel,exon==2)
novel2.1<-filter(novel2,V11== " exon_number 1")
novel2.1<-novel2.1[,c("V1","V5","V9","V10","V11")]

novel2.2<-filter(novel2,V11== " exon_number 2")
novel2.2<-novel2.2[,c("V1","V4","V9","V10","V11")]

novel2.janno<-full_join(novel2.1,novel2.2,by=c("V9","V10","V1"))%>%rename("chrom"="V1","start"="V5","end"="V4")
novel2.janno.final<-novel2.janno[,c("chrom","start","end","V9","V10")]

##3 exons##
novel3<-filter(novel,exon==3)
novel3.1<-filter(novel3,V11== " exon_number 2")
novel3.1<-novel3.1[,c("V1","V5","V9","V10","V11")]

novel3.2<-filter(novel3,V11== " exon_number 3")
novel3.2<-novel3.2[,c("V1","V4","V9","V10","V11")]

novel3.janno<-full_join(novel3.1,novel3.2,by=c("V9","V10","V1"))%>%rename("chrom"="V1","start"="V5","end"="V4")
novel3.janno.final<-novel3.janno[,c("chrom","start","end","V9","V10")]

##final##
novel.janno.final<-rbind(novel2.janno.final,novel3.janno.final)
novel.janno.final<-unique(novel.janno.final[,1:4])

writein(novel.janno.final,TRUE,FALSE)

