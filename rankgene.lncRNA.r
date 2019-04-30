#########2019-04-29##########################################
#########################QC BMI1 data########################
##rank genes/lncRNAs (known or novel) by read counts#########
#############################################################

library("biomaRt")
library(dplyr)
library(tidyr)
library(stringr)

BMI1<-read.delim("BMI1_id_counts.txt",stringsAsFactors = FALSE,skip = 1,header = TRUE)
##calculate rpkm##
total1<-sum(BMI1$BMI1_1_sorted.bam) #13606540
total2<-sum(BMI1$BMI1_2_sorted.bam) #13733602

BMI1$rpkm1<-BMI1$BMI1_1_sorted.bam*10^9/total1/BMI1$Length
BMI1$rpkm2<-BMI1$BMI1_2_sorted.bam*10^9/total2/BMI1$Length

BMI1<-BMI1%>%
arrange(desc(rpkm1),desc(rpkm1))

##rank##
BMI1$rank1<-rank(-BMI1$rpkm1)
BMI1$rank2<-rank(-BMI1$rpkm2)

##choose top 200s##
BMI1top200<-filter(BMI1,rank1<=200|rank2<=200)

#####to gene symbol###
BMI1known<-BMI1top200 %>%
  filter(str_detect(Geneid, "ENSG"))
BMI1known<-BMI1known%>%
  mutate(gene_name=str_split_fixed(BMI1known$Geneid,"[.]",n=2)[,1])


ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
genesym<-getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=BMI1known$gene_name, mart=ensembl)
#BMI1name<-BMI1known[BMI1known$gene_name%in%genesym$ensembl_gene_id,]
BMI1namefinal<-full_join(BMI1known,genesym,by=c("gene_name"="ensembl_gene_id"))
##known gene and lncRNAs##
BMI1namefinalsub<-BMI1namefinal%>%rename("gene"="hgnc_symbol") 

##novel lncRNA##
BMI1lncRNA<-BMI1top200 %>%
  filter(str_detect(Geneid, "^G."))%>%
  mutate(gene_name=Geneid)%>%
  mutate(gene=Geneid)

##combine top200 gene and lncRNA##
BMI1all<-rbind(BMI1namefinalsub,BMI1lncRNA)

geneid_top200<-as.data.frame(BMI1all$Geneid)
  
write.table(geneid_top200,"geneid_top200.hg38.txt",col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")


##top200 gtf##
##read in rank200 gtf file##
rand200gtf<-read.delim("hg38_genelnRNA_sorted_top200.gtf",header = FALSE,stringsAsFactors = FALSE)
rand200gtf<-rand200gtf%>%separate(V9,c(paste0("V",9:26)),";")

rand200gtfS<-filter(rand200gtf,V2=="StringTie")
rand200gtfS<-rand200gtfS%>%filter(V3=="transcript")
rand200gtfS<-rand200gtfS[,1:10]
rand200gtfS$gene_id<-str_replace_all(rand200gtfS$V9, fixed("gene_id "), "")%>%
  str_replace_all(fixed(" "), "")
rand200gtfS<-rand200gtfS[,c(1,4,5,7,9,11)]

rand200gtfH<-filter(rand200gtf,V2!="StringTie")
rand200gtfH<-filter(rand200gtfH,V3=="gene")
rand200gtfH$gene_id<-str_replace_all(rand200gtfH$V10, fixed("gene_id "), "")%>%
  str_replace_all(fixed(" "), "")
rand200gtfH<-rand200gtfH[,-c(21:26)]
rand200gtfH<-rand200gtfH[,c(1,4,5,7,11,21)]
rand200gtfH<-rand200gtfH%>%rename("V9"="V11")

##combine##
rank200_allgtf<-rbind(rand200gtfS,rand200gtfH)

rank200_allgtf.rpkm<-full_join(rank200_allgtf,BMI1all,by=c("gene_id"="Geneid"))

rank200_allgtf.rpkm<-unique(rank200_allgtf.rpkm)

##output for top 200 genes and lncRNAs##
write.table(rank200_allgtf.rpkm,"rank200_genelncRNA.gtf.rpkm.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")


##check the Geneid in BMI1namefinalsub is lncRNA or genes##
####gtf file find lncRNAs,produced from lncRNAs.py with genecode.v21.annotation.gtf##
gtf<-read.delim("lncRNAs_YR.gtf",header = FALSE,stringsAsFactors = FALSE)
gtf<-gtf%>%separate(V9,c(paste0("V",9:26)),";")
gtfvector<-unlist(gtf)
gtf_geneid<-str_subset(gtfvector,"ENSG")
##lncRNAs geneid##
gtf_geneid_u<-str_replace_all(gtf_geneid,"gene_id ","")
gtf_geneid_u<-as.data.frame(gtf_geneid_u)
gtf_geneid_u$gtf_geneid_u<-as.character(gtf_geneid_u$gtf_geneid_u)
gtf_geneid_u<-unique(gtf_geneid_u)

for (i in 1:nrow(gtf_geneid_u)){
  gtf_geneid_u$gene_id[i]<-strsplit(gtf_geneid_u$gtf_geneid_u[i],'[.]')[[1]][1]
}
##remove leading whitespace##
gtf_geneid_u$gene_id<-str_replace_all(gtf_geneid_u$gene_id, fixed(" "), "")

##novel lncrna##
#novel<-c("G.5257","G.5260","G.7790")

##rank200_allgtf.rpkm lncRNAs##
#for (i in 1:nrow(rank200_allgtf.rpkm)){
#  rank200_allgtf.rpkm$gene_type[i]<-ifelse(rank200_allgtf.rpkm$gene_name[i]%in%gtf_geneid_u$gene_id|rank200_allgtf.rpkm$gene_name[i]%in%novel,"lncRNA","protein_coding gene")
#}

##join##
##no lncRNAS##only three novel lncRNA and known genes##
rank200_allgtf.rpkm.type<-left_join(rank200_allgtf.rpkm,gtf_geneid_u,by=c("gene_name"="gene_id"))


#########with gene_type###############
##top200 gtf##
##read in rank200 gtf file##
rand200gtf<-read.delim("hg38_genelnRNA_sorted_top200.gtf",header = FALSE,stringsAsFactors = FALSE)
rand200gtf<-rand200gtf%>%separate(V9,c(paste0("V",9:26)),";")

rand200gtfS<-filter(rand200gtf,V2=="StringTie")
rand200gtfS<-rand200gtfS%>%filter(V3=="transcript")
#rand200gtfS<-rand200gtfS[,1:10]
rand200gtfS$gene_id<-str_replace_all(rand200gtfS$V9, fixed("gene_id "), "")%>%
  str_replace_all(fixed(" "), "")

#rand200gtfS<-rand200gtfS[,c(1,4,5,7,9,11)]

rand200gtfH<-filter(rand200gtf,V2!="StringTie")
rand200gtfH<-filter(rand200gtfH,V3=="gene")
rand200gtfH$gene_id<-str_replace_all(rand200gtfH$V10, fixed("gene_id "), "")%>%
  str_replace_all(fixed(" "), "")
#rand200gtfH<-rand200gtfH[,-c(21:26)]
#rand200gtfH<-rand200gtfH[,c(1,4,5,7,11,21)]
#rand200gtfH<-rand200gtfH%>%rename("V9"="V11")

##combine##
rank200_allgtf<-rbind(rand200gtfS,rand200gtfH)

rank200_allgtf.rpkm<-full_join(rank200_allgtf,BMI1all,by=c("gene_id"="Geneid"))

rank200_allgtf.rpkm<-unique(rank200_allgtf.rpkm)

##output for top 200 genes and lncRNAs##
write.table(rank200_allgtf.rpkm,"rank200_genelncRNA.gtf.rpkm.type.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")

