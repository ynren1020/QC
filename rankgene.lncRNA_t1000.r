#########2019-04-29##########################################
#########################QC BMI1 data########################
##rank genes/lncRNAs (known or novel) by read counts#########
##choose top 1000############################################
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

##choose top 1000s##
BMI1top1000<-filter(BMI1,rank1<=1000|rank2<=1000)

#####to gene symbol###
BMI1known<-BMI1top1000 %>%
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
BMI1lncRNA<-BMI1top1000 %>%
  filter(str_detect(Geneid, "^G."))%>%
  mutate(gene_name=Geneid)%>%
  mutate(gene=Geneid)

##combine top1000 gene and lncRNA##
BMI1all<-rbind(BMI1namefinalsub,BMI1lncRNA)

geneid_top1000<-as.data.frame(BMI1all$Geneid)

write.table(geneid_top1000,"geneid_top1000.hg38.txt",col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")





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

##rank1000_allgtf.rpkm lncRNAs##
#for (i in 1:nrow(rank1000_allgtf.rpkm)){
#  rank1000_allgtf.rpkm$gene_type[i]<-ifelse(rank1000_allgtf.rpkm$gene_name[i]%in%gtf_geneid_u$gene_id|rank1000_allgtf.rpkm$gene_name[i]%in%novel,"lncRNA","protein_coding gene")
#}

##join##
##no lncRNAS##only three novel lncRNA and known genes##
rank1000_allgtf.rpkm.type<-left_join(rank1000_allgtf.rpkm,gtf_geneid_u,by=c("gene_name"="gene_id"))


#########with gene_type###############
##top1000 gtf##
##read in rank1000 gtf file##
rand1000gtf<-read.delim("hg38_genelnRNA_sorted_top1000.gtf",header = FALSE,stringsAsFactors = FALSE)
rand1000gtf<-rand1000gtf%>%separate(V9,c(paste0("V",9:26)),";")

rand1000gtfS<-filter(rand1000gtf,V2=="StringTie")
rand1000gtfS<-rand1000gtfS%>%filter(V3=="transcript")
#rand1000gtfS<-rand1000gtfS[,1:10]
rand1000gtfS$gene_id<-str_replace_all(rand1000gtfS$V9, fixed("gene_id "), "")%>%
  str_replace_all(fixed(" "), "")

#rand1000gtfS<-rand1000gtfS[,c(1,4,5,7,9,11)]

rand1000gtfH<-filter(rand1000gtf,V2!="StringTie")
rand1000gtfH<-filter(rand1000gtfH,V3=="gene")
rand1000gtfH$gene_id<-str_replace_all(rand1000gtfH$V10, fixed("gene_id "), "")%>%
  str_replace_all(fixed(" "), "")
#rand1000gtfH<-rand1000gtfH[,-c(21:26)]
#rand1000gtfH<-rand1000gtfH[,c(1,4,5,7,11,21)]
#rand1000gtfH<-rand1000gtfH%>%rename("V9"="V11")

##combine##
rank1000_allgtf<-rbind(rand1000gtfS,rand1000gtfH)

rank1000_allgtf.rpkm<-full_join(rank1000_allgtf,BMI1all,by=c("gene_id"="Geneid"))

rank1000_allgtf.rpkm<-unique(rank1000_allgtf.rpkm)

##output for top 1000 genes and lncRNAs##
##1174 protein coding##

write.table(rank1000_allgtf.rpkm,"rank1000_genelncRNA.gtf.rpkm.type.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")

##lncRNA types##
lncRNAs = c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding",
  "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA",
  "bidirectional_promoter_lncrna")

##choose top 20 genes and 17novel and known lncRNAS######################
##processed_transcript(DANCR);
##lincRNA(NORAD,ILF3-DT,NCBP2-AS2,ENSG00000265519,ENSG00000254531) or (LINC00657, ILF3-AS1, NCBP2-AS2,CTD-3157E16.1,FLJ20021)
##antisense(CTD-2396E7.11,MID1IP1-AS1,GATA2-AS1,AC113189.5)
##G.(G.5257,G.5260,G.7790,G.5441,G.4473,G.7490,G.8242, G.1149)

##genes##
rank1000_allgtf.rpkm.top20gene<-rank1000_allgtf.rpkm[rank1000_allgtf.rpkm$rank1<=22&rank1000_allgtf.rpkm$rank2<=22,]
top20genename<-rank1000_allgtf.rpkm.top20gene$V11
top20genename<-top20genename%>%str_replace_all(fixed(" gene_name "), "")
write.table(rank1000_allgtf.rpkm.top20gene,"rank1000_allgtf.rpkm.top20gene.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")
write.table(top20genename,"top20genename.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")


