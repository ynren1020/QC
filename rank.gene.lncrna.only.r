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

##lncRNAs and protein coding genes#############################################
##check the Geneid in BMI1 is lncRNA or not##
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
gtf_geneid_u.trim<-str_trim(gtf_geneid_u$gtf_geneid_u) #remove whitespace

##check the Geneid in BMI1 is protein coding gene or not##
####gtf file find protein coding gene,produced from proteincoding.py with genecode.v21.annotation.gtf##
coding<-read.delim("protein_coding_YR_hg38.gtf",header = FALSE,stringsAsFactors = FALSE)
coding<-coding%>%separate(V9,c(paste0("V",9:26)),";")
codingvector<-unlist(coding)
coding_geneid<-str_subset(codingvector,"ENSG")
##protein coding genes' geneid##
coding_geneid_u<-str_replace_all(coding_geneid,"gene_id ","")
coding_geneid_u<-as.data.frame(coding_geneid_u)
coding_geneid_u$coding_geneid_u<-as.character(coding_geneid_u$coding_geneid_u)
coding_geneid_u<-unique(coding_geneid_u)
coding_geneid_u.trim<-str_trim(coding_geneid_u$coding_geneid_u)
###subset BMI1 as ENSG or G.##
BMI1known<-BMI1 %>%
  filter(str_detect(Geneid, "ENSG")) ##choose ENSG geneid
BMI1knownsub<-BMI1known[BMI1known$Geneid%in%coding_geneid_u.trim|BMI1known$Geneid%in%gtf_geneid_u.trim,] ##choose genes or lncRNAs

##novel lncRNA##
BMI1lncRNA<-BMI1 %>%
  filter(str_detect(Geneid, "^G."))
  
##combine known gene and lncrna with novel lncrna##
BMI1.gene.lncrna<-rbind(BMI1knownsub,BMI1lncRNA)

##rank within genes and lncRNAs##
BMI1.gene.lncrna$rank1<-rank(-BMI1.gene.lncrna$rpkm1)
BMI1.gene.lncrna$rank2<-rank(-BMI1.gene.lncrna$rpkm2)

##choose top 1000s##
BMI1top1000<-filter(BMI1.gene.lncrna,rank1<=1000&rank2<=1000)

#####ENSG to gene symbol###
BMI1top1000sym<-BMI1top1000 %>%
  filter(str_detect(Geneid, "ENSG"))
BMI1top1000sym<-BMI1top1000sym%>%
  mutate(gene_name=str_split_fixed(BMI1top1000sym$Geneid,"[.]",n=2)[,1])

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
genesym<-getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=BMI1top1000sym$gene_name, mart=ensembl)
#BMI1name<-BMI1known[BMI1known$gene_name%in%genesym$ensembl_gene_id,]
BMI1namefinal<-full_join(BMI1top1000sym,genesym,by=c("gene_name"="ensembl_gene_id"))
##known gene and lncRNAs##
BMI1namefinalsub<-BMI1namefinal%>%rename("gene"="hgnc_symbol") 

##add genetype annotation##
for (i in 1:nrow(BMI1namefinalsub)){
  BMI1namefinalsub$gene_type[i]<-ifelse(BMI1namefinalsub$Geneid[i]%in%gtf_geneid_u.trim,"lncRNAs","protein_coding genes")
}

BMI1top1000novel<-BMI1top1000 %>%
  filter(str_detect(Geneid, "^G"))%>%
  mutate(gene_name=Geneid)%>%
  mutate(gene=Geneid)%>%
  mutate(gene_type="novel lncRNAs")

##combine known genes/lncRNAs and novel lncRNA##
BMI1all<-rbind(BMI1namefinalsub,BMI1top1000novel)
##only get first 720 in rank1 and rank2##FINALLY USE THIS ONE FOR FURTHUR ANALYSIS##
BMI1all<-filter(BMI1all,rank1<=720&rank2<=720) #7 lncRNAs(2 known),541protein coding genes
write.table(BMI1all,"BMI1_genelncRNA_updated_720.txt",col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
##
geneid_top1000<-as.data.frame(BMI1all$Geneid)
write.table(geneid_top1000,"geneid_top1000.hg38.txt",col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")

##top20genes##
genes20<-filter(BMI1all,rank1<=23&gene_type=="protein_coding genes")
write.table(genes20,"genes20_allcol.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")
genes20fortest<-as.data.frame(genes20$gene)
write.table(genes20fortest,"genes20.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep="\t")






