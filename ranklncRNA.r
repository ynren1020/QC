#########2019-04-24##########################################
#########################QC BMI1 data########################
##rank lncRNAs (known or novel) by read counts###############
#############################################################

library("biomaRt")
library(dplyr)
library(tidyr)
library(stringr)

BMI1<-read.delim("BMI1_id_counts.txt",stringsAsFactors = FALSE,skip = 1,header = TRUE)
BMI1known<-BMI1 %>%
  filter(str_detect(Geneid, "ENSG"))
BMI1known<-BMI1known%>%
  mutate(gene_name=str_split_fixed(BMI1known$Geneid,"[.]",n=2)[,1])


ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
genesym<-getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=BMI1known$gene_name, mart=ensembl)
BMI1name<-BMI1known[BMI1known$gene_name%in%genesym$ensembl_gene_id,]
  

BMI1namefinal<-full_join(BMI1name,genesym,by=c("gene_name"="ensembl_gene_id"))
BMI1namefinalsub<-BMI1namefinal[,c(10,7,8,6,1)]
BMI1namefinalsub<-BMI1namefinalsub%>%rename("gene"="hgnc_symbol")

##output geneid to find which one is lncRNA##
geneid_hg38<-as.data.frame(BMI1namefinalsub$Geneid)
write.table(geneid_hg38,"geneid_hg38.txt",col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")

##novel lncRNA##
BMI1lncRNA<-BMI1 %>%
  filter(str_detect(Geneid, "^G."))%>%
  mutate(gene=Geneid)%>%
  select(gene,BMI1_1_sorted.bam:BMI1_2_sorted.bam,Length,Geneid)



####gtf file find lncRNAs,produced from lncRNAs.py with genecode.v21.annotation.gtf##
gtf<-read.delim("lncRNAs.gtf",header = FALSE,stringsAsFactors = FALSE)
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

##find matched lncRNAs with BMI1namefinalsub##
BMI1namefinalsub$Geneid<-as.character(BMI1namefinalsub$Geneid)
for (i in 1:nrow(BMI1namefinalsub)){
  BMI1namefinalsub$gene_id[i]<-strsplit(BMI1namefinalsub$Geneid[i],'[.]')[[1]][1]
}

#BMI1namefinalsub.matched<-BMI1namefinalsub[BMI1namefinalsub$Geneid%in%gtf_geneid_u,]
BMI1namefinalsub.matched<-full_join(gtf_geneid_u,BMI1namefinalsub,by="gene_id")
BMI1namefinalsub.matched<-na.omit(BMI1namefinalsub.matched)

BMI1namefinalsub.matched<-BMI1namefinalsub.matched[,3:7]

##combine##
BMI1final<-rbind(BMI1namefinalsub.matched,BMI1lncRNA)%>%
  arrange(desc(BMI1_1_sorted.bam),desc(BMI1_2_sorted.bam))
##RPKM##
BMI1final<-BMI1final%>%
  mutate(BMI1_1_rpkm=BMI1_1_sorted.bam/Length*100,BMI1_2_rpkm=BMI1_2_sorted.bam/Length*100)%>%
  arrange(desc(BMI1_1_rpkm),desc(BMI1_2_rpkm))

##rank##
BMI1final$rank1<-rank(-BMI1final$BMI1_1_rpkm)
BMI1final$rank2<-rank(-BMI1final$BMI1_2_rpkm)
##output##
write.table(BMI1final,"BMI1final_selected_all.txt",col.names = TRUE,row.names = FALSE,quote = FALSE,sep="\t")
##only geneid##
all_geneid<-BMI1final$Geneid
all_geneid<-as.data.frame(all_geneid)
write.table(all_geneid,"all_geneid.txt",col.names = FALSE,row.names = FALSE,quote = FALSE,sep="\t")

##rank <=50##
BMI1finalrank50<-filter(BMI1final,rank1<=50|rank2<=50)
write.table(BMI1finalrank50,"BMI1final_selected_rank50.txt",col.names = TRUE,row.names = FALSE,quote = FALSE,sep="\t")
##geneid of first 50s##
rank50_geneid<-as.data.frame(BMI1finalrank50$Geneid)
write.table(rank50_geneid,"rank50_geneid.txt",col.names = FALSE,row.names = FALSE,quote = FALSE,sep="\t")

##read in rank50 gtf file##
rand50gtf<-read.delim("rank50_geneid_lncRNA.hg38.gtf",header = FALSE,stringsAsFactors = FALSE)
rand50gtf<-rand50gtf%>%separate(V9,c(paste0("V",9:26)),";")

rand50gtfS<-filter(rand50gtf,V2=="StringTie")
rand50gtfS<-filter(rand50gtfS,V3=="transcript")%>%select(V1:V10)
rand50gtfS$gene_id<-str_replace_all(rand50gtfS$V9, fixed("gene_id "), "")%>%
  str_replace_all(fixed(" "), "")
rand50gtfS<-rand50gtfS[,c(1,4,5,7,9,11)]

rand50gtfH<-filter(rand50gtf,V2!="StringTie")
rand50gtfH<-filter(rand50gtfH,V3=="gene")
rand50gtfH$gene_id<-str_replace_all(rand50gtfH$V10, fixed("gene_id "), "")%>%
  str_replace_all(fixed(" "), "")
rand50gtfH<-rand50gtfH%>%select(-c(V21:V26))
rand50gtfH<-rand50gtfH[,c(1,4,5,7,11,21)]%>%rename("V9"="V11")

##combine##
rank50_allgtf<-rbind(rand50gtfS,rand50gtfH)

rank50_allgtf.rpkm<-full_join(rank50_allgtf,BMI1finalrank50,by=c("gene_id"="Geneid"))

rank50_allgtf.rpkm<-unique(rank50_allgtf.rpkm)

##output##
write.table(rank50_allgtf.rpkm,"rank50_allgtf.rpkm.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")





                    