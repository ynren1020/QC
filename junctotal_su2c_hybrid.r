##################################2019-04-26#############################################
##top expressed known lncRNA DANCR and MALAT1 in BMI1 data###############################
##summarize their first junction:exon1end-exon2start numbers,normalize by total juncs####
##TCGA(PRAD:normal and tumor),SU2C metastic cancer#######################################
##boxplot,data:/home/tywang/Projects/Exitron/SU2C/Hybrid_Selection#######################
##calculate total junction for each janno file###########################################
#########################################################################################

#library(tidyverse)
#library(ggpubr)
#library(dplyr)
#library(tidyr)

args <- commandArgs(TRUE)
hybrid_janno<-read.delim(args[1],header = TRUE,stringsAsFactors = TRUE)
total<-sum(hybrid_janno$score)
hybrid_total<-data.frame(sample=args[1],total=total)
write.table(hybrid_total,paste0(args[1],".total.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE,sep="\t")