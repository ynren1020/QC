###########################2019-05-08##################################
##lncRNAs expression in prad tumor samples correlation with survival###
##rpkm_prad.t.join.txt; TCGA-CDR.info.txt##############################
#######################################################################

library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(stringr)

input1<-"rpkm_prad.t.join.txt"
input2<-"TCGA-CDR.info.txt"

df.rpkm<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
df.surv<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

lncrnas<-c("ENSG00000260032","ENSG00000226950","ENSG00000267100","ENSG00000275234","ENSG00000244300",
           "ENSG00000254531","ENSG00000261770","ENSG00000259330","ENSG00000253352","ENSG00000230733",
           "ENSG00000235280","ENSG00000261373","ENSG00000273344","ENSG00000255198")
#lncrnas<-c("","","","","","","","")
lncrnas<-as.data.frame(lncrnas)
lncrnas$lncrnas2<-lncrnas$lncrnas<-as.character(lncrnas$lncrnas)
df.rpkm.sub<-full_join(df.rpkm,lncrnas,by=c("Geneid"="lncrnas"))
df.rpkm.sub<-df.rpkm.sub[!is.na(df.rpkm.sub$lncrnas2),] 
df.rpkm.sub<-na.omit(df.rpkm.sub)
df.rpkm.sub<-df.rpkm.sub[,c(2:492,495)]
##transpose##
df.rpkm.sub.long<-as.data.frame(t(df.rpkm.sub))
##first row as column names##
colnames(df.rpkm.sub.long) <- as.character(unlist(df.rpkm.sub.long[492,]))
df.rpkm.sub.long <- df.rpkm.sub.long[-492, ]
##all factor columns to numeric##
df.rpkm.sub.long[] <- lapply(df.rpkm.sub.long, function(x) as.numeric(as.character(x)))
#df.rpkm.sub.long$sample<-rownames(df.rpkm.sub.long)

##quantiles of all columns##
quants <- c(0.05,0.10,0.20,0.25,0.4,0.50,0.75)
df.rpkm.sub.long.quantile<-apply( df.rpkm.sub.long[1:12] , 2 , quantile , probs = quants , na.rm = TRUE )
df.rpkm.sub.long.quantile<-as.data.frame(df.rpkm.sub.long.quantile)
#ENSG00000230733.2 ENSG00000235280.2 ENSG00000244300.2 ENSG00000254531.1 ENSG00000255198.4 ENSG00000259330.1 ENSG00000260032.1
#5%          0.8927452          1.867047          2.207025          13.29672          3.085909          8.552635          29.13948
#10%         1.1205594          2.228413          2.614374          19.99545          3.737030          9.858740          33.18951
#20%         1.4648650          2.636121          3.212582          28.13676          4.931422         11.896790          38.59204
#25%         1.6520230          2.890785          3.484603          30.63338          5.370368         12.541152          40.36347
#40%         2.0876681          3.485143          4.335797          41.24361          7.218073         14.073285          44.89203
#50%         2.3491391          3.905147          4.788013          48.15356          9.159649         14.974625          47.50317
#75%         3.2679360          5.262363          6.027805          70.34682         15.192654         18.092180          56.08138
#ENSG00000261373.1 ENSG00000261770.1 ENSG00000267100.1 ENSG00000273344.1 ENSG00000275234.1
#5%           1.600179        0.05609604          6.747811          2.808585          10.73961
#10%          1.972882        0.08237338          7.702494          3.279569          12.61466
#20%          2.651968        0.12134930          9.084642          3.726964          14.51006
#25%          2.955809        0.13519245          9.670472          3.885310          15.25765
#40%          3.914801        0.19015030         11.453198          4.362691          17.11611
#50%          4.512871        0.23223920         12.197686          4.732943          18.40405
#75%          6.594325        0.42769795         15.410895          5.950982          22.69855

##based on median separate patient into two groups##
#for (i in 1:nrow(df.rpkm.sub.long)){
#  df.rpkm.sub.long$ENSG00000230733.2_status[i]<-ifelse(df.rpkm.sub.long$ENSG00000230733.2[i]>=median(df.rpkm.sub.long$ENSG00000230733.2),"high","low")
#  df.rpkm.sub.long$ENSG00000235280.2_status[i]<-ifelse(df.rpkm.sub.long$ENSG00000235280.2[i]>=median(df.rpkm.sub.long$ENSG00000235280.2),"high","low")
#}
#function to assgin status to genes based on median expression##
status<-function(x){
  y<-ifelse(x>=median(x),"high","low")
}
#df.rpkm.sub.long$ENSG00000230733.2_status<-status(df.rpkm.sub.long$ENSG00000230733.2)
df.rpkm.sub.long.status<-apply(df.rpkm.sub.long,2,status)
df.rpkm.sub.long.status<-as.data.frame(df.rpkm.sub.long.status)
df.rpkm.sub.long.status$sample<-rownames(df.rpkm.sub.long.status)
df.rpkm.sub.long.status$sample_id<-str_replace_all(df.rpkm.sub.long.status$sample,'[.]',"-")

##prad survival##
df.surv.prad<-df.surv[df.surv$cohort=="PRAD",]
prad.surv.status<-full_join(df.surv.prad,df.rpkm.sub.long.status,by=c("barcode"="sample_id"))
prad.surv.status<-prad.surv.status[!is.na(prad.surv.status$sample),]

#[1] "patient"           "barcode"           "file_id"           "cohort"            "OS"                "OS_time"          
#[7] "DSS"               "DSS_time"          "DFI"               "DFI_time"          "PFI"               "PFI_time"         
#[13] "ENSG00000230733.2" "ENSG00000235280.2" "ENSG00000244300.2" "ENSG00000254531.1" "ENSG00000255198.4" "ENSG00000259330.1"
#[19] "ENSG00000260032.1" "ENSG00000261373.1" "ENSG00000261770.1" "ENSG00000267100.1" "ENSG00000273344.1" "ENSG00000275234.1"
#[25] "sample"  
#[1] "patient"       "barcode"       "file_id"       "cohort"        "OS"            "OS_time"       "DSS"          
#[8] "DSS_time"      "DFI"           "DFI_time"      "PFI"           "PFI_time"      "DANCR"         "AC092171.4"   
#[15] "MCF2L-AS1"     "GATA2-AS1"     "TUG1"          "FLJ20021"      "SNHG9"         "INAFM2"        "LINC00657"    
#[22] "VPS9D1-AS1"    "CTC-459F4.1"   "ILF3-AS1"      "PAXIP1-AS1"    "CTD-2396E7.11" "sample" 

##survival analysis##
lncrnas_names<-c("DANCR","AC092171.4","MCF2L-AS1","GATA2-AS1","TUG1","FLJ20021","SNHG9","INAFM2", "LINC00657" ,"VPS9D1-AS1","CTC-459F4.1","ILF3-AS1","PAXIP1-AS1","CTD-2396E7.11")
##disease specific##
for (i in 13:26){
  status<-prad.surv.status[,i]

fit1<-survfit(Surv(DSS_time,DSS)~status,data=prad.surv.status)
res1<-ggsurvplot(fit1,data=prad.surv.status,
           xlab = "Days",
           ylab = "Disease Specific Survival Probability (%)",
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata",
           legend.labs = c(paste0(names(prad.surv.status[i]),"-high"),paste0(names(prad.surv.status[i]),"-low"))
           )
ggsave(paste0(names(prad.surv.status[i]),".DS.png"), plot = print(res1), width = 8, height = 8, dpi = 500)

##overall survival##
fit2<-survfit(Surv(OS_time,OS)~status,data=prad.surv.status)
res2<-ggsurvplot(fit2,data=prad.surv.status,
                 xlab = "Days",
                 ylab = "Overall Survival Probability (%)",
                 conf.int=TRUE,
                 pval = TRUE,
                 fun="pct",
                 risk.table = TRUE,
                 size=1,
                 linetype = "strata",
                 legend.labs = c(paste0(names(prad.surv.status[i]),"-high"),paste0(names(prad.surv.status[i]),"-low")))
ggsave(paste0(names(prad.surv.status[i]),".OS.png"), plot = print(res2), width = 8, height = 8, dpi = 500)

##disease free##
fit3<-survfit(Surv(DFI_time,DFI)~status,data=prad.surv.status)
res3<-ggsurvplot(fit3,data=prad.surv.status,
           xlab = "Days",
           ylab = "Disease Free Survival Probability (%)",
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata",
           legend.labs = c(paste0(names(prad.surv.status[i]),"-high"),paste0(names(prad.surv.status[i]),"-low")))

ggsave(paste0(names(prad.surv.status[i]),".DF.png"), plot = print(res3), width = 8, height = 8, dpi = 500)
#AC092171.4,MCF2L-AS1,GATA2-AS1,
##progress free##
fit4<-survfit(Surv(PFI_time,PFI)~status,data=prad.surv.status)
res4<-ggsurvplot(fit4,data=prad.surv.status,
          xlab = "Days",
          ylab = "Progression Free Survival Probability (%)",
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata",
           legend.labs = c(paste0(names(prad.surv.status[i]),"-high"),paste0(names(prad.surv.status[i]),"-low")))

ggsave(paste0(names(prad.surv.status[i]),".PF.png"), plot = print(res4), width = 8, height = 8, dpi = 500)
}
