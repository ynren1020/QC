#########################2019-09-11########################
##TUG1 and SNHG9 expression in prad and su2c ##############
##survival of TUG1, ref to SNHG9 (tumor expression#########
###########################################################
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(stringr)

input1<-"rpkm_prad.t.join.txt"
input2<-"TCGA-CDR.info.txt"

df.rpkm<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
df.surv<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

#lncrnas<-c("ENSG00000253352","ENSG00000226950","ENSG00000267100","ENSG00000275234","ENSG00000244300",
#           "ENSG00000254531","ENSG00000261770","ENSG00000259330","ENSG00000253352","ENSG00000230733",
#           "ENSG00000235280","ENSG00000261373","ENSG00000273344","ENSG00000255198")
#lncrnas<-c("","","","","","","","")
#lncrnas<-as.data.frame(lncrnas)
#lncrnas$lncrnas2<-lncrnas$lncrnas<-as.character(lncrnas$lncrnas)
df.rpkm.sub<-df.rpkm[df.rpkm$Geneid=="ENSG00000253352",]

df.rpkm.sub<-df.rpkm.sub[,c(2:492,495)]
##transpose##
rownames(df.rpkm.sub)<-df.rpkm.sub$V2
df.rpkm.sub$V2<-NULL
df.rpkm.sub.long<-as.data.frame(t(df.rpkm.sub))

##quantiles of all columns##
quants <- c(0.05,0.10,0.20,0.25,0.4,0.50,0.75)
df.rpkm.sub.long.quantile<-apply( df.rpkm.sub.long[1] , 2 , quantile , probs = quants , na.rm = TRUE )
df.rpkm.sub.long.quantile<-as.data.frame(df.rpkm.sub.long.quantile)
#        TUG1
#5%  10.06895
#10% 11.48459
#20% 13.72061
#25% 14.35893
#40% 16.27833
#50% 17.27569
#75% 20.62186

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
#[8] "DSS_time"      "DFI"           "DFI_time"      "PFI"           "PFI_time"      "AC092171.4"    "MCF2L-AS1"    
#[15] "GATA2-AS1"     "FLJ20021"      "SNHG9"         "INAFM2"        "LINC00657"     "VPS9D1-AS1"    "CTC-459F4.1"  
#[22] "ILF3-AS1"      "PAXIP1-AS1"    "CTD-2396E7.11" "sample"

##survival analysis##
#lncrnas_names<-c("AC092171.4","MCF2L-AS1","GATA2-AS1","FLJ20021","SNHG9","INAFM2", "LINC00657" ,"VPS9D1-AS1","CTC-459F4.1","ILF3-AS1","PAXIP1-AS1","CTD-2396E7.11")
##disease specific##


fit1<-survfit(Surv(DSS_time,DSS)~TUG1,data=prad.surv.status)
#summary(fit1)
res<-ggsurvplot(fit1,data=prad.surv.status,
                xlab = "Days",
                ylab = "Disease Specific Survival Probability (%)",
                conf.int=TRUE,
                pval = TRUE,
                fun="pct",
                risk.table = TRUE,
                size=1,
                linetype = "strata"
)
ggsave("TUG1.DS.pdf", plot = print(res), width = 8, height = 8, dpi = 500)



##overall survival##
fit2<-survfit(Surv(OS_time,OS)~TUG1,data=prad.surv.status)
summary(fit2)
ggsurvplot(fit2,data=prad.surv.status,
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata")

ggsave("TUG1.OS.pdf", plot = print(res), width = 8, height = 8, dpi = 500)

##disease free##
fit3<-survfit(Surv(DFI_time,DFI)~TUG1,data=prad.surv.status)
summary(fit3)
ggsurvplot(fit3,data=prad.surv.status,
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata")
ggsave("TUG1.DF.pdf", plot = print(res), width = 8, height = 8, dpi = 500)
#AC092171.4,MCF2L-AS1,GATA2-AS1,
##progress free##
fit4<-survfit(Surv(PFI_time,PFI)~TUG1,data=prad.surv.status)
summary(fit4)
ggsurvplot(fit4,data=prad.surv.status,
           conf.int=TRUE,
           pval = TRUE,
           fun="pct",
           risk.table = TRUE,
           size=1,
           linetype = "strata")
ggsave("TUG1.PF.pdf", plot = print(res), width = 8, height = 8, dpi = 500)





