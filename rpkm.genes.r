#####################2019-04-30###########################
##su2c.featurecount.txt,PRAD.txt,PRAD.normal.txt##########
##top20genes rpkm ########################################
##########################################################

##function to read in files##
read<-function(input,header){
  temp<-read.delim(input,header = header,stringsAsFactors = FALSE)
  return(temp)
}

rpkm_su2c<-read("su2c.featurecount.txt") #[1] 62757    59
rpkm_prad.t<-read("PRAD.txt") #60488   492
rpkm_prad.n<-read("PRAD.normal.txt") #60488    53
mapname<-read("mapping_names.dedup.corrected.txt",FALSE)

##strsplit##
for (i in 1:nrow(rpkm_prad.t)){
  rpkm_prad.t$Geneid[i]<-strsplit(rpkm_prad.t$GeneID[i],"[.]")[[1]][1]
}

for (i in 1:nrow(rpkm_prad.n)){
  rpkm_prad.n$Geneid[i]<-strsplit(rpkm_prad.n$GeneID[i],"[.]")[[1]][1]
}

for (i in 1:nrow(mapname)){
  mapname$Geneid[i]<-strsplit(mapname$V1[i],"[.]")[[1]][1]
}

##left_join and then calculated rpkm##
rpkm_su2c.join<-left_join(rpkm_su2c,mapname,by="Geneid")
rpkm_prad.n.join<-left_join(rpkm_prad.n,mapname,by="Geneid")
rpkm_prad.t.join<-left_join(rpkm_prad.t,mapname,by="Geneid")

#calrpkm<-function(x,y){x*1000000000/(sum(x)*y)} #x is the column of read counts, y is the length
#t1<-rpkm_su2c.join$X01115153.TA2.sorted.bam/sum(rpkm_su2c.join$X01115153.TA2.sorted.bam)/rpkm_su2c.join$V3*10^9
for (i in 2:59){
  rpkm_su2c.join[,i]<-rpkm_su2c.join[,i]/sum(rpkm_su2c.join[,i])/rpkm_su2c.join$V3*1000000000
}

for (i in 2:53){
  rpkm_prad.n.join[,i]<-rpkm_prad.n.join[,i]/sum(rpkm_prad.n.join[,i])/rpkm_prad.n.join$V3*1000000000
}

for (i in 2:492){
  rpkm_prad.t.join[,i]<-rpkm_prad.t.join[,i]/sum(rpkm_prad.t.join[,i])/rpkm_prad.t.join$V3*1000000000
}

##use > deparse(substitute(rpkm_prad.t.join))   ##dataframe's name changed to character
#[1] "rpkm_prad.t.join"
writein<-function(input,colname,rowname){
  write.table(input,paste0(deparse(substitute(input)),".txt"),quote=FALSE,col.names = colname,row.names = rowname,sep="\t")
}

writein(rpkm_su2c.join,TRUE,FALSE)
writein(rpkm_prad.n.join,TRUE,FALSE)
writein(rpkm_prad.t.join,TRUE,FALSE)

##subset dataframe##
genes<-read("genes20_allcol.txt")
#ENSG00000226950	DANCR
#ENSG00000260032	NORAD
lncrnas<-c("ENSG00000226950","ENSG00000260032")

subfunc<-function(input){
  temp<-input[input$Geneid%in%genes$gene_name|input$Geneid%in%lncrnas,]
  return(temp)
}
  
rpkm_su2c.sub<-subfunc(rpkm_su2c.join) #22 64
rpkm_prad.t.sub<-subfunc(rpkm_prad.t.join) ##22 498
rpkm_prad.n.sub<-subfunc(rpkm_prad.n.join) ##22 59

##wide to long##
su2c_wide<-rpkm_su2c.sub[,c(2:59,61)]
su2c_long<-as.data.frame(t(su2c_wide))
su2c_long<-su2c_long[c(nrow(su2c_long),1:nrow(su2c_long)-1),]
colnames(su2c_long)<-as.character(unlist(su2c_long[1,]))
su2c_long<-su2c_long[-1,]
su2c_long[] <- lapply(su2c_long, function(x) as.numeric(as.character(x)))
su2c_long$type<-"Metastatic (N=58)"


prad.n_wide<-rpkm_prad.n.sub[,c(2:53,56)]
prad.n_long<-as.data.frame(t(prad.n_wide))
prad.n_long<-prad.n_long[c(nrow(prad.n_long),1:nrow(prad.n_long)-1),]
colnames(prad.n_long)<-as.character(unlist(prad.n_long[1,]))
prad.n_long<-prad.n_long[-1,]
prad.n_long[] <- lapply(prad.n_long, function(x) as.numeric(as.character(x)))
prad.n_long$type<-"Normal (N=52)"


prad.t_wide<-rpkm_prad.t.sub[,c(2:492,495)]
prad.t_long<-as.data.frame(t(prad.t_wide))
prad.t_long<-prad.t_long[c(nrow(prad.t_long),1:nrow(prad.t_long)-1),]
colnames(prad.t_long)<-as.character(unlist(prad.t_long[1,]))
prad.t_long<-prad.t_long[-1,]
prad.t_long[] <- lapply(prad.t_long, function(x) as.numeric(as.character(x)))
prad.t_long$type<-"Tumor (N=491)"


##rbind all three for anova ##
df4test<-rbind(su2c_long,prad.n_long,prad.t_long)
df4test$type<-factor(df4test$type,levels=c("Metastatic (N=58)","Normal (N=52)","Tumor (N=491)"))

##parametric anova##

for (i in 1:22)
{
  formula <- paste(colnames(df4test)[i], " ~ type", sep="")
  
  p <- summary(aov(as.formula(formula), data=df4test))[[1]][["Pr(>F)"]][1]
  
  print(paste(formula, ": p=", p, sep=""))
}

#[1] "ENO1 ~ type: p=6.26820942746804e-36"
#[1] "EPRS ~ type: p=1.88629937784393e-24"
#[1] "HSPD1 ~ type: p=9.37473146317397e-29"
#[1] "EIF4G1 ~ type: p=8.7000377865538e-46"
#[1] "DANCR ~ type: p=8.7264883920902e-10"
#[1] "HSP90AB1 ~ type: p=8.2279326641772e-38"
#[1] "EEF1A1 ~ type: p=0.00122208462144051"
#[1] "ACTB ~ type: p=1.49874973756669e-17"
#[1] "GARS ~ type: p=3.42537344789294e-46"
#[1] "PABPC1 ~ type: p=1.65879278063285e-23"
#[1] "PRUNE2 ~ type: p=1.59073785316368e-07"
#[1] "HSPA8 ~ type: p=7.8464352928622e-20"
#[1] "ATP5B ~ type: p=2.79283137582914e-11"
#[1] "CKB ~ type: p=0.0210306882758719"
#[1] "AARS ~ type: p=8.59966783180506e-10"
#[1] "ACTG1 ~ type: p=1.38401717492706e-08"
#[1] "FASN ~ type: p=1.18611218096267e-31"
#[1] "LINC00657 ~ type: p=0.114243611511153"
#[1] "EEF1A2 ~ type: p=3.52406511942726e-14"
#[1] "EEF2 ~ type: p=3.24551316398549e-11"
#[1] "TRIM28 ~ type: p=3.36022102319783e-42"
#[1] "RPL3 ~ type: p=4.77134632644341e-17"

##use this one output##
for (i in 1:22)
{
  formula <- paste(colnames(df4test)[i], " ~ type", sep="")
  
  p <- summary(res.aov<-aov(as.formula(formula), data=df4test))[[1]][["Pr(>F)"]][1]
  
  p.pair<-TukeyHSD(res.aov)
  print(unlist(strsplit(paste(paste(colnames(df4test)[i],p, sep=":"),paste(p.pair$type[,4],collapse=":"),sep=":",collapse =NULL),":")))
  #print(p.pair$type[,4])
}

#p.avova Normal (N=52)-Metastatic (N=58) Tumor (N=491)-Metastatic (N=58)     Tumor (N=491)-Normal (N=52) 


##average expression per gene per group##
df4test.mean<-df4test%>%group_by(type)%>%
  summarize_all("mean")

df4test.mean.t<-as.data.frame(t(df4test.mean))


df4test.sd<-df4test%>%group_by(type)%>%
  summarize_all("sd")

#type   ENO1  EPRS HSPD1 EIF4G1 DANCR HSP90AB1 EEF1A1  ACTB  GARS PABPC1 PRUNE2 HSPA8 ATP5B   CKB  AARS ACTG1  FASN LINC00657 EEF1A2  EEF2
#<fct> <dbl> <dbl> <dbl>  <dbl> <dbl>    <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl>  <dbl> <dbl>
#  1 Meta… 168.  12.5   71.6  45.1  40.2     222.    541.  405.  9.30  285.    9.76 146.  120.  164.  25.6   440. 526.       34.0  151.   634.
#2 Norm…  40.7  3.17  14.3   8.11  8.29     69.0   194.  411.  2.02   61.1   7.49  46.1  28.7  92.8  4.98  148.  40.4      11.3   12.3  157.
#3 Tumo…  43.5  5.93  24.1  14.0  13.2      93.0   216.  222.  4.76  133.    5.39  64.2  50.2 117.   8.85  182.  99.9      12.8   69.2  265.
