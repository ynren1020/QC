########################2019-04-30##################################
##hyb_dancr.join.total (su2c)
##boxplot for normal,tumor and su2c_hybrid##########################

library(dplyr)
library(tidyr)
library(ggpubr)

##su2c_hybrid##
janno_su2c<-read.delim("su2c_hyb_dancr.join.total.txt",header = TRUE,stringsAsFactors = FALSE)
janno_su2c$type<-"Metastatic (N=179)"
janno_su2c<-janno_su2c%>%rename("total"="V2")

##prad,prostate cancer##
prad<-read.delim("prad.metainfo.txt",header = TRUE,stringsAsFactors = FALSE)
janno_prad<-read.delim("janno.dancr.malat1.prad.txt",header = FALSE,stringsAsFactors = FALSE)

for (i in 1:nrow(janno_prad)){
  janno_prad$sample[i]<-strsplit(janno_prad$V1[i],"[.]")[[1]][1]
}

janno_prad.join<-full_join(janno_prad,prad,by=c("sample"="FILE_ID"))%>%
  rename("gene"="V2","score"="V3","total"="V4","scoreN"="V5","type"="TUMOR_TYPE")%>%
  select(sample,gene,score,total,scoreN,type)

##label##
for (i in 1:nrow(janno_prad.join)){
  janno_prad.join$type[i]<-ifelse(janno_prad.join$type[i]=="Primary Tumor","Primary Tumor (N=499)","Solid Tissue Normal (N=52)")
}

janno_prad.join.dancr<-filter(janno_prad.join,gene=="DANCR")
janno_prad.join.dancr$gene<-NULL
janno_prad.join.dancr$scoreNlog<-log10(janno_prad.join.dancr$scoreN)

##combine su2c and prad##
df<-rbind(janno_su2c,janno_prad.join.dancr)
##use log10scale to make the distribution more easier to see#
#df$scoreNlog<-log10(df$scoreN)
#primary tummor 499;solid tissue normal 52,su2c 59
##plot##
p <- ggboxplot(df, x = "type", y = "scoreNlog",
               color = "type", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = "type")+
               rotate_x_text(90)
p
# Change the plot orientation: horizontal
#ggpar(p, orientation = "horiz",xlab="Metastatic/n (N=59)")
# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list(c("Solid Tissue Normal (N=52)", "Primary Tumor (N=499)"), c("Primary Tumor (N=499)", "Metastatic (N=179)"), c("Solid Tissue Normal (N=52)", "Metastatic (N=179)"))
p1<-p + stat_compare_means(comparisons = my_comparisons,method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 4)                   # Add global p-value
p1
#ggpar(p1,yscale = "log2")
# Violin plots with box plots inside
# :::::::::::::::::::::::::::::::::::::::::::::::::::
# Change fill color by groups: dose
# add boxplot with white fill color
ggviolin(df, x = "type", y = "scoreN", fill = "type",
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"),xlab = FALSE,ylab="Normalized Junction Numbers")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 250)+                                      # Add global the p-value 
  rotate_x_text(90)
