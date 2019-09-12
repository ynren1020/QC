###############################2019-09-11####################
##boxplot of TUG1 and SNHG9 expression prad and su2c#########
#################rpkm_prad.n.join.txt
#################"rpkm_prad.t.join.txt"
#################rpkm_su2c.join.txt
#############################################################

library(tidyverse)
library(ggpubr)

input1<-"rpkm_prad.n.join.txt"
input2<-"rpkm_prad.t.join.txt"
input3<-"rpkm_su2c.join.txt"
prad.n<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
prad.t<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
su2c<-read.delim(input3,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

##TUG1
prad.n.sub<-prad.n[prad.n$V2=="TUG1",]
prad.t.sub<-prad.t[prad.t$V2=="TUG1",]
su2c.sub<-su2c[su2c$V2=="TUG1",]

prad.n.sub<-na.omit(prad.n.sub)
prad.t.sub<-na.omit(prad.t.sub)
su2c.sub<-na.omit(su2c.sub)

##transpose##
rownames(prad.n.sub)<-prad.n.sub$V2
rownames(prad.t.sub)<-prad.t.sub$V2
rownames(su2c.sub)<-su2c.sub$V2
prad.n.sub<-prad.n.sub[,2:(ncol(prad.n.sub)-6)]
prad.t.sub<-prad.t.sub[,2:(ncol(prad.t.sub)-6)]
su2c.sub<-su2c.sub[,2:(ncol(su2c.sub)-5)]

prad.n.subT<-as.data.frame(t(prad.n.sub))
prad.t.subT<-as.data.frame(t(prad.t.sub))
su2c.subT<-as.data.frame(t(su2c.sub))


prad.n.subT$group<-"Normal"
prad.t.subT$group<-"Localized"
su2c.subT$group<-"Metastasis"

##rbind##
all.join<-rbind(prad.n.subT,prad.t.subT,su2c.subT)
all.join$group<-factor(all.join$group,levels = c("Normal","Localized","Metastasis"))
all.join$logTUG1<-log10(all.join$TUG1)
##plot##
p <- ggboxplot(all.join, x = "group", y = "logTUG1",xlab = FALSE,ylab="TUG1 expression (log10FPKM)",
               color = "group", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter")

p
# Change the plot orientation: horizontal
#ggpar(p, orientation = "horiz",xlab="Metastatic/n (N=59)")
# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list( c("Normal", "Localized"), c("Localized", "Metastasis"), c("Metastasis", "Normal") )
p1<-p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 2.5)+rotate_x_text(90)+                  # Add global p-value
  scale_x_discrete(labels=c("Normal" = "Normal\n (N=52)", "Localized" = "Localized\n (N=491)",
                            "Metastasis" = "Metastasis\n (N=58)"))+theme(legend.position = "none")
p1

ggsave("TUG1_TCGASU2C.boxplot.log10.pdf",width = 8,height = 8)




