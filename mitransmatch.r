####################2019-05-03#################
##lncrna.mitrans.match.annotated.eck.gtf#######
##all lncRNAs match to mitrans (hg19)##########
###############################################
library(dplyr)
library(tidyr)
library(stringr)

##choose matched lncrna(with mitranscriptome)##
milncrna<-read.delim("lncrna.mitrans.match.annotated.eck.gtf",header = FALSE,stringsAsFactors = FALSE)
lncrna.rank<-read.delim("BMI1final_selected_all.txt",header=TRUE,stringsAsFactors = FALSE)

milncrna<-separate(milncrna,"V9",paste0("V",9:16),";")
milncrna$V10<-str_replace_all(milncrna$V10," gene_id ","")

mi.rank<-full_join(lncrna.rank,milncrna,by=c("Geneid"="V10"))
mi.rank<-na.omit(mi.rank)

mi.rank$rank1<-as.integer(mi.rank$rank1)
mi.rank$rank2<-as.integer(mi.rank$rank2)

match.lncrna<-unique(mi.rank$Geneid) #5809
for (i in 1:length(match.lncrna)){
  match.lncrna[i]<-strsplit(match.lncrna[i],'[.]')[[1]][1]
}

##calculate those lncrnas rpkm in su2c and prad##
rpkm_su2c.join<-read.delim("rpkm_su2c.join.txt",header = TRUE,stringsAsFactors = FALSE)
rpkm_prad.n.join<-read.delim("rpkm_prad.n.join.txt",header = TRUE,stringsAsFactors = FALSE)
rpkm_prad.t.join<-read.delim("rpkm_prad.t.join.txt",header = TRUE,stringsAsFactors = FALSE)

subfunc<-function(input){
  temp<-input[input$Geneid%in%match.lncrna,]
  return(temp)
}

rpkm_su2c.sub<-subfunc(rpkm_su2c.join) #5708 64
rpkm_prad.t.sub<-subfunc(rpkm_prad.t.join) ##5806 498
rpkm_prad.n.sub<-subfunc(rpkm_prad.n.join) ##5806 59

rpkm_su2c.sub$type<-"su2c"
rpkm_prad.t.sub$type<-"tumor"
rpkm_prad.n.sub$type<-"normal"


two.join<-full_join(rpkm_prad.n.sub,rpkm_prad.t.sub,by="Geneid")
three.join<-full_join(two.join,rpkm_su2c.sub,by="Geneid") 
three.join<-na.omit(three.join) #5705 622

##separeate again##
rpkm_prad.n.sub.s<-three.join[,1:60]
rpkm_prad.t.sub.s<-three.join[,61:558]
rpkm_su2c.sub.s<-three.join[,559:622]



##wide to long##
su2c_wide<-rpkm_su2c.sub.s[,1:59]
su2c_long<-as.data.frame(t(su2c_wide))
su2c_long<-su2c_long[c(nrow(su2c_long),1:nrow(su2c_long)-1),]
colnames(su2c_long)<-as.character(unlist(su2c_long[1,]))
su2c_long<-su2c_long[-1,]
su2c_long[] <- lapply(su2c_long, function(x) as.numeric(as.character(x)))
su2c_long$type<-"Metastatic (N=58)"


prad.n_wide<-rpkm_prad.n.sub.s[,c(2:53,55)]
prad.n_long<-as.data.frame(t(prad.n_wide))
prad.n_long<-prad.n_long[c(nrow(prad.n_long),1:nrow(prad.n_long)-1),]
colnames(prad.n_long)<-as.character(unlist(prad.n_long[1,]))
prad.n_long<-prad.n_long[-1,]
prad.n_long[] <- lapply(prad.n_long, function(x) as.numeric(as.character(x)))
prad.n_long$type<-"Normal (N=52)"


prad.t_wide<-rpkm_prad.t.sub.s[,c(2:493)]
prad.t_long<-as.data.frame(t(prad.t_wide))
prad.t_long<-prad.t_long[c(nrow(prad.t_long),1:nrow(prad.t_long)-1),]
colnames(prad.t_long)<-as.character(unlist(prad.t_long[1,]))
prad.t_long<-prad.t_long[-1,]
prad.t_long[] <- lapply(prad.t_long, function(x) as.numeric(as.character(x)))
prad.t_long$type<-"Tumor (N=491)"

##rbind all three for anova ##
df4test<-rbind(su2c_long,prad.n_long,prad.t_long)
df4test$type<-factor(df4test$type,levels=c("Metastatic (N=58)","Normal (N=52)","Tumor (N=491)")) #601 5706


##average expression per gene per group##
df4test.mean<-df4test%>%group_by(type)%>%
  summarize_all("mean")

df4test.mean.t<-as.data.frame(t(df4test.mean))
df4test.mean.t<-rename(df4test.mean.t,"Metastatic (N=58)"="V1","Normal (N=52)"="V2","Tumor (N=491)"="V3")
df4test.mean.t<-df4test.mean.t[-1,]
df4test.mean.t$Geneid<-rownames(df4test.mean.t)

##join with lncrna.rank##
df4test.mean.t.join<-full_join(df4test.mean.t,lncrna.rank,by="Geneid")

df4test.mean.t.join<-na.omit(df4test.mean.t.join)
df4test.mean.t.join$rank1<-as.integer(df4test.mean.t.join$rank1)
df4test.mean.t.join$rank2<-as.integer(df4test.mean.t.join$rank2)
df4test.mean.t.join<-df4test.mean.t.join[,c(4:5,1:3,6:12)]

write.table(df4test.mean.t.join,"mitrans.lncrna.rank.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")


df4test.sd<-df4test%>%group_by(type)%>%
  summarize_all("sd")

