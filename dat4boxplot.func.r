############################2019-09-11##################
##data of expression for boxplot #######################
########################################################

dat4box<-function(input,genes,col_1,col_2,groups){
  ##read data##
  tmp<-read.delim(input,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
  ##choose target gene##
  tmp.sub<-dplyr::filter(tmp,V2==genes)
  #tmp.sub<-tmp[tmp$col==genes,col_1:(ncol(tmp)-col_2)]
  tmp.sub<-na.omit(tmp.sub)
  ##transpose##
  rownames(tmp.sub)<-tmp.sub$V2
  tmp.sub<-tmp.sub[,col_1:(ncol(tmp.sub)-col_2)]
  tmp.subT<-as.data.frame(t(tmp.sub))
  ##create group factor variable##
  tmp.subT$group<-groups
  return(tmp.subT)
}

normal<-dat4box("rpkm_prad.n.join.txt","TUG1",2,6,"Normal")
tumor<-dat4box("rpkm_prad.t.join.txt","TUG1",2,6,"Localized")
mets<-dat4box("rpkm_su2c.join.txt","TUG1",2,5,"Mets")



