# March 30, 2022 
# Input xpore DRACH results of CvsE and CvsM comparison 
# RNAmod ppt results are from this script ---


ce <- data.table::fread("./ont/xpore/CvsE_woSM_majority_direction_kmer_diffmod.table_DRACH.txt")
cm <- data.table::fread("./ont/xpore/CvsM_woSM_majority_direction_kmer_diffmod.table_DRACH.txt")
dim(ce) # 170780 125350
dim(cm) # 157624 115265

# create a index for finding unique gene:position:kmer 
ce$index <- paste0(ce$id,";",ce$position,";",ce$kmer)
cm$index <- paste0(cm$id,";",cm$position,";",cm$kmer)
length(unique(c(ce$index,cm$index))) # 163037
# overlap between 
dim(ce)[1]+dim(cm)[1]-length(unique(c(ce$index,cm$index))) # 77578

# sig (p.value < 0.001)
ce.sig <- ce[ce$pval_KO_vs_WT < 0.001, ] # 2758
cm.sig <- cm[cm$pval_KO_vs_WT < 0.001, ] # 3903

# sig (p.value < 0.01)
#ce.sig <- ce[ce$pval_KO_vs_WT < 0.01, ] # 
#cm.sig <- cm[cm$pval_KO_vs_WT < 0.01, ] # 

# full_join
ce.cm.sig <- dplyr::full_join(ce.sig, cm.sig, by = "index") # 6212
# remove na
ce.cm.sig.sub <- ce.cm.sig[(!is.na(ce.cm.sig$id.y))&(!is.na(ce.cm.sig$id.x)),] # 449

# multiply diff rate of CvsE(.x) and CvsM(.y) then order in decreasing order
ce.cm.sig.sub$diff_mod_rate_XY <- ce.cm.sig.sub$diff_mod_rate_KO_vs_WT.x*ce.cm.sig.sub$diff_mod_rate_KO_vs_WT.y
ce.cm.sig.sub <- ce.cm.sig.sub[order(-ce.cm.sig.sub$diff_mod_rate_XY),]

# remove columns "m6a_3.x" "m6a_2.x" "m6a_4.x" 
ce.cm.sig.sub.order <- ce.cm.sig.sub[,c(27,54,1:21,28:48)]
names(ce.cm.sig.sub.order) <-stringr::str_replace_all(names(ce.cm.sig.sub.order),"[.]x",".CvsE")
names(ce.cm.sig.sub.order) <-stringr::str_replace_all(names(ce.cm.sig.sub.order),"[.]y",".CvsM")
names(ce.cm.sig.sub.order)[2] <- "diff_mod_rate_CvsE*CvsM"

data.table::fwrite(ce.cm.sig.sub.order,"./ont/xpore/DRACH_motif_overlap_CvsE.CvsM_sig.txt",
                   quote = FALSE, sep = "\t")

# overlapped genes for GSEA ---
ce.cm.sig.sub.order <- cbind(ce.cm.sig.sub.order, read.table(text= as.character(ce.cm.sig.sub.order$id.CvsE),sep="|", fill=TRUE, na.strings=''))
ce.cm.sig.sub.order.sub <- ce.cm.sig.sub.order[,c(50,2)]
data.table::fwrite(ce.cm.sig.sub.order.sub,"./ont/xpore/DRACH_motif_overlap_CvsE.CvsM_sig4GSEA.txt",
                   quote = FALSE, sep = "\t")

# overlap by index ---
overlap.index <- intersect(ce.sig$index, cm.sig$index) # 449
ce.sig.overlap <- ce.sig[ce.sig$index%in%overlap.index, ]
cm.sig.overlap <- cm.sig[cm.sig$index%in%overlap.index, ]

# unique index --
length(unique(ce.cm.sig.sub.order$index)) #449
# kmer freq 
table(ce.cm.sig.sub.order$kmer.CvsE)
kmer.freq <- data.frame(table(ce.cm.sig.sub.order$kmer.CvsE))
# kmer mean diff mod rate of X*Y
library(dplyr)
kmer.mean <- ce.cm.sig.sub.order[,c(2,5)]

kmer.mean.sum <- kmer.mean %>% group_by(kmer.CvsE)%>%summarise(diff_mod_rate_Mean = mean(`diff_mod_rate_CvsE*CvsM`))

# join kmer freq and mean 
kmer.freq.mean <- dplyr::full_join(kmer.freq,kmer.mean.sum,by=c("Var1"="kmer.CvsE"))



# test sites (p.value < 0.001) overlap between CvsE and CvsM 
# refer to https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html
# hypergeometric test ---
# also see this https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
#https://stackoverflow.com/questions/8382806/hypergeometric-test-phyper
# an example 
# sum(dhyper(t:b, a, n - a, b)) # b<=a
# sum(dhyper(253:297,13448,13448+297-253-13448,297))

# use all unique index as n and dim(ce.sig)[1] < dim(cm.sig)[1]
sum(dhyper(length(overlap.index):dim(ce.sig)[1],dim(cm.sig)[1],length(unique(c(ce$index,cm$index)))-dim(cm.sig)[1],dim(ce.sig)[1])) # pvalue = 4.725653e-231
# fisher 
# matrix(c(n - union(A,B), setdiff(A,B), setdiff(B,A), intersect(A,B)), nrow=2)
fisher.test(matrix(c(length(unique(c(ce$index,cm$index)))-dim(ce.cm.sig)[1], dim(cm.sig)[1]-length(overlap.index),dim(ce.sig)[1]-length(overlap.index),length(overlap.index)), nrow=2), alternative="greater")

#Fisher's Exact Test for Count Data

#data:  
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is greater than 1
#95 percent confidence interval:
# 8.057418      Inf
#sample estimates:
#odds ratio 
#  8.828761 


# venn diagram for overlap index between ce.more and cm.more --
# Load library
library(VennDiagram)

# Generate 3 sets of 200 words
CvsE_sigIndex <- ce.sig$index
CvsM_sigIndex <- cm.sig$index


# Chart
venn.diagram(
  x = list(CvsM_sigIndex,CvsE_sigIndex),
  category.names = c("CvsM_sigIndex" , "CvsE_sigIndex "),
  imagetype = "tiff",
  filename = './ont/xpore/CvsEsigAndCvsMsig_overlap.tiff',
  output=TRUE
)

# color 
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(CvsM_sigIndex, CvsE_sigIndex),
  category.names = c("CvsM" , "CvsE"),
  filename = './ont/xpore/CvsEsigAndCvsMsig_overlap.tiff',
  output=TRUE,
  total.population = length(unique(c(ce$index,cm$index))),
  hyper.test = TRUE,
  lower.tail = FALSE,
  
  # Output features
  imagetype="tiff" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[-3],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-10, 10),
  cat.dist = c(0.035, 0.035),
  cat.fontfamily = "sans"
)







# add transcript, gene_id, gene_name, length, type(position)e.g. retained_intron
# try an easy way 
ce.more <- cbind(ce, read.table(text= as.character(ce$id),sep="|", fill=TRUE, na.strings=''))
cm.more <- cbind(cm, read.table(text= as.character(cm$id),sep="|", fill=TRUE, na.strings=''))

# remove columns 
ce.more$V3 <- ce.more$V4 <- ce.more$V9 <- cm.more$V3 <- cm.more$V4 <- cm.more$V9 <- NULL

# keep p.value < 0.001 and diff rate < -0.5
ce.more.dmr <- ce.more[ce.more$diff_mod_rate_KO_vs_WT <= -0.5 & ce.more$pval_KO_vs_WT < 0.001, ] # 506
cm.more.dmr <- cm.more[cm.more$diff_mod_rate_KO_vs_WT <= -0.5 & cm.more$pval_KO_vs_WT < 0.001, ] # 699

# unqiue gene 
length(unique(c(ce.more.dmr$V6, cm.more.dmr$V6))) # 1013
length(unique(ce.more.dmr$V6)) #449
length(unique(cm.more.dmr$V6)) #614

# sub for gene and diff mod rate 
ce.dmr.gene <- ce.more.dmr[,c(31,4)]
cm.dmr.gene <- cm.more.dmr[,c(31,4)]

# !duplicate
ce.dmr.gene <- ce.dmr.gene[!duplicated(ce.dmr.gene$V6), ]
cm.dmr.gene <- cm.dmr.gene[!duplicated(cm.dmr.gene$V6), ]

write.table(ce.dmr.gene, "./ont/xpore/ce.dmr.gene4GSEA.txt",col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)
write.table(cm.dmr.gene, "./ont/xpore/cm.dmr.gene4GSEA.txt",col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)


# visualize DRACH motif modification by volcano plot --------------------------
library(EnhancedVolcano)
ce <- data.table::fread("./ont/xpore/CvsE_woSM_majority_direction_kmer_diffmod.table_DRACH.txt")
cm <- data.table::fread("./ont/xpore/CvsM_woSM_majority_direction_kmer_diffmod.table_DRACH.txt")

ce$index <- paste0(ce$id,";",ce$position,";",ce$kmer)
cm$index <- paste0(cm$id,";",cm$position,";",cm$kmer)

ce$pval_KO_vs_WT <- ifelse(ce$pval_KO_vs_WT==0,10e-40,ce$pval_KO_vs_WT)
cm$pval_KO_vs_WT <- ifelse(cm$pval_KO_vs_WT==0,10e-40,cm$pval_KO_vs_WT)
# gene name + kmer + position ---
ce.more <- cbind(ce, read.table(text= as.character(ce$id),sep="|", fill=TRUE, na.strings=''))
cm.more <- cbind(cm, read.table(text= as.character(cm$id),sep="|", fill=TRUE, na.strings=''))

ce.more$label <- paste0(ce.more$V5,":",ce.more$kmer,":",ce.more$position)
cm.more$label <- paste0(cm.more$V5,":",cm.more$kmer,":",cm.more$position)

pdf("./ont/xpore/Figure_DMR_shEZH2.vs.Ctrl.pdf",width=10,height=10)
EnhancedVolcano(ce.more,
                lab = ce.more$label,
                selectLab = ce.more$label[abs(ce.more$diff_mod_rate_KO_vs_WT)>0.5&ce.more$pval_KO_vs_WT< 0.00001],
                x = 'diff_mod_rate_KO_vs_WT',
                y = 'pval_KO_vs_WT',
                pCutoff = 10e-50,
                FCcutoff = 0.1,
                xlab = bquote(~DMR(shEZH2.vs.Control)),
                xlim = c(-2,2)
)

dev.off()
# CvsM---
pdf("./ont/xpore/Figure_DMR_shMETTLE3.vs.Ctrl.pdf",width=10,height=10)
EnhancedVolcano(cm.more,
                lab = cm.more$label,
                selectLab = cm.more$label[abs(cm.more$diff_mod_rate_KO_vs_WT)>0.5&cm.more$pval_KO_vs_WT< 0.00001],
                x = 'diff_mod_rate_KO_vs_WT',
                y = 'pval_KO_vs_WT',
                pCutoff = 10e-50,
                FCcutoff = 0.1,
                xlab = bquote(~DMR(shMETTLE3.vs.Control)),
                xlim = c(-2,2)
)

dev.off()


# write output in excel 
library(xlsx)

# DMR > 0.5 or DMR < -0.5 and also p.value < 0.001

ce.more.KOvsWT.neg <- ce.more[ce.more$diff_mod_rate_KO_vs_WT < -0.5&ce.more$pval_KO_vs_WT<0.001, ]
ce.more.KOvsWT.pos <- ce.more[ce.more$diff_mod_rate_KO_vs_WT > 0.5&ce.more$pval_KO_vs_WT<0.001, ]

cm.more.KOvsWT.neg <- cm.more[cm.more$diff_mod_rate_KO_vs_WT < -0.5&cm.more$pval_KO_vs_WT<0.001, ]
cm.more.KOvsWT.pos <- cm.more[cm.more$diff_mod_rate_KO_vs_WT > 0.5&cm.more$pval_KO_vs_WT<0.001, ]

my_path <- "./ont/xpore/"
write.xlsx2(ce.more.KOvsWT.neg, paste0(my_path, "DRACH_modification_DMR.xlsx"), row.names = FALSE, sheetName = "shEZH2vsCtrl_negativeModRate") 
write.xlsx2(ce.more.KOvsWT.pos, paste0(my_path, "DRACH_modification_DMR.xlsx"), row.names = FALSE, sheetName = "shEZH2vsCtrl_positiveModRate", append = TRUE)                       # Append to second sheet
write.xlsx2(cm.more.KOvsWT.neg, paste0(my_path, "DRACH_modification_DMR.xlsx"), row.names = FALSE, sheetName = "shMETTLE3vsCtrl_negativeModRate", append = TRUE)                       # Append to third sheet
write.xlsx2(cm.more.KOvsWT.pos, paste0(my_path, "DRACH_modification_DMR.xlsx"), row.names = FALSE, sheetName = "shMETTLE3vsCtrl_positiveModRate", append = TRUE)                       # Append to third sheet


# April 05 #
# find the global A modification level for control siEZH2 siMETTLE 

ce <- data.table::fread("./ont/xpore/CvsE_woSM_majority_direction_kmer_diffmod.table_NNANN.sub.txt")
ce <- data.table::fread("./ont/xpore/CvsE_woSM_majority_direction_kmer_diffmod.table_DRACH.txt")

cm <- data.table::fread("./ont/xpore/CvsM_woSM_majority_direction_kmer_diffmod.table_NNANN.sub.txt")
cm <- data.table::fread("./ont/xpore/CvsM_woSM_majority_direction_kmer_diffmod.table_DRACH.txt")

ce.mod <- ce[,c(1:3,7:10)]
cm.mod <- cm[,c(1:3,7:10)]

# calculate mean mod rate ---
ce.mod$KO_mean <- rowSums(ce.mod[,c(4:5)],na.rm = TRUE)
ce.mod$WT_mean <- rowSums(ce.mod[,c(6:7)],na.rm = TRUE)

cm.mod$KO_mean <- rowSums(cm.mod[,c(4:5)],na.rm = TRUE)
cm.mod$WT_mean <- rowSums(cm.mod[,c(6:7)],na.rm = TRUE)

# wide to long for boxplot 
library(tidyr)
ce.mod.long <- gather(ce.mod, Condition, ModificationRate, KO_mean:WT_mean)
cm.mod.long <- gather(cm.mod, Condition, ModificationRate, KO_mean:WT_mean)

#ce.mod.long$Condition <- factor(ce.mod.long$Condition, levels = c("KO_mean","WT_mean"))
# boxplot ---
library(ggpubr)

p <- ggboxplot(ce.mod.long, x = "Condition", y = "ModificationRate",
               fill = "Condition", palette =c("#00AFBB", "#FC4E07"))
p + stat_compare_means(label.y = 2.5, 
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"))   # opposite compared to wilcox.test only   

ggsave("./ont/xpore/CvsE_global.mod.rate.DRACH.boxplot.pdf")
ggsave("./ont/xpore/CvsE_global.mod.rate.NNANN.boxplot.pdf")


# and the one-sided alternative "greater" is that x is shifted to the right of y
wilcox.test(ce.mod$KO_mean,ce.mod$WT_mean, alternative = "less")

#Wilcoxon rank sum test with continuity correction
#data:  ce.mod$KO_mean and ce.mod$WT_mean
#W = 7964386210, p-value = 1.223e-09
#alternative hypothesis: true location shift is greater than 0

# and the one-sided alternative "less" is that x is shifted to the left of y
wilcox.test(cm.mod$KO_mean,cm.mod$WT_mean, alternative = "less")
#Wilcoxon rank sum test with continuity correction
#data:  cm.mod$KO_mean and cm.mod$WT_mean
#W = 6255819092, p-value < 2.2e-16
#alternative hypothesis: true location shift is less than 0


# density plot for CvsE and CvsM samples
# Density plot with mean lines and marginal rug
# :::::::::::::::::::::::::::::::::::::::::::::::::::
# Change outline and fill colors by groups ("sex")
# Use custom palette
p2 <- ggdensity(cm.mod.long, x = "ModificationRate",
          add = "mean", rug = TRUE,
          color = "Condition", fill = "Condition",
          palette = c("#00AFBB", "#FC4E07"))

ggsave("./ont/xpore/CvsM_global.mod.rate.DRACH.density.pdf")
ggsave("./ont/xpore/CvsM_global.mod.rate.NNANN.density.pdf")



# April 29, 2022 Check if CvsE and CvsM has some correlation ---
ce <- data.table::fread("./ont/xpore/CvsE_woSM_majority_direction_kmer_diffmod.table_DRACH.txt")
cm <- data.table::fread("./ont/xpore/CvsM_woSM_majority_direction_kmer_diffmod.table_DRACH.txt")

# subset ---
ce.sub <- ce[,1:4]
cm.sub <- cm[,1:4]

# join two data set by first three columns ---
ce.cm.sub <- dplyr::full_join(ce.sub, cm.sub, by = c("id","position","kmer"))
ce.cm.sub <- na.omit(ce.cm.sub)
names(ce.cm.sub)[c(4,5)] <- c("diff_mod_rate_CvsE","diff_mod_rate_CvsM")

plot(ce.cm.sub$diff_mod_rate_CvsE, ce.cm.sub$diff_mod_rate_CvsM)
# corrleation 
ggscatter(ce.cm.sub, x = "diff_mod_rate_CvsE", y = "diff_mod_rate_CvsM",
          add = "reg.line", # Add regression line
          color = "#bdbdbd", 
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = -0.85, label.y = 0.6)  # Add correlation coefficient
# R = 0.57, p < 2.2e-16
ggsave("./ont/xpore/CvsE_CvsM_corplot.pdf")


# Xpore developer Yuk suggested to modify nanopolish eventalign output, modify the id column to only contain transcript id
# try this and see if xpore dataprep succeed
m.event <- data.table::fread("./ont/xpore/M2_eventalign_reads_index2_head.tsv")

contig <- NULL
for (i in 1:nrow(m.event)){
  contig[i] <- strsplit(m.event$contig[i],"[|]")[[1]][1]
}

m.event$contig <- contig

data.table::fwrite(m.event, "./ont/xpore/M2_eventalign_reads_index2_head.tsv")



