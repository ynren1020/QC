# July 22, 2022 Qingshu PARP project 
# compare other groups with W groups
# /panfs/home/yang4414/renxx275/project_yren/QC/qingshu
# besides limma, the other method used wilcox test for differential gene expression analysis

# load package ---
library("edgeR")
#library("sva")
library("R.utils")
library("limma")
library(gplots)
library(dplyr)

# read data ---
dat <- read.delim("./qingshu/PARP_qingshu.txt",
                  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
names(dat)[7:18]
#[1] "W1_sorted.bam" "W2_sorted.bam" "W3_sorted.bam" "A1_sorted.bam" "A2_sorted.bam" "A3_sorted.bam" "E1_sorted.bam"
#[8] "E2_sorted.bam" "E3_sorted.bam" "F1_sorted.bam" "F2_sorted.bam" "F3_sorted.bam"

group <- c(rep(1,3), rep(2,3), rep(3,3),
           rep(4,3))
group <- as.factor(group)  

# count only
countdata <- dat[,-(1:6)]
rownames(countdata) <- dat[,1]
colnames(countdata)

# create DGE object
y<-DGEList(countdata, group = group)
y$genes <- data.frame(Length=dat$Length)

# cpm 
myCPM <- cpm(y$counts)
head(myCPM)
plot(myCPM[,1],countdata[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 1.7 CPM 
abline(h=10,v=1.7)

# Which values in myCPM are greater than 1.7? count > 10
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
table(rowSums(thresh))

# we would like to keep genes that have at least 15 TRUES in each row of thresh
keep <- rowSums(thresh) >= 3
# Subset the rows of countdata to keep the more highly expressed genes
y <- y[keep,keep.lib.sizes=FALSE]
summary(keep)
#Mode   FALSE    TRUE 
#logical   43065   16320 

##unsupervised clustering##
plotMDS(y)

# create design matrix 
design <- model.matrix(~0+group);
rownames(design) <- colnames(y)
design
#group1 group2 group3 group4
#W1_sorted.bam      1      0      0      0
#W2_sorted.bam      1      0      0      0
#W3_sorted.bam      1      0      0      0
#A1_sorted.bam      0      1      0      0
#A2_sorted.bam      0      1      0      0
#A3_sorted.bam      0      1      0      0
#E1_sorted.bam      0      0      1      0
#E2_sorted.bam      0      0      1      0
#E3_sorted.bam      0      0      1      0
#F1_sorted.bam      0      0      0      1
#F2_sorted.bam      0      0      0      1
#F3_sorted.bam      0      0      0      1

# normalization by the library sizes
y <- calcNormFactors(y,method="TMM")

tmm <- cpm(y$counts)

# output rpkm --
logrpkm<-rpkm(y$counts,gene.length = y$genes$Length,log = TRUE)
#logrpkm <- as.data.frame(logrpkm)
#logrpkm$gene <- rownames(logrpkm)
write.table(logrpkm,file="./qingshu/PARP_WAEF_logrpkm.txt",row.names=T,quote=F,sep="\t")


# normalise the read counts with 'voom' function
v <- voom(y,design,plot = TRUE)
# extract the normalised read counts
counts.voom <- v$E

# fit linear model for each gene given a series of libraries
fit <- lmFit(v, design)

# construct the contrast matrix corresponding to specified contrasts of a set of parameters
# AvsW
matrix.2vs1 <- makeContrasts(group2-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.2vs1 <- contrasts.fit(fit, matrix.2vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.2vs1 <- eBayes(fit.2vs1)
options(digits=3)

summary(decideTests(fit.2vs1, p.value=0.05,lfc=0.5))
#group2 - group1
#Down                143
#NotSig           16145
#Up                   32

num = length(fit.2vs1$genes$Length)
degs.2vs1 <- topTable(fit.2vs1, coef="group2 - group1", confint=TRUE, number = num)
degs.2vs1$gene <- rownames(degs.2vs1)
degs.2vs1 <- degs.2vs1[order(-degs.2vs1$logFC),]
degs.2vs1 <- degs.2vs1[,c(10,1:9)]
write.table(degs.2vs1, file=paste0('./qingshu/PARP_AvsW','.degs_lowcutoff.txt'), sep='\t',row.names = TRUE, quote = FALSE)

# gsea ---
avsw.rnk <- degs.2vs1[,c(1,3)]
write.table(avsw.rnk, file=paste0('./qingshu/PARP_AvsW','.gsea.txt'), sep='\t',row.names = FALSE, quote = FALSE,col.names = FALSE)

# logFC *abs(log(p.adj))
avsw.rnk <- data.frame(gene_name=rownames(degs.2vs1),metric=degs.2vs1$logFC*abs(log(degs.2vs1$adj.P.Val+0.001)))
avsw.rnk.order <- avsw.rnk[order(-avsw.rnk$metric),]
write.table(avsw.rnk.order, file=paste0('./qingshu/PARP_AvsW','.gsea.logFCandPval.txt'), sep='\t',row.names = FALSE, quote = FALSE,col.names = FALSE)

# second compare----------
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
# EvsW

matrix.3vs1 <- makeContrasts(group3-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.3vs1 <- contrasts.fit(fit, matrix.3vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.3vs1 <- eBayes(fit.3vs1)
options(digits=3)

summary(decideTests(fit.3vs1, p.value=0.05,lfc=0.5))
#group3 - group1
#Down               811
#NotSig           15232
#Up                  277

num = length(fit.3vs1$genes$Length)
degs.3vs1 <- topTable(fit.3vs1, coef="group3 - group1", confint=TRUE, number = num)
degs.3vs1$gene <- rownames(degs.3vs1)
degs.3vs1 <- degs.3vs1[order(-degs.3vs1$logFC),]
degs.3vs1 <- degs.3vs1[,c(10,1:9)]
write.table(degs.3vs1, file=paste0('./qingshu/PARP_EvsW','.degs_lowcutoff.txt'), sep='\t',row.names = TRUE, quote = FALSE)

# gsea ---
avsw.rnk <- degs.3vs1[,c(1,3)]
write.table(avsw.rnk, file=paste0('./qingshu/PARP_EvsW','.gsea.txt'), sep='\t',row.names = FALSE, quote = FALSE,col.names = FALSE)

# logFC *abs(log(p.adj))
avsw.rnk <- data.frame(gene_name=rownames(degs.3vs1),metric=degs.3vs1$logFC*abs(log(degs.3vs1$adj.P.Val+0.001)))
avsw.rnk.order <- avsw.rnk[order(-avsw.rnk$metric),]
write.table(avsw.rnk.order, file=paste0('./qingshu/PARP_EvsW','.gsea.logFCandPval.txt'), sep='\t',row.names = FALSE, quote = FALSE,col.names = FALSE)

# compare FvsW ---------------------------------------------------------------------------------------------------------
matrix.4vs1 <- makeContrasts(group4-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.4vs1 <- contrasts.fit(fit, matrix.4vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.4vs1 <- eBayes(fit.4vs1)
options(digits=3)

summary(decideTests(fit.4vs1, p.value=0.05,lfc=0.5))
#group2 - group1
#Down                59
#NotSig           16296
#Up                   36

num = length(fit.4vs1$genes$Length)
degs.4vs1 <- topTable(fit.4vs1, coef="group4 - group1", confint=TRUE, number = num)
degs.4vs1$gene <- rownames(degs.4vs1)
degs.4vs1 <- degs.4vs1[order(-degs.4vs1$logFC),]
degs.4vs1 <- degs.4vs1[,c(10,1:9)]
write.table(degs.4vs1, file=paste0('./qingshu/PARP_FvsW','.degs_lowcutoff.txt'), sep='\t',row.names = TRUE, quote = FALSE)

# gsea ---
avsw.rnk <- degs.4vs1[,c(1,3)]
write.table(avsw.rnk, file=paste0('./qingshu/PARP_FvsW','.gsea.txt'), sep='\t',row.names = FALSE, quote = FALSE,col.names = FALSE)

# logFC *abs(log(p.adj))
avsw.rnk <- data.frame(gene_name=rownames(degs.4vs1),metric=degs.4vs1$logFC*abs(log(degs.4vs1$adj.P.Val+0.001)))
avsw.rnk.order <- avsw.rnk[order(-avsw.rnk$metric),]
write.table(avsw.rnk.order, file=paste0('./qingshu/PARP_FvsW','.gsea.logFCandPval.txt'), sep='\t',row.names = FALSE, quote = FALSE,col.names = FALSE)


# similarity of ordered gene list 
library(OrderedList)

EvsWandAwsW <- compareLists(degs.3vs1$gene, degs.2vs1$gene, mapping = NULL, two.sided = TRUE, B = 1000, alphas = NULL, min.weight = 1e-5, invar.q = 0.5)
 
EvsWandFwsW <- compareLists(degs.3vs1$gene, degs.4vs1$gene, mapping = NULL, two.sided = TRUE, B = 1000, alphas = NULL, min.weight = 1e-5, invar.q = 0.5)


# p.val.adj < 0.05, logFC > 1 
AvsWandEvsW.comm <- intersect(degs.3vs1$gene[degs.3vs1$logFC>1&degs.3vs1$adj.P.Val<0.05],
                              degs.2vs1$gene[degs.2vs1$logFC>1&degs.2vs1$adj.P.Val<0.05])

AvsWandEvsW.comm
#[1] "SDCCAG8"       "RP11-341D18.8" "MT-TM"         "PRKXP1"   

FvsWandEvsW.comm <- intersect(degs.4vs1$gene[degs.4vs1$logFC>1&degs.4vs1$adj.P.Val<0.05],
                              degs.3vs1$gene[degs.3vs1$logFC>1&degs.3vs1$adj.P.Val<0.05])
FvsWandEvsW.comm
#[1] "RP11-341D18.8" "KB-1440D3.16"  "KCTD9P2"       "SDCCAG8" 

# heatmap by using differntial genes between EvsW (3vs1)
library(ComplexHeatmap)
library(circlize)

col_fun2 = colorRamp2(c(-3.0, 0, 3.0), c("blue", "white", "red"))

df = data.frame(type = c(rep("W", 3), rep("A", 3), rep("E", 3), 
                         rep("F", 3)))

ha = HeatmapAnnotation(df = df,col = list(type = c("W" =  "#2b8cbe", "A" = "#fa9fb5",
                                                   "E" = "#ffeda0", "F" = "#a1d99b")),
                       simple_anno_size =unit(3, "mm"))
# use DEG from EvsW 
logrpkm.sub <- logrpkm[rownames(logrpkm)%in%degs.3vs1$gene[degs.3vs1$adj.P.Val<0.05&abs(degs.3vs1$logFC)>0.5],]

# use DEG from AvsW
logrpkm.avsw <- logrpkm[rownames(logrpkm)%in%degs.2vs1$gene[degs.2vs1$adj.P.Val<0.05&abs(degs.2vs1$logFC)>0.5],]


# to zscore ---
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# EvsW
logrpkm.z <- as.matrix(apply(logrpkm.sub, 1, cal_z_score))
# AvsW 
logrpkm.z <- as.matrix(apply(logrpkm.avsw, 1, cal_z_score))
dim(logrpkm.z) #175 12

logrpkm.z.b<-logrpkm.z
logrpkm.z.b[logrpkm.z.b>3]<-3
logrpkm.z.b[logrpkm.z.b< -3]<- -3

# transpose 
logrpkm.z.bT<-as.data.frame(t(logrpkm.z.b))
colnames(logrpkm.z.bT) <- c(paste0("W",1:3),paste0("A",1:3),paste0("E",1:3),paste0("F",1:3))

# heatmap

pdf("./qingshu/FigureEvsW.DEGs.heatmap.pdf") # July27,2022
Heatmap(as.matrix(logrpkm.z.bT),
        col = col_fun2,
        name = "Z-Score", 
        cluster_rows=TRUE,
        cluster_columns =TRUE,
        cluster_column_slices = TRUE,
        top_annotation = ha,
        #split = genes.bind$Cytoband,
        #row_labels = arm$arm,
        row_names_gp = gpar(fontsize = 6,fontface = "bold"),
        show_column_names = TRUE,
        show_row_names = FALSE,
        show_row_dend = FALSE,
        width = unit(12, "cm"), 
        height = unit(8, "cm")
)

dev.off()

# AvsW DEGs ---
pdf("./qingshu/FigureAvsW.DEGs.heatmap.pdf") # July27,2022
Heatmap(as.matrix(logrpkm.z.bT),
        col = col_fun2,
        name = "Z-Score", 
        cluster_rows=TRUE,
        cluster_columns =TRUE,
        cluster_column_slices = TRUE,
        top_annotation = ha,
        #split = genes.bind$Cytoband,
        #row_labels = arm$arm,
        row_names_gp = gpar(fontsize = 6,fontface = "bold"),
        show_column_names = TRUE,
        show_row_names = FALSE,
        show_row_dend = FALSE,
        width = unit(6, "cm"), 
        height = unit(8, "cm")
)

dev.off()

# Venn diagram for DEGs --------------------
# up.regulated --- 
avsw.up <- rownames(degs.2vs1)[degs.2vs1$logFC>0.5&degs.2vs1$adj.P.Val<=0.05] # 32
fvsw.up <- rownames(degs.4vs1)[degs.4vs1$logFC>0.5&degs.4vs1$adj.P.Val<=0.05] # 36

# down.regulated ---
avsw.down <- rownames(degs.2vs1)[degs.2vs1$logFC< -0.5&degs.2vs1$adj.P.Val<=0.05] # 141
fvsw.down <- rownames(degs.4vs1)[degs.4vs1$logFC< -0.5&degs.4vs1$adj.P.Val<=0.05] # 59

# library
library(VennDiagram)
#Make the plot
venn.diagram(
  x = list(
    avsw.up,
    fvsw.up,
    avsw.down,
    fvsw.down
  ),
  category.names = c("AvsW.Up (32)" , "FvsW.Up (36)" , "AvsW.Down (141)", "FvsW.Down (59)"),
  filename = './qingshu/Figure.AvsWandFvsW.DEGsoverlap.venn.png',
  output = TRUE ,
  imagetype="tiff" ,
  height = 500 , 
  width = 500 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c('#c51b7d', '#4d9221', '#fde0ef', '#e6f5d0'),
  fill = c(adjustcolor("#c51b7d",alpha.f=0.3), adjustcolor('#4d9221',alpha.f=0.3), adjustcolor('#fde0ef',alpha.f=0.3), adjustcolor('#e6f5d0',alpha.f=0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-190, 190, 180, -180),
  cat.dist = c(-0.20, -0.20, -0.085, -0.085),
  cat.fontfamily = "sans",
  cat.col = c('#c51b7d', '#4d9221', '#fde0ef', '#e6f5d0')
  #rotation = 1
)


