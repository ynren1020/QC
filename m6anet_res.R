# April 13, 2022 analyze m6A results from m6anet method 
# July 21, 2022 replot 

c1 <- data.table::fread("./ont/m6anet/C1.data.result.csv.gz")
c2 <- data.table::fread("./ont/m6anet/C2.data.result.csv.gz")
e1 <- data.table::fread("./ont/m6anet/E1.data.result.csv.gz")
e2 <- data.table::fread("./ont/m6anet/E2.data.result.csv.gz")
m1 <- data.table::fread("./ont/m6anet/M1.data.result.csv.gz")
m2 <- data.table::fread("./ont/m6anet/M2.data.result.csv.gz")

# replication
c1$rep <- "C1"
c2$rep <- "C2"
e1$rep <- "E1"
e2$rep <- "E2"
m1$rep <- "M1"
m2$rep <- "M2"

# condition 
c1$condition <- "Control"
c2$condition <- "Control"
e1$condition <- "shEZH2"
e2$condition <- "shEZH2"
m1$condition <- "shMETTLE3"
m2$condition <- "shMETTLE3"

# rbind all 
mod <- data.table::rbindlist(list(c1,c2,e1,e2,m1,m2))
mod$condition <- factor(mod$condition, levels=c("Control","shEZH2","shMETTLE3"))
# boxplot to check the global methylation level across samples ---
library(ggpubr)
stat.test <- compare_means(probability_modified ~ condition, data= mod, method = "wilcox.test",p.adjust.method = "fdr") #, ref.group = "Control"
stat.test2 <- compare_means(probability_modified ~ rep, data= mod, method = "wilcox.test",p.adjust.method = "bonferroni")

  
stat.test2
# opposite compared to wilcox.test only   
# A tibble: 3 × 8
#.y.                  group1  group2            p    p.adj p.format p.signif method  
#<chr>                <chr>   <chr>         <dbl>    <dbl> <chr>    <chr>    <chr>   
#1 probability_modified Control shEZH2    0         0        < 2e-16  ****     Wilcoxon
#2 probability_modified Control shMETTLE3 0         0        < 2e-16  ****     Wilcoxon
#3 probability_modified shEZH2  shMETTLE3 0.0000159 0.000016 1.6e-05  ****     Wilcoxon

#my_comparisons <- list( c("Control", "shEZH2"), c("shEZH2", "shMETTLE3"), c("Control", "shMETTLE3") )
p <- ggboxplot(mod, x = "condition", y = "probability_modified",
               fill = "condition", 
               palette =c("#a50f15", "#fb6a4a","#fcbba1"),
               xlab = "",
               ylab = "Probability of m6A",
               width = 0.3) +
theme(legend.position="none")

p + stat_pvalue_manual(
  stat.test, 
  tip.length = 0.02,
  step.increase = 0.075,
  y.position = 1.1,
  label = "p.signif"
  #position = position_dodge(0.2)
)

ggsave("./ont/m6anet/global.mod.probablity.boxplot.pdf")

#july 21 2022
ggsave("./ont/figure/m6anet.global.mod.probablity.boxplot.svg",width=5,height=5)


# find common transcript position across samples -------------------------------
c1 <- data.table::fread("./ont/m6anet/C1.data.result.csv.gz")
c2 <- data.table::fread("./ont/m6anet/C2.data.result.csv.gz")
e1 <- data.table::fread("./ont/m6anet/E1.data.result.csv.gz")
e2 <- data.table::fread("./ont/m6anet/E2.data.result.csv.gz")
m1 <- data.table::fread("./ont/m6anet/M1.data.result.csv.gz")
m2 <- data.table::fread("./ont/m6anet/M2.data.result.csv.gz")


control <- dplyr::full_join(c1,c2,by=c("transcript_id","transcript_position"))
control <- na.omit(control)

shE <- dplyr::full_join(e1,e2,by=c("transcript_id","transcript_position"))
shE <- na.omit(shE)

shM <- dplyr::full_join(m1,m2,by=c("transcript_id","transcript_position"))
shM <- na.omit(shM)

mod2 <- purrr::reduce(list(control,shE,shM), dplyr::full_join, by = c("transcript_id","transcript_position"))

# rename columns ---
names(mod2)[c(3:14)] <- c(paste0(c("n_reads","probability_modified"),".C1"),
                          paste0(c("n_reads","probability_modified"),".C2"),
                          paste0(c("n_reads","probability_modified"),".E1"),
                          paste0(c("n_reads","probability_modified"),".E2"),
                          paste0(c("n_reads","probability_modified"),".M1"),
                          paste0(c("n_reads","probability_modified"),".M2"))

# mean within condition
mod2.mean <- data.frame(transcript_id = mod2$transcript_id,
                        transcript_position = mod2$transcript_position,
                        probability.control = rowMeans(mod2[,c(4,6)],na.rm=TRUE),
                        probability.shEZH2 = rowMeans(mod2[,c(8,10)],na.rm=TRUE),
                        probability.shMETTLE3 = rowMeans(mod2[,c(12,14)],na.rm=TRUE))

mean(mod2.mean$probability.control,na.rm=TRUE) # 0.2892524
mean(mod2.mean$probability.shEZH2,na.rm=TRUE) # 0.2626877
mean(mod2.mean$probability.shMETTLE3,na.rm=TRUE) # 0.2556337

# test ---
wilcox.test(mod2.mean$probability.control, mod2.mean$probability.shEZH2) # p-value < 2.2e-16
wilcox.test(mod2.mean$probability.control, mod2.mean$probability.shMETTLE3) # p-value < 2.2e-16
wilcox.test(mod2.mean$probability.shMETTLE3, mod2.mean$probability.shEZH2) # p-value = 0.001795

# calculate the mean difference ---
mod2.mean.diff <- mod2.mean
mod2.mean.diff$diff_ctrl.shE <- mod2.mean.diff$probability.control - mod2.mean.diff$probability.shEZH2
mod2.mean.diff$diff_ctrl.shM <- mod2.mean.diff$probability.control - mod2.mean.diff$probability.shMETTLE3
plot(mod2.mean.diff$diff_ctrl.shE,mod2.mean.diff$diff_ctrl.shM)
cor.test(mod2.mean.diff$diff_ctrl.shE,mod2.mean.diff$diff_ctrl.shM)
#Pearson's product-moment correlation
#data:  mod2.mean.diff$diff_ctrl.shE and mod2.mean.diff$diff_ctrl.shM
#t = 267.54, df = 139862, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.5783480 0.5852814
#sample estimates:
#      cor 
#0.5818253 

# scatter plot ---
ggscatter(mod2.mean.diff, x = "diff_ctrl.shE", y = "diff_ctrl.shM",
          add = "reg.line", # Add regression line
          color = "#bdbdbd", 
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 0.3, label.y = 0.5)  # Add correlation coefficient

ggsave("./ont/m6anet/ctrl.shEZH2_ctrl.shMETTLE2_corplot.pdf")


# gather, wide to long to do boxplot ---
mod2.mean.long <- tidyr::gather(mod2.mean, condition, probability, -c(transcript_id,transcript_position))
mod2.mean.long$condition <- stringr::str_replace_all(mod2.mean.long$condition,"probability.","")
mod2.mean.long$condition <- factor(mod2.mean.long$condition, levels = c("control","shEZH2","shMETTLE3"))

p <- ggboxplot(mod2.mean.long, x = "condition", y = "probability",
               fill = "condition", palette =c("#a50f15", "#fb6a4a","#fcbba1"))

my_comparisons <- list( c("control", "shEZH2"), c("control", "shMETTLE3"), c("shEZH2", "shMETTLE3") )
p + stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 1.5)

ggsave("./ont/m6anet/global.mod.probablity.boxplot.commonTranscript.pdf")

####################annotate with gtf ##################################################################
transcript.gtf <- data.table::fread("gencode.v38.annotation.transcript.gtf")
transcript.gtf$transcript_id <- stringr::str_sub(transcript.gtf$transcript_id,1,15)
mod2.mean.long.gtf <- dplyr::left_join(mod2.mean.long, transcript.gtf, by = "transcript_id")

names(mod2.mean.long.gtf)[5] <- "seqnames"
mod2.mean.long.gtf <- as.data.table(mod2.mean.long.gtf)
# overlap with supp 12 table ---
results <- data.table::fread("./ont/Yi/SupplementaryTable12_hg38.txt")
results$group_name <- NULL
results$group <- NULL

setkey(mod2.mean.long.gtf,seqnames, start, end)
overlaps.res.m6anet <- foverlaps(results,mod2.mean.long.gtf,type="within",nomatch=NULL) 

# boxplot 
p <- ggboxplot(overlaps.res.m6anet, x = "condition", y = "probability",
               fill = "condition", palette =c("#a50f15", "#fb6a4a","#fcbba1"))

my_comparisons <- list( c("control", "shEZH2"), c("control", "shMETTLE3"), c("shEZH2", "shMETTLE3") )
p + stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 1.5)

ggsave("./ont/m6anet/global.mod.probablity.boxplot.overlapSupp12.pdf")



