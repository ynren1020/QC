# April 21 2022 ---
# use coordiante (ranges) to find overlap between xpore results and Yi results ---
library(data.table)
library(dplyr)

# Yi results ---
results <- data.table::fread("./ont/Yi/SupplementaryTable12_hg38.txt")
results$group_name <- NULL
results$group <- NULL
# Xpore results ---
# read in gtf v38 transcript annotation file ---
transcript.gtf <- data.table::fread("gencode.v38.annotation.transcript.gtf")
# check the results of CvsE from Xpore ------
ce <- data.table::fread("./ont/xpore/CvsE_woSM_majority_direction_kmer_diffmod.table_NNANN.sub.txt") # 1645175      10
#ce.wona <- na.omit(ce) # 1052220      10

ce.more <- cbind(ce, read.table(text= as.character(ce$id),sep="|")) # 1052220      18
ce.more$V9 <- ce.more$V3 <- ce.more$V4 <- NULL

# remove .* from ENST***.* for join ---
ce.more$V1 <- stringr::str_sub(ce.more$V1,1,15)
transcript.gtf$transcript_id <- stringr::str_sub(transcript.gtf$transcript_id,1,15)

# join with transcript.gtf
ce.more.gtf <- dplyr::left_join(ce.more,transcript.gtf, by = c("V1"="transcript_id"))

# find overlap between ce.more.gtf and results ---
# as x
# 1
results.sub <- data.table(seqnames=results$seqnames,start=results$start,end=results$end,name=results$name) #results[,c(1,2,3,6)]
# 2
results
# as y
# 1
ce.more.sub <- data.table(seqnames=ce.more.gtf$chrom,start=ce.more.gtf$start,end=ce.more.gtf$end,id=ce.more.gtf$id)
# 2
names(ce.more.gtf)[17] <- "seqnames"
# 1
setkey(ce.more.sub, seqnames,start, end)
# 2
setkey(ce.more.gtf, seqnames,start, end)
# 1
overlaps.res <- foverlaps(results.sub,ce.more.sub,type="within",nomatch=NULL) # 
# 2
overlaps.res <- foverlaps(results,ce.more.gtf,type="within",nomatch=NULL) # 
#names(overlaps.res)[c(12,13,14,15)] <- c("transcript_id","gene_id","transcript_name","gene_name")
#overlaps.res$V7 <- NULL
#overlaps.res <- overlaps.res[,c(2:17,1,18:20,32,21:31)]
overlaps.res$gene_type <- NULL
overlaps.res$transcript_type <- NULL
overlaps.res$V5 <- overlaps.res$V6 <- overlaps.res$V7 <- overlaps.res$V8 <- NULL

# significant ones in xpore ---
overlaps.res.sig <- overlaps.res[overlaps.res$pval_KO_vs_WT<0.05,]
# require two replicates ---
overlaps.res.sig.wona <- na.omit(overlaps.res.sig)
names(overlaps.res.sig.wona)[c(21:32)] <- paste0("supp12_",names(overlaps.res.sig.wona)[c(21:32)] )
overlaps.res.sig.wona <- overlaps.res.sig.wona[order(overlaps.res.sig.wona$diff_mod_rate_KO_vs_WT),]

data.table::fwrite(overlaps.res.sig.wona,"./ont/Yi/ce.more.supp12.v38.overlap.txt",quote = FALSE,
                   sep = "\t")

# p.value < 0.001
overlaps.res.sig.wona.sub <- overlaps.res.sig.wona[overlaps.res.sig.wona$pval_KO_vs_WT<0.001, ]
data.table::fwrite(overlaps.res.sig.wona.sub,"./ont/Yi/ce.more.supp12.v38.overlap.sub.txt",quote = FALSE,
                   sep = "\t")
# - strand use end - position to find the modification site
# + strand use start + position to find the modification site 
length(unique(overlaps.res.sig.wona.sub$supp12_gene)) #1683
length(unique(overlaps.res.sig$gene)) #3300
length(unique(results$gene)) #6001

###########################################################################
# after transcript overlapping, find the region overlap, exon, start_codon, stop_codon, UTR
# x
overlaps.res.sig.wona
# y
type.gtf <- data.table::fread("gencode.v38.annotation.type.gtf")
names(type.gtf)[c(1,4,5)] <- c("seqnames","supp12_i.start","supp12_i.end")
 
setkey(type.gtf, seqnames,supp12_i.start, supp12_i.end)
overlaps.res.type <- foverlaps(overlaps.res.sig.wona,type.gtf,type="within",nomatch=NULL) # 

overlaps.res.type$source <- overlaps.res.type$V6 <- overlaps.res.type$V8 <- overlaps.res.type$gene_type <- NULL
overlaps.res.type$gene_id2 <- overlaps.res.type$transcript_id2 <- NULL




