# overlap results between m6anet and supp 12 ---


c1 <- data.table::fread("./ont/m6anet/C1.data.result.csv.gz")

# read in gtf v38 transcript annotation file ---
transcript.gtf <- data.table::fread("gencode.v38.annotation.transcript.gtf")
transcript.gtf$transcript_id <- stringr::str_sub(transcript.gtf$transcript_id,1,15)

# join together ---
c1.transcript.gtf <- dplyr::left_join(c1, transcript.gtf, by = "transcript_id")

c1.transcript.gtf$type <- NULL
c1.transcript.gtf$gene_id <- c1.transcript.gtf$gene_type <- NULL
names(c1.transcript.gtf)[5] <- "seqnames"
# overlap with supp 12 table ---
results <- data.table::fread("./ont/Yi/SupplementaryTable12_hg38.txt")
results$group_name <- NULL
results$group <- NULL

setkey(c1.transcript.gtf,seqnames, start, end)
overlaps.res.m6a <- foverlaps(results,c1.transcript.gtf,type="within",nomatch=NULL) 

# order ---
overlaps.res.m6a.order <- overlaps.res.m6a %>% group_by(transcript_id) %>% arrange(desc(probability_modified))

write.table(overlaps.res.m6a.order,"./ont/m6anet/C1_Yi_overlap.txt",quote = FALSE,
            sep = "\t", row.names = FALSE)


