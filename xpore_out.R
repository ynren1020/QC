# filtering xpore diffmod results by selecting NNANN pattern kmers

# for linux 
args <- commandArgs(trailingOnly = TRUE)
#args[1] <- "./ont/xpore/majority_direction_kmer_diffmod.table.top100"
xpore <- data.table::fread(args[1])

# find NNANN pattern (middle position is A)
xpore$m6a <- stringr::str_sub(xpore$kmer,3,3)
# keep NNANN pattern
xpore.sub <- xpore[xpore$m6a=="A",]
# keep p.val < 0.05 
xpore.sub.sig <- xpore.sub[xpore.sub$pval_KO_vs_WT<0.05,]

# increase diff mod rate 
xpore.sub.sig <- xpore.sub.sig[order(xpore.sub.sig$diff_mod_rate_KO_vs_WT),]

# output
#args[2] <- "./ont/xpore/majority_direction_kmer_diffmod.table.top100.NNANN.sig"
data.table::fwrite(xpore.sub.sig, args[2],
                   quote = FALSE, sep = "\t")


################################################################################
# Find DRACH motif      --------------------------------------------------------
# [ATG][AG]AC[ATC]    --------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args[1] <- "./ont/xpore/majority_direction_kmer_diffmod.table.top100"
xpore <- data.table::fread(args[1])

# find NNANN pattern (middle position is A)
xpore$m6a_3 <- stringr::str_sub(xpore$kmer,3,3)
xpore$m6a_2 <- stringr::str_sub(xpore$kmer,2,2)
xpore$m6a_4 <- stringr::str_sub(xpore$kmer,4,4)
xpore$m6a_1 <- stringr::str_sub(xpore$kmer,1,1)
xpore$m6a_5 <- stringr::str_sub(xpore$kmer,5,5)
# keep NNANN pattern
xpore.sub <- xpore[xpore$m6a_3=="A"&xpore$m6a_4=="C"&(xpore$m6a_2%in%c("A","G"))&(xpore$m6a_1%in%c("A","G","T"))&(xpore$m6a_5%in%c("A","C","T")),]
# keep p.val < 0.05 
#xpore.sub.sig <- xpore.sub[xpore.sub$pval_KO_vs_WT<0.05,]

# increase diff mod rate 
xpore.sub <- xpore.sub[order(xpore.sub$diff_mod_rate_KO_vs_WT),]

# output
#args[2] <- "./ont/xpore/majority_direction_kmer_diffmod.table.top100.NNANN.sig"
data.table::fwrite(xpore.sub, args[2],
                   quote = FALSE, sep = "\t")

