# March 14, 2022 ---
# change guppy_output summary_.txt filename for tombo 
# refer to this issue posted on tombo github
# https://github.com/nanoporetech/tombo/issues/371

c1 <- read.delim("./ont/sequencing_summary_1.txt",header = TRUE, check.names = FALSE,
                 stringsAsFactors = FALSE)

c1$filename <- paste0(c1$read_id,".fast5")

write.table(c1, "./ont/sequencing_summary_C1_4tombo.txt",col.names = TRUE,row.names = FALSE,
            sep = "\t", quote = FALSE)

# script for linux 
args <- commandArgs()

c1 <- read.delim(args[1],header = TRUE, check.names = FALSE,
                 stringsAsFactors = FALSE)

c1$filename <- paste0(c1$read_id,".fast5")

write.table(c1, args[2],col.names = TRUE,row.names = FALSE,
            sep = "\t", quote = FALSE)
