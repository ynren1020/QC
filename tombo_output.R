# March 22, 2022 convert wig or bedgraph output of tombo to txt file ---

# p-value output of CvsE and CvsM 

ce <- read.delim("./ont/tombo/tombo_results.statistic.plus.gene_name.CvsE.txt", header = FALSE,
                 stringsAsFactors = FALSE, sep=" ")
cm <- read.delim("./ont/tombo/tombo_results.statistic.plus.gene_name.CvsM.txt", header = FALSE,
                 stringsAsFactors = FALSE, sep=" ")

overlap <- intersect(ce$V2,cm$V2)