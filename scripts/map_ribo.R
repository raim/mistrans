
## fasta
seq <- segmenTools::readFASTA("NM_000016.fa")

## ribo density from HDPR
rib <- unlist(strsplit(readLines("NM_000016.txt"), "_"))

length(unlist(strsplit(seq[[1]]$seq,"")))
