library(segmenTools)
options(stringsAsFactors=FALSE)

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")

feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")
tpmap.file <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
saap.file <- file.path(proj.path,"processedData","unique_saap_blast.tsv")
bp.ids <- file.path(proj.path,"processedData","all_proteins.tsv")

fig.path <- file.path(proj.path,"figures","saap_mapping")
dir.create(fig.path)

ftyp <- "png"
if ( !interactive() ) ftyp <- "pdf"

## READ SAAP BLAST RESULTS
saap <- read.delim(saap.file, header=FALSE)
colnames(saap) <- c("SAAP","protein","identity", "mismatches",
                  "alen", "qlen", "slen", "sstart", "send", "e", "bitscore")

## TODO: sort for best hits and remove duplicates as in get_protein_match.R
## TODO: fuse this into get_protein_match.R or map_peptides3.R
## SORT & FILTER
## sort by
## 1) identity,
## 2) -mismatches (perhaps redundant?),
## 3) whether protein is the annotated MANE, and
## 4) number of annotations.
## 5) remove duplicated after sorting!
saap <- saap[order(saap$identity,
               -saap$mismatches),]

## remove duplicated, which should appear below better hits!
saap <- saap[!duplicated(paste0(saap$SAAP,saap$gene)),]

plotdev(file.path(fig.path,"saap_blast"),
        res=300, width=7, height=3.5, type=ftyp)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(saap$identity)
barplot(table(saap$mismatches), xlab="number of mismatches")
dev.off()

pm <- saap[which(saap$mismatches==0),]
