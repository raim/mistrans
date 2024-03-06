
## load blast of base peptides against ensembl+mutationes, and
## find best matching protein for each unique BP

library(segmenTools)
options(stringsAsFactors=FALSE)

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")

feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")
tpmap.file <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
bp.file <- file.path(proj.path,"processedData","unique_bp_blast.tsv")
bp.ids <- file.path(proj.path,"processedData","all_proteins.tsv")

out.file <- file.path(proj.path,"processedData","bp_mapped.tsv")

## LOAD BLAST RESULTS
bp <- read.delim(bp.file, header=FALSE)
colnames(bp) <- c("BP","protein","identity", "mismatches",
                  "alen", "qlen", "slen", "sstart", "send", "e", "bitscore")
bp$ensembl <- sub("_.*", "", bp$protein)

## LOAD FEATURE FILE
features <- read.delim(feature.file)
## use number of GO annotations as a filter criterium: best annotated
gol <- unlist(lapply(strsplit(features$GO, ";"),length))

## load protein-transcript map
tpmap <- read.delim(file=tpmap.file, header=FALSE, row.names=2)
## reverse map transcript-protein
pampt <- matrix(rownames(tpmap), ncol=1)
rownames(pampt) <- tpmap[,1]

## map to genes
## takes long - store as file and reload
index.file <- file.path(proj.path,"processedData","protein_index.txt")
if ( file.exists(index.file) ) {
    gtab <- read.delim(index.file, header=FALSE, row.names=1)
    gidx <- unlist(gtab)
    names(gidx) <- rownames(gtab)
} else {
    uprt <- unique(bp$ensembl)
    gidx <- sapply(uprt, function(x) grep(x, features$proteins))
    gidx <- unlist(gidx)
    write.table(cbind(names(gidx), gidx), file=index.file, sep="\t",
                col.names=FALSE, row.names=FALSE, quote=FALSE)
}

bp$gene <- features$ID[gidx[bp$ensembl]]
bp$transcript <- tpmap[bp$ensembl,1]

## FUSE MANE and canonical annotation
mains <- data.frame(MANE=features$MANE[gidx[bp$ensembl]],
                    canonical=features$canonical[gidx[bp$ensembl]])
mains[mains==""] <- NA
tmp <- apply(mains,1, unique)
tmp <- lapply(tmp, function(x) x[!is.na(x)])
tmp <- lapply(tmp, function(x) {if (length(x)==0) x<-NA;x})

bp$MANE <- unlist(tmp)
bp$MANE.protein <- pampt[match(bp$MANE,rownames(pampt)),1]

bp$numGO <- gol[gidx[bp$ensembl]]


## filter - only those where we found a gene
all <- bp ## store to analyze missing below
bp <- bp[!is.na(bp$gene),] # proteins annotated to a scaffold

## filter too many mismatches
bp <- bp[bp$mismatches<3,]

## SORT & FILTER
## sort by
## 1) identity,
## 2) -mismatches (perhaps redundant?),
## 3) whether protein is the annotated MANE, and
## 4) number of annotations.
## 5) remove duplicated after sorting!
bp <- bp[order(bp$identity,
               -bp$mismatches,
               !is.na(bp$MANE.protein),
               bp$protein==bp$MANE.protein,
               bp$numGO, decreasing=TRUE),]

## remove duplicated, which should appear below better hits!
bp <- bp[!duplicated(paste0(bp$BP,bp$gene)),]

## first try to find best match MANE, then handle missing
bad <- which(bp$mismatches >0 | bp$MANE=="" | bp$identity<100)

god <- bp[-bad,] # MANE, 100% identity and 0 mismatches
bad <- bp[ bad,]
bad <- bad[!bad$BP%in%god$BP,]

## of multiple that survived best filter,
## take the one (i) with most GO annotations,
## (ii) the longest, and if still multiple are remaining,
## take (iii) the first

god <- god[order(god$numGO, decreasing=TRUE),]
god <- god[!duplicated(god$BP),]

## missing: sort by identity, then MAN, then numGO
bad <- bad[order(bad$identity, bad$MANE, -bad$numGO, decreasing=TRUE),]
bad <- bad[!duplicated(bad$BP),]

## rest:
rest <- all[!all$BP%in%c(god$BP,bad$BP),]
rest <- rest[order(rest$identity, rest$MANE, -rest$numGO, decreasing=TRUE),]
rest <- rest[!duplicated(rest$BP),]

best <- rbind(cbind(god,match="good"),
              cbind(bad, match="bad"),
              cbind(rest, match="wrong"))

## replace tagged proteins with their original long name
##lids <- read.delim(bp.ids, header=FALSE, row.names=2, sep=" ")
##haveid <- which(best$protein%in%rownames(lids))
##best$protein[haveid] <- lids[best$protein[haveid],1]

## TAG CATEGORIES
type <- features$type[gidx[best$ensembl]]
best$IG <- FALSE
best$IG[grep("^IG",type)] <- TRUE

desc <- features$description[gidx[best$ensembl]]
best$albumin <- FALSE
best$albumin[grep("albumin",desc)] <- TRUE
best$globin <- FALSE
best$globin[grep("globin",desc)] <- TRUE

## load GO terms
## extracellular region: GO:0005576
## |- extracellular space: GO:0005615
gots <- features$GO[gidx[best$ensembl]]
best$extracellular <- FALSE
best$extracellular[grep("GO:0005576", gots)] <- TRUE

## exclude tag  - TODO: any K/R                 
best$exclude <- best$globin|best$albumin|best$IG


write.table(best, file=out.file, sep="\t",
            quote=FALSE, na="", row.names=FALSE)


