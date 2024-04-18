
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

## map proteins to feature table!
plst <- strsplit(features$proteins,";")
gidx <- rep(1:nrow(features), lengths(plst))
names(gidx) <- unlist(plst)

## genes and transcripts for each protein
bp$gene <- features$ID[gidx[bp$ensembl]]
bp$transcript <- tpmap[bp$ensembl,1]

## FUSE MANE and canonical annotation
## NOTE: NO CONFLICTS noted between MANE and canonical
## TODO: stats on how many MANE vs. canonical
## NA: scaffold proteins w/o match in genome feature table
## "": no MANE or canonical annotation
if ( FALSE ) {
    mains <- data.frame(MANE=features$MANE[gidx[bp$ensembl]],
                        canonical=features$canonical[gidx[bp$ensembl]])
    mains[mains==""] <- NA 
    tmp <- apply(mains,1, unique) # merge unique MANE/canonical
    tmp <- lapply(tmp, function(x) x[!is.na(x)]) # replace non-found by NAs
    tmp <- lapply(tmp, function(x) {if (length(x)==0) x<-NA;x})
    
    if ( any(lengths(tmp)!=1) )
    stop("check MANE/canonical mapping")
}
## ONLY USE MANE!!
tmp <- features$MANE[gidx[bp$ensembl]]
tmp[tmp==""] <- NA
bp$MANE <- features$MANE[gidx[bp$ensembl]] ##unlist(tmp)
bp$MANE.protein <- pampt[match(bp$MANE,rownames(pampt)),1]

bp$numGO <- gol[gidx[bp$ensembl]]


## filter - only those where we found a gene, others
## w/o gene are encoded on scaffolds, TODO: check all.
all <- bp ## store to analyze missing below

ina <- is.na(bp$gene)
cat(paste("removing", sum(ina), "blast hits without a mapped gene!\n"))
bp <- bp[!ina,] # proteins annotated to a scaffold

## filter too many mismatches
mm <- bp$mismatches>=3
cat(paste("removing", sum(mm), "blast hits with >= 3 mismatches\n"))
bp <- bp[!mm,]

## SORT & FILTER
## sort by
## 1) identity,
## 2) -mismatches (perhaps redundant?),
## 3) whether protein is the annotated MANE, and
##    prefer if the blasted protein IS the MANE,
## 4) protein length - TAKE SHORTEST,
## 5) number of annotations.
## -> remove duplicated after sorting!
bp <- bp[order(bp$identity,
               -bp$mismatches,
               !is.na(bp$MANE.protein),
               bp$protein==bp$MANE.protein,
               -bp$slen,
               bp$numGO, decreasing=TRUE),]

## .. and remove duplicated, which should appear below better hits!!

## TODO: instead use a "the winner takes it all" approach and
## and assign AAS to the protein with most overall matches.

## TODO: analyze this step a bit, e.g
## * we retain BCL2L2-PABN1, a read-through
## fusion protein with only 1 AAS; could be the single gene!
##bp[bp$BP=="SIYVGNVDYGATAEELEAHFHGCGSVNR",]
## * all proteasome results are lost  if Q:G in PSMA1 is mapped to
## the SHORTER
## PSMA1 hits lost by new preference for smaller proteins:
## AQSELAAHQK (Q->A and Q->G) and NQYDNDVTVWSPQGR (V->N), the former
## matches ENSP00000457299, a "novel protein" that does have a canonical
## transcript annotation, but no other information.

## REMOVE DUPLICATE BP
dupb <- duplicated(paste0(bp$BP))
cat(paste("removing", sum(dupb), "multiple matches\n"))
bp <- bp[!dupb,]

## TAG QUALITY
bad <- bp$mismatches >0 | bp$MANE=="" | bp$identity<100
bp$match <- ifelse(bad, "bad", "good")
bp <- bp[order(bp$match, decreasing=TRUE),]

## rest:
rest <- all[!all$BP%in%c(bp$BP),]
rest <- rest[order(rest$identity, rest$MANE, -rest$numGO, decreasing=TRUE),]
rest <- rest[!duplicated(rest$BP),]

## FUSE
best <- rbind(bp,
              cbind(rest, match="wrong"))

## replace tagged proteins with their original long name
##lids <- read.delim(bp.ids, header=FALSE, row.names=2, sep=" ")
##haveid <- which(best$protein%in%rownames(lids))
##best$protein[haveid] <- lids[best$protein[haveid],1]

### TAG CATEGORIES

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

## ADD GENE NAME
best$name <- features$name[gidx[best$ensembl]]

if ( !interactive() ) {
    write.table(best, file=out.file, sep="\t",
                quote=FALSE, na="", row.names=FALSE)
}

