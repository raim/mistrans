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
colnames(saap) <- c("SAAP","protein","identity", "mismatches", "gaps",
                    "alen", "qlen", "qstart", "qend",
                    "slen", "sstart", "send", "e", "bitscore")

## simply count perfect matches
## alen>=qlen-1, gaps=0, identity=100
best <- saap[saap$alen >= saap$qlen-1 &
             saap$gaps == 0 ,]
## sort by mismatches to get best hit, then remove duplicates
best <- best[order(best$mismatches),]
tid <- best$SAAP[which(duplicated(best$SAAP))[2]]
best[best$SAAP==tid,]

## perfect matches: full length SAAP match (qstart=1, qlen=qend,
## no gaps, mismatches and 100% identity)
perfect <- saap[((saap$qstart==1 & saap$qend==saap$qlen)
    & saap$gaps==0) &
    (saap$mismatches==0 & saap$identity==100),]

## sort by SAAP
perfect <- perfect[order(perfect$SAAP),]


cat(paste("finding perfect matches for", sum(!duplicated(perfect$SAAP)),
          "SAAPs in the full set, with substitutions\n"))

## exlude those with a substitution tag
pcore <- perfect[grep("_", perfect$protein, invert=TRUE),]

cat(paste("finding perfect matches for", sum(!duplicated(pcore$SAAP)),
          "SAAPs in the ensembl proteome\n"))

## length distribution

brks <- range(saap$qlen)
brks <- brks[1]:brks[2]
hist(saap$qlen[!duplicated(saap$SAAP)], breaks=brks,
     xlab="SAAP length", main=NA)
hist(perfect$qlen[!duplicated(perfect$SAAP)], add=TRUE, col=2,
     breaks=brks)
legend("topright", c("all SAAP", "perfect match"), col=c(8,2),
       lty=1, lwd=5)

## plot RAAS
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")

tmtf <- read.delim(tmt.file)

pkeep <- perfect[perfect$SAAP%in%tmtf$SAAP[tmtf$Keep.SAAP],]
pkore <- pkeep[grep("_", pkeep$protein, invert=TRUE),]

cat(paste("finding perfect matches for", sum(!duplicated(pkeep$SAAP)),
          "SAAPs that are not tagged to be removed\n"))
cat(paste("finding perfect matches for", sum(!duplicated(pkore$SAAP)),
          "SAAPs that are not tagged to be removed and have no mutation tage\n"))

## TESTING SOME CASES
## BP:   LGFYGLDESDLDKVFHLPTTTFIGGQESALPLR, matches ENSP00000222673 (MANE)
## SAAP: GGFYGLDESDLDKVFHLPTTTFIGGQESALPLR, matches ENSP00000388183,ENSP00000414662
## both SAAP hits are proteoforms

## BP:   DTEEEDFHVDQVTTVK, ENSP00000376802, SERPINA1
## SAAP: DTEEEDFHVDQATTVK, ENSP00000486067, SERPINA1 description:serpin family A member 1, but annotated as scaffold!

## BP:   QGVEDAFYTLVR ENSP00000309845, HRAS, P01112
## SAAP: QGVDDAFYTLVR ENSP00000308495, KRAS, P01116

hist(tmtf$RAAS[tmtf$Keep.SAAP])
hist(tmtf$RAAS[tmtf$Keep.SAAP & tmtf$SAAP %in% unique(perfect$SAAP)],
     col=2, add=TRUE)

nsaap <- sum(unique(tmtf$SAAP[tmtf$Keep.SAAP])%in%unique(perfect$SAAP))
nraas <- length(tmtf$RAAS[tmtf$Keep.SAAP & tmtf$SAAP %in% unique(perfect$SAAP)])


plotdev(file.path(fig.path,"saap_blast_raas"),
        res=300, width=7, height=3.5, type=ftyp)
par(mfcol=c(1,1), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(tmtf$RAAS[tmtf$Keep.SAAP], freq=FALSE, ylim=c(0,.4),
     xlab=expression(log[10](RAAS)), main=NA)
hist(tmtf$RAAS[tmtf$Keep.SAAP & tmtf$SAAP %in% unique(perfect$SAAP)],
     col="#ff000077", border=2, add=TRUE, freq=FALSE)
legend("topright", c("all SAAP",
                     paste0("perfect match: ", nraas)),
       col=c(8,2), lty=1, lwd=5)
dev.off()

table(tmtf$Keep.SAAP[tmtf$SAAP %in% unique(perfect$SAAP)])
table(tmtf$AAS[tmtf$Keep.SAAP & tmtf$SAAP %in% unique(perfect$SAAP)])

## TODO: sort for best hits and remove duplicates as in get_protein_match.R
## TODO: fuse this into get_protein_match.R or map_peptides3.R
## SORT & FILTER
## sort by
## 1) identity,
## 2) -mismatches (perhaps redundant?),
## x) whether protein is the annotated MANE, and
## x) number of annotations.
## x) remove duplicated after sorting!
saap <- saap[order(saap$identity, saap$alen, 
                   -saap$gaps, -saap$mismatches, decreasing=TRUE),]

## check sorting for first entries
tid <- saap$SAAP[which(duplicated(saap$SAAP))[1]]

head(saap[saap$SAAP==tid,])

## filter all where alignment length equals query (SAAP) length,
## w/o gaps, then count mismatches



## remove duplicated, which should appear below best hits!
saap <- saap[!duplicated(paste0(saap$SAAP,saap$gene)),]

plotdev(file.path(fig.path,"saap_blast"),
        res=300, width=7, height=3.5, type=ftyp)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(saap$identity)
bp <- table(saap$mismatches)
barplot(bp, xlab="number of mismatches")
dev.off()

pm <- saap[which(saap$mismatches==0),]
