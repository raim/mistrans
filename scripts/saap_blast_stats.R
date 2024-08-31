library(segmenTools)
options(stringsAsFactors=FALSE)

## TEST SAAP BLASTED AGAINST ALL ENSEMBL PROTEINS
## incl. patient-specific mutations

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

## remove SAAP that are tagged to be removed
pkeep <- perfect[perfect$SAAP%in%tmtf$SAAP[tmtf$Keep.SAAP],]
## subset of SAAP that have a mutation tag
pkore <- pkeep[grep("_", pkeep$protein, invert=TRUE),]

cat(paste("finding perfect matches for", sum(!duplicated(pkeep$SAAP)),
          "SAAPs that are not tagged to be removed\n"))
cat(paste("finding perfect matches for", sum(!duplicated(pkore$SAAP)),
          "SAAPs that are not tagged to be removed and have no mutation tag\n"))

## TESTING SOME CASES
## BP:   LGFYGLDESDLDKVFHLPTTTFIGGQESALPLR, matches ENSP00000222673 (MANE)
## SAAP: GGFYGLDESDLDKVFHLPTTTFIGGQESALPLR, matches ENSP00000388183,ENSP00000414662
## both SAAP hits are proteoforms, with N-term difference
##ENST00000222673 ENSP00000222673
##ENST00000447398 ENSP00000388183 cds_3   7       44647444        44647534        +
##ENST00000444676 ENSP00000414662 cds_3   7       44647657        44647759        +
## confirm only difference by diff these files:

##grep ENST00000222673 ~/data/mammary/originalData/transcript_coordinates.tsv |cut -f 2,3,4,5 > data/ENST00000222673.txt
##grep ENST00000447398 ~/data/mammary/originalData/transcript_coordinates.tsv |cut -f 2,3,4,5 > data/ENST00000447398.txt
##grep ENST00000444676 ~/data/mammary/originalData/transcript_coordinates.tsv |cut -f 2,3,4,5 > data/ENST00000444676.txt

##raim > diff data/ENST00000222673.txt data/ENST00000447398.txt 
##3c3,4
##< 7     44647657        44647759
##---
##> 7     44647444        44647534
##> 7     44656318        44656362
##
##raim > diff data/ENST00000222673.txt data/ENST00000444676.txt 
##3a4
##> 7     44656318        44656362
##raim > diff data/ENST00000447398.txt data/ENST00000444676.txt 
##3c3
##< 7     44647444        44647534
##---
##> 7     44647657        44647759



## BP:   DTEEEDFHVDQVTTVK, ENSP00000376802, SERPINA1
## SAAP: DTEEEDFHVDQATTVK, ENSP00000486067, SERPINA1 description:serpin family A member 1, but annotated as scaffold!

## BP:   QGVEDAFYTLVR ENSP00000309845, HRAS, P01112
## SAAP: QGVDDAFYTLVR ENSP00000308495, KRAS, P01116

## count SAAP
nsaap <- sum(unique(tmtf$SAAP[tmtf$Keep.SAAP])%in%unique(perfect$SAAP))
## NOTE: count contains Inf values
nraas <- length(tmtf$RAAS[tmtf$Keep.SAAP & tmtf$SAAP %in% unique(perfect$SAAP)])
asaap <- length(unique(tmtf$SAAP[tmtf$Keep.SAAP]))

allraas <- tmtf$RAAS[tmtf$Keep.SAAP]
mtcraas <- tmtf$RAAS[tmtf$Keep.SAAP & tmtf$SAAP %in% unique(perfect$SAAP)]

## handle infinite RAAS

mxraas <- max(allraas[is.finite(allraas)])
mnraas <- min(allraas[is.finite(allraas)])

allraas[allraas==Inf] <-mxraas +1
allraas[allraas==-Inf] <- mnraas -1

mtcraas[mtcraas==Inf] <- mxraas +1
mtcraas[mtcraas==-Inf] <- mnraas -1

brks <- -7:5

plotdev(file.path(fig.path,"saap_blast_raas_total"),
        res=300, width=4, height=3.5, type=ftyp)
par(mfcol=c(1,1), mai=c(.5,.5,.1,.5), mgp=c(1.3,.3,0), tcl=-.25)
hist(allraas, xlab=expression(log[10](RAAS)), main=NA, breaks=brks)
hist(mtcraas, col=2, add=TRUE, breaks=brks)
legend("topleft", c(paste0("all SAAP: ", asaap),
                     paste0("isoform match: ", nsaap)),
       col=c(8,2), pch=15, bty="n")
dev.off()


plotdev(file.path(fig.path,"saap_blast_raas"),
        res=300, width=4, height=3.5, type=ftyp)
par(mfcol=c(1,1), mai=c(.5,.5,.1,.5), mgp=c(1.3,.3,0), tcl=-.25)
hist(allraas, freq=FALSE, ylim=c(0,.4),
     col="#00000099", border=NA,
     xlab=expression(log[10](RAAS)), main=NA, breaks=brks)
par(new=TRUE)
hist(mtcraas, breaks=brks, ylim=c(0,110),
     col="#ff000077", border=2, xlab=NA, ylab=NA, main=NA, axes=FALSE)
mtext("Frequency", 4, 1.3, col=2)
axis(4, col=2, col.axis=2)
legend("topleft", c(paste0("all SAAP: ", asaap, ", n=",sum(!is.na(allraas))),
                     paste0("isoform match: ", nsaap, ", n=",nraas)),
       col=c("#00000099","#ff000077"), pch=15, bty="n")
dev.off()

## AAS TYPES
table(tmtf$AAS[tmtf$Keep.SAAP & tmtf$SAAP %in% unique(pkeep$SAAP)])


## TODO:
## count ensembl vs. ensembl+mutation match
## non-mutation examples
