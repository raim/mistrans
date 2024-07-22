
### ANALYZE SEQUENCE CONTEXT around AAS

## FOCUS ON CLEAVAGE SITE ENRICHMENT (the KRAQ)
## AA and AAS enrichment along peptides in AAS set and all peptides,
## as supplied by Shiri Tsour.

## TODO: repeat diAA positional enrichment from get_aas_context.R

library(stringr)
library(DiffLogo)
library(Biostrings)
library(segmenTools)
AAS <- sort(unique(GENETIC_CODE))
AAT <- AAS[AAS!="*"]
diAAT <- sort(paste(AAT, rep(AAT,each=length(AAT)),sep=""))

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")


sfig.path <- file.path(fig.path,"saap")
dir.create(sfig.path, showWarnings=FALSE)

### remove duplicate sites!
## TODO: is there a cleaner way to do this? here we randomly loose AAS types
cdat <- bdat[!duplicated(paste(bdat$ensembl,bdat$pos, bdat$fromto)),]


### WOULD SAAP MAP TO OTHER BP?
library(stringdist)
## rows are BP and columns are SAAP
mtx <- stringdistmatrix(unique(bdat$BP), unique(bdat$SAAP),
                        method = "hamming", useNames = "strings")
mtx <- mtx==1 # distance 1?
mtb <-  apply(mtx, 2, sum)

cat(paste(sum(mtb>1), "SAAP map to more than one BP\n"))

ids <- names(which.max(mtb))

ids <- names(which(mtb>1))[1]
bps <- names(which(mtx[,ids]))

cat(paste(paste(c(ids, bps), collapse="\n"),"\n"))

##15080 AEGPEVDVNLPK AEGPEVDVTLPK    PDAC
##17325 GEGPEVDVTLPK GEGPEVDVSLPK    PDAC

pnms[unique(bdat[bdat$BP%in%bps,"protein"])]

ids <- "ADQCYEDIR"
bpr <- "ADQCYEDVR" # actual paired BP
bpm <- "NDQCYEDIR" # missed BP with equal distance
bps <- names(which(mtx[,ids]))

tmtf[tmtf$BP%in%c(bps),c("BP","SAAP", "Dataset")]

## SAAPs that map to more than 1 BP
msaaps <- names(which(mtb>1))

mdat <- bdat[bdat$SAAP%in%msaaps,]
mlst <- split(mdat, mdat$SAAP)

ml <- mlst[[2]]
unique(ml$SAAP)
unique(ml$SAAP)

## TODO: investigate all:
## ## * mutated proteins? would blast have reported the mutated protein?
## * just look at Healthy,
## * ignore tissues and take a median intensity for each BP and SAAP over
##   all tissues they are observed in,
## * calculate RAAS for all possible (hamming distance==1) BP:SAAP pairs, 
##   make a scatter plot showing each SAAP that has >1 matching BP.

tmth <- tmtf[tmtf$Dataset=="Healthy",]

bpa <- split(tmth$BP.abundance, tmth$BP)
spa <- split(tmth$SAAP.abundance, tmth$SAAP)
bpa <- unlist(lapply(bpa,  median))
spa <- unlist(lapply(spa,  median))

## rows are BP and columns are SAAP
mtx <- stringdistmatrix(names(bpa), names(spa),
                        method = "hamming", useNames = "strings")
mtx <- mtx==1 # distance 1?
mtb <-  apply(mtx, 2, sum)

table(mtb)

spm <- names(which(mtb>1))


## matrix of SAAP/BP pairs
sbm <- sapply(spm, function(x) {
    bps <- names(which(mtx[,x]))
    cbind(SAAP=rep(x, length(bps)),
          BP=bps)})
sbm <- as.data.frame(do.call(rbind, sbm))

sbm$identified <- paste(sbm$BP, sbm$SAAP) %in% paste(tmtf$BP, tmtf$SAAP)
sbm$BP.abundance <- bpa[sbm$BP]
sbm$SAAP.abundance <- spa[sbm$SAAP]
sbm$RAAS <- log10(sbm$SAAP.abundance/sbm$BP.abundance)
sbl <- split(sbm, sbm$SAAP)

sbrl <- lapply(sbl, function(x) {
    n <- sum(!x$identified)
    m <- data.frame(x=rep(x$RAAS[x$identified], n),
                    y=x$RAAS[!x$identified],
                    SAAP=rep(x$SAAP[x$identified], n),
                    BPi=rep(x$BP[x$identified], n),
                    BPe=x$BP[!x$identified])
})
sbr <- do.call(rbind, sbrl)

plotdev(file.path(sfig.path,paste0("saap_multiple_bp")),
        height=2.5, width=2.5, res=300, type=ftyp)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.05)
plotCor(sbr$x, sbr$y, density=FALSE, xlab=expression(log[10](RAAS[identified])),
        ylab=expression(log[10](RAAS["alternative BP"])), line.methods="",
        title=TRUE, cor.legend=FALSE)
abline(a=0, b=1)
dev.off()
