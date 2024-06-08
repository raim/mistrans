
library(Biostrings) # for blosum62
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

dfig.path <- file.path(fig.path,"domains")
dir.create(dfig.path, showWarnings=FALSE)

p.show <- p.txt #1e-3

## load additional data
## pfam clans and readable names
clans.file <- file.path(mam.path,"originalData", "pfam", "Pfam-A.clans.tsv.gz")
pfams <- read.delim(clans.file, row.names=1 , header=FALSE)
clans <- pfams[,1:2]
clans <- clans[!duplicated(clans[,1]),]
clans <- clans[clans[,1]!="",]
rownames(clans) <- clans[,1]


## tag unique - TODO: in init?
tmtf$BP.SAAP <- paste0(tmtf$BP,".",tmtf$SAAP)

## add "none" class for AAS outside pfam/clans
tmtf$clan.ebi[tmtf$clan.ebi==""] <- "No_pfam" 
tmtf$pfam.ebi[tmtf$pfam.ebi==""] <- "No_pfam" 
tmtf$clan[tmtf$clan==""] <- "No_pfam" 
tmtf$pfam[tmtf$pfam==""] <- "No_pfam" 


## TODO: replace by clan/pfam names

## PFAM CLANS, EBI

plst <- strsplit(tmtf$clan.ebi,";")
dnms <- unique(unlist(plst))

## generate boolean matrix of pfam associations
pmat <- matrix(FALSE, ncol=length(dnms), nrow=nrow(tmtf))
colnames(pmat) <- dnms
for ( i in 1:length(plst) ) 
    pmat[i,plst[[i]]] <- TRUE

## replace column by name
dnms <- clans[sub("\\..*","",colnames(pmat)),2]
dnms[is.na(dnms)] <- colnames(pmat)[is.na(dnms)]
colnames(pmat) <- dnms

ovd <- raasProfile(x=tmtf, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value="RAAS", 
                   col.srt=uds,
                   use.test=use.test, 
                   xlab=xl.raas, verb=0)

## cut significant
ovc <- sortOverlaps(ovd, p.min=p.txt, cut=TRUE)
omai <- c(.8,1.5,.6,.6)
omai[2] <- .1*max(nchar(rownames(ovc$p.value)))

plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("clans_ebi_",SETID,"")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             llab="InterPro", rlab="CLAN",
             ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)#, plot.all=TRUE)

## PFAMS, EBI

plst <- strsplit(tmtf$pfam.ebi,";")
dnms <- unique(unlist(plst))

## generate boolean matrix of pfam associations
pmat <- matrix(FALSE, ncol=length(dnms), nrow=nrow(tmtf))
colnames(pmat) <- dnms
for ( i in 1:length(plst) ) 
    pmat[i,plst[[i]]] <- TRUE

## replace column by name
dnms <- pfams[sub("\\..*","",colnames(pmat)),3]
dnms[is.na(dnms)] <- colnames(pmat)[is.na(dnms)]
colnames(pmat) <- dnms

ovd <- raasProfile(x=tmtf, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value="RAAS", 
                   col.srt=uds,
                   use.test=use.test, 
                   xlab=xl.raas, verb=0)

ovc <- sortOverlaps(ovd, p.min=p.txt, cut=TRUE)
omai <- c(.8,1.5,.6,.6)
omai[2] <- .1*max(nchar(rownames(ovc$p.value)))

plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("pfams_ebi_",SETID,"")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             llab="InterPro", rlab="PFAM",
             ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)#, plot.all=TRUE)



## PFAM CLANS, own

plst <- strsplit(tmtf$clan,";")
dnms <- unique(unlist(plst))

## generate boolean matrix of pfam associations
pmat <- matrix(FALSE, ncol=length(dnms), nrow=nrow(tmtf))
colnames(pmat) <- dnms
for ( i in 1:length(plst) ) 
    pmat[i,plst[[i]]] <- TRUE

## replace column by name
dnms <- clans[sub("\\..*","",colnames(pmat)),2]
dnms[is.na(dnms)] <- colnames(pmat)[is.na(dnms)]
colnames(pmat) <- dnms

ovd <- raasProfile(x=tmtf, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value="RAAS", 
                   col.srt=uds,
                   use.test=use.test, 
                   xlab=xl.raas, verb=0)

ovc <- sortOverlaps(ovd, p.min=p.txt, cut=TRUE) #, sign=1)
omai <- c(.8,1.5,.6,.6)
omai[2] <- .1*max(nchar(rownames(ovc$p.value)))

plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("clans_",SETID,"")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             llab="hmmer", rlab="CLAN", ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)#, plot.all=TRUE)

## only high RAAS
ovc <- sortOverlaps(ovd, p.min=p.txt, cut=TRUE, sign=1)
omai <- c(.8,1.5,.6,.6)
omai[2] <- .1*max(nchar(rownames(ovc$p.value)))

plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("clans_",SETID,"_highRAAS")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             llab="hmmer", rlab="CLAN", ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, plot.all=FALSE)

## filter and sort
ovc <- sortOverlaps(ovd, p.min=p.dot, cut=TRUE)#, sign=1)
frequent <- names(which(ovc$num.query[,1]>10))
ovc <- sortOverlaps(ovc, srt=frequent, cut=TRUE)
ovc <- sortOverlaps(ovc, p.min=p.txt, cut=FALSE)#, sign=1)

omai <- c(.05,1.5,.5,.5)
##CMAIL <- 1.54 ## commonly used between motif and domain figures, defined hered
##              ## to fit x-axis label
omai[2] <- CMAIL # .11*max(nchar(rownames(ovc$p.value)))

plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("clans_",SETID,"_frequent")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             ##llab="hmmer", rlab="CLAN",
             ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             bg=NA, tot.cex=.8,
             gcols=gcols, plot.all=TRUE, ffam=FONT)

## PFAMS, own

plst <- strsplit(tmtf$pfam,";")
dnms <- unique(unlist(plst))

## generate boolean matrix of pfam associations
pmat <- matrix(FALSE, ncol=length(dnms), nrow=nrow(tmtf))
colnames(pmat) <- dnms
for ( i in 1:length(plst) ) 
    pmat[i,plst[[i]]] <- TRUE

## replace column by name
dnms <- pfams[sub("\\..*","",colnames(pmat)),3]
dnms[is.na(dnms)] <- colnames(pmat)[is.na(dnms)]
colnames(pmat) <- dnms

ovd <- raasProfile(x=tmtf, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value="RAAS", 
                   col.srt=uds,
                   use.test=use.test, 
                   xlab=xl.raas, verb=0)

ovc <- sortOverlaps(ovd, p.min=p.txt, cut=TRUE)
omai <- c(.8,1.5,.6,.6)
omai[2] <- .1*max(nchar(rownames(ovc$p.value)))

plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("pfams_",SETID,"")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             llab="hmmer", rlab="PFAM", ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)#, plot.all=TRUE)


## TODO: repeat for unique SAAP/BP!
