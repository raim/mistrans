
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

## additional data
goslim.file  <- file.path(mam.path,"processedData","goslim.tsv")


dfig.path <- file.path(fig.path,"domains")
dir.create(dfig.path, showWarnings=FALSE)

p.show <- p.txt #1e-3

## TIGHT PLOT FOR MAIN
p.tgt <- 1e-15 # tighter cutoff for GO analysis
p.go <- 1e-15 # tighter cutoff for GO analysis

fmin <- 20 # minimal number of available RAAS measurements

## load additional data
## pfam clans and readable names
clans.file <- file.path(mam.path,"originalData", "pfam", "Pfam-A.clans.tsv.gz")
pfams <- read.delim(clans.file, row.names=1 , header=FALSE)
clans <- pfams[,1:2]
clans <- clans[!duplicated(clans[,1]),]
clans <- clans[clans[,1]!="",]
rownames(clans) <- clans[,1]


## TODO: do this in and merge with raasprofiles3_init.R

## tag unique 
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

cle.ovl <- raasProfile(x=tmtf, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value="RAAS", 
                   col.srt=uds,
                   use.test=use.test, 
                   xlab=xl.raas, verb=0)

## cut significant
ovc <- sortOverlaps(cle.ovl, p.min=p.txt, cut=TRUE)
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

pfe.ovl <- raasProfile(x=tmtf, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value="RAAS", 
                   col.srt=uds,
                   use.test=use.test, 
                   xlab=xl.raas, verb=0)

ovc <- sortOverlaps(pfe.ovl, p.min=p.txt, cut=TRUE)
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

cl.ovl <- raasProfile(x=tmtf, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value="RAAS", 
                   col.srt=uds,
                   use.test=use.test, 
                   xlab=xl.raas, verb=0)

ovc <- sortOverlaps(cl.ovl, p.min=p.txt, cut=TRUE) #, sign=1)
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
ovc <- sortOverlaps(cl.ovl, p.min=p.txt, cut=TRUE, sign=1)
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
ovc <- sortOverlaps(cl.ovl, p.min=p.dot, cut=TRUE)#, sign=1)
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

pf.ovl <- raasProfile(x=tmtf, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value="RAAS", 
                   col.srt=uds,
                   use.test=use.test, 
                   xlab=xl.raas, verb=0)

ovc <- sortOverlaps(pf.ovl, p.min=p.txt, cut=TRUE)
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

ovc <- sortOverlaps(pf.ovl, p.min=p.dot, cut=TRUE)#, sign=1)
frequent <- names(which(ovc$num.query[,1]>10))
ovc <- sortOverlaps(ovc, srt=frequent, cut=TRUE)
ovc <- sortOverlaps(ovc, p.min=p.txt, cut=FALSE)#, sign=1)

omai <- c(.05,1.5,.5,.5)
##CMAIL <- 1.54 ## commonly used between motif and domain figures, defined hered
##              ## to fit x-axis label
omai[2] <- CMAIL # .11*max(nchar(rownames(ovc$p.value)))

plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("pfams_",SETID,"_frequent")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             ##llab="hmmer", rlab="CLAN",
             ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             bg=NA, tot.cex=.8,
             gcols=gcols, plot.all=TRUE, ffam=FONT)


### PROTEIN RAAS PROFILE

pr.ovl <- raasProfile(x=tmtf, id="SAAP", 
                   rows="name", cols="Dataset",
                   col.srt=uds,
                    bg=TRUE, value="RAAS", 
                   use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                   verb=0)

ovc <- sortOverlaps(pr.ovl, axis=2, p.min=p.min, cut=TRUE)
omai <- c(.8,1,.6,.6)
plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("proteins_",SETID,"")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam=FONT)#, plot.all=TRUE)

## TODO: repeat for unique SAAP/BP!


### GO RAAS Profiles

## get GOslim table for use with clusterAnnotation
got <- parseAnnotationList(genes[,c("ID","GOslim")]) 
## replace GO IDs by terms
terms <- read.delim(goslim.file)
trms <- terms[,2]
names(trms) <- terms[,1]
colnames(got) <- trms[colnames(got)]

## shorten names
colnames(got) <- sub("plasma membrane", "PM",
                     sub("binding", "bnd.",
                         sub("templated", "templ.",
                             sub("regulation", "reg.",
                                 sub("transcription", "transcr.",
                                     sub("localization","local.",
                                         colnames(got)))))))

## REDUCE AGAIN to use only available proteins with RAAS
## TODO: define proper background set!
local <- FALSE
if ( local ) {
    ##gen.cls <- gen.cls[gidx,]
    ##got <- got[gidx,]
}

go.ovl <- raasProfile(x=tmtf, id="SAAP", 
                    rows=got[tmtf$gene,], cols="Dataset",
                    bg=TRUE, value="RAAS", 
                    col.srt=uds,
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=0)

ovc <- sortOverlaps(go.ovl, axis=2, p.min=p.min, sign=1, cut=TRUE)
omai <- c(.8,3.5,.6,.6)
plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("go_",SETID,"")),
             mai=omai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             p.min=p.min, p.txt=p.txt,
             ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="sans")#, plot.all=TRUE)


## collect all, and do common plot
## pf.ovl, cl.ovl, pr.ovl, go.ovl.

## cut at minimal p-value
cids <- c("pfe","cle","pf","cl","pr","go")
names(cids) <- cids
cnms <- cids
cnms["cle"] <- "CLAN/EBI"
cnms["pfe"] <- "PFAM/EBI"
cnms["cl"] <- "CLAN"
cnms["pf"] <- "PFAM"
cnms["pr"] <- "protein"
cnms["go"] <- "GOslim"

## TODO: why 102 in PSMA1 and in Proteasome PFAM,
## although there are others PSMA4/5/6/... with the same fold,
## by accident with the Q:G outside this domain?

source("~/programs/segmenTools/R/clusterTools.R")
omai <- c(.8,CMAIL,.5,.5)
for ( cid in cids ) {

    ovf <- get(paste0(cid, ".ovl"))
    did <- paste0("type_", cid, "_")

    ## p.txt, both 
    ovc <- sortOverlaps(ovf, p.min=p.txt, cut=TRUE)
    fname <- file.path(dfig.path,paste0(did,SETID,"_ptxt"))

    lmai <- omai
    lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
    plotProfiles(ovc, fname=fname,
                 mai=lmai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=p.dot,
                 p.min=p.min, p.txt=p.txt,
                 rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT, bg=NA)#, plot.all=TRUE)

    
    ## p.min, both & frequent
    ovc <- sortOverlaps(ovf, p.min=p.min, cut=TRUE)
    fname <- file.path(dfig.path,paste0(did,SETID,"_pmin_frequent"))
    ## .. and FREQUENTLY MEASURED
    frequent <- names(which(ovc$num.query[,1]> fmin))
    ovc <- sortOverlaps(ovc, srt=frequent, cut=TRUE)
      
    lmai <- omai
    lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
    plotProfiles(ovc, fname=fname,
                 mai=lmai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=p.dot,
                 p.min=p.min, p.txt=p.txt,
                 rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT, bg=NA)#, plot.all=TRUE)

    ## p.min, both 
    ovc <- sortOverlaps(ovf, p.min=p.min, cut=TRUE)
    fname <- file.path(dfig.path,paste0(did,SETID,"_pmin"))
    
    lmai <- omai
    lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
    plotProfiles(ovc, fname=fname,
                 mai=lmai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=p.dot,
                 p.min=p.min, p.txt=p.txt,
                 rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT, bg=NA)#, plot.all=TRUE)

    ## p.min, high RAAS
    ovc <- sortOverlaps(ovf, p.min=p.min, cut=TRUE, sign=1)
    fname <- file.path(dfig.path,paste0(did,SETID,"_pmin_high"))
    
    lmai <- omai
    lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
    plotProfiles(ovc, fname=fname,
                 mai=lmai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=p.dot,
                 p.min=p.min, p.txt=p.txt,
                 rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT, bg=NA)#, plot.all=TRUE)

    ## p.tgt, both 
    ovc <- sortOverlaps(ovf, p.min=p.tgt, cut=TRUE)
    fname <- file.path(dfig.path,paste0(did,SETID,"_ptgt"))
    
    lmai <- omai
    lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
    plotProfiles(ovc, fname=fname,
                 mai=lmai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=p.dot,
                 p.min=p.min, p.txt=p.txt,
                 rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT, bg=NA)#, plot.all=TRUE)

    ## p.tgt, high RAAS
    ovc <- sortOverlaps(ovf, p.min=p.tgt, cut=TRUE, sign=1)
    fname <- file.path(dfig.path,paste0(did,SETID,"_ptgt_high"))
    
    lmai <- omai
    lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
    plotProfiles(ovc, fname=fname,
                 mai=lmai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=p.dot,
                 p.min=p.min, p.txt=p.txt,
                 rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT, bg=NA)#, plot.all=TRUE)
    
  
    
    ## MOST RESTRICTIVE

    ## ptgt, high RAAS
    ovc <- sortOverlaps(ovf, p.min=p.tgt, cut=TRUE, sign=1)
    fname <- file.path(dfig.path,paste0(did,SETID,"_frequent"))
    ## .. and FREQUENTLY MEASURED
    frequent <- names(which(ovc$num.query[,1]> fmin))
    ovc <- sortOverlaps(ovc, srt=frequent, cut=TRUE)
    ## and present in all data sets
    ## sort rest by lower cutoff p-value
    ovc <- sortOverlaps(ovc, p.min=p.txt, cut=FALSE)#, sign=1)
        
    lmai <- omai
    lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
    plotProfiles(ovc, fname=fname,
                 mai=lmai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=p.dot,
                 p.min=p.min, p.txt=p.txt,
                 rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT)#, plot.all=TRUE)

    ## TIGHTEST FOR MAIN: ptgt, high RAAS

    p.tight <- p.min
    if ( cid=="go" )
        p.tight <- p.go
    
    ovc <- sortOverlaps(ovf, p.min=p.tight, cut=TRUE, sign=1)
    fname <- file.path(dfig.path,paste0(did,SETID,"_all"))
    ## .. and FREQUENTLY MEASURED
    all <- rownames(ovc$count)[apply(ovc$count, 1, function(x) all(x>0))]
    ovc <- sortOverlaps(ovc, srt=all, cut=TRUE)
    ## and present in all data sets
    ## sort rest by lower cutoff p-value
    ovc <- sortOverlaps(ovc, p.min=p.txt, cut=FALSE)#, sign=1)
        
    lmai <- omai
    lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
    plotProfiles(ovc, fname=fname,
                 mai=lmai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=p.dot,
                 p.min=p.min, p.txt=p.txt,
                 rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT)#, plot.all=TRUE)
    
}
