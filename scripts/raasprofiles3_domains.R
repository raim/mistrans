
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



## only use unique BP/SAAP per Dataset
do.unique <- FALSE # TRUE # 
dfig.path <- file.path(fig.path,"domains")
dir.create(dfig.path, showWarnings=FALSE)

tmtu <- tmtf
value <- "RAAS"

## unique BP/SAAP/Dataset
if ( do.unique ) {
    ## REDUCE TMT SET TO UNIQUE BP/SAAP/Dataset
    ## calculated in _init.R
    tmtu <- tmtf[!duplicated(tmtf$unique),]
    dfig.path <- file.path(dfig.path,"unique")
    dir.create(dfig.path, showWarnings=FALSE)
    value <- "RAAS.median"
} 

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
tmtu$BP.SAAP <- paste0(tmtu$BP,".",tmtu$SAAP)
## add "none" class for AAS outside pfam/clans
tmtu$clan.ebi[tmtu$clan.ebi==""] <- "No_pfam" 
tmtu$pfam.ebi[tmtu$pfam.ebi==""] <- "No_pfam" 
tmtu$clan[tmtu$clan==""] <- "No_pfam" 
tmtu$pfam[tmtu$pfam==""] <- "No_pfam" 

### GENERATE MEDIAN RAAS TABLES



## PFAM CLANS, EBI

plst <- strsplit(tmtu$clan.ebi,";")
dnms <- unique(unlist(plst))

## generate boolean matrix of pfam associations
pmat <- matrix(FALSE, ncol=length(dnms), nrow=nrow(tmtu))
colnames(pmat) <- dnms
for ( i in 1:length(plst) ) 
    pmat[i,plst[[i]]] <- TRUE

## replace column by name
dnms <- clans[sub("\\..*","",colnames(pmat)),2]
dnms[is.na(dnms)] <- colnames(pmat)[is.na(dnms)]
colnames(pmat) <- dnms

cle.ovl <- raasProfile(x=tmtu, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value=value, 
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

plst <- strsplit(tmtu$pfam.ebi,";")
dnms <- unique(unlist(plst))

## generate boolean matrix of pfam associations
pmat <- matrix(FALSE, ncol=length(dnms), nrow=nrow(tmtu))
colnames(pmat) <- dnms
for ( i in 1:length(plst) ) 
    pmat[i,plst[[i]]] <- TRUE

## replace column by name
dnms <- pfams[sub("\\..*","",colnames(pmat)),3]
dnms[is.na(dnms)] <- colnames(pmat)[is.na(dnms)]
colnames(pmat) <- dnms

pfe.ovl <- raasProfile(x=tmtu, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value=value, 
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

plst <- strsplit(tmtu$clan,";")
dnms <- unique(unlist(plst))

## generate boolean matrix of pfam associations
pmat <- matrix(FALSE, ncol=length(dnms), nrow=nrow(tmtu))
colnames(pmat) <- dnms
for ( i in 1:length(plst) ) 
    pmat[i,plst[[i]]] <- TRUE

## replace column by name
dnms <- clans[sub("\\..*","",colnames(pmat)),2]
dnms[is.na(dnms)] <- colnames(pmat)[is.na(dnms)]
colnames(pmat) <- dnms

cl.ovl <- raasProfile(x=tmtu, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value=value, 
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

plst <- strsplit(tmtu$pfam,";")
dnms <- unique(unlist(plst))

## generate boolean matrix of pfam associations
pmat <- matrix(FALSE, ncol=length(dnms), nrow=nrow(tmtu))
colnames(pmat) <- dnms
for ( i in 1:length(plst) ) 
    pmat[i,plst[[i]]] <- TRUE

## replace column by name
dnms <- pfams[sub("\\..*","",colnames(pmat)),3]
dnms[is.na(dnms)] <- colnames(pmat)[is.na(dnms)]
colnames(pmat) <- dnms

pf.ovl <- raasProfile(x=tmtu, id="BP.SAAP", 
                   rows=pmat, cols="Dataset",
                   bg=TRUE, value=value, 
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

pr.ovl <- raasProfile(x=tmtu, id="SAAP", 
                   rows="name", cols="Dataset",
                   col.srt=uds,
                    bg=TRUE, value=value, 
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

go.ovl <- raasProfile(x=tmtu, id="BP.SAAP", 
                    rows=got[tmtu$gene,], cols="Dataset",
                    bg=TRUE, value=value, 
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
    if ( nrow(ovc$p.value)>0 ) {        
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
    }
    
    
    ## p.min, both & frequent
    ovc <- sortOverlaps(ovf, p.min=p.min, cut=TRUE)
    fname <- file.path(dfig.path,paste0(did,SETID,"_pmin_frequent"))
    ## .. and FREQUENTLY MEASURED
    frequent <- names(which(ovc$num.query[,1]> fmin))
    if ( nrow(ovc$p.value)>0 ) {        
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
    }

    ## p.min, both 
    ovc <- sortOverlaps(ovf, p.min=p.min, cut=TRUE)
    if ( nrow(ovc$p.value)>0 ) {        
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
    }
    
    ## p.min, high RAAS
    ovc <- sortOverlaps(ovf, p.min=p.min, cut=TRUE, sign=1)
    if ( nrow(ovc$p.value)>0 ) {        

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
    }

    ## p.tgt, both 
    ovc <- sortOverlaps(ovf, p.min=p.tgt, cut=TRUE)
    if ( nrow(ovc$p.value)>0 ) {        
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
    }
    
    ## p.tgt, high RAAS
    ovc <- sortOverlaps(ovf, p.min=p.tgt, cut=TRUE, sign=1)
    if ( nrow(ovc$p.value)>0 ) {        
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
    }
  
    
    ## MOST RESTRICTIVE

    ## ptgt, high RAAS
    ovc <- sortOverlaps(ovf, p.min=p.tgt, cut=TRUE, sign=1)
    fname <- file.path(dfig.path,paste0(did,SETID,"_frequent"))
    ## .. and FREQUENTLY MEASURED
    frequent <- names(which(ovc$num.query[,1]> fmin))
    ovc <- sortOverlaps(ovc, srt=frequent, cut=TRUE)
    if ( nrow(ovc$p.value)>0 ) {        
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
    if ( nrow(ovc$p.value)>0 ) {
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
}


## TODO: use tightest GOslim result and dissect into proteins,
## do flowCharts?

ovc <- sortOverlaps(go.ovl, p.min=p.tight, cut=TRUE, sign=1)
## .. and FREQUENTLY MEASURED
all <- rownames(ovc$count)[apply(ovc$count, 1, function(x) all(x>0))]
ovc <- sortOverlaps(ovc, srt=all, cut=TRUE)
ovc <- sortOverlaps(ovc, p.min=p.txt, cut=FALSE)#, sign=1)

## count unique proteins for each enriched GO
gol <- list()
for ( i in 1:nrow(ovc$p.value) ) 
    gol[[i]] <- rownames(got)[which(got[,rownames(ovc$p.value)[i]])]
## filter those in TMT set
gol <- lapply(gol, function(x) unique(x[x%in%tmtu$gene]))
names(gol) <- rownames(ovc$p.value)
## protein names
gon <- lapply(gol, function(x) tmtu$name[tmtu$gene%in%x])
## BP
gob <- lapply(gol, function(x) unique(tmtu$BP.SAAP[tmtu$gene%in%x]))

## PLOT GO WITH ADDITIONAL BP and protein count columns
lmai <- omai
lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
lmai[4] <- 1.5
nw <- ncol(ovc$p.value)*.2 + lmai[2] + lmai[4]
nh <- nrow(ovc$p.value)*.2 + lmai[1] + lmai[3]
fname <- file.path(dfig.path,paste0("type_go_",SETID,"_all_manual_dotplot"))

plotdev(fname,  height=nh, width=nw, res=300, type=ftyp)
par(mai=lmai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
dotprofile(ovc, value="median", xlab=NA, ylab=NA,
           dot.sze=dot.sze, p.dot=p.dot,
           vcols=acols, vbrks=abrks, axis=1:2, show.total=TRUE)
axis(4, at=length(gol):1, labels=format(lengths(gol), big.mark=",", trim=TRUE),
     line=5, las=2, cex.axis=.8, col=1)
axis(4, at=length(gob):1, labels=format(lengths(gob), big.mark=",", trim=TRUE),
     line=2.5, las=2, cex.axis=.8, col=1)
text(x=c(9.5,11.5, 14), y=0.25, srt=90,
     labels=c("#RAAS","#BP/SAAP","#proteins"), xpd=TRUE, pos=2, cex=.8)
#text(x=11, length(gol):1, labels=lengths(gol),xpd=TRUE, cex=.8, pos=4)
dev.off()



## PLOT GENES FROM GO ENRICHMENT CATEGORIES 
prc.ovl <- pr.ovl #sortOverlaps(pr.ovl, axis=1, srt=uds[uds!="Healthy"])
for ( i in 1:length(gon) ) {

    goid <- names(gon)[i]

    ## get list of proteins that contribute to GO class
    gopr <- sortOverlaps(prc.ovl, axis=2, srt=unique(unlist(gon[[i]])))
    ## filter with a very mild p-value cutoff, at high RAAS
    gopr <- sortOverlaps(gopr, p.min=1e-1, sign=1, cut=TRUE)
    ## sort and cut by number of RAAS values,
    nsrt <- names(sort(gopr$num.query[gopr$num.query[,1]>1,1]))
    gopr <- sortOverlaps(pr.ovl, axis=2, srt=nsrt)
    
    fname <- file.path(dfig.path,paste0("goslim_",SETID,"_",
                                        gsub(" ","_",goid)))
    omai <- c(.8,1.1,.5,.5)
    nw <- ncol(gopr$p.value)*.2 + omai[2] + omai[4]
    nh <- nrow(gopr$p.value)*.2 + omai[1] + omai[3]

    plotdev(fname,  height=nh, width=nw, res=300, type=ftyp)
    par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="sans")
    dotprofile(gopr, value="median", vcols=acols, vbrks=abrks, p.dot=p.dot,
               axis=1:2,  show.total=TRUE, xlab=NA, ylab=NA)
    figlabel(goid, pos="bottomleft", font=2, cex=1)
    dev.off()
}

## TODO: for each dot, calculate % of signal that comes from a single BP

if ( interactive() ) {
    clusterFlow(hdat[,c("iupred3.bins","MMSeq2.bins")],
                srt=list(levels(iupred3.bins),
                         rev(levels(MMSeq2.bins))))



    cr <- plotCor(lengths(gol), ovc$num.query[names(gol),1], density=FALSE)#,
    ##              xlim=c(1e2,8e3), ylim=c(1e2,8e3), legpos="bottomright")
    abline(a=0, b=1)
    legend("right", paste0("slope=", round(cr$tls$beta,1)))
    unique(tmtf$name[tmtf$gene%in%gol[["protein catabolic process"]]])
}
## TODO: flow chart from ovc$p.value <-> genes
## TODO: dissect high RAAS go to high RAAS proteins
