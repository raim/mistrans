
### ANALYZE SEQUENCE CONTEXT around AAS

## 1: all -> just do hypergeo enrichment tests,
##    as in previous get_aas_context.R
## 2: analyze position of AAS in peptide.

## Motivated by above analysis, select subsets of BP/sequences,
## and do a DiffLogo and a RAAS profile.

## 1: M or W in -2/-1/+1/+2
## 2: higher than threshold E/D
## 3: AAS in peptide position 1-3.

## ... and finally loop over all of certain type:
## 1: from AA,
## 2: to AA,
## 3: from/to AA,
## 4: from/to by AA prop.


library(DiffLogo)
library(Biostrings)
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


mfig.path <- file.path(fig.path,"motifs")
dir.create(mfig.path, showWarnings=FALSE)

## NOTE: 1e-6 doesn't seem to be worth it, little gain.
## but further experience required
mp.min <- 1e-6 #p.min
mp.dot <- 1e-6 #p.dot
mp.txt <- 1e-3 #p.txt


#### GENERAL SEQUENCE LOGOS

### AA MATRIX
aam <- do.call(rbind, strsplit(bdat$AA,""))
nc <- (ncol(aam)-1)/2
colnames(aam) <- -nc:nc

## AA MATRIX WITH AAS
aas <- aam
aas[,"0"] <- bdat$to

## HYPERGEO TESTS - tigther context
ovl <- aaProfile(aam[,as.character(-7:7)], abc=AAT)
ovl <- sortOverlaps(ovl, p.min=mp.txt)
plotdev(file.path(mfig.path,paste0("AAS_overlap")),
        height=5, width=5, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
plotOverlaps(ovl, p.min=mp.min, p.txt=mp.txt, show.total=TRUE,
             xlab="Distance from AAS")
dev.off()

plotdev(file.path(mfig.path,paste0("AAS_dotplot")),
        height=5, width=5, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
dotprofile(ovl, value="ratio", vcols=ttcols,
           xlab="Distance from AAS",
           p.dot=mp.dot, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), ylab="amino acid", axis=1:2, show.total=TRUE)
dev.off()

## DEFINE CLASSES
## TODO 20240531
## * to [AG] w/o position filter,
## * test RAAS for several KRAQ definitions vs. Q->G,
## * rm *** for obvious positions,
## * test dotplot w lower p-value cutoff,
## * use -2 to 2 in M/C/W context, label MxM, etc.,
## * ->G context.

do.all.motifs <-  FALSE #TRUE # for interactive use to get pngs

aatight <- apply(aam[,as.character(-2:2)], 1, paste, collapse="")

kraq <- rep(FALSE, nrow(bdat))
kraq[grep("[KR]AQ",aatight)] <- TRUE
krgq <- rep(FALSE, nrow(bdat))
krgq[grep("[KR]GQ",aatight)] <- TRUE
krxq <- rep(FALSE, nrow(bdat))
krxq[grep("[KR][A-Z]Q",aatight)] <- TRUE

CTXT <- as.character(c(-2,-1,1,2))
CCxCC <- apply(aam[,CTXT], 1, function(x) any(x%in%c("C")))
MMxMM <- apply(aam[,CTXT], 1, function(x) any(x%in%c("M")))
WWxWW <- apply(aam[,CTXT], 1, function(x) any(x%in%c("W")))
GGxGG <- apply(aam[,CTXT], 1, function(x) any(x%in%c("G")))
CTXTT <- as.character(c(-1,1))
CxC <- apply(aam[,CTXTT], 1, function(x) any(x%in%c("C"))) 
MxM <- apply(aam[,CTXTT], 1, function(x) any(x%in%c("M")))
WxW <- apply(aam[,CTXTT], 1, function(x) any(x%in%c("W")))
GxG <- apply(aam[,CTXTT], 1, function(x) any(x%in%c("G")))

Wbranch <- WWxWW & bdat$from%in%c("I","L","V")
Mphos <- MMxMM & bdat$from%in%c("S","T")

classes <- cbind(
    kr=apply(aam[,as.character(-3:-1)], 1,
             function(x) any(x%in%c("K","R"))),
    KRAQ=kraq,
    KRGQ=krgq,
    KRXQ=krxq,
    CCxCC=CCxCC,
    MMxMM=MMxMM,
    WWxWW=WWxWW,
    GGxGG=GGxGG,
    CxC=CxC,
    MxM=MxM,
    WxW=WxW,
    GxG=GxG,
    Wbranch=Wbranch,
    Mphos=Mphos,
    Acidic= apply(aam[,as.character(c(-6:-1,1:6))], 1,
               function(x) sum(x%in%c("E","D")))>3,
    short=bdat$len <= 1000  & bdat$median > -.5,
    long=bdat$len > 1000  & bdat$median > -.5,
    high=bdat$median > -1,
    low=bdat$median <= -1,
    A=bdat$from=="A",
    A1= bdat$from=="A"& bdat$site%in%1,
    AG1= bdat$fromto%in%c("A:G","G:A")& bdat$site%in%1,
    P=bdat$from=="P",
    to_AG1=bdat$to%in%c("G","A") & bdat$site%in%1,
    to_AG23=bdat$to%in%c("G","A") & bdat$site%in%2:3,
    to_AG13=bdat$to%in%c("G","A") & bdat$site%in%1:3,
    to_AG=bdat$to%in%c("G","A"),
    Q=bdat$from=="Q",
    QA=bdat$fromto=="Q:A",
    QG=bdat$fromto=="Q:G",
    VQ=bdat$fromto=="V:Q",
    TV=bdat$fromto=="T:V",
    branched=bdat$from%in%c("I","L","V"),
    Nterm1=bdat$site%in%1,
    Nterm2=bdat$site%in%2)
selected <- colnames(classes)


fts <- unique(bdat$fromto)
ftcls <- matrix(NA, ncol=length(fts), nrow=nrow(aam))
colnames(ftcls) <-  paste0("fromto_",fts)
for ( i in seq_along(fts) ) 
    ftcls[,i] <- bdat$fromto==fts[i]
if ( !interactive() | do.all.motifs ) classes <- cbind(classes, ftcls)

frm <- unique(bdat$from)
ftcls <- matrix(NA, ncol=length(frm), nrow=nrow(aam))
colnames(ftcls) <- paste0("from_",frm)
for ( i in seq_along(frm) ) 
    ftcls[,i] <- bdat$from==frm[i]
if ( !interactive() | do.all.motifs  ) classes <- cbind(classes, ftcls)

frm <- unique(bdat$to)
ftcls <- matrix(NA, ncol=length(frm), nrow=nrow(aam))
colnames(ftcls) <- paste0("to_",frm)
for ( i in seq_along(frm) ) 
    ftcls[,i] <- bdat$to==frm[i]
if ( !interactive() | do.all.motifs  ) classes <- cbind(classes, ftcls)

cols <- as.character(c(-3:3))
rngs <- rep(list(cols), ncol(classes))
names(rngs) <- colnames(classes)

## custom ranges
rngs$QG <- as.character(-4:1)
rngs$QA <- as.character(-3:1)
rngs$KRAQ <- as.character(-4:4)
rngs$Acidic <- rngs$Long <- rngs$High <-as.character(-10:10)
rngs$MMxMM <- rngs$WWxWW <- rngs$CCxCC <-as.character(-2:2)
rngs$MxM <- rngs$WxW <- rngs$CxC <-as.character(-1:1)
rngs$to_G13 <- as.character(-3:0)

psig <- 10^-c(3,5,10)

## use our internal AA colors
ASN2 <- ASN
ASN2$cols <- aa.cols[ASN2$chars]
ASN2$cols["V"] <- "#2eb774"

log.path <- file.path(mfig.path, "logos")
dir.create(log.path, showWarnings=FALSE)
tmp.path <- file.path(mfig.path, "selected")
dir.create(tmp.path, showWarnings=FALSE)
##pwm <- getPFM(aam)
##seqLogo(pwm, stackHeight=sumProbabilities, alphabet=ASN2)

for ( i in 1:ncol(classes) ) {

    id <- colnames(classes)[i]
    lb <- id # sub(".*_", "", id)
    filt <- classes[,id]
    cols <- rngs[[id]]
    axlab <- as.character(cols)
    axlab[axlab=="0"] <- ""
 
    if ( sum(filt)<2 ) next

    tmp.path <- log.path
    if ( id %in% selected )
        tmp.path <- file.path(mfig.path, "selected")
    
    if ( interactive() )
        boxplot(bdat$median ~ filt)

    ## position weight matrices and diffLogo
    pfm1 <- getPFM(aam[ filt,cols,drop=FALSE])
    pfm2 <- getPFM(aam[!filt,cols,drop=FALSE])
    n1 <- sum(filt)
    n2 <- sum(!filt)
    dfob <- createDiffLogoObject(pfm1, pfm2, alphabet=ASN2)
    dfop <- enrichDiffLogoObjectWithPvalues(dfob, n1=n1, n2=n2)

    if ( any(dfop$pvals[cols!=0]<min(psig)) | id %in% selected ) {

        mmai <- c(.5,.5,.1,.1)
        wd <- .4*length(cols) + .6
        ht <- 3
        mn <- dfob$ylim.negMax
        mx <- dfob$ylim.posMax
        plotdev(file.path(tmp.path,paste0("logos_", id)),
                height=ht, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=mn, ymax=mx)
        axis(1, at=1:length(cols), labels=axlab)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels="AAS", las=2)
        figlabel(paste0(lb,": ",sum(filt)), pos="topright", font=2)
        diffLogo_addPvals(dfop, ymin=mx)
        dev.off()
        ## again w/o lower part
        plotdev(file.path(tmp.path,paste0("logos_", id, "_top")),
                height=ht/2, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=0, ymax=mx)
        axis(1, at=1:length(cols), labels=axlab)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels="AAS", las=2)
        figlabel(paste0(lb,": ",sum(filt)), pos="topright", font=2)
        diffLogo_addPvals(dfop, ymin=mx)
        dev.off()
    }

    ## position weight matrices and diffLogo
    pfm1 <- getPFM(aas[ filt,cols,drop=FALSE])
    pfm2 <- getPFM(aas[!filt,cols,drop=FALSE])
    n1 <- sum(filt)
    n2 <- sum(!filt)
    dfob <- createDiffLogoObject(pfm1, pfm2, alphabet=ASN2)
    dfop <- enrichDiffLogoObjectWithPvalues(dfob, n1=n1, n2=n2)

    if ( any(dfop$pvals[cols!=0]<min(psig))  | id %in% selected ) {

        mmai <- c(.5,.5,.1,.1)
        wd <- .4*length(cols) + .6
        ht <- 3
        mn <- dfob$ylim.negMax
        mx <- dfob$ylim.posMax
        plotdev(file.path(tmp.path,paste0("logos_", id,"_S")),
                height=ht, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=mn, ymax=mx)
        axis(1, at=1:length(cols), labels=axlab)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels=expression(INC), las=2)
        figlabel(paste0(lb,": ",sum(filt)), pos="topright", font=2)
        diffLogo_addPvals(dfop, ymin=mx)
        dev.off()
    }
    
}

plotdev(file.path(mfig.path,"protein_location_G"), ftyp,
        width=3, height=3, res=300)
hist(bdat$rpos[bdat$to=="G"])
dev.off()

## TODO:
## * look at BP with multiple sites,

### RAAS dot profiles for EACH selection:
## * spatial selections, as in _kraq.R;
##   all ->G at position 1-3 (high RAAS), and
##   PNEM-involving substitutions towards the center protein (AAS enrichments),
## * simply by cleavage dimer (BP enrichments),
## * all BP with [KR]AG in pos.1, or [KR]xQL in pos 2,  independent of AAS
##   (BP enrichments),
## * KRAQ (high RAAS),
## * W/M context (current motif figures).

## TODO: define motifs in bdat and map to tmtf
## MAP protein level info to TMT Data

tmtm <- tmtf

idx <- match(paste(tmtm$BP, tmtm$SAAP), paste(bdat$BP, bdat$SAAP))
ina <- which(is.na(idx))
if ( length(ina)>0 ) {
    cat(paste("TODO:", length(ina), "missing from unique saap file.\n"))
    tmtm <- tmtm[-ina,]
    idx <- idx[-ina]    
}

table(classes[,"WWxWW"],classes[,"QG"])

### NOTE: top-down filtering emphasizes last!!
### TODO: allow overlapping classes here!!
###       * use listProfile? and generate ovw structure
##          from multiple listProfile matrices?  

motclass <- rep("n.a.", nrow(classes))
##motclass[classes[,"QG"]] <- "QG"
##motclass[classes[,"QA"]] <- "QA"
##motclass[classes[,"KRXQ"]] <- "[KR]XQ"
##motclass[classes[,"KRGQ"]] <- "[KR]GQ"
##motclass[classes[,"KRAQ"]] <- "[KR]AQ"
motclass[classes[,"CCxCC"]] <- "CCxCC"
motclass[classes[,"MMxMM"]] <- "MMxMM"
motclass[classes[,"WWxWW"]] <- "WWxWW" ## NOTE: overrules more frequent MMxMM
##motclass[classes[,"GGxGG"]] <- "GGxGG" 
motclass[classes[,"CxC"]] <- "CxC"
motclass[classes[,"MxM"]] <- "MxM"
motclass[classes[,"WxW"]] <- "WxW" ## NOTE: overrules more frequent MMxMM
##motclass[classes[,"GxG"]] <- "GxG"
##motclass[classes[,"Mphos"]] <- "Mphos"
##motclass[classes[,"Wbranch"]] <- "Wbranch"
##motclass[classes[,"TV"]] <- "TV"
##
##motclass[classes[,"A1"]] <- "A1>"
##motclass[classes[,"AG1"]] <- "A>G1"
##motclass[classes[,"to_G1"]] <- ">G1"
motclass[classes[,"to_AG13"]] <- "KRAQ"

msrt <- c(
    ##"A1>",
    ##">G1",
    ##"A>G1",
    ##">[GA]123",
    ##"[KR]XQ",
    ##"[KR]AQ",
    ##"QG",
    ##"[KR]GQ",
    ##"QA",
    "KRAQ",
    "CCxCC",
    "CxC",
    ##"Mphos","Wbranch",
    "MMxMM",
    "MxM",
    ##"TV",
    "WWxWW",
    "WxW" )


### separate plots for motifs
dot.path <- file.path(mfig.path, "dotplots")
dir.create(dot.path)
for ( j in 1:ncol(classes) ) {

    mid <- colnames(classes)[j]
    mclass <- rep("n.a.", nrow(classes))
    mclass[classes[,j]] <- mid

    tmtm$motifs <- mclass[idx]
    
    ovw <- raasProfile(x=tmtm, id="SAAP", 
                   rows="motifs", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=mid,
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0)

    plotProfiles(ovw,
                 fname=file.path(dot.path,paste0("motifs_",SETID,"_",mid)),
                 mai=c(0,.7,0,.6), ttcols=ttcols, value="median",
                 p.min=mp.min, p.txt=mp.txt,
                 dot.sze=dot.sze, p.dot=mp.dot,
                 rlab=LAB,  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam="monospace")
 
}

## AA/codon/structure mapping
tmtm$motifs <- motclass[idx]

ovw <- raasProfile(x=tmtm, id="SAAP", 
                   rows="motifs", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=msrt,
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("motifs_",SETID,"_")))

plotProfiles(ovw, fname=file.path(mfig.path,paste0("motifs_",SETID)),
             mai=c(.8,.7,.6,.6), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace")

## first AA in BP
bp1 <- unlist(lapply(strsplit(bdat$BP,""), function(x) x[1]))
bp1[bp1==bdat$from & bdat$site==1] <- "AAS"
##bp1[bdat$fromto=="Q:G"] <- "QG" ##- NOTE: this makes A at 1st insign.
##bp1[bdat$fromto=="Q:A"] <- "QA" ##- NOTE: this makes A at 1st insign.
## AA/codon/structure mapping
tmtm$bp1 <- bp1[idx]

ovw <- raasProfile(x=tmtm, id="SAAP", 
                   rows="bp1", cols="Dataset",
                   bg=TRUE, value="RAAS", ##row.srt=motsrt,
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("bp1_",SETID,"_")))

plotProfiles(ovw, fname=file.path(mfig.path,paste0("bp1_",SETID)),
             mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="AA at position 1 in BP", mtxt.line=1.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)
ovwp <- sortOverlaps(ovw, p.min=mp.txt, cut=TRUE)
plotProfiles(ovwp, fname=file.path(mfig.path,paste0("bp1_",SETID,"_cut")),
             mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="AA at pos. 1 in BP", mtxt.line=1.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)

## second AA in BP
bp2 <- unlist(lapply(strsplit(bdat$BP,""), function(x) x[2]))
bp2[bp2==bdat$from & bdat$site==2] <- "AAS"

## AA/codon/structure mapping
tmtm$bp2 <- bp2[idx]

ovw <- raasProfile(x=tmtm, id="SAAP", 
                   rows="bp2", cols="Dataset",
                   bg=TRUE, value="RAAS", ##row.srt=motsrt,
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("bp2_",SETID,"_")))
plotProfiles(ovw, fname=file.path(mfig.path,paste0("bp2_",SETID)),
             mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="AA at position 2 in BP", mtxt.line=1.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace")
ovwp <- sortOverlaps(ovw, p.min=mp.txt, cut=TRUE)
plotProfiles(ovwp, fname=file.path(mfig.path,paste0("bp2_",SETID,"_cut")),
             mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="AA at pos. 2 in BP", mtxt.line=1.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)
## second AA in BP
bp3 <- unlist(lapply(strsplit(bdat$BP,""), function(x) x[3]))
bp3[bp3==bdat$from & bdat$site==3] <- "AAS"

## AA/codon/structure mapping
tmtm$bp3 <- bp3[idx]

ovw <- raasProfile(x=tmtm, id="SAAP", 
                   rows="bp3", cols="Dataset",
                   bg=TRUE, value="RAAS", ##row.srt=motsrt,
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("bp3_",SETID,"_")))
plotProfiles(ovw, fname=file.path(mfig.path,paste0("bp3_",SETID)),
             mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="AA at position 3 in BP", mtxt.line=1.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace")
ovwp <- sortOverlaps(ovw, p.min=mp.txt, cut=TRUE)
plotProfiles(ovwp, fname=file.path(mfig.path,paste0("bp3_",SETID,"_cut")),
             mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="AA at pos. 3 in BP", mtxt.line=1.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)

