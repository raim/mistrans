
### ANALYZE SEQUENCE CONTEXT around AAS

## AA ENRICHMENT AROUND AAS SITES


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

seq.path <- file.path(out.path,"motifs")
dir.create(seq.path, showWarnings=FALSE)
mfig.path <- file.path(fig.path,"motifs")
dir.create(mfig.path, showWarnings=FALSE)

## NOTE: 1e-6 doesn't seem to be worth it, little gain.
## but further experience required
mp.min <- p.min
mp.dot <- p.dot
mp.txt <- p.txt

## motif p-value cutoffs *, **, ***
psig <- 10^-c(3,5,10)

#### GENERAL SEQUENCE LOGOS

### AA MATRIX
aam <- do.call(rbind, strsplit(bdat$AA,""))
nc <- (ncol(aam)-1)/2
colnames(aam) <- -nc:nc
rownames(aam) <- paste0(bdat$BP,"_", bdat$SAAP)

## AA MATRIX WITH AAS
aas <- aam
aas[,"0"] <- bdat$to

## TODO: AAS_overlap for RAAS BINS

### HYPERGEO TESTS - tigther context

omai <- c(.5,.5,.6,.6)

ovl <- aaProfile(aam[,as.character(-7:7)], abc=AAT)
ovl <- sortOverlaps(ovl, p.min=p.txt, sign=1)

nw <- ncol(ovl$p.value)*.35 + omai[2] + omai[4]
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("AAS_overlap")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS", ylab="Encoded AA")
##figlabel("all", pos="bottomleft")
dev.off()

plotdev(file.path(mfig.path,paste0("AAS_dotplot")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
dotprofile(ovl, value="ratio", vcols=ttcols,
           xlab="Distance from AAS", ylab="Encoded AA",
           p.dot=p.dot, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), axis=1:2, show.total=TRUE)
figlabel("all", pos="bottomleft")
dev.off()

ovc <- sortOverlaps(ovl, p.min=p.txt, sign=1, cut=TRUE)
nh <- nrow(ovc$p.value)*.5 + omai[1] + omai[3]
plotdev(file.path(mfig.path,paste0("AAS_overlap_cut")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS", ylab="Encoded AA")
figlabel("all", pos="bottomleft")
dev.off()

## HIGH RAAS log10(RAAS)>-1

source("~/programs/segmenTools/R/clusterTools.R")

ovl <- aaProfile(aam[bdat$median> -1,as.character(-7:7)], abc=AAT)
ovl <- sortOverlaps(ovl, p.min=p.txt, sign=1, cut=TRUE)

nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("AAS_overlap_lRAAS-1")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS")
figlabel("RAAS > 0.1", pos="bottomleft")
dev.off()

plotdev(file.path(mfig.path,paste0("AAS_dotplot_lRAAS-1")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
dotprofile(ovl, value="ratio", vcols=ttcols,
           xlab="Distance from AAS",
           p.dot=p.dot, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), ylab="amino acid", axis=1:2, show.total=TRUE)
figlabel("RAAS > 0.1", pos="bottomleft")
dev.off()

## MAX RAAS  log10(RAAS)>0.5
ovl <- aaProfile(aam[bdat$median>0,as.character(-7:7)], abc=AAT)
ovl <- sortOverlaps(ovl, p.min=p.txt, sign=1, cut=TRUE)
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("AAS_overlap_lRAAS-.5")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS")
figlabel("RAAS > 0.32", pos="bottomleft")
dev.off()

plotdev(file.path(mfig.path,paste0("AAS_dotplot_lRAAS-.5")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
dotprofile(ovl, value="ratio", vcols=ttcols,
           xlab="Distance from AAS",
           p.dot=p.dot, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), ylab="amino acid", axis=1:2, show.total=TRUE)
figlabel("RAAS > 0.32", pos="bottomleft")
dev.off()


### DEFINE CLASSES

## TODO 20240531
## * to [AG] w/o position filter,
## * test RAAS for several KRAQ definitions vs. Q->G,
## * rm *** for obvious positions,
## * test dotplot w lower p-value cutoff,
## * use -2 to 2 in M/C/W context, label MxM, etc.,
## * ->G context.

do.all.motifs <-  FALSE # TRUE # 

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

## subclasses

Wbranch <- WWxWW & bdat$from%in%c("I","L","V")
Mphos <- MMxMM & bdat$from%in%c("S","T")

CCxPP <- apply(aam[,c("-1","-2")], 1, function(x) any(x%in%c("C"))) &
    apply(aam[,c("1","2")], 1, function(x) any(x%in%c("P")))

## A in position 1 of BP
bpA1 <- unlist(lapply(strsplit(bdat$AA, ""), function(x) x[1]=="A"))
bpQ2 <- unlist(lapply(strsplit(bdat$AA, ""), function(x) x[2]=="Q"))

aaall <- apply(aam, 1, paste, collapse="")
kraqa <- rep(FALSE, nrow(bdat))
kraqa[grep("[KR]AQ",aaall)] <- TRUE

## Mid-Peptide (non N/C-terminal)
## TODO: CHECK THIS SELECTION

MIDPEP <- 3
nt <- bdat$site
ct <- nchar(bdat$BP)-bdat$site +1
## fuse central AAS
nt[nt>MIDPEP] <- ct[ct>MIDPEP] <- MIDPEP+1
at <- paste0("N",nt)
at[ct<nt] <- paste0("C",ct[ct<nt])

enmpc <- c("E","N","M","P","C")
ENM <- bdat$from%in%enmpc | bdat$to%in%enmpc 


classes <- cbind(
    kr=apply(aam[,as.character(-3:-1)], 1,
             function(x) any(x%in%c("K","R"))),
    frQ=bdat$from=="Q",
    QG=bdat$fromto=="Q:G",
    QA=bdat$fromto=="Q:A",
    KRAQ=kraq,
    KRGQ=krgq,
    KRxQ=krxq,
    KRAQ_all=kraqa,
    KRAQ_o=kraq & ! bdat$from%in%c("Q"),
    bpA1=bpA1,
    bpQ2=bpQ2,
    bpAQ=bpQ2 & bpA1,
    toAG1 =bdat$to%in%c("G","A") & bdat$site%in%1,
    toAG23=bdat$to%in%c("G","A") & bdat$site%in%2:3,
    toAG13=bdat$to%in%c("G","A") & bdat$site%in%1:3,
    toAG  =bdat$to%in%c("G","A"),
    disord.  =bdat$iupred3 > .6,
    disord._high  =bdat$iupred3 > .6 & bdat$median > -1,
    binding  =bdat$DisoRDPbind > .6,
    binding_high  =bdat$DisoRDPbind > .6 & bdat$median > -1,
    noncons.  =bdat$MMSeq2 < 1,
    noncons._high  =bdat$MMSeq2 < 1 & bdat$median > -1,
    N1=bdat$site%in%1,
    N2=bdat$site%in%2,
    N3=bdat$site%in%3,
    N13=bdat$site%in%1:3,
    frA=bdat$from=="A",
    frA1= bdat$from=="A"& bdat$site%in%1,
    AG1= bdat$fromto%in%c("A:G","G:A")& bdat$site%in%1,

    center= at==paste0("N",MIDPEP),
    ENM= ENM, 
    ENM.center= at==paste0("N",MIDPEP) & ENM, 

    CCxCC=CCxCC,
    CxC=CxC,
    CCxPP=CCxPP,
    MMxMM=MMxMM,
    MxM=MxM,
    "MM[ST]MM"=Mphos,
    WWxWW=WWxWW,
    WxW=WxW,
    "WW[ILV]WW"=Wbranch,
    "[ILV]"=bdat$from%in%c("I","L","V"),
    GGxGG=GGxGG,
    GxG=GxG,
    frW=bdat$from=="W",
    
    QNrich= apply(aam[,as.character(c(-6:-1,1:6))], 1,
               function(x) sum(x%in%c("Q","N")))>2,
    Acidic= apply(aam[,as.character(c(-6:-1,1:6))], 1,
               function(x) sum(x%in%c("E","D")))>2,
    short=bdat$len <= 1000,
    long=bdat$len > 1000,
    longlow=bdat$len > 1000 & bdat$median <= -.5,
    longhigh=bdat$len > 1000 & bdat$median > -.5,
    longlow1=bdat$len > 1000 & bdat$median <= -1,
    longhigh1=bdat$len > 1000 & bdat$median > -1,
    max=bdat$median > 0,
    high=bdat$median > -1,
    low=bdat$median < -3)

## NA values in filters
classes[is.na(classes)] <- FALSE

## store manual selection before (otpionally)
## expanding to full scan of AAS types
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

cols <- as.character(c(-5:5))

## custom ranges
rngs <- rep(list(cols), ncol(classes))
names(rngs) <- colnames(classes)

rngs$QG <- as.character(-4:1)
rngs$QA <- as.character(-3:1)
rngs$KRAQ <- as.character(-4:4)

### LONG RANGE

rngs$Acidic <- rngs$long <- rngs$QNrich <- rngs$High <-
    rngs$longhigh <- rngs$longlow <- rngs$longhigh1 <- rngs$longlow1 <-
        rngs$disord. <-rngs$disord._high <-
            rngs$noncons. <-  rngs$binding <-
                rngs$noncons._high <-  rngs$binding_high <- as.character(-10:10)

rngs$MMxMM <- rngs$WWxWW <- rngs$CCxCC <- rngs$CCxPP <-as.character(-2:2)
rngs$MxM <- rngs$WxW <- rngs$CxC <- as.character(-1:1)
rngs$toAG13 <- as.character(-3:0)
rngs$toAG <- as.character(-3:1)

## suppress p-value indicators
opval <- rep(list(""), ncol(classes))
names(opval) <- colnames(classes)
opval$MMxMM <- opval$WWxWW <- opval$CCxCC <- opval$CCxPP <-
    as.character(c(-2,-1,1,2))
opval$MxM <- opval$WxW <- opval$CxC <- as.character(c(-1,1))
## suppress on incorporated side only
npval <- opval
npval[grep("toAG",names(npval))] <- "0"

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

    

    ## PWM ENCODED
    ## position weight matrices and diffLogo
    pfm1 <- getPFM(aam[ filt,cols,drop=FALSE])
    pfm2 <- getPFM(aam[!filt,cols,drop=FALSE])

    n1 <- sum(filt)
    n2 <- sum(!filt)
    dfob <- createDiffLogoObject(pfm1, pfm2, alphabet=ASN2)
    dfop <- enrichDiffLogoObjectWithPvalues(dfob, n1=n1, n2=n2)

    ## suppress p.vals where ALL AA are equal
    ## -> obvious selection bias
    setpo <- apply(pfm1, 2, function(x) all(x%in%c(0.0,1.0)))
    ## manual pval suppression
    setpo[cols%in%opval[[id]]] <- TRUE
    dfop$pvals[setpo] <- 1

    ## PWM INCORPORATED
    ## position weight matrices and diffLogo
    pfn1 <- getPFM(aas[ filt,cols,drop=FALSE])
    pfn2 <- getPFM(aas[!filt,cols,drop=FALSE])
    n1 <- sum(filt)
    n2 <- sum(!filt)

    dfnb <- createDiffLogoObject(pfn1, pfn2, alphabet=ASN2)
    dfnp <- enrichDiffLogoObjectWithPvalues(dfnb, n1=n1, n2=n2)

    ## suppress p.vals where ALL AA are equal
    ## -> obvious selection bias
    setpn <- apply(pfm1, 2, function(x) all(x%in%c(0.0,1.0)))
    ## manual pval suppression
    setpn[cols%in%npval[[id]]] <- TRUE
    dfnp$pvals[setpn] <- 1

    ## write out sequences for more detailed analysis
    ## TODO: currently not used
    if ( FALSE ) {
        motseq <- aam[ filt,]
        refseq <- aam[!filt,]
        fname <- file.path(seq.path,paste0(id,".fa")) 
        if ( file.exists(fname) ) unlink(fname)
        for ( j in 1:nrow(motseq) ) {
            osq <- gsub("-","",paste0(motseq[j,], collapse=""))
            cat(paste0(">", rownames(motseq)[j],"\n", osq, "\n"),
                file=fname, append=TRUE)
        }
    }

    ## write out POSITION WEIGHT MATRICES
    ## for genome-wide scan
    fname <- file.path(seq.path,paste0(id,"pwm.tsv")) 
    colnames(pfm1) <- cols
    write.table(pfm1, file=fname, sep="\t", quote=FALSE, row.names=TRUE)
        
    ## PLOT LOGOS

    ## ENCODED
    if ( any(dfop$pvals[cols!=0]<min(psig)) | id %in% selected ) {

        mmai <- c(.5,.5,.25,.1)
        wd <- .4*length(cols) + mmai[2] + mmai[4]
        ht <- 2.5
        mn <- dfob$ylim.negMax
        mx <- dfob$ylim.posMax
        plotdev(file.path(tmp.path,paste0("logos_", id,"_encoded")),
                height=ht, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=0, ymax=mx)
        axis(2)
        axis(1, at=1:length(cols), labels=axlab)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels="AAS", las=2)
        figlabel(paste0(lb,": ",sum(filt)), pos="topright", font=2)
        diffLogo_addPvals(dfop, ymin=mx, levels=psig)
        ##mtext("encoded", 4,-.25,adj=.05)
        dev.off()
    }

    ## INCORPORATED
    if ( any(dfnp$pvals[cols!=0]<min(psig))  | id %in% selected ) {

        mmai <- c(.5,.5,.25,.1)
        wd <- .4*length(cols) + mmai[2] + mmai[4]
        ht <- 2.5
        mn <- dfnb$ylim.negMax
        mx <- dfnb$ylim.posMax
        plotdev(file.path(tmp.path,paste0("logos_", id,"_incorporated")),
                height=ht, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfnb, sparse=TRUE, ymin=0, ymax=mx)
        axis(1, at=1:length(cols), labels=axlab)
        axis(2)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels=expression(INC), las=2)
        figlabel(paste0(lb,": ",sum(filt)), pos="topright", font=2)
        diffLogo_addPvals(dfnp, ymin=mx, levels=psig)
        ##mtext("incorporated", 4,-.25,adj=.05)
        dev.off()
    }
    
    ## TIGHT LOGO PLOT FOR MAIN
    if ( (any(dfop$pvals[cols!=0]<min(psig))  |
          any(dfnp$pvals[cols!=0]<min(psig))) |
         id %in% selected ) {
        
        mno <- dfob$ylim.negMax
        mxo <- dfob$ylim.posMax
        mnn <- dfnb$ylim.negMax
        mxn <- dfnb$ylim.posMax
        
        ##mxo <- mxn <- max(mxo, mxn)
        
        ## figure heights
        ht <- 3
        omai <- c(.15,.25,.25,.05)
        nmai <- c(.25,.25,.15,.05)

        ## heights
        htn <- ht/2 + nmai[1] + nmai[3]
        hto <- ht/2 + omai[1] + omai[3]
        
        ## adaptive width
        wd <- .4*length(cols) + nmai[2] + nmai[4]
        
        dfop$ylab <- "JS divergence"
        plotdev(file.path(tmp.path,paste0("logos_", id, "")),
                height=hto+htn , width=wd, res=300, type=ftyp)
        layout(t(t(1:2)), heights=c(hto, htn))

        ## encoded
        par(mai=omai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=0, ymax=mxo)
        axis(1, at=1:length(cols), labels=FALSE)
        axis(2, cex.axis=1.2)
        figlabel(paste0(lb,": ",sum(filt),"  "),
                 pos="topright", font=2, cex=1.2)
        diffLogo_addPvals(dfop, ymin=mxo, levels=psig)
        ##mtext("encoded", 4,-.25,adj=.05)
        
        ## incorporated
        par(mai=nmai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfnb, sparse=TRUE, ymin=0, ymax=mxn)
        axis(1, at=1:length(cols), labels=axlab, cex.axis=1.2)
        axis(2, cex.axis=1.2)
        ##mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels="AAS", las=1, cex.axis=1.2)
        ##text(par("usr")[2], mxn*.75, "incorporated", pos=2)
        diffLogo_addPvals(dfnp, ymin=mxn, levels=psig)
        ##mtext("incorporated", 4,-.25,adj=.05)
        dev.off()
  
   
    }
}

## y-axis label for tight motif plot
## figure heights
ht <- 3
omai <- c(.15,.25,.25,.15)
nmai <- c(.25,.25,.15,.15)
ymai <- c(.25,0,.25,0)
## heights
htn <- ht/2 + nmai[1] + nmai[3]
hto <- ht/2 + omai[1] + omai[3]
 
plotdev(file.path(mfig.path,paste0("logos_ylab")),
        height=hto+htn , width=.6, res=300, type=ftyp)
par(mai=ymai)
plot(1,1, axes=FALSE, col=NA, col.axis=NA)
text(1.28,1,label="JS divergence", xpd=TRUE,srt=90, cex=1.3)
text(0.775,1.25,label="encoded", xpd=TRUE,srt=90, cex=1.7)
text(0.8,.775,label="incorporated", xpd=TRUE,srt=90, cex=1.7)
dev.off()


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

## TODO: is this filtering necessary??
TIDX <- match(paste(tmtm$BP, tmtm$SAAP), paste(bdat$BP, bdat$SAAP))
ina <- which(is.na(TIDX))
if ( length(ina)>0 ) {
    cat(paste("TODO:", length(ina), "missing from unique saap file.\n"))
    tmtm <- tmtm[-ina,]
    TIDX <- TIDX[-ina]    
}

table(classes[,"WWxWW"],classes[,"QG"])




### CALCULATE ALL MOTIF RAAS PROFILES
## TODO: implement overlapping classes in raasProfile
plot.all.motifs <- FALSE
dot.path <- file.path(mfig.path, "dotplots")
dir.create(dot.path)
covw <- list()
for ( j in 1:length(selected) ) {
    
    mid <- selected[j]
    mclass <- rep("n.a.", nrow(classes))
    mclass[classes[,j]] <- mid

    tmtm$motifs <- mclass[TIDX]
    
    ovw <- raasProfile(x=tmtm, id="SAAP", 
                       rows="motifs", cols="Dataset",
                       bg=TRUE, value="RAAS", row.srt=mid,
                       col.srt=uds,
                       use.test=use.test, do.plots=FALSE,
                       xlab=xl.raas,
                       verb=0)
    covw[[mid]] <- ovw

    if ( plot.all.motifs ) {
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
}

### MERGE PROFILES
ovm <- mergeProfiles(covw)

source("~/work/mistrans/scripts/saap_utils.R")

## use RAAS profile with matrix input (new 20240604)
tcls <- classes[TIDX,]
ovm2 <- raasProfile(x=tmtm, id="SAAP", 
                    rows=tcls, cols="Dataset",
                    bg=TRUE, value="RAAS", 
                    col.srt=uds,
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=0)
## TODO: test new RAAS profiler with matrix input
## and remove above loop if all is good.
if ( all(ovm$p.value==ovm2$.p.value))
    warning("new profile ok, rm old code")


## sort
ovm <- sortOverlaps(ovm, axis=2, srt=selected)

amai <- c(0.8,1.2,0.6,.6)

plotProfiles(ovm,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"_all")),
             mai=amai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace")


plotProfiles(ovm,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_pval10")),
             mai=amai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=1e-10,
             llab="p10",  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace",plot.all=FALSE)
plotProfiles(ovm,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_pval6")),
             mai=amai, ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=1e-6,
             llab="p6",  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace",plot.all=FALSE)
plotProfiles(ovm,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_pval100")),
             mai=amai, ttcols=ttcols, value="median",
             p.min=1e-100, p.txt=1e-50,
             dot.sze=dot.sze, p.dot=1e-100,
             llab="p100",  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace",plot.all=FALSE)
plotProfiles(ovm,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_pval30")),
             mai=amai, ttcols=ttcols, value="median",
             p.min=1e-30, p.txt=1e-15,
             dot.sze=dot.sze, p.dot=1e-30,
             llab="p30",  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace",plot.all=FALSE)

plotProfiles(ovm,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_vcols")),
             mai=amai, ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             llab="vcols",  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols, ffam="monospace")

## SELECTED MOTIS FOR MAIN
msrt <- c("QG",
          "toAG13",
          "CCxCC",
          "MMxMM",
          "WWxWW")

ovs <- sortOverlaps(ovm, axis=2, srt=msrt)
plotProfiles(ovs,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"")),
             mai=c(0.05,CMAIL,0.05,.6), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace")

ssrt <- c("disord.",
          "binding",
          "noncons.",
          "long")
ovs <- sortOverlaps(ovm, axis=2, srt=ssrt)
plotProfiles(ovs,
             fname=file.path(mfig.path,paste0("structure_",SETID,"")),
             mai=c(0.05,CMAIL,0.05,.6), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace")


### ANALYZE OVERLAP MOTIFS vs STRUCTURE


CLS <- classes[,c(ssrt, msrt)]

## CORRELATION OF MEASURES
image_matrix(cor(classes), col=ttcols, breaks=seq(-1,1,length=length(ttcols)+1),
             axis=1:2)


ovl <- clusterCluster(CLS[,"disord."], CLS[,"noncons."])
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt)

rbins <- cut(bdat$median, breaks=c(-6,-4,-2,-1,0,3))
rsrt <- levels(rbins)
rbins <- as.character(rbins)
names(rbins) <- rownames(CLS)
rbins[is.na(rbins)] <- "na."
##
ovl <- clusterAnnotation(cls=rbins, data=CLS, cls.srt=rsrt)
omai <- c(.75,1,.6,.6)
nw <- ncol(ovl$p.value)*.35 + omai[2] + omai[4]
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_raas")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, xlab="", ylab="",
             show.total=TRUE)
figlabel("median RAAS",pos="bottomright")
dev.off()


## for all non-conserved bins
csrt <- as.character(CLS[,"noncons."])
csrt[is.na(csrt)] <- "na."

rbins <- cut(bdat$median, breaks=c(-6,-4,-2,-1,0,3))
rsrt <- levels(rbins)
rbins <- as.character(rbins)
names(rbins) <- rownames(CLS)
rbins[is.na(rbins)] <- "na."

ovl <- clusterAnnotation(cls=bdat$MMSeq2.bins, data=CLS,
                         cls.srt=levels(MMSeq2.bins))
omai <- c(.75,1,.6,.6)
nw <- ncol(ovl$p.value)*.35 + omai[2] + omai[4]
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_conservation")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, xlab="", ylab="",
             show.total=TRUE)
mtext("conservation",1, 2.5)
dev.off()



## ALTERNATIVE RAAS PROFILES
ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="iupred3.bins", cols="MMSeq2.bins",
                    row.srt=rev(levels(iupred3.bins)),
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value="RAAS", 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,.5,.6,.6)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_conservation_disorder_raas")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE)
polygon(y=c(1, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("disorder",2, 1.3)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()

ovls <- raasProfile(x=tmtt, id="SAAP", 
                    rows=tcls[,c(msrt,ssrt)], cols="MMSeq2.bins",
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value="RAAS", 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,1,.6,.6)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]
plotdev(file.path(mfig.path,paste0("classes_conservation_motifs_raas")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, axis=2,
           xlab=NA, ylab=NA, show.total=TRUE)
mtext("conservation",1, 1.3)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
dev.off()

ovls <- raasProfile(x=tmtt, id="SAAP", 
                    rows=tcls[,c(msrt,ssrt)], cols="iupred3.bins",
                    col.srt=levels(iupred3.bins),
                    bg=FALSE, value="RAAS", 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,1,.6,.6)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]
plotdev(file.path(mfig.path,paste0("classes_disorder_motifs_raas")),
        height=nh, width=nw, res=300)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="monospace")
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, axis=2,
           xlab=NA, ylab=NA, show.total=TRUE)
polygon(x=c(1, ncol(ovls$p.value), ncol(ovls$p.value)),
            y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("disorder",1, 1.3)
dev.off()

## SELECTED MOTIFS FOR P-VAL COMPARISON
msrt <- c("toAG13",
          "CCxCC",
          "MMxMM",
          "WWxWW",
          "QNrich",
          "disord.",
          "Acidic",
          "long"
          )

ovs <- sortOverlaps(ovm, axis=2, srt=msrt)
plotProfiles(ovs,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"_selected")),
             mai=c(0.8,.9,0.6,.6), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace")

##


### POSITION WISE RAAS

## TODO: add to motifs?

## first AA in BP
bp1 <- unlist(lapply(strsplit(bdat$BP,""), function(x) x[1]))
bp1[bp1==bdat$from & bdat$site==1] <- "AAS"
##bp1[bdat$fromto=="Q:G"] <- "QG" ##- NOTE: this makes A at 1st insign.
##bp1[bdat$fromto=="Q:A"] <- "QA" ##- NOTE: this makes A at 1st insign.
## AA/codon/structure mapping
tmtm$bp1 <- bp1[TIDX]

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
tmtm$bp2 <- bp2[TIDX]

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
tmtm$bp3 <- bp3[TIDX]

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

