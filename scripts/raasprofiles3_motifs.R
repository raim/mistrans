
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

## remove identical sequences before motif plots
remove.duplicated.sequences <- TRUE

#### GENERAL SEQUENCE LOGOS

### remove duplicate sites!
## TODO: is there a cleaner way to do this? here we randomly loose
## 
cdat <- bdat[!duplicated(paste(bdat$ensembl,bdat$pos, bdat$fromto)),]


### AA MATRIX
aam <- do.call(rbind, strsplit(cdat$AA,""))
nc <- (ncol(aam)-1)/2
colnames(aam) <- -nc:nc
rownames(aam) <- paste0(cdat$BP,"_", cdat$SAAP)

## AA MATRIX WITH AAS
aas <- aam
aas[,"0"] <- cdat$to

## TODO: AAS_overlap for RAAS BINS

### HYPERGEO TESTS - tigther context

omai <- c(.5,.5,.6,.6)

ovl <- aaProfile(aam[,as.character(-7:7)], abc=AAT)
ovl <- sortOverlaps(ovl, p.min=p.txt, sign=1)

nw <- ncol(ovl$p.value)*.35 + omai[2] + omai[4]
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("AAS_overlap")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS", ylab="Encoded AA")
##figlabel("all", pos="bottomleft")
dev.off()

plotdev(file.path(mfig.path,paste0("AAS_dotplot")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovl, value="ratio", vcols=ttcols,
           xlab="Distance from AAS", ylab="Encoded AA",
           p.dot=p.dot, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), axis=1:2, show.total=TRUE)
figlabel("all", pos="bottomleft")
dev.off()

ovc <- sortOverlaps(ovl, p.min=p.txt, sign=1, cut=TRUE)
nh <- nrow(ovc$p.value)*.5 + omai[1] + omai[3]
plotdev(file.path(mfig.path,paste0("AAS_overlap_cut")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS", ylab="Encoded AA")
figlabel("all", pos="bottomleft")
dev.off()

## HIGH RAAS log10(RAAS)>-1

source("~/programs/segmenTools/R/clusterTools.R")

ovl <- aaProfile(aam[cdat$median> -1,as.character(-7:7)], abc=AAT)
ovl <- sortOverlaps(ovl, p.min=p.txt, sign=1, cut=TRUE)

nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("AAS_overlap_lRAAS-1")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS")
figlabel("RAAS > 0.1", pos="bottomleft")
dev.off()

plotdev(file.path(mfig.path,paste0("AAS_dotplot_lRAAS-1")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovl, value="ratio", vcols=ttcols,
           xlab="Distance from AAS",
           p.dot=p.dot, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), ylab="amino acid", axis=1:2, show.total=TRUE)
figlabel("RAAS > 0.1", pos="bottomleft")
dev.off()

## MAX RAAS  log10(RAAS)>0.5
ovl <- aaProfile(aam[cdat$median>0,as.character(-7:7)], abc=AAT)
ovl <- sortOverlaps(ovl, p.min=p.txt, sign=1, cut=TRUE)
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("AAS_overlap_lRAAS-.5")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS")
figlabel("RAAS > 0.32", pos="bottomleft")
dev.off()

plotdev(file.path(mfig.path,paste0("AAS_dotplot_lRAAS-.5")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
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

kraq <- rep(FALSE, nrow(cdat))
kraq[grep("[KR]AQ",aatight)] <- TRUE
krgq <- rep(FALSE, nrow(cdat))
krgq[grep("[KR]GQ",aatight)] <- TRUE
krxq <- rep(FALSE, nrow(cdat))
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

Wbranch <- WWxWW & cdat$from%in%c("I","L","V")
Mphos <- MMxMM & cdat$from%in%c("S","T")

xP <- aam[,"1"] == "P"


CCxPP <- apply(aam[,c("-1","-2")], 1, function(x) any(x%in%c("C"))) &
    apply(aam[,c("1","2")], 1, function(x) any(x%in%c("P")))

## A in position 1 of BP
bpA1 <- unlist(lapply(strsplit(cdat$AA, ""), function(x) x[1]=="A"))
bpQ2 <- unlist(lapply(strsplit(cdat$AA, ""), function(x) x[2]=="Q"))

aaall <- apply(aam, 1, paste, collapse="")
kraqa <- rep(FALSE, nrow(cdat))
kraqa[grep("[KR]AQ",aaall)] <- TRUE

## Mid-Peptide (non N/C-terminal)
## TODO: CHECK THIS SELECTION

MIDPEP <- 3
nt <- cdat$site
ct <- nchar(cdat$BP)-cdat$site +1
## fuse central AAS
nt[nt>MIDPEP] <- ct[ct>MIDPEP] <- MIDPEP+1
at <- paste0("N",nt)
at[ct<nt] <- paste0("C",ct[ct<nt])

enmpc <- c("E","N","M","P","C")
ENM <- cdat$from%in%enmpc | cdat$to%in%enmpc 


classes <- cbind(
    "Q:x"=cdat$from=="Q",
    "Q:G"=cdat$fromto=="Q:G",
    "Q:A"=cdat$fromto=="Q:A",
    KRAQ_d=kraq,
    KRGQ=krgq,
    KRxQ=krxq,
    KRAQ_all=kraqa,
    KRAQ_o=kraq & ! cdat$from%in%c("Q"),
    bpA1=bpA1,
    bpQ2=bpQ2,
    bpAQ=bpQ2 & bpA1,
    "[KR]n-3"=apply(aam[,as.character(-3:-1)], 1,
             function(x) any(x%in%c("K","R"))),
    n1=cdat$site%in%1,
    n2=cdat$site%in%2,
    n3=cdat$site%in%3,
    n13=cdat$site%in%1:3,
    "A:x"=cdat$from=="A",
    "x:Gn1"= cdat$to=="G"& cdat$site%in%1,
    "A:xn1"= cdat$from=="A"& cdat$site%in%1,
    "[AG]:[AG]n1"= cdat$fromto%in%c("A:G","G:A")& cdat$site%in%1,
    "x:[AG]n1"=cdat$to%in%c("G","A") & cdat$site%in%1,
    "x:[AG]n2"=cdat$to%in%c("G","A") & cdat$site%in%2:3,
    "KRAQ"=cdat$to%in%c("G","A") & cdat$site%in%1:3,
    "x:[AG]"  =cdat$to%in%c("G","A"),

    center= at==paste0("N",MIDPEP),
    ENM= ENM, 
    ENM.center= at==paste0("N",MIDPEP) & ENM, 

    CCxCC=CCxCC,
    CxC=CxC,
    CCxPP=CCxPP,
    xP=xP,
    MMxMM=MMxMM,
    MMxMM_PDAC=MMxMM & cdat$Datasets=="PDAC",
    MMxMM_others=MMxMM & cdat$Datasets!="PDAC",
    "MM[T:V]MM"=MMxMM  &  cdat$fromto=="T:V",
    MxM=MxM,
    "MM[ST]MM"=Mphos,
    WWxWW=WWxWW,
    WxW=WxW,
    "WW[ILV]WW"=Wbranch,
    "[ILV]:x"=cdat$from%in%c("I","L","V"),
    GGxGG=GGxGG,
    GxG=GxG,
    "W:x"=cdat$from=="W",

    "T:V"=cdat$fromto=="T:V",
    "T:V_PDAC"=cdat$fromto=="T:V" & cdat$Datasets=="PDAC",
    "T:V_others"=cdat$fromto=="T:V" & cdat$Datasets!="PDAC",

    QNrich= apply(aam[,as.character(c(-6:-1,1:6))], 1,
               function(x) sum(x%in%c("Q","N")))>2,
    Acidic= apply(aam[,as.character(c(-6:-1,1:6))], 1,
               function(x) sum(x%in%c("E","D")))>2,

    disord.  =cdat$iupred3 > .6,
    disord._high  =cdat$iupred3 > .6 & cdat$median > -1,
    binding  =cdat$DisoRDPbind > .6,
    binding_high  =cdat$DisoRDPbind > .6 & cdat$median > -1,
    noncons.  =cdat$MMSeq2 < 1,
    noncons._high  =cdat$MMSeq2 < 1 & cdat$median > -1,

    short=cdat$len <= 1000,
    short_high=cdat$len <= 1000 & cdat$median > -1,
    long=cdat$len > 1000,
    longlow=cdat$len > 1000 & cdat$median <= -.5,
    longhigh=cdat$len > 1000 & cdat$median > -.5,
    longlow1=cdat$len > 1000 & cdat$median <= -1,
    longhigh1=cdat$len > 1000 & cdat$median > -1,
    max=cdat$median > 0,
    high=cdat$median > -1,
    low=cdat$median < -3)

## NA values in filters
classes[is.na(classes)] <- FALSE

## store manual selection before (otpionally)
## expanding to full scan of AAS types
selected <- colnames(classes)


fts <- unique(cdat$fromto)
ftcls <- matrix(NA, ncol=length(fts), nrow=nrow(aam))
colnames(ftcls) <-  paste0("fromto_",fts)
for ( i in seq_along(fts) ) 
    ftcls[,i] <- cdat$fromto==fts[i]
if ( !interactive() | do.all.motifs ) classes <- cbind(classes, ftcls)

frm <- unique(cdat$from)
ftcls <- matrix(NA, ncol=length(frm), nrow=nrow(aam))
colnames(ftcls) <- paste0("from_",frm)
for ( i in seq_along(frm) ) 
    ftcls[,i] <- cdat$from==frm[i]
if ( !interactive() | do.all.motifs  ) classes <- cbind(classes, ftcls)

frm <- unique(cdat$to)
ftcls <- matrix(NA, ncol=length(frm), nrow=nrow(aam))
colnames(ftcls) <- paste0("to_",frm)
for ( i in seq_along(frm) ) 
    ftcls[,i] <- cdat$to==frm[i]
if ( !interactive() | do.all.motifs  ) classes <- cbind(classes, ftcls)

cols <- as.character(c(-3:3))

## custom ranges
rngs <- rep(list(cols), ncol(classes))
names(rngs) <- colnames(classes)

rngs$"Q:G" <- as.character(-4:1)
rngs$"Q:A" <- as.character(-3:1)
rngs$KRAQ_d <- as.character(-4:4)

### LONG RANGE

rngs$Acidic <- rngs$long <- rngs$QNrich <- rngs$High <-
    rngs$longhigh <- rngs$longlow <- rngs$longhigh1 <- rngs$longlow1 <-
        rngs$disord. <-rngs$disord._high <-
            rngs$noncons. <-  rngs$binding <-
                rngs$noncons._high <-  rngs$binding_high <- as.character(-10:10)

rngs$MMxMM <- rngs$WWxWW <- rngs$CCxCC <- rngs$CCxPP <-as.character(-2:2)
rngs$"MM[T:V]MM" <- rngs$"MMxMM_PDAC" <- rngs$"MMxMM_others" <- as.character(-2:2)
rngs$"T:V" <- rngs$"T:V_PDAC" <- rngs$"T:V_others" <- as.character(-2:2)
rngs$MxM <- rngs$WxW <- rngs$CxC <- as.character(-1:1)
rngs$"KRAQ" <- as.character(-3:0)
rngs$"x:[AG]" <- as.character(-3:1)

## suppress p-value indicators
opval <- rep(list(""), ncol(classes))
names(opval) <- colnames(classes)
opval$MMxMM <- opval$WWxWW <- opval$CCxCC <- opval$CCxPP <-
    as.character(c(-2,-1,1,2))
opval$MxM <- opval$WxW <- opval$CxC <- as.character(c(-1,1))

## suppress on incorporated side only
npval <- opval
npval[grep("x:\\[AG\\]",names(npval))] <- "0"
npval["KRAQ"] <- "0"

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

### MANUAL INVESTIGATION OF T->V

if ( interactive() ) {

    sum(classes[,"T:V"]) # 215 unique BP/SAAP with T:V
    sum(classes[,"T:V"]& classes [,"MMxMM"]) # 128
    table(cdat$Datasets[classes[,"T:V"] & classes [,"MMxMM"]])

    ovl <- clusterCluster(classes[,"T:V"], classes[,"MMxMM"])
    plotOverlaps(ovl, ylab="T:V", xlab="MMxMM")

    mclasses <- rep("other", nrow(classes))
    mclasses[classes[,"MMxMM"]] <- "MMxMM"
    mclasses[classes[,"T:V"]] <- "T:V"
    mclasses[classes[,"MMxMM"] & classes[,"T:V"]] <- "MMxMM & T:V"
    msrt <- c("MMxMM","T:V","MMxMM & T:V","other")


    dsets <- cdat$Datasets
    ##dsets[dsets=="Healthy"] <- site$Tissues[dsets=="Healthy"]
    nsets <- stringr::str_count(dsets, ";")+1
    ##dsets[nsets > 1] <- ">1"
    ##dsets[nsets > 2] <- ">2"
    dsets[nsets > 3] <- ">3"
    ds.srt <- unique(dsets)
    ds.srt <- ds.srt[order(nchar(ds.srt))]

    ##ovl <- clusterCluster(mclasses, cdat$Datasets, cl1.srt=msrt)

    ovl <- clusterAnnotation(cls=dsets,
                             data=classes[,c("MMxMM",
                                        "T:V")], cls.srt=ds.srt)

    ds.cut <- colnames(ovl$overlap)[apply(ovl$overlap,2,sum)>0]
    ovl <- sortOverlaps(ovl, axis=1, srt=ds.cut)
    
    omai <- c(2,1.5,.6,.6)
    nw <- ncol(ovl$p.value)*.35 + omai[2] + omai[4]
    nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

    plotdev(file.path(mfig.path,paste0("MMxMM_dissection")),
            height=nh, width=nw, res=300,  type=ftyp)
    par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
    plotOverlaps(ovl, ylab=NA, xlab=NA, show.total=TRUE)
    figlabel("count of unique BP/SAAP in Datasets", pos="bottom", font=2, cex=1.2)
    dev.off()
}


i=which(colnames(classes)=="fromto_I:P")
i=which(colnames(classes)=="KRAQ")

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
        boxplot(cdat$median ~ filt)

    

    ## PWM ENCODED
    ## position weight matrices and diffLogo

    aam1 <- aam[ filt,,drop=FALSE]
    aam2 <- aam[!filt,,drop=FALSE]

    if ( remove.duplicated.sequences ) {
        aam1 <- aam1[!duplicated(apply(aam1,1,paste,collapse="")),]
        aam2 <- aam2[!duplicated(apply(aam2,1,paste,collapse="")),]
    }

    ## filter tight context
    aam1 <- aam1[,cols]
    aam2 <- aam2[,cols]

    pfm1 <- getPFM(aam1)
    pfm2 <- getPFM(aam2)

    n1 <- nrow(aam1)
    n2 <- nrow(aam2)
    
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
    aas1 <- aas[ filt,,drop=FALSE]
    aas2 <- aas[!filt,,drop=FALSE]

    if ( remove.duplicated.sequences ) {
        aas1 <- aas1[!duplicated(apply(aas1,1,paste,collapse="")),]
        aas2 <- aas2[!duplicated(apply(aas2,1,paste,collapse="")),]
    }

    ## filter tight context
    aas1 <- aas1[,cols]
    aas2 <- aas2[,cols]

    pfn1 <- getPFM(aas1)
    pfn2 <- getPFM(aas2)

    n1 <- nrow(aas1)
    n2 <- nrow(aas2)
  
    dfnb <- createDiffLogoObject(pfn1, pfn2, alphabet=ASN2)
    dfnp <- enrichDiffLogoObjectWithPvalues(dfnb, n1=n1, n2=n2)

    ## suppress p.vals where ALL AA are equal
    ## -> obvious selection bias
    setpn <- apply(pfn1, 2, function(x) all(x%in%c(0.0,1.0)))
    ## manual pval suppression
    setpn[cols%in%npval[[id]]] <- TRUE
    dfnp$pvals[setpn] <- 1


    ## write out POSITION WEIGHT MATRICES
    ## for genome-wide scan
    fname <- file.path(seq.path,paste0(id,"pwm.tsv")) 
    colnames(pfm1) <- cols
    write.table(pfm1, file=fname, sep="\t", quote=FALSE, row.names=TRUE)
        
    ## PLOT LOGOS

    if ( length(grep(":",lb)) ) {
        ft <- unlist(strsplit(lb,":"))
        ft[ft=="x"] <- ""
           lb <- as.expression(bquote(.(ft[1]) %->% .(ft[2])))
    } else lb <- as.expression(bquote(.(lb)))
   

    ## ENCODED
    if ( any(dfop$pvals<min(psig)) | id %in% selected ) {

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
        figlabel(paste0("n=",nrow(aam1)), pos="topright", font=2)
        figlabel(lb, pos="topleft", cex=1.3, family=FONT)
        diffLogo_addPvals(dfop, ymin=mx, levels=psig)
        ##mtext("encoded", 4,-.25,adj=.05)
        dev.off()
    }

    ## INCORPORATED
    if ( any(dfnp$pvals<min(psig))  | id %in% selected ) {

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
        figlabel(paste0("n=",nrow(aas1)), pos="topright", font=2)
        figlabel(lb, pos="topleft", cex=1.3, family=FONT)
        diffLogo_addPvals(dfnp, ymin=mx, levels=psig)
        ##mtext("incorporated", 4,-.25,adj=.05)
        dev.off()
    }
    
    ## TIGHT LOGO PLOT FOR MAIN
    if ( (any(dfop$pvals<min(psig))  |
          any(dfnp$pvals<min(psig))) |
         id %in% selected ) {
        
        mno <- dfob$ylim.negMax
        mxo <- dfob$ylim.posMax
        mnn <- dfnb$ylim.negMax
        mxn <- dfnb$ylim.posMax
        
        ##mxo <- mxn <- max(mxo, mxn)
        
        ## figure heights
        ht <- 3
        omai <- c(.15,.3,.25,.05)
        nmai <- c(.25,.3,.15,.05)

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
        ## TODO: custom axis
        xat <- pretty(par("usr")[3:4])
        xlb <- c("0",xat[length(xat)-2])
        xat <- c(0,xat[length(xat)-2])
        axis(2, cex.axis=1.2, at=xat, labels=xlb)
        axis(2, at=pretty(par("usr")[3:4]), labels=FALSE)
        axis(1, at=1:length(cols), labels=FALSE)
        figlabel(paste0("n=",nrow(aam1),"  "),
                 pos="topright", font=1, cex=1.3)
        figlabel(lb, pos="topleft", cex=1.3, family=FONT, font=2)
        diffLogo_addPvals(dfop, ymin=mxo, levels=psig)
        ##mtext("encoded", 4,-.25,adj=.05)
        
        ## incorporated
        par(mai=nmai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfnb, sparse=TRUE, ymin=0, ymax=mxn)
        xat <- pretty(par("usr")[3:4])
        xlb <- c("0",xat[length(xat)-1])
        xat <- c(0,xat[length(xat)-1])
        axis(2, cex.axis=1.2, at=xat, labels=xlb)
        axis(2, at=pretty(par("usr")[3:4]), labels=FALSE)
        ##mtext("JS divergence", 2, 1.3)
        axis(1, at=1:length(cols), labels=axlab, cex.axis=1.2)
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
        height=hto+htn , width=.7, res=300, type=ftyp)
par(mai=ymai)
plot(1,1, axes=FALSE, col=NA, col.axis=NA)
text(1.25,1,label="JS divergence", xpd=TRUE, srt=90, cex=1.5)
rect(xleft=.6, ybottom=1.02, xright=1.05, ytop=par("usr")[2],
     xpd=TRUE, col="darkgray")
text(0.8,1.22,label="encoded", xpd=TRUE,srt=90, cex=1.6,
     col="white", font=2)
rect(xleft=.6, ybottom=par("usr")[3],
     xright=1.05, ytop=1.02, xpd=TRUE, col="darkgray")
text(0.82,.79,label="incorporated", xpd=TRUE,srt=90, cex=1.6, 
     col="white", font=2)
dev.off()



plotdev(file.path(mfig.path,"protein_location_G"), 
        width=3, height=3, res=300, type=ftyp)
hist(cdat$rpos[cdat$to=="G"])
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

## TODO: define motifs in cdat and map to tmtf
## MAP protein level info to TMT Data

do.unique <- FALSE ## ONLY FOR TESTING

tmtm <- tmtf
value <- "RAAS"

if ( do.unique ) {
    tmtm <- tmtf[!duplicated(tmtf$unique), ]
    value <- "RAAS.median"

    ## NOTE: no Dataset - only works for non Dataset profiles
    ##tmtm <- cdat
    ##value <- "median"
}
    
### MAP MOTIF TABLE w unique protein sites TO FULL RAAS TABLE
## THIS is used below for RAAS profiles
TIDX <- match(paste(tmtm$ensembl, tmtm$pos, tmtm$fromto),
              paste(cdat$ensembl, cdat$pos, cdat$fromto))

## is anyone missing
ina <- which(is.na(TIDX))
if ( length(ina)>0 ) 
    stop("TODO:", length(ina), "missing from unique site table.\n")



### CALCULATE ALL MOTIF RAAS PROFILES
## NOTE: obsolete since raasProfile can now handle class tables
## instead of a single "rows" vector in input data, but keeping
## as a reminder to move mergeOverlaps to segmenTools
if  ( FALSE ) {
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
                           bg=TRUE, value=value, row.srt=mid,
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
                     gcols=gcols, ffam=FONT)
        }
    }

### MERGE PROFILES - TODO: integrate this in segmenTools
    ovm <- mergeProfiles(covw)
}

## use RAAS profile with matrix input (new 20240604)
if ( "Dataset"%in%colnames(tmtm) ) {
    
    tcls <- classes[TIDX,]
    ovm <- raasProfile(x=tmtm, id="SAAP", 
                       rows=tcls, cols="Dataset",
                       bg=TRUE, value=value, 
                       col.srt=uds,
                       use.test=use.test, do.plots=FALSE,
                       xlab=xl.raas,
                       verb=0)
    
    
    ## sort
    ovm <- sortOverlaps(ovm, axis=2, srt=selected)
    
    amai <- c(0.8,1.4,0.6,.6)
    
    plotProfiles(ovm,
                 fname=file.path(mfig.path,paste0("motifs_",SETID,"_all")),
                 mai=amai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=mp.dot,
                 rlab=LAB,  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT)
    
    
    plotProfiles(ovm,
                 fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_pval10")),
                 mai=amai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=1e-10,
                 llab="p10",  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT,plot.all=FALSE)
    plotProfiles(ovm,
                 fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_pval6")),
                 mai=amai, ttcols=ttcols, value="median",
                 dot.sze=dot.sze, p.dot=1e-6,
                 llab="p6",  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT,plot.all=FALSE)
    plotProfiles(ovm,
                 fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_pval100")),
                 mai=amai, ttcols=ttcols, value="median",
                 p.min=1e-100, p.txt=1e-50,
                 dot.sze=dot.sze, p.dot=1e-100,
                 llab="p100",  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT,plot.all=FALSE)
    plotProfiles(ovm,
                 fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_pval30")),
                 mai=amai, ttcols=ttcols, value="median",
                 p.min=1e-30, p.txt=1e-15,
                 dot.sze=dot.sze, p.dot=1e-30,
                 llab="p30",  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, ffam=FONT,plot.all=FALSE)
    
    plotProfiles(ovm,
                 fname=file.path(mfig.path,paste0("motifs_",SETID,"_all_vcols")),
                 mai=amai, ttcols=ttcols, value="median",
                 p.min=mp.min, p.txt=mp.txt,
                 dot.sze=dot.sze, p.dot=mp.dot,
                 llab="vcols",  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols, ffam=FONT)


    ## SELECTED MOTIS FOR MAIN
    msrt <- c(##"Q:G",
        "KRAQ",
        "CCxCC",
        "MMxMM",
        "WWxWW")
    
    ovs <- sortOverlaps(ovm, axis=2, srt=msrt)
    plotProfiles(ovs,
                 fname=file.path(mfig.path,paste0("motifs_",SETID,"")),
                 mai=c(0.75,CMAIL,0.05,.5), ttcols=ttcols, value="median",
                 p.min=mp.min, p.txt=mp.txt,
                 dot.sze=dot.sze, p.dot=mp.dot,
                 rlab=LAB,  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 bg=NA, 
                 gcols=gcols, ffam=FONT)
    
    ssrt <- c("disord.",
              "binding",
              "noncons.",
              "long")
    ovs <- sortOverlaps(ovm, axis=2, srt=ssrt)
    plotProfiles(ovs,
                 fname=file.path(mfig.path,paste0("structure_",SETID,"")),
                 mai=c(0.05,CMAIL,0.05,.5), ttcols=ttcols, value="median",
                 p.min=mp.min, p.txt=mp.txt,
                 dot.sze=dot.sze, p.dot=mp.dot,
                 rlab=LAB,  ftyp=ftyp,
                 mtxt="", mtxt.line=2.3,
                 vcols=acols, vbrks=abrks,
                 bg=NA, 
                 gcols=gcols, ffam=FONT)
}


### ANALYZE OVERLAP MOTIFS vs STRUCTURE


CLS <- classes[,c(ssrt, msrt)]

## CORRELATION OF MEASURES
if ( interactive() ) {
    image_matrix(cor(classes), col=ttcols,
                 breaks=seq(-1,1,length=length(ttcols)+1), axis=1:2)

    ovl <- clusterCluster(CLS[,"disord."], CLS[,"noncons."])
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt)
}

rbins <- cut(cdat$median, breaks=c(-6,-4,-2,-1,0,3), include.lowest = TRUE)
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
        height=nh, width=nw, res=300,  type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, xlab="", ylab="",
             show.total=TRUE)
figlabel("median RAAS",pos="bottomright")
dev.off()


## for all non-conserved bins
csrt <- as.character(CLS[,"noncons."])
csrt[is.na(csrt)] <- "na."

rbins <- cut(cdat$median, breaks=c(-6,-4,-2,-1,0,3), include.lowest = TRUE)
rsrt <- levels(rbins)
rbins <- as.character(rbins)
names(rbins) <- rownames(CLS)
rbins[is.na(rbins)] <- "na."

ovl <- clusterAnnotation(cls=cdat$MMSeq2.bins, data=CLS,
                         cls.srt=levels(MMSeq2.bins))
omai <- c(.75,1,.6,.6)
nw <- ncol(ovl$p.value)*.35 + omai[2] + omai[4]
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_conservation_motifs_overlap")),
        height=nh, width=nw, res=300,  type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, xlab="", ylab="",
             show.total=TRUE)
mtext("conservation",1, 2.5)
dev.off()



### ALTERNATIVE RAAS PROFILES

## SELECTED MOTIFS
rows <- c(msrt,ssrt)
axex <- rep("",length(rows))
names(axex) <- rows
for ( i in seq_along(rows) ) {
    if ( length(grep(":",rows[i])) ) {
        ft <- unlist(strsplit(rows[i],":"))
        ft[ft=="x"] <- ""
        axex[i] <- as.expression(bquote(.(ft[1]) %->% .(ft[2])))
    } else axex[i] <- as.expression(bquote(.(rows[i])))
}


ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="loglen.bins", cols="Dataset",
                    row.srt=c(rev(levels(loglen.bins)),"na"),
                    col.srt=uds,
                    bg=TRUE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.8,CMAIL,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_loglen_raas")),
        height=nh, width=nw, res=300, type=ftyp, bg=NA)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8, axis=1)
polygon(y=c(1, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("protein length",2, 1.3)
dev.off()

ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="iupred3.bins", cols="loglen.bins",
                    row.srt=c(rev(levels(iupred3.bins)),"na"),
                    col.srt=c(levels(loglen.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,CMAIL,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_loglen_disorder_raas")),
        height=nh, width=nw, res=300, type=ftyp, bg=NA)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8)
polygon(y=c(2, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("disorder",2, 1.3)
polygon(x=c(1, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("protein length",1, 1.3)
dev.off()

ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="MMSeq2.bins", cols="loglen.bins",
                    row.srt=c(rev(levels(MMSeq2.bins)),"na"),
                    col.srt=c(levels(loglen.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,CMAIL,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_loglen_conservation_raas")),
        height=nh, width=nw, res=300, type=ftyp, bg=NA)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8)
polygon(y=c(2, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("conservation",2, 1.3)
polygon(x=c(1, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("protein length",1, 1.3)
dev.off()

ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="iupred3.bins", cols="MMSeq2.bins",
                    row.srt=c(rev(levels(iupred3.bins)),"na"),
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,CMAIL,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_conservation_disorder_raas")),
        height=nh, width=nw, res=300, type=ftyp, bg=NA)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8)
polygon(y=c(2, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("disorder",2, 1.3)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()

## hypergeo overlap in unique BP/SAAP
if ( sum(duplicated(paste(cdat$BP,cdat$SAAP)))>0 )
    stop("non-unique BP/SAAP tuple at critical point")
ovl <- clusterCluster(cdat$iupred3.bins, cdat$MMSeq2.bins,
                      cl1.srt=c(rev(levels(iupred3.bins)), "na"),
                      cl2.srt=c("na", levels(MMSeq2.bins)),
                      alternative="two.sided")
omai <- c(.5,.5,.5,.5)
nw <- ncol(ovl$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]
plotdev(file.path(mfig.path,paste0("classes_conservation_disorder_overlap")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, xlab=NA, ylab=NA, axis=NA,
             show.total=TRUE, text.cex=.5)
polygon(y=c(2, nrow(ovl$p.value), nrow(ovl$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("disorder",2, 1.5)
polygon(x=c(2, ncol(ovl$p.value), ncol(ovl$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()

## TODO: move this to structure and only do on unique sites and median RAAAS
## full distributions for conserved and disordred


library(vioplot)

plotdev(file.path(mfig.path,
                  paste0("classes_conservation_raas_violin")),
        height=2, width=3, res=300, type=ftyp)
par(mfrow=c(1,1), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, family=FONT)
levs <- c("na",levels(MMSeq2.bins))
vioplot(cdat$RAAS ~ factor(cdat$MMSeq2.bins, levels=levs),
        ylab=xl.raas, xlab="conservation", axes=FALSE, 
        col.axis="#ffffff", col.ticks="#ffffff")

## determine position of axis triangle
## reused below for all violin plots
yds <- diff(par("usr")[3:4])/10
ydf <- par("usr")[3] - abs(yds)
ypos <- c(ydf, ydf-(yds/2), ydf+(yds/2))

axis(2)
axis(1, at=1, label="na", las=2)
polygon(x=c(2, length(levs), length(levs)),
        y=ypos, xpd=TRUE, col="#aaaaaa", border=1)
dev.off()

plotdev(file.path(mfig.path,
                  paste0("classes_disorder_raas_violin")),
        height=2, width=3, res=300, type=ftyp)
par(mfrow=c(1,1), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, family=FONT)
levs <- c("na",levels(iupred3.bins))
vp <- vioplot(cdat$RAAS ~ factor(cdat$iupred3.bins, levels=levs),
        ylab=xl.raas, xlab="disorder", axes=FALSE, col.axis="#ffffff")
polygon(x=c(2, length(levs), length(levs)),
        y=ypos, xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
axis(2)
dev.off()

plotdev(file.path(mfig.path,
                  paste0("classes_loglen_raas_violin")),
        height=2, width=3, res=300, type=ftyp)
par(mfrow=c(1,1), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, family=FONT)
fct <- factor(cdat$loglen.bins)
levs <- c(levels(fct))
vp <- vioplot(cdat$RAAS ~ fct,
              ylab=xl.raas, xlab=expression(log[10](length)),
              axes=FALSE, col.axis="#ffffff")
polygon(x=c(1, length(levs), length(levs)),
        y=ypos, xpd=TRUE, col="#aaaaaa", border=1)
##axis(1, at=1, label="na", las=2)
axis(2)
dev.off()

## TODO: move this to structure?

stnd <- function(x,na.rm=TRUE) {
    (x-min(x, na.rm=na.rm))/(max(x, na.rm=na.rm) - min(x, na.rm=na.rm)) 
}

plotdev(file.path(mfig.path,paste0("peptides_conservation+disorder_raas")),
        height=3, width=3, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotCor(cdat$RAAS, stnd(-cdat$MMSeq2)+stnd(cdat$iupred3), xlab=xl.raas,
        ylab="disorder + -conservation", legpos="topright")
dev.off()

plotdev(file.path(mfig.path,paste0("peptides_conservationXdisorder_raas")),
        height=3, width=3, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotCor(cdat$RAAS, stnd(-cdat$MMSeq2)*stnd(cdat$iupred3), xlab=xl.raas,
        ylab="disorder x -conservation", legpos="topright")
dev.off()
plotdev(file.path(mfig.path,paste0("peptides_conservation_raas")),
        height=3, width=3, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotCor(cdat$RAAS, stnd(-cdat$MMSeq2), xlab=xl.raas,
        ylab="-conservation", legpos="topright")
dev.off()
plotdev(file.path(mfig.path,paste0("peptides_disorder_raas")),
        height=3, width=3, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotCor(cdat$RAAS, stnd(cdat$iupred3), xlab=xl.raas,
        ylab="disorder", legpos="topright")
dev.off()
plotdev(file.path(mfig.path,paste0("peptides_conservation_disorder")),
        height=3, width=3, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotCor(stnd(-cdat$MMSeq2), stnd(cdat$iupred3), xlab="conservation",
        ylab="disorder")
dev.off()

## TODO: test "marginal distribution" thingy, cons. vs. disordered
if ( interactive() ) {
    plotCor(cdat$RAAS, cdat$iupred3)
    plotCor(cdat$RAAS, cdat$MMSeq2)
    plotCor(cdat$RAAS*cdat$iupred3, cdat$RAAS*cdat$MMSeq2,
            xlab=expression(log[10](RAAS)~x~disorder),
            ylab=expression(log[10](RAAS)~x~conservation))
}

## MOTIFS VS CONS
ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows=tcls[,rows], cols="MMSeq2.bins",
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,1.1,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]
plotdev(file.path(mfig.path,paste0("classes_conservation_motifs_raas")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, axis=NA,
           xlab=NA, ylab=NA, show.total=TRUE)
mtext("conservation",1, 1.3)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
axis(2, length(axex):1, labels=axex, las=2, family=FONT,
     font=2, cex.axis=1.2)
dev.off()

ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows=tcls[,rows], cols="iupred3.bins",
                    col.srt=c("na",levels(iupred3.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,1.1,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]
plotdev(file.path(mfig.path,paste0("classes_disorder_motifs_raas")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, axis=NA,
           xlab=NA, ylab=NA, show.total=TRUE)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
            y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("disorder",1, 1.3)
axis(2, length(axex):1, labels=axex, las=2, family=FONT,
     font=2, cex.axis=1.2)
dev.off()


ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows=tcls[,rows], cols="loglen.bins",
                    col.srt=c(levels(loglen.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,1.1,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]
plotdev(file.path(mfig.path,paste0("classes_loglen_motifs_raas")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",vbrks=abrks, vcols=acols, axis=NA,
           xlab=NA, ylab=NA, show.total=TRUE)
polygon(x=c(1, ncol(ovls$p.value), ncol(ovls$p.value)),
            y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
##axis(1, at=1, label="na", las=2)
mtext(expression(log[10](length)),1, 1.3)
axis(2, length(axex):1, labels=axex, las=2, family=FONT,
     font=2, cex.axis=1.2)
dev.off()

### MUTUAL INFORMATION: conservation v disorder
## JOINT DISTRIBUTION PROBABILITY


source("~/work/mistrans/scripts/saap_utils.R")
diso <- raasProfile(x=tmtm, id="SAAP", 
                    rows="iupred3.bins", #cols="MMSeq2.bins",
                    row.srt=c(rev(levels(iupred3.bins)),"na"),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
plotdev(file.path(mfig.path,paste0("classes_disorder")),
        height=1.4, width=.8, res=300, type=ftyp, bg="white")
par(mai=c(.1,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(diso, vcols=acols, vbrks=abrks, value="median",
           axis=NA, xlab=NA, ylab=NA)
polygon(y=c(2, nrow(diso$p.value), nrow(diso$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("disorder",2, 1.5)
axis(2, at=1, label="na", las=2)
dev.off()

cons <- raasProfile(x=tmtm, id="SAAP", 
                    rows="MMSeq2.bins",
                    row.srt=c(rev(levels(MMSeq2.bins)),"na"),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
plotdev(file.path(mfig.path,paste0("classes_conservation")),
        height=1.4, width=.8, res=300, type=ftyp, bg="white")
par(mai=c(.1,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(cons, vcols=acols, vbrks=abrks, value="median",
           axis=NA, xlab=NA, ylab=NA)
polygon(y=c(2, nrow(cons$p.value), nrow(cons$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("conservation",2, 1.5)
axis(2, at=1, label="na", las=2)
dev.off()

## joint distribution, only one-sided p.value for high RAAS
disop <- diso$p.value[,1]
consp <- cons$p.value[,1]
disop[diso$sign<0] <- 1
consp[cons$sign<0] <- 1

jpv <- t(t(disop)) %*% t(consp)
jpv <- jpv[c(rev(levels(iupred3.bins)),"na"),]
jpv <- jpv[,c("na",levels(MMSeq2.bins))]

dison <- diso$p.value[,1]
consn <- cons$p.value[,1]
dison[diso$sign>0] <- 1
consn[cons$sign>0] <- 1

jpn <- t(t(dison)) %*% t(consn)
jpn <- jpn[c(rev(levels(iupred3.bins)),"na"),]
jpn <- jpn[,c("na",levels(MMSeq2.bins))]


pj.min <- 1e-100

## manual plot
jps <- jpv 
jps[jps<pj.min] <- pj.min
cols <- grDevices::gray(seq(1,0,length.out=100))
brks <- seq(0,-log10(pj.min),length.out=101)

omai <- c(.5,.5,.1,.1)
nw <- ncol(jps)*.2 + omai[2] + omai[4]
nh <- nrow(jps)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_conservation_disorder_joint")),
        height=nh, width=nw, res=300, type=ftyp, bg=NA)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
image_matrix(-log10(jps), col=cols, breaks=brks, axis=NA,
             xlab="", ylab="")
polygon(y=c(2, nrow(jps), nrow(jps)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
mtext("disorder",2, 1.5)
polygon(x=c(2, ncol(jps), ncol(jps)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()

## as dotplot

disocons <- raasProfile(x=tmtm, id="SAAP", 
                    rows="iupred3.bins", cols="MMSeq2.bins",
                    row.srt=c(rev(levels(iupred3.bins)),"na"),
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)

omai <- c(.5,.5,.5,.5)
nw <- ncol(disocons$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(disocons$p.value)*.2 + omai[1] + omai[3]


## high RAAS
tmp <- disocons$p.value
disocons$p.value <- jpv[rownames(tmp),colnames(tmp)]


plotdev(file.path(mfig.path,paste0("classes_conservation_disorder_joint_raas")),
        height=nh, width=nw, res=300, type=ftyp, bg="white")
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(disocons, value="median",vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8,
           p.dot=pj.min)
polygon(y=c(2, nrow(disocons$p.value), nrow(disocons$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("disorder",2, 1.3)
polygon(x=c(2, ncol(disocons$p.value), ncol(disocons$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()

## low RAAS
tmp <- disocons$p.value
disocons$p.value <- jpn[rownames(tmp),colnames(tmp)]


plotdev(file.path(mfig.path,paste0("classes_conservation_disorder_joint_raas_low")),
        height=nh, width=nw, res=300, type=ftyp, bg="white")
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(disocons, value="median",vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8,
           p.dot=pj.min)
polygon(y=c(2, nrow(disocons$p.value), nrow(disocons$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("disorder",2, 1.3)
polygon(x=c(2, ncol(disocons$p.value), ncol(disocons$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()


### POSITION WISE RAAS

## TODO: add to motifs?

## first AA in BP
bp1 <- unlist(lapply(strsplit(cdat$BP,""), function(x) x[1]))
bp1[bp1==cdat$from & cdat$site==1] <- "AAS"
##bp1[cdat$fromto=="Q:G"] <- "QG" ##- NOTE: this makes A at 1st insign.
##bp1[cdat$fromto=="Q:A"] <- "QA" ##- NOTE: this makes A at 1st insign.
## AA/codon/structure mapping
tmtm$bp1 <- bp1[TIDX]

ovw <- raasProfile(x=tmtm, id="SAAP", 
                   rows="bp1", cols="Dataset",
                   bg=TRUE, value=value, ##row.srt=motsrt,
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
bp2 <- unlist(lapply(strsplit(cdat$BP,""), function(x) x[2]))
bp2[bp2==cdat$from & cdat$site==2] <- "AAS"

## AA/codon/structure mapping
tmtm$bp2 <- bp2[TIDX]

ovw <- raasProfile(x=tmtm, id="SAAP", 
                   rows="bp2", cols="Dataset",
                   bg=TRUE, value=value, ##row.srt=motsrt,
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
             gcols=gcols, ffam=FONT)
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
bp3 <- unlist(lapply(strsplit(cdat$BP,""), function(x) x[3]))
bp3[bp3==cdat$from & cdat$site==3] <- "AAS"

## AA/codon/structure mapping
tmtm$bp3 <- bp3[TIDX]

ovw <- raasProfile(x=tmtm, id="SAAP", 
                   rows="bp3", cols="Dataset",
                   bg=TRUE, value=value, ##row.srt=motsrt,
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
             gcols=gcols, ffam=FONT)
ovwp <- sortOverlaps(ovw, p.min=mp.txt, cut=TRUE)
plotProfiles(ovwp, fname=file.path(mfig.path,paste0("bp3_",SETID,"_cut")),
             mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
             p.min=mp.min, p.txt=mp.txt,
             dot.sze=dot.sze, p.dot=mp.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="AA at pos. 3 in BP", mtxt.line=1.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)

