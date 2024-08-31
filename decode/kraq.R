
### ANALYZE SEQUENCE CONTEXT along BASE PEPTIDES
### and around splice sites

## TODO: splice sites, resolve 0: 5' or 3'?
SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source(file.path(SRC.PATH,"raas_init.R"))

## additionally required data files
mam.path <- file.path(Sys.getenv("MAMDATA")) 
if ( mam.path=="" ) # author's local path
    mam.path <- "/home/raim/data/mammary"
cdspos <- file.path(mam.path,"processedData","protein_cds_coordinates.tsv")

## local output path
kfig.path <- file.path(fig.path,"kraq")
dir.create(kfig.path, showWarnings=FALSE)

### remove duplicate sites!
cdat <- bdat[!duplicated(paste(bdat$ensembl,bdat$pos, bdat$fromto)),]

### POSITION OF AAS IN PEPTIDES

## Categorize by N+i and C-i
mx <- 9

nt <- cdat$site
ct <- nchar(cdat$BP)-cdat$site +1

## fuse central AAS
nt[nt>mx] <- ct[ct>mx] <- mx+1

## take smaller
at <- paste0("N",nt)
at[ct<nt] <- paste0("C",ct[ct<nt])

## re-name central AAS
at[at==paste0("N",mx+1)] <- paste0(">",mx)
at.srt <- c(
    paste0("N",1:mx),
    paste0(">",mx),
    paste0("C",mx:1))

asite.bins <- at
asite.srt <- at.srt

## calculate overlap enrichment by cumulative hypergeometric
## distribution tests
cls <- clusterCluster(cdat$fromto, asite.bins, cl2.srt=asite.srt)

## sort by only enriched
cls <- sortOverlaps(cls, p.min=1e-3)
clc <- sortOverlaps(cls, p.min=p.txt, cut=TRUE, sign=1)

plotdev(file.path(kfig.path,paste0("peptides_AAS_AAStype_tight")),
        height=.2*nrow(clc$p.value)+1.1, width=.25*nrow(clc$p.value)+1.25,
        res=300, type=ftyp)
par(mai=c(.5,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=2)
axex <- ftlabels(rownames(clc$p.value))
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("AAS type", 2, 2.7)
mtext("position of AAS in peptide", 1, 1.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()


### SPLICE SITE TEST

## TODO: distinguish splice site distance from N/C-terminal ends

## NOTE: unique sites for splice site distance distributions
## but using full sites X AAS for hypergeo test below
usite <- csite ## unique chrom.sites


plotCor(log10(abs(usite$ssd)), usite$RAAS.median)
plotCor(log10(abs(usite$ssd)[usite$fromto=="Q:G"]),
        usite$RAAS.median[usite$fromto=="Q:G"])


hist(usite$exl, breaks=seq(0,2e4,1), xlim=c(0,500))
hist(usite$ssd, breaks=seq(-2e4,2e4,1), xlim=c(-500,500))

## RELATIVE POS. in EXON
rpos <- abs(usite$ssd/usite$exl)


## distance to closest splice site vs. exon length
plotCor(log10(abs(usite$ssd)), log10(usite$exl),
        xlab=expression(log[10]("distance from s.s.")),
        ylab=expression(log[10]("exon length")))

plotCor(usite$ssd, usite$exl,
        xlab=expression("distance from s.s."),
        ylab=expression("exon length"))
## zoom in
ssd <- usite$ssd
exl <- usite$exl
ssd[ssd>  1000]<-  1000
ssd[ssd< -1000]<- -1000
exl[exl>  1000]<- 1000
exl[exl< -1000]<- -1000

## absolute positions
if ( interactive() ) {
    hist(exl,breaks=100, xlab="exon length") 
    hist(ssd,breaks=100, xlab="distance from s.s.")
}

## load all exons
cpos <- read.delim(cdspos, header=FALSE, sep=" ")
cpos <- cpos[cpos[,1] %in% genes$MANE.protein,]
colnames(cpos) <- c("ID","chr","start","end","strand")

axl <- cpos$end-cpos$start+1
axl[axl>1000] <- 1000
if ( interactive() ) {
    hist(axl,breaks=100, xlab="exon length", freq=FALSE)
    hist(exl,breaks=100, add=TRUE, border=2, col=NA, freq=FALSE) 
}


## relative position
plotdev(file.path(kfig.path,paste0("splicesites_relativepos_hist")),
        height=3.5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(rpos, breaks=20,
     xlab="relative distance of AAS to closest splice site", main=NA,
     ylim=c(0,800))
text(0,800, "no spatial bias of AAS within exons",pos=4)
arrows(x0=.1, x1=.4, y0=750, code=3, length=.05)
text(.5, 200, "closer to 5'/3' end\nthan to splice site", pos=4)
arrows(x0=.6, x1=.9, y0=150, code=3, length=.05)
dev.off()

plotCor(ssd, exl,
        xlab=expression("distance from s.s."),
        ylab=expression("exon length"))

plotCor(rpos, log10(usite$exl),
        xlab=expression("relative position in exon"),
        ylab=expression(log[10]("exon length")))

plotCor(rpos, usite$RAAS.median,
        xlab=expression("relative position in exon"),
        ylab=xl.raas)

hist(usite$exl, breaks=c(0:max(usite$exl, na.rm=TRUE)),
     xlim=c(0,1000))

col4 = c("#000000", "#DF536B", "#61D04F", "#2297E6",
         "#28E2E5", "#CD0BBC", "#EEC21F", "#9E9E9E")

ssd <- usite$ssd
## pos 1: 3' exon
## TODO: is ssdst==0 5' or 3'
ssd[ssd<=0] <- ssd[ssd<=0]-1
mxl <- 100
#ssd[ssd>  mxl] <-  mxl
#ssd[ssd< -mxl] <- -mxl
plotdev(file.path(kfig.path,paste0("splicesites_absolutepos_hist")),
        height=3.5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.25), mgp=c(1.3,.3,0), tcl=-.25)
#hist(usite$exl, breaks=2000,
#     xlim=1.5*c(-mxl,mxl),col="#aaaaaa", border=NA,axes=FALSE,
#     xlab=NA, ylab=NA, main=NA)
#shadowtext(x=150, y=350, labels="exon\nlengths", font=2, pos=4, col="#aaaaaa")
##abline(v=0)
#par(new=TRUE)
hist(ssd, breaks=seq(-1e5,1e5,1), xlim=c(-mxl,mxl),
     xlab="AAS (2nd cod.) distance from closest splice site", main=NA, col="#DF536B",border=NA)
##shadowtext(x=-100, y=50, labels="AAS", font=2, col="#ff0000")
axis(3, labels=FALSE)
text(40,80, labels="pref. splicing\nbetween codons", pos=4, font=2, xpd=TRUE)
dev.off()

## smooth and including 0
ssz <- usite$ssd 
##ssz[ssz<0] <-ssz[ssz<0]+1

plotdev(file.path(kfig.path,paste0("splicesites_absolutepos_hist_smooth")),
        height=3.5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.25), mgp=c(1.3,.3,0), tcl=-.25)
hist(ssz, breaks=seq(-1e5+2,1e5,12), xlim=c(-mxl,mxl), axes=FALSE,
     xlab=NA, ylab=NA, main=NA, col="#99999999", border=NA)
axis(4, col="#999999", col.axis="#999999")
par(new=TRUE)
hist(ssz, breaks=seq(-1e5+2,1e5,3), xlim=c(-mxl,mxl),
     xlab="AAS (2nd cod.) distance from closest splice site",
     main=NA, col="#DF536B",border=NA)
axis(3, labels=FALSE)
legend("topright", c("12 nt", " 3 nt"), col=c("#99999999", "#DF536B"),
       pch=15, pt.cex=1.5, bty="n") 
dev.off()
plotdev(file.path(kfig.path,paste0("splicesites_absolutepos_hist_annotated")),
        height=3.5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.25), mgp=c(1.3,.3,0), tcl=-.25)
hist(ssz, breaks=seq(-1e5+1,1e5,12), xlim=c(-mxl,mxl), axes=FALSE,
     xlab=NA, ylab=NA, main=NA, col="#99999999", border=NA)
axis(4, col="#999999", col.axis="#999999")
par(new=TRUE)
hist(ssz, breaks=seq(-1e5+1,1e5,3), xlim=c(-mxl,mxl),
     xlab="AAS (2nd cod.) distance from closest splice site",
     main=NA, col="#DF536B",border=NA)
axis(3, labels=FALSE)
legend("topright", c("12 nt", " 3 nt"), col=c("#99999999", "#DF536B"),
       pch=15, pt.cex=1.5, bty="n") 
##text(50,150, labels="with 0, and\nover 3 nt", pos=4, font=2)
text(-50,150, labels="N-terminal\nisoforms?", pos=2, font=2,
     col="#2297E6", xpd=TRUE)
arrows(x0=-50, x1=-15, y0=157, length=.05, lwd=4, col="#2297E6")
dev.off()

plotdev(file.path(kfig.path,paste0("splicesites_exonlength_hist")),
        height=3.5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(usite$exl, breaks=2000,
     xlim=c(0,3.5*mxl),
     col="#aaaaaa", border=NA,
     xlab="exon length", main=NA)
shadowtext(x=200, y=350, labels="exon\nlengths", font=2, pos=4, col="#aaaaaa")
##abline(v=0)
axis(3, labels=FALSE)
dev.off()

## ZOOM IN
## OVERLAP ENRICHMENT BY CUMULATIVE HYPERGEOMETRIC
## DISTRIBUTION TESTS

## NOTE: here using site X AAS
usite <- site ## incl. all distinct AAStypes, w. redundant chrom.sites

ssd <- usite$ssd

## pos 1: 3' exon, increase by 1, 0 is BETWEEN splice site
## TODO: is ssdst==0 5' or 3'
ssd[ssd<=0] <- ssd[ssd<=0]-1

mxl <- 6
ssd[ssd>  mxl] <- mxl
ssd[ssd< -mxl] <- -mxl
ssd[ssd== mxl] <- ">"
ssd[ssd==-mxl] <- "<"
ssd.srt <- c("<",(-mxl+1):(mxl-1), ">")
ssd.srt <- ssd.srt[ssd.srt%in%unique(ssd)]

## overlap type/exon distance


cls <- clusterCluster(usite$fromto, ssd, cl2.srt=ssd.srt)


## liberal: cut at 1e-2 (sort at 1e-3)
clc <- sortOverlaps(cls, p.min=1e-2, cut=TRUE, sign=1)
clc <- sortOverlaps(clc, p.min=1e-3)

plotdev(file.path(kfig.path,paste0("splicesites_AAStype")),
        height=.2*nrow(clc$p.value)+1.25, width=.25*ncol(clc$p.value)+1.5,
        res=300, type=ftyp)
par(mai=c(.75,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=1,
     cex.axis=.8)
axex <- ftlabels(rownames(clc$p.value))
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("AAS type", 2, 2.7)
mtext("pos. of AAS (2nd codon pos.)\nrelative to closest s.s.", 1, 2.5)
dev.off()

## TIGHTER, cut at 1e-3
clc <- sortOverlaps(cls, p.min=p.txt, cut=TRUE, sign=1)

plotdev(file.path(kfig.path,paste0("splicesites_AAStype_tight")),
        height=.2*nrow(clc$p.value)+1.25, width=.25*ncol(clc$p.value)+1.5,
        res=300, type=ftyp)
par(mai=c(.75,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=1,
     cex.axis=.8)
axex <- ftlabels(rownames(clc$p.value))
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("AAS type", 2, 2.7)
mtext("pos. of AAS (2nd codon pos.)\nrelative to closest s.s.", 1, 2.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

## NOTE: WHY IS OGDH gene, L->G MISSING, should be at -2
## ssd is 2, should be -2; strand mix-up?
## check PSMA1: ssd at -8, should be +8
## -> strand-mix up confirmed,fix in map_peptides.!
if ( FALSE ) {

    saap="GGFYGLDESDLDKVFHLPTTTFIGGQESALPLR"
    csite[csite$SAAP==saap,]
    csite[csite$name=="PSMA1" & csite$fromto=="Q:G;Q:A",c("SAAP","ssd","pos")]
}


tog <- grep(":G",rownames(cls$p.value), value=TRUE)
clc <- sortOverlaps(cls, axis=2, srt=tog, cut=TRUE)

plotdev(file.path(kfig.path,paste0("splicesites_AAStype_tight_toG")),
        height=.2*nrow(clc$p.value)+1.25, width=.25*ncol(clc$p.value)+1.5,
        res=300, type=ftyp)
par(mai=c(.75,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=1,
     cex.axis=.8)
axex <- ftlabels(rownames(clc$p.value))
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("AAS type", 2, 2.7)
mtext("pos. of AAS (2nd codon pos.)\nrelative to closest s.s.", 1, 2.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

toa <- grep(":A",rownames(cls$p.value), value=TRUE)
clc <- sortOverlaps(cls, axis=2, srt=toa, cut=TRUE)

plotdev(file.path(kfig.path,paste0("splicesites_AAStype_tight_toA")),
        height=.2*nrow(clc$p.value)+1.25, width=.25*ncol(clc$p.value)+1.5,
        res=300, type=ftyp)
par(mai=c(.75,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=1,
     cex.axis=.8)
axex <- ftlabels(rownames(clc$p.value))
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("AAS type", 2, 2.7)
mtext("pos. of AAS (2nd codon pos.)\nrelative to closest s.s.", 1, 2.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

## BY INCORPORATED vs ENCODED AA
## add to/from to site table
usite$to <-
    sapply(strsplit(usite$fromto, ":"), "[[", 2)
usite$from <-
    sapply(strsplit(usite$fromto, ":"), "[[", 1)

## by INCORPORATED
##usite$toAG <- usite$to
##usite$toAG[usite$to%in%c("A","G")] <- "AG"

cls <- clusterCluster(usite$to, ssd, cl2.srt=ssd.srt)
## sort by only enriched
clc <- sortOverlaps(cls, p.min=1e-2)

plotdev(file.path(kfig.path,paste0("splicesites_AAStype_to")),
        height=.2*nrow(clc$p.value)+1.25, width=.25*ncol(clc$p.value)+1.5,
        res=300, type=ftyp)
par(mai=c(.75,.5,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=1,
     cex.axis=.8)
axex <- rownames(clc$p.value)
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("Incorporated AA", 2, 1.3)
mtext("pos. of AAS (2nd codon pos.)\nrelative to closest s.s.", 1, 2.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

clc <- sortOverlaps(cls, p.min=1e-3, cut=TRUE, sign=1)
plotdev(file.path(kfig.path,paste0("splicesites_AAStype_to_tight")),
        height=.2*nrow(clc$p.value)+1.25, width=.25*ncol(clc$p.value)+1.5,
        res=300, type=ftyp)
par(mai=c(.75,.5,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=1,
     cex.axis=.8)
axex <- rownames(clc$p.value)
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("Incorp. AA", 2, 1.3)
mtext("pos. of AAS (2nd codon pos.)\nrelative to closest s.s.", 1, 2.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()


## by ENCODED
cls <- clusterCluster(usite$from, ssd, cl2.srt=ssd.srt)
## sort by only enriched
clc <- sortOverlaps(cls, p.min=1e-2)

plotdev(file.path(kfig.path,paste0("splicesites_AAStype_from")),
        height=.2*nrow(clc$p.value)+1.25, width=.25*ncol(clc$p.value)+1.5,
        res=300, type=ftyp)
par(mai=c(.75,.5,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=1,
     cex.axis=.8)
axex <- rownames(clc$p.value)
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("Encoded AA", 2, 1.3)
mtext("pos. of AAS (2nd codon pos.)\nrelative to closest s.s.", 1, 2.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

clc <- sortOverlaps(cls, p.min=1e-3, cut=TRUE, sign=1)
plotdev(file.path(kfig.path,paste0("splicesites_AAStype_from_tight")),
        height=.2*nrow(clc$p.value)+1.25, width=.25*ncol(clc$p.value)+1.5,
        res=300, type=ftyp)
par(mai=c(.75,.5,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=1,
     cex.axis=.8)
axex <- rownames(clc$p.value)
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("Encoded AA", 2, 1.3)
mtext("pos. of AAS (2nd codon pos.)\nrelative to closest s.s.", 1, 2.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

## D: GAT,GAC
## W: TGG

##

tmtm <- tmtf
rpos <- abs(tmtm$ssd/tmtm$exl)

ssd <- tmtf$ssd
ssd[ssd>  1e3] <-  1e3
ssd[ssd< -1e3] <- -1e3
hist(ssd, breaks=seq(-1.1e3,1.1e3,3), xlim=c(-300,300),
     xlab="distance from closest s.s.",)

ssd.bins <- cut(rpos, breaks=seq(0,.5,length.out=6),
                include.lowest = TRUE)
#levels(ssd.bins) <- c(levels(ssd.bins),"NA")
#ssd.bins[is.na(ssd.bins)] <- "NA"
tmtm$ssd.bins <- ssd.bins

tmtm$all <- "all"
ovw  <- raasProfile(x=tmtm, id="SAAP", value="RAAS",
                    delog=TRUE, rows="ssd.bins", cols="Dataset",
                    bg=TRUE, col.srt=uds,
                    use.test=use.test, do.plots=FALSE, verb=0)

mai <- c(.75,1,.5,.5)
nr <- nrow(ovw$p.value)
nc <- ncol(ovw$p.value)
nh <- nr *fh + mai[1] + mai[3]
nw <- nc *fw + mai[2] + mai[4]

plotdev(file.path(kfig.path,paste0("splicesites_relativepos_raas")),
        type=ftyp, res=300, width=nw,height=nh)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(x=ovw, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=1:2, xlab=NA, ylab=NA, show.total=TRUE)
mtext("relative pos. of s.s.", 2, 3.5)
dev.off()

