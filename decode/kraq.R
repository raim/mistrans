
### ANALYZE SEQUENCE CONTEXT along BASE PEPTIDES

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
#if ( !exists("bdat") )
    source(file.path(SRC.PATH,"raas_init.R"))

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

usite <- site ## incl. all distinct AAStypes, w. redundant chrom.sites


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
hist(exl,breaks=100, xlab="exon length") 
hist(ssd,breaks=100, xlab="distance from s.s.")

## load all exons
cdspos <- file.path(mam.path,"processedData","protein_cds_coordinates.tsv")
cpos <- read.delim(cdspos, header=FALSE, sep=" ")
cpos <- cpos[cpos[,1] %in% genes$MANE.protein,]
colnames(cpos) <- c("ID","chr","start","end","strand")

axl <- cpos$end-cpos$start+1
axl[axl>1000] <- 1000
hist(axl,breaks=100, xlab="exon length", freq=FALSE)
hist(exl,breaks=100, add=TRUE, border=2, col=NA, freq=FALSE) 



## relative position
plotdev(file.path(kfig.path,paste0("splicesites_relativepos_hist")),
        height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(rpos, breaks=10, xlab="relative position in exon", main=NA)
dev.off()

plotCor(ssd, exl,
        xlab=expression("distance from s.s."),
        ylab=expression("exon length"))

plotCor(rpos, log10(usite$exl),
        xlab=expression("relative position in exon"),
        ylab=expression(log[10]("exon length")))

plotCor(rpos, usite$RAAS,
        xlab=expression("relative position in exon"),
        ylab=xl.raas)

hist(usite$exl, breaks=c(0:max(usite$exl, na.rm=TRUE)),
     xlim=c(-1000,1000))
hist(usite$ssd, breaks=min(usite$ssd, na.rm=TRUE):max(usite$ssd, na.rm=TRUE),
     xlab="distance from closest s.s.", add=TRUE, border="#ff000099")


ssd <- usite$ssd
mxl <- 200
ssd[ssd>  mxl] <-  mxl
ssd[ssd< -mxl] <- -mxl
hist(ssd, breaks=seq(-1.1e3,1.1e3,3), xlim=c(-mxl,mxl),
     xlab="distance from closest s.s.",)


## ZOOM IN
## calculate overlap enrichment by cumulative hypergeometric
## distribution tests
ssd <- usite$ssd
mxl <- 6
ssd[ssd>  mxl] <- mxl
ssd[ssd< -mxl] <- -mxl
ssd[ssd== mxl] <- ">"
ssd[ssd==-mxl] <- "<"
ssd.srt <- c("<",(-mxl+1):(mxl-1), ">")

## overlap type/exon distance
usite$to <-
    sapply(strsplit(usite$fromto, ":"), "[[", 2)
usite$from <-
    sapply(strsplit(usite$fromto, ":"), "[[", 1)

cls <- clusterCluster(usite$fromto, ssd, cl2.srt=ssd.srt)



## sort by only enriched
cls <- sortOverlaps(cls, p.min=1e-3)
clc <- sortOverlaps(cls, p.min=1e-3, cut=TRUE, sign=1)

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
#levels(ssd.bins) <- c("NA", levels(ssd.bins))
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

plotdev(file.path(kfig.path,paste0("splicesite_relativepos_raas")),
        type=ftyp, res=300, width=nw,height=nh)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(x=ovw, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=1:2, xlab=NA, ylab=NA, show.total=TRUE)
mtext("relative pos. of s.s.", 2, 3.5)
dev.off()
