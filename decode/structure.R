
## BP/SAAP-LEVEL CORRELATION OF CONSERVATION AND DISORDER

SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source(file.path(SRC.PATH, "raas_init.R"))

## local output path
sfig.path <- file.path(fig.path,"structure")
dir.create(sfig.path, showWarnings=FALSE)

### START ANALYSIS

tmtm <- tmtf

## add protein intensity

## via leading razor
hist(log10(tmtm$razor.intensity))
rint.bins <- cut(log10(tmtm$razor.intensity),
                 breaks=seq(6,14,length.out=6), include.lowest = TRUE)
levels(rint.bins) <- c("NA", levels(rint.bins))
rint.bins[is.na(rint.bins)] <- "NA"

tmtm$rint.bins <- rint.bins


## mapping back from total ensembl proteins
hist(log10(tmtm$protein.intensity))


pint.bins <- cut(log10(tmtm$protein.intensity),
                 breaks=seq(7,14,length.out=7), include.lowest = TRUE)
levels(pint.bins) <- c("NA", levels(pint.bins))
pint.bins[is.na(pint.bins)] <- "NA"

tmtm$pint.bins <- pint.bins

### RAAS dotplot: INTENSITY - CONSERVATION - DISORDER
## TODO: move this to structure.R

## NOTE: tmtf level leading razor protein intensity, not our remapping
ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="rint.bins", cols="MMSeq2.bins",
                    row.srt=c(rev(levels(rint.bins))),
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,.5,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(sfig.path,paste0("classes_conservation_intensity_raas_razor")),
        height=nh, width=nw, res=300, type=ftyp, bg="white")
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",
           p.dot=p.dot, dot.sze=dot.sze, vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8)#, axis=2)
polygon(y=c(2, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("razor intensity",2, 1.3)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()

## NOTE: ensembl protein level intensity
ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="pint.bins", cols="MMSeq2.bins",
                    row.srt=c(rev(levels(pint.bins))),
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,.5,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(sfig.path,paste0("classes_conservation_intensity_raas_protein")),
        height=nh, width=nw, res=300, type=ftyp, bg="white")
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",
           p.dot=p.dot, dot.sze=dot.sze, vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8)#, axis=2)
polygon(y=c(2, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("protein intensity",2, 1.3)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()

## NOTE: ensembl protein level intensity
ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="pint.bins", cols="iupred3.bins",
                    row.srt=c(rev(levels(pint.bins))),
                    col.srt=c("na", levels(iupred3.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,.5,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(sfig.path,paste0("classes_disorder_intensity_raas_protein")),
        height=nh, width=nw, res=300, type=ftyp, bg="white")
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",
           p.dot=p.dot, dot.sze=dot.sze, vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8)#, axis=2)
polygon(y=c(2, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("protein intensity",2, 1.3)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("disorder",1, 1.3)
dev.off()


ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="iupred3.bins", cols="MMSeq2.bins",
                    row.srt=c(rev(levels(iupred3.bins)),"white"),
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,.5,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(sfig.path,paste0("classes_conservation_disorder_raas")),
        height=nh, width=nw, res=300, type=ftyp, bg="white")
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",
           p.dot=p.dot, dot.sze=dot.sze, vbrks=abrks, vcols=acols, 
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





## UNIQUE SITE CORRELATION

fname <- file.path(sfig.path,paste0("structure_cor_MMSeq2_RAAS_sites"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(asite$MMSeq2, asite$RAAS.median,
        xlab="conservation, MMSeq2",
        ylab=xl.raas, legpos="topright")
dev.off()


fname <- file.path(sfig.path,paste0("structure_cor_iupred3_RAAS_sites"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(asite$iupred3, asite$RAAS.median,
        xlab="disordered score, iupred3",
        ylab=xl.raas, legpos="bottomleft")
dev.off()

## TODO: protein level measure, compare on protein level
fname <- file.path(sfig.path,paste0("structure_cor_intensity_RAAS_sites"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(log10(asite$protein.intensity), asite$RAAS.median,
        xlab="median protein intensity",
        ylab=xl.raas, legpos="bottomleft")
dev.off()

