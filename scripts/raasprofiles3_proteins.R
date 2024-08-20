
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

pfig.path <- file.path(fig.path,"proteins")
dir.create(pfig.path, showWarnings=FALSE)

## file of windows with >2 AAS
window.file <- file.path(out.path, "aas_windows.tsv")

corW <- corH <- 2.5
pmai <- c(.5,.5,.25,.25)
pmpg <- c(1.3,.3,0)

## maximal SAAP/peptide, etc.
max.saap <- 12

## axis labels
xl.hlfm <- expression(protein~"half-life"/h)
xl.hlf <- expression(protein~"half-life"/h)

## overrule specific y-axis label
xl.prota <- xl.raas



### START ANALYSIS




### BP/SAAP/PROTEIN RAAS TABLES
## TODO: handle proteins/BP/SAAP missing from current tmt file!

## bp median raas
bpl <- split(tmtf$RAAS, tmtf$BP)
bpstat <- listProfile(bpl, y=tmtf$RAAS, use.test=use.test, min=3)
## add data
bpstat$bplen <- nchar(rownames(bpstat))
## count distinct SAAP per BP
bps <- unlist(lengths(lapply(split(bdat$SAAP, bdat$BP), unique)))
bpstat$nsaap <- bps[rownames(bpstat)]

## add median intensities per BP+SAAP
bsints <- tmtf$BP.abundance+tmtf$SAAP.abundance
bpint <- unlist(lapply(split(bsints, tmtf$BP), median, na.rm=TRUE))

bpstat$intensity <- bpint[rownames(bpstat)]

## protein median raas per protein w/o site-specific median first
ptl <- split(tmtf$RAAS, tmtf$ensembl)
ptstat <- listProfile(ptl, y=tmtf$RAAS, use.test=use.test, min=3)

## count distinct BP per  protein
pts <- unlist(lengths(lapply(split(bdat$BP, bdat$ensembl), unique)))
ptstat$nbp <- pts[rownames(ptstat)]
## count distinct SAAP per protein
pts <- unlist(lengths(lapply(split(bdat$SAAP, bdat$ensembl), unique)))
ptstat$nsaap <- pts[rownames(ptstat)]

## rank proteins by RAAS
ptstat$rank <- rank(ptstat$median)

## order site matrix by protein RAAS rank (for hotspot plot)
site$rank <- ptstat[site$ensembl,"rank"]
site <-site[order(site$rank, site$pos),]

if ( interactive() ) {
    ## test alternative measures of consistent RAAS
    ## CV:
    dense2d(site$RAAS.cv, -log10(site$RAAS.p.value),
            xlab="CV", ylab=expression(-log10(p)))
    dense2d(site$RAAS.cv, site$RAAS.sd, xlab="CV", ylab="SD")
    hist(site$RAAS.sd)
}






### COLLECT PROTEIN DATA

#### TODO: why do we have n=747 at exon but n=713 for
#### the log intensity plot on exon/intron?

## PROTEIN ABUNDANCE
ptstat$intensity <- pint[rownames(ptstat)]

## PROTEIN LENGTH
ptstat$length <- plen[rownames(ptstat)]

## THERMOSTABILITY @Savitski2014
ptstat$Tmelt <- pmlt[match(pnms[rownames(ptstat)], names(pmlt))]

## delta melting without ATP
## TODO: analyze this better
ptstat$DeltaTmelt <- pdmlt[match(pnms[rownames(ptstat)], names(pdmlt))]

## PROTEIN HALF-LIFE @Mathieson2018
ptstat$halflife <- phlv[match(pnms[rownames(ptstat)], names(phlv))]

## protstab2 - predicted 
ptstat$ProtStab2 <- protstab$Human_predict_Tm[match(rownames(ptstat),
                                                    protstab$ensembl)]

## GET WHOLE PROTEIN MEAN IUPRED3 SCORE (from saap_peptides4.tsv
iu3 <- split(hdat$iupred3.protein, hdat$ensembl)
## QC: all protein level
table(unlist(lapply(iu3, function(x) length(unique))))
iu3 <- unlist(lapply(iu3, unique))

ptstat$iupred3 <- iu3[rownames(ptstat)]

## p20 core proteasome in vitro targets
ptstat$p20_lgfc <- p20.lg2fc[rownames(ptstat)]
ptstat$p20_pval <- p20.pval[rownames(ptstat)]


## write-out collected protein results
write.table(cbind(protein=rownames(ptstat), ptstat),
            file=file.path(pfig.path, "proteins_raas.tsv"), sep="\t",
            na="", row.names=FALSE, quote=FALSE)

### HOTSPOTS

## TODO:
## * RAAS profiles by sliding window

pns <- ptstat$n
pns[pns>15] <- paste0(">", max.saap)
paas <- table(pns)[c(as.character(1:15),paste0(">", max.saap))]
plotdev(file.path(pfig.path,"hotspots_RAAS_per_protein"), type=ftyp,
        width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp <- barplot(paas, xlab="#RAAS per protein",
              ylab=NA, axes=FALSE, axisnames = FALSE)
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
axis(2)
mtext("# of proteins", 2, 1.6)
legend("topright",
       paste(sum(ptstat$n>1), "proteins with >1 RAAS"), bty="n")
dev.off()

pns <- ptstat$nbp
pns[pns>15] <- paste0(">", max.saap)
paas <- table(pns)[c(as.character(1:15),paste0(">", max.saap))]
plotdev(file.path(pfig.path,"hotspots_BP_per_protein"), type=ftyp,
        width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp <- barplot(paas, xlab="BP per protein",
              ylab=NA, axes=FALSE, axisnames = FALSE)
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
axis(2)
mtext("# of proteins", 2, 1.6)
legend("topright",
       paste(sum(ptstat$nbp>1,na.rm=TRUE), "proteins with >1 BP"), bty="n")
dev.off()

pns <- ptstat$nsaap
pns[pns>15] <- paste0(">", max.saap)
paas <- table(pns)[c(as.character(1:15),paste0(">", max.saap))]
plotdev(file.path(pfig.path,"hotspots_SAAP_per_protein"), type=ftyp,
        width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp <- barplot(paas, xlab="SAAP per protein",
              ylab=NA, axes=FALSE, axisnames = FALSE)
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
axis(2)
mtext("# of proteins", 2, 1.6)
legend("topright",
       paste(sum(ptstat$nsaap>1,na.rm=TRUE), "proteins with >1 SAAP"), bty="n")
dev.off()

plotdev(file.path(pfig.path,"hotspots_RAAS_per_protein_length"), type=ftyp,
        width=corW, height=corH, res=200)
par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
plotCor(log10(ptstat$length), log10(ptstat$n), axes=FALSE,
        ylab="#RAAS per protein",  xlab="protein length")
for ( ax in 1:2 ) {
    axis(ax, at=1:10, labels=10^(1:10))
    axis(ax, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125,
         labels=FALSE)
}
dev.off()


## RAAS DISTRIBUTION of UNIQUE SITES

site[which.max(site$RAAS.n),]
pid <- site[which.max(site$RAAS.n),"ensembl"]
pos <- site[which.max(site$RAAS.n),"pos"]
tidx <- which(tmtf$ensembl==pid & tmtf$pos==pos)

plotdev(file.path(pfig.path,"hotspots_example_site_RAAS"), type=ftyp,
        width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.1),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
hist(tmtf$RAAS[tidx], xlab=xl.raas, main=NA)
mtext(paste0(pnms[pid], ", ", unique(tmtf[tidx,"from"]),pos), 3,0, cex=.7)
dev.off()

for ( i in 1:30 ) {
    sid <- rownames(site)[order(site$RAAS.n,decreasing=TRUE)[i]]


    pid <- site[sid,"ensembl"]
    pos <- site[sid,"pos"]
    tidx <- which(tmtf$ensembl==pid & tmtf$pos==pos)
    plotdev(file.path(pfig.path,paste0("hotspots_example_site_RAAS_",i)),
            type=ftyp, width=3, height=3, res=200)
    par(mai=c(0.5,.5,.2,.1),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    hist(tmtf$RAAS[tidx], xlab=xl.raas, main=NA)
    mtext(paste0(pnms[pid], ", ", unique(tmtf[tidx,"from"]),pos), 3,0, cex=1)

    ##legend("right",  pnms[pid], bty="n")
    from <- unique(tmtf[tidx,"from"])
    to <- paste0(names(sort(table(tmtf[tidx,"to"]))),collapse=",")

    legend("topleft",  legend=c(bquote(.(from) %->% .(to)),
                                paste0("#RAAS=",site[sid,"RAAS.n"])), bty="n",
           seg.len=0, x.intersp=0)
    dev.off()
}

## RAAS DISTRIBUTIONS of SINGLE BP


pns <- bpstat$nsaap
pns[pns>max.saap] <- paste0(">",max.saap)
paas <- table(pns)[c(as.character(1:max.saap),paste0(">", max.saap))]

bid <- rownames(bpstat)[which(bpstat$nsaap>80)]
plotdev(file.path(pfig.path,"hotspots_example_RAAS"), type=ftyp,
        width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.1),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
hist(tmtf$RAAS[tmtf$BP==bid],
     xlab=xl.raas, main=NA)
mtext(bid, 3,0, cex=.7)
legend("topright",  unique(bdat$name[bdat$BP==bid]), bty="n")
dev.off()

for ( i in 1:30 ) {
    bid <- rownames(bpstat)[order(bpstat$nsaap,decreasing=TRUE)[i]]
    plotdev(file.path(pfig.path,paste0("hotspots_example_RAAS_",i)), type=ftyp,
        width=3, height=3, res=200)
    par(mai=c(0.5,.5,.15,.1),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    y <- tmtf$RAAS[tmtf$BP==bid]
    y[y>  0] <- .5
    y[y< -4] <- -4.5
    hist(y,
         xlab=xl.raas, main=NA, breaks=seq(-10,10,.5), xlim=c(-4.5,.5))
    mtext(bid, 3,0, cex=.7)
    legend("topright",  unique(bdat$name[bdat$BP==bid]), bty="n")
    legend("topleft",  c(paste0("#SAAP=",bpstat[bid,"nsaap"]),
                         paste0("#RAAS=",bpstat[bid,"n"])), bty="n",
           seg.len=0, x.intersp=0)
    dev.off()
}

plotdev(file.path(pfig.path,"hotspots_SAAP_per_peptide"), type=ftyp,
        width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp <- barplot(paas, xlab="SAAPs / base peptide",
              ylab=NA, axes=FALSE, axisnames = FALSE)
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
axis(2)
mtext("Number of base peptides", 2, 1.6)
legend("topright",
       paste(sum(bpstat$nsaap>1,nr.rm=TRUE), "base peptides\nwith >1 SAAP"), bty="n")
dev.off()

plotdev(file.path(pfig.path,paste0("hotspots_SAAP_per_peptide_log")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp <- barplot(paas, xlab="SAAPs / base peptide",ylab=NA,
              axes=FALSE, axisnames = FALSE, log="y", ylim=c(1,max(paas)))
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
options(scipen=999)
axis(2, at=1)
axis(2, at=10^(1:10))#, las=1)
axis(2, at=c(1:10,(1:10)*10,(1:10)*100, (1:10)*1000), tcl=-.125,
     labels=FALSE)
options(scipen=0)
mtext("Number of base peptides", 2, 1.3)
legend("topright",
       paste(sum(bpstat$nsaap>1,nr.rm=TRUE), "base peptides\nwith >1 SAAP"), bty="n")
dev.off()


pns <- bpstat$n
pns[pns>max.saap] <- paste0(">", max.saap)
paas <- table(pns)[c(as.character(1:max.saap),paste0(">", max.saap))]

plotdev(file.path(pfig.path,"hotspots_RAAS_per_peptide"), type=ftyp,
        width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp <- barplot(paas, xlab="#RAAS per base peptide",
              ylab=NA, axes=FALSE, axisnames = FALSE)
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
axis(2)
mtext("# of peptides", 2, 1.6)
legend("topright",
       paste(sum(bpstat$n>1,na.rm=TRUE), "peptides with >1 RAAS"), bty="n")
dev.off()

## TODO: fix
##plotdev(file.path(pfig.path,"hotspots_RAAS_per_peptide_length"), type=ftyp,
##        width=3, height=3, res=200)
##par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
##bp <- barplot(paas, xlab="#RAAS per base peptide",
##              ylab=NA, axes=FALSE, axisnames = FALSE)
##axis(1, at=bp, labels=names(paas), las=2, cex=.8)
##axis(2)
##mtext("# of peptides", 2, 1.6)
##legend("topright",
##       paste(sum(ubpn>1), "peptides with >1 AAS"), bty="n")
##dev.off()

## TODO: BP dotprofile


if ( interactive() )
    plot(log(ptstat$n), log(ptstat$length))



## list sites per protein

## filter here for only proteins>1, but not used!
multip <- split(site$pos, site$ensembl)
multip <- names(lengths(multip)[lengths(multip)> 0 ])
sites <- site[site$ensembl%in%multip,]

cpos <- cumsum(sites$pos) ## TODO: fix this!!!

### TODO: get correct cumulative position
## get cumulative position of each AAS for nice plot
sites$plen <- plen[sites$ensembl]
sites$clen[!duplicated(sites$ensembl)] <-
    cumsum(sites$plen[!duplicated(sites$ensembl)])
## fill up
nl <- NA
for ( i in 1:nrow(sites) ) {
    if ( !is.na(sites$clen[i]))
        nl <- sites$clen[i]
    else sites$clen[i] <- nl
}

## AAS and RAAS along all concatenated proteins

intv <- findInterval(sites$RAAS.median, vec=vbrks)
## raise to min/max
intv[intv==0] <- 1
intv[intv>length(vcols)] <- length(vcols)
raas.col <- vcols[intv]
plotdev(file.path(pfig.path,"hotspots_problem"), type=ftyp,
        width=10, height=3, res=200)
layout(t(t(1:3)), heights=c(.1,.4,.6))
par(mai=c(0.1,.5,.05,.1),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
## TODO: number of RAAS
par(mai=c(0,.5,.05,.1))
dense2d(cpos, rep(1, length(cpos)), axes=FALSE, ylab=NA, xlab=NA)
par(mai=c(0.1,.5,.01,.1))
plot(x=cpos, y=sites$RAAS.n, col=raas.col, pch=19, cex=.5,
     ylab="#RAAS per site",
     axes=FALSE, log="y", xpd=TRUE)
axis(2)
axis(3, at=cpos[sites$n==1], col=2, col.axis=2, labels=FALSE, tcl=.25)
par(mai=c(.35,.5,0.1,.1))
plot(cpos, sites$RAAS.median, type="h", col=raas.col, xpd=TRUE,
     xlab="pos. concatenated proteins with >1 AAS, sorted by protein RAAS",
     ylab=xl.site)
points(cpos, sites$RAAS.median, col=raas.col, pch=19, cex=.2, xpd=TRUE)
## TODO: protein lines
axis(3, at=cpos[sites$n==1], col=2, col.axis=2, labels=FALSE, tcl=.5)
pois <- c("AHNAK")
for ( poi in pois ) {
    gidx<- which(pnms[sites$ensembl]==poi)[1]
    text(cpos[gidx], 2, labels=poi, pos=4)
    arrows(x0=cpos[gidx], y0=2, y1=1, length=.05)
}

dev.off()


### MEDIAN RAAS FOR EACH UNIQUE PROTEIN POSITION



## PROTEIN MEDIAN RAAS 

plotdev(file.path(pfig.path,paste0("proteins_volcano_all")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=pmpg, tcl=-.25)
res <- volcano(ptstat, value="median",
               p.txt=12, v.txt=c(-2,-1), cut=50, mid=0,
               ids=pnms, xlab=xl.prota)
mtext("protein RAAS, median of all RAAS", 3,0)
dev.off()


#### COMPARE PROTEIN RAAS TO VARIOUS PROTEIN LEVEL MEASURES
## TODO: collect those values for all proteins e.g. in halflives,
## and load here


### PROTEIN INTENSITIES

plotdev(file.path(pfig.path,paste0("protein_intensities_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$intensity), ptstat$median, ylab=xl.raas,
        xlab=expression(log[10](intensity)), title=TRUE, cor.legend=FALSE)
dev.off()

plotdev(file.path(pfig.path,paste0("protein_intensities_length")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$intensity), log10(ptstat$length), ylab="protein length",
        xlab=expression(log[10](intensity)), title=TRUE, cor.legend=FALSE,
        axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
dev.off()

plotdev(file.path(pfig.path,paste0("protein_intensities_halflives")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$intensity), log10(ptstat$halflife),
        ylab=xl.hlfm, 
        xlab=expression(log[10](intensity)), title=TRUE, cor.legend=FALSE,
        axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
dev.off()

plotdev(file.path(pfig.path,paste0("protein_intensities_Tmelt")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$intensity), ptstat$Tmelt,
        ylab=expression(protein~melting~point/"°C"),
        xlab=expression(log[10](intensity)), title=TRUE, cor.legend=FALSE,
        axes=FALSE)
axis(1)
axis(2)
dev.off()

plotdev(file.path(pfig.path,paste0("protein_intensities_nsites")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$intensity), log10(ptstat$n),
        ylab="#RAAS per protein",
        xlab=expression(log[10](median~intensity)),
        title=TRUE, cor.legend=FALSE,
        axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
dev.off()

## NOTE: for hotspots
plotdev(file.path(pfig.path,paste0("peptide_intensities_nsites")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(bpstat$intensity), log10(bpstat$n),
        ylab="#RAAS per base peptide",
        xlab=expression(log[10](median~intensity["BP+SAAP"])),
        title=TRUE, cor.legend=FALSE,
        axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
dev.off()

plotdev(file.path(pfig.path,paste0("peptide_intensities_RAAS")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(bpstat$intensity), bpstat$median, 
        ylab=expression(log[10](median~RAAS)),
        xlab=expression(log[10](median~intensity["BP+SAAP"])),
        title=TRUE, cor.legend=FALSE,
        axes=FALSE)
axis(1)
axis(2)
dev.off()
plotdev(file.path(pfig.path,paste0("peptide_RAAS_nsites")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(bpstat$median, log10(bpstat$n), 
        xlab="peptide log10(median RAAS)",
        ylab="#RAAS per base peptide",
        title=TRUE, cor.legend=FALSE,
        axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
dev.off()


### PROTEIN LENGTH

plotdev(file.path(pfig.path,paste0("protein_lengths_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$length), ptstat$median,
        ylim=range(ptstat$median), axes=FALSE,
        ylab=xl.prota, xlab="protein length", legpos="bottomright",
        legbg="#ffffff77", title=TRUE, cor.legend=FALSE)
axis(2)
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
##box()
dev.off()

### AAS Density

plotdev(file.path(pfig.path,paste0("protein_nsites_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$n), ptstat$median, ylim=range(ptstat$median),
        ylab=xl.prota, xlab=expression("AAS per protein"),
        axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(2)
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
dev.off()

plotdev(file.path(pfig.path,paste0("protein_nsites_density_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$n/ptstat$length), ptstat$median, ylim=range(ptstat$median),
        ylab=xl.prota, xlab=expression(log[10](AAS/AA)),
        axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(2)
axis(1, at=-10:10, labels=10^(-10:10))
axis(1, at=log10(rep(1:10, 7) * 10^rep(-4:2, each=10)), tcl=-.125, labels=FALSE)
dev.off()


## PROTEIN MELTING TEMPERATURE


## melting point all
plotdev(file.path(pfig.path,paste0("protein_Tmelt_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(ptstat$Tmelt, ptstat$median, ylim=range(ptstat$median),
        ylab=xl.prota, xlab=expression(protein~melting~point/"°C"),
        axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(1)
axis(2)
##box()
dev.off()

## Diff melting point with ATP

## melting point difference with ATP
plotdev(file.path(pfig.path,paste0("protein_TmeltDelta_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(ptstat$DeltaTmelt, ptstat$median, ylim=range(ptstat$median),
        ylab=xl.prota,
        xlab=expression(melting~point~difference~Delta*T[ATP-vehicle]/"°C"),
        axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(1)
axis(2)
##box()
dev.off()


## PROTEIN THERMOSTABILITY 

## predicted stability 
plotdev(file.path(pfig.path,paste0("protein_protstab2_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(ptstat$ProtStab2, ptstat$median,
        ylim=range(ptstat$median),
        ylab=xl.prota, xlab="ProtStab2", axes=FALSE,
        title=TRUE, cor.legend=FALSE)
axis(1)
axis(2)
##box()
dev.off()

## TODO: thermostability and ATP/GTP from @Sridharan2019:
## Proteome-wide solubility and thermal stability profiling reveals
## distinct regulatory roles for ATP
## https://www.nature.com/articles/s41467-019-09107-y

## PROTEIN HALF-LIVES, @Mathieson2018

## halflives all
plotdev(file.path(pfig.path,paste0("protein_halflives_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$halflife), ptstat$median, ylim=range(ptstat$median),
        ylab=xl.prota, xlab=xl.hlfm, axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(2)
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
##box()
dev.off()

## loop through cell types
hlvd <- readxl::read_xlsx(math18.file)
for ( ctype in c("Bcells", "NK", "hepatocytes", "monocytes", "neurons",
                 "Mouse Neurons") ) {

    cidx <- grep("half_life", colnames(hlvd), value=TRUE)
    cidx <- cidx[grep(ctype, cidx, ignore.case=TRUE)]
    hlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
    names(hlv) <- unlist(hlvd[,1])
    

    plotdev(file.path(pfig.path,paste0("protein_halflives_all_",ctype)),
            type=ftyp, res=300, width=corW,height=corH)
    idx <- match(pnms[rownames(ptstat)], names(hlv))
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plotCor(ptstat$median, log10(hlv[idx]), signif = 1,
        xlab=xl.prota, ylab=xl.hlf, axes=FALSE)
    axis(1)
    axis(2, at=1:10, labels=10^(1:10))
    figlabel(ctype, pos="bottomleft")
    dev.off()
## halflives site
    plotdev(file.path(pfig.path,paste0("protein_halflives_lengths_",ctype)),
            type=ftyp, res=300, width=corW,height=corH)
    idx <- match(pnms[names(plen)], names(hlv))
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plotCor(log10(plen), log10(hlv[idx]),
            xlab="protein length", ylab=xl.hlf, axes=FALSE)
    axis(1, at=1:10, labels=10^(1:10))
    axis(2, at=1:10, labels=10^(1:10))
    ## TODO: log tick marks
    axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)),
         tcl=-.125, labels=FALSE)
    axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)),
         tcl=-.125, labels=FALSE)
    ##box()
    figlabel(ctype, pos="bottomleft")
    dev.off()
}

### PROTEIN RAAS vs. 20S Targets - @Pepelnjak2024
## TODO: p20 data is on peptide level with many duplicated prots
    

plotdev(file.path(pfig.path,paste0("protein_p20_lg2fc_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(ptstat$p20_lgfc, ptstat$median, ylim=range(ptstat$median),ylab=xl.prota,
        xlab=xl.20s, signif=2, round=2, legpos="bottomright",
        title=TRUE, cor.legend=FALSE)
dev.off()

plotdev(file.path(pfig.path,paste0("protein_p20_pval_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(ptstat$p20_pval, ptstat$median, ylim=range(ptstat$median),ylab=xl.prota,
        xlab=xl.20p, signif=2, round=2, title=TRUE, cor.legend=FALSE)
dev.off()


### PROTEIN RAAS vs. IUPRED


plotdev(file.path(pfig.path,paste0("protein_iupred3_all")),
        type=ftyp, res=300, width=corW, height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(ptstat$iupred3, ptstat$median, ylim=range(ptstat$median),
        ylab=xl.prota, xlab="disordered score, IUPred3",
        axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(1)
axis(2)
dev.off()


### POSITION IN PROTEIN

if ( FALSE ) {
    ## TODO: where was rpos added? makes no sense, since proteins can have
    ## multiple rpos
    plotdev(file.path(pfig.path,paste0("protein_rposition_all")),
            type=ftyp, res=300, width=corW,height=corH)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plotCor(ptstat$rpos, ptstat$median, ylim=range(ptstat$median),
        ylab=xl.prota, xlab=expression(log[10](AAS/AA)),
        axes=FALSE, title=TRUE, cor.legend=FALSE)
    axis(2)
    axis(1, at=-10:10, labels=10^(-10:10))
    axis(1, at=log10(rep(1:10, 7) * 10^rep(-4:2, each=10)),
         tcl=-.125, labels=FALSE)
    dev.off()
}

plotdev(file.path(pfig.path,paste0("sites_position_hist")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
hist(asite$rpos,
     xlab="relative position in proteins", main=NA)
dev.off()
plotdev(file.path(pfig.path,paste0("sites_position_raas")),
            type=ftyp, res=300, width=corW,height=corH)
par(mai=c(.5,.5,.25,.1), mgp=pmpg, tcl=-.25)
plotCor(asite$rpos, asite$RAAS.median, ylab=xl.raas,
        xlab="relative position in proteins",
        cor.legend=FALSE, title=TRUE, cex=.5)
dev.off()

### PAIRWISE DISTANCES of HIGH RAAS SITES

## nslavov at slack, 20240629:
## Thus, we focus on inferences that are less affected by those
## factors and by missing data. For example:

## * We can report the # of protein loci having more than 2 AAS within
##   10 amino acid segments. This of course underestimates the
##   prevalence of AAS (as the rest of our results) though it helps
##   describe the phenomena.
## * If we detect multiple low RAAS AAS from a segment but no high
##   RAAS AAS, we can conclude that likely (though not certainly) this
##   segment did not have high ratio AAS.


## TODO:
## * why so many small distances BETWEEN BP,
## * investigate overlapping BP; define group IDs
table(stringr::str_count(site$BP, ";"))
# 396 sites with more than one BP



bsite <- asite
bsite$RAAS <- bsite$RAAS.median
bsite$BP <- bsite$all.BP
filters <- list(all=1:nrow(bsite),
                high=which(bsite$RAAS.median> -1))
fnms <- c(all="all", high="RAAS> 0.1")

## TODO: asite-based
## * for asite-based distances, define BP groups,
## * asite-based: why within-BP at very long distances?

for ( j in length(filters):1 ) {

    nm <- names(filters)[j]
    rdat <- bsite[filters[[j]],]
    
    dsts <- matrix(NA, nrow=nrow(rdat), ncol=nrow(rdat))
    same <- dsts ## same basepeptide filter
    for ( k in 1:nrow(rdat) ) {

        dst <- abs(rdat$pos[k] - rdat$pos)
        dst[rdat$ensembl!=rdat$ensembl[k]] <- Inf
        dsts[k,] <- dst
        same[k,] <- rdat$BP[k] == rdat$BP
        
    }
    ## remove 0 distance
    diag(dsts) <-  NA
    dsts[which(dsts==0)] <- NA
    ## remove redundant
    dsts[upper.tri(dsts)] <- NA

    ## TODO: plot distances separately
    ## for same and distinct base peptides
    dsts.same <- dsts.other <- dsts
    dsts.other[same] <- NA
    dsts.same[!same] <- NA

    ## max. distance
    mxd <- max(c(dsts[is.finite(dsts)]))
    ## max. length BP
    bpmx <- max(nchar(rdat$BP))
    
    plotdev(file.path(pfig.path,paste0("AAS_distances_hist_",nm)),
            type=ftyp, res=300, width=3, height=3)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    hst <- hist(c(dsts), xlab="pairwise distances between AAS", breaks=0:mxd, 
         main=NA)
    if ( nm!="all" ) {
        par(xpd=TRUE)
        figlabel(fnms[nm], pos="topright", region="plot", xpd=TRUE)
    }
    dev.off()

    ## plots using above hists
    if ( interactive() )
        plot(hst$mids+.5, log10(hst$counts), type="l", xlim=c(0,500),
             xlab="pairwise distances between AAS",
             ylab="count")
    
 
    plotdev(file.path(pfig.path,paste0("AAS_distances_hist_zoom1_",nm)),
            type=ftyp, res=300, width=3, height=3)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    hist(c(dsts), xlab="pairwise distances between AAS",
         breaks=0:mxd, xlim=c(0,500), main=NA)
    text(x=bpmx, par("usr")[4]*.9, label="max. length BP", col=2, pos=4)
    
    hst.other <- hist(c(dsts.other), breaks=c(0:mxd), border=4, add=TRUE)
    hst.same <- hist(c(dsts.same), breaks=c(0:mxd), border=2, add=TRUE)

    abline(v=bpmx, col=2)
    if ( nm!="all" ) {
        par(xpd=TRUE)
        figlabel(fnms[nm], pos="topright", region="plot", xpd=TRUE)
    }
    dev.off()

    plotdev(file.path(pfig.path,paste0("AAS_distances_zoom1_",nm)),
            type=ftyp, res=300, width=3, height=3)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plot(hst$mids+.5, hst$counts, type="l", xlim=c(0,500),
         xlab="pairwise distances between AAS",
         ylab="count", col=NA)
    polygon(x=c(hst$mids+.5, rev(hst.same$mids)+.5),
            y=c(hst$counts,rev(hst$counts-hst.same$counts)),
            col="#ff0000", border=NA)
    polygon(x=c(hst.other$mids+.5, hst.other$mids[1]+.5),
            y=c(hst.other$counts,0),
            col="#0000ff77", border=NA)
    abline(v=bpmx, col=2)
    if ( nm!="all" ) 
        figlabel(fnms[nm], pos="topright")
    legend("topright", c("within BP","between BP"), lty=1,
           col=c(2,4))
    dev.off()

    plotdev(file.path(pfig.path,paste0("AAS_distances_zoom2_",nm)),
            type=ftyp, res=300, width=3, height=3)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plot(hst$mids+.5, hst$counts, type="l", xlim=c(0,100),
         xlab="pairwise distances between AAS",
         ylab="count", col=NA)

    polygon(x=c(hst$mids+.5, rev(hst.same$mids)+.5),
            y=c(hst$counts,rev(hst$counts-hst.same$counts)),
            col="#ff0000", border=NA)
    polygon(x=c(hst.other$mids+.5, hst.other$mids[1]+.5),
            y=c(hst.other$counts,0), col="#0000ff77", border=NA)
    if ( nm=="all" ) {
        hist(nchar(rdat$BP), breaks=0:bpmx, add=TRUE, border="yellow")
        legend("topright", c("within BP","between BP", "BP lengths"), lty=1,
               col=c(2,4,7))
    } 
    if ( nm!="all" )  {
        figlabel(fnms[nm], pos="topright")
        legend("topright", c("within BP","between BP"), lty=1,
               col=c(2,4))
    }
    dev.off()

    plotdev(file.path(pfig.path,paste0("AAS_distances_zoom3_",nm)),
            type=ftyp, res=300, width=3, height=3)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plot(hst$mids+.5, hst$counts, type="l", xlim=c(0,50),
         xlab="pairwise distances between AAS",
         ylab="count", col=NA)
    polygon(x=c(hst$mids+.5, rev(hst.same$mids)+.5),
            y=c(hst$counts,rev(hst$counts-hst.same$counts)),
            col="#ff0000", border=NA)
    polygon(x=c(hst.other$mids+.5, hst.other$mids[1]+.5),
            y=c(hst.other$counts,0), col="#0000ff77", border=NA)
    if ( nm=="all" ) {
        hist(nchar(rdat$BP), breaks=0:bpmx, add=TRUE, border="yellow")
        legend("topright", c("within-BP","between BP", "BP lengths"), lty=1,
               col=c(2,4,7))
    } 
    if ( nm!="all" )  {
        figlabel(fnms[nm], pos="topright")
        legend("topright", c("within-BP","between BP"), lty=1,
               col=c(2,4))
    }
    ##abline(v=c(3.7,10), lty=2)
    dev.off()
}


### AUTOCORRELATIONS 

## TODO: do we need to clean data for this?

## * concatenated proteins
## * concatenated BP.

N <- 150
tosc <- N/1:(N/2) # periods for fourier of ACF
freq <- c(0, 1/tosc)
trng <- which(tosc>2 & tosc<12)
     
## autocorrelation of of site median RAAS with concatenated protein sequences

## * initiate vector for concatenated proteins with median RAAS,
## * set actual RAAS values,
## * calculate auto-correlation.

## NOTE: asite, using really unique sites

osite <-asite[order(site$ensembl, site$pos),]
ositel <- split(osite, osite$ensembl)
##ositel <- ositel[unlist(lapply(ositel, nrow))>5]


## RAAS

## expand each protein to full sequence vector
## initialized to NA 
araas <- lapply(ositel, function(x) {
    aras <- rep(NA, plen[unique(x$ensembl)])
    aras[x$pos] <- x$RAAS.median # TODO: max(abs(median-RAAS))
    aras
})
araas <- unlist(araas)

## MANUAL autocorrelation due to NA
## TODO: different handling of NA in R's acf
acor <- rep(NA, N)
for ( k in 1:N ) {
    rng <- 1:(length(araas)-N)
    acor[k] <- cor(araas[rng], araas[rng+k], use="pairwise.complete")
}

## FOURIER of ACF
## TODO: fourier with correct period
afft <- fft(acor)[1:(N/2+1)]
dc <- afft[1]/N 
aamp <- abs(afft[1:length(afft)])/N


## compare to R's acf
if ( interactive() ) {

    bcor <- acf(araas, lag.max=N, na.action = na.pass)
    ##abline(v=10, col=2)
    
    plot(bcor$acf[1:N+1], type="l")
    lines(acor, col=2)
    plotCor(acor, bcor$acf[1:N+1])
    

}

## RAAS initialized to global median

## expand each protein to full sequence vector
## initialized to global median RAAS
global.raas <-  median(tmtf$RAAS, na.rm=TRUE)
init.raas <- global.raas # 
rraas <- lapply(ositel, function(x) {
    aras <- rep(init.raas, plen[unique(x$ensembl)])
    aras[x$pos] <- x$RAAS.median # TODO: max(abs(median-RAAS))
    aras
})
rraas <- unlist(rraas)

racf <- acf(rraas, lag.max=N, plot=FALSE)
rcor <- racf$acf[1:N+1]

## FOURIER of ACF
## TODO: fourier with correct period
rfft <- fft(rcor)[1:(N/2+1)]
dc <- rfft[1]/N 
ramp <- abs(rfft[1:length(rfft)])/N


## MEASUREMENT FREQUENCY

## expand each protein to full sequence vector
## initialized to 0 and add RAAS counts
nraas <- lapply(ositel, function(x) {
    aras <- rep(0, plen[unique(x$ensembl)])
    aras[x$pos] <- x$RAAS.n
    aras
})
nraas <- unlist(nraas)
nacf <- acf(nraas, lag.max=N, plot=FALSE)
ncor <- nacf$acf[1:N+1]

## FOURIER of ACF
## TODO: fourier with correct period
nfft <- fft(ncor)[1:(N/2+1)]
dc <- nfft[1]/N # TODO: compare with mean cor
namp <- abs(nfft[1:length(nfft)])/N


### NOTE: R's acf is equivalent to manual pearson correlation
### but behaves differently with NAs
if ( FALSE ) {

    ## calculate autocorrelation manually
    mcor <- rep(NA, N)
    for ( k in 1:N ) {
        rng <- 1:(length(nraas)-N)
        mcor[k] <- cor(nraas[rng], nraas[rng+k])
    }
    ##abline(v=10, col=2)
    plot(ncor, type="l")
    lines(mcor, col=2)
    plotCor(ncor, mcor)
}



## PLOT AUTO-CORRELATION AND DFT
plotdev(file.path(pfig.path,paste0("autocorrelation_protein_RAAS")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, acor, type="l", xlab="lag k", ylab=expression(C[RAAS]))
points(1:N, acor, cex=.25)
abline(h=aamp[1], col=2)
abline(h=0, col=1)
figlabel(paste("NA sites: NA"), pos="topright", region="plot") 
dev.off()

plotdev(file.path(pfig.path,paste0("autocorrelation_protein_RAAS_DFT_full")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(freq, aamp, type="h", 
     xlab="AA frequency", ylab="amplitude", axes=FALSE)
axis(1);axis(2)
lines(0, aamp[1], type="h",col=2)
text(0, aamp[1], "mean: DC component", pos=4, col=2, xpd=TRUE)
lines(1/tosc[1], aamp[2], type="h",col=4)
text(1/tosc[1], aamp[2], "trend: 1 cycle", pos=4, col=4)
lines(1/tosc[N/2], aamp[N/2+1], type="h",col=3)
text(1/tosc[N/2], aamp[N/2+1], "noise: Nyquist freq.",
     srt=90, pos=4, col=3, xpd=TRUE)
dev.off()

plotdev(file.path(pfig.path,paste0("autocorrelation_protein_RAAS_DFT")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
mxi <- which.max(aamp[trng+1])
mxa <- max(aamp[trng+1])
mxp <- tosc[trng][mxi]
plot(tosc[trng], aamp[trng+1], type="h", 
     ylim=c(0,mxa),
     xlab="AA period", ylab="amplitude", axes=FALSE)
axis(1);axis(2)
points(mxp, mxa, col=2)
text(mxp, mxa, labels=round(mxp,2), font=2, col=2, xpd=TRUE, pos=4)
figlabel(paste("NA sites: NA"), pos="topright", region="plot") 
dev.off()

## RAAS autocor, init. to global median RAAS
plotdev(file.path(pfig.path,paste0("autocorrelation_protein_globalRAAS")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, rcor, type="l", xlab="lag k", ylab=expression(C[RAAS]))
points(1:N, rcor, cex=.25)
abline(h=ramp[1], col=2)
abline(h=0, col=1)
figlabel(paste("NA sites:", round(init.raas,1)), pos="topright", region="plot") 
dev.off()

plotdev(file.path(pfig.path,paste0("autocorrelation_protein_globalRAAS_DFT_full")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(freq, ramp, type="h", 
     xlab="AA frequency", ylab="amplitude", axes=FALSE)
axis(1);axis(2)
lines(0, ramp[1], type="h",col=2)
text(0, ramp[1], "mean: DC component", pos=4, col=2, xpd=TRUE)
lines(1/tosc[1], ramp[2], type="h",col=4)
text(1/tosc[1], ramp[2], "trend: 1 cycle", pos=4, col=4)
lines(1/tosc[N/2], ramp[N/2+1], type="h",col=3)
text(1/tosc[N/2], ramp[N/2+1], "noise: Nyquist freq.",
     srt=90, pos=4, col=3, xpd=TRUE)
dev.off()

plotdev(file.path(pfig.path,paste0("autocorrelation_protein_globalRAAS_DFT")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
mxi <- which.max(ramp[trng+1])
mxa <- max(ramp[trng+1])
mxp <- tosc[trng][mxi]
plot(tosc[trng], ramp[trng+1], type="h", 
     ylim=c(0,mxa),
     xlab="AA period", ylab="amplitude", axes=FALSE)
axis(1);axis(2)
points(mxp, mxa, col=2)
text(mxp, mxa, labels=round(mxp,2), font=2, col=2, xpd=TRUE, pos=4)
figlabel(paste("NA sites:", round(init.raas,1)), pos="topright", region="plot") 
dev.off()

## AAS measurement frequency
plotdev(file.path(pfig.path,paste0("autocorrelation_protein_frequency")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, ncor, type="l", xlab="lag k", ylab=expression(C[frequency]))
points(1:N, ncor, cex=.25)
abline(h=namp[1], col=2)
dev.off()

plotdev(file.path(pfig.path,paste0("autocorrelation_protein_frequency_DFT_full")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(freq, namp, type="h", 
     xlab="AA frequency", ylab="amplitude", axes=FALSE)
axis(1);axis(2)
lines(0, namp[1], type="h",col=2)
text(0, namp[1], "mean: DC component", pos=4, col=2, xpd=TRUE)
lines(1/tosc[1], namp[2], type="h",col=4)
text(1/tosc[1], namp[2], "trend: 1 cycle", pos=4, col=4)
lines(1/tosc[N/2], namp[N/2+1], type="h",col=3)
text(1/tosc[N/2], namp[N/2+1], "noise: Nyquist freq.",
     srt=90, pos=4, col=3, xpd=TRUE)
dev.off()

plotdev(file.path(pfig.path,paste0("autocorrelation_protein_frequency_DFT")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
mxi <- which.max(namp[trng+1])
mxa <- max(namp[trng+1])
mxp <- tosc[trng][mxi]
plot(tosc[trng], namp[trng+1], type="h", 
     ylim=c(0,mxa),
     xlab="AA period", ylab="amplitude", axes=FALSE)
axis(1);axis(2)
points(mxp, mxa, col=2)
text(mxp, mxa, labels=round(mxp,2), font=2, col=2, xpd=TRUE, pos=4)
dev.off()

plotdev(file.path(pfig.path,
                  paste0("autocorrelation_protein_frequency_highRAAS")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, rcor, type="l", xlab="lag k", ylab=expression(C[frequency]))
points(1:N, rcor, cex=.25)
figlabel(expression(RAAS > 0.1), pos="topright", region="plot")
dev.off()



### COLLECT WINDOWS

## count windows in maxdst distance and plot length distribution
maxdst <- 10

## ALL WINDOWS - no RAAS filter

osite <-site[order(site$ensembl, site$pos),]

wsite <- split(osite, osite$ensembl) # split by protein

## use function in saap_utils to get windows
wsite <- lapply(wsite, get.windows, maxdst=maxdst)
wsite <- do.call(rbind, wsite)
rownames(wsite) <- rownames(osite)

wsite$windowID <- paste(wsite$ensembl, wsite$window, sep="-")

wcnt <- table(wsite$windowID)
wtxt <- wcnt
wtxt[wtxt>max.saap] <- paste0(">", max.saap)
wtxt <- table(wtxt)[c(as.character(1:max.saap),paste0(">", max.saap))]

plotdev(file.path(pfig.path,paste0("windows_AAS")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.15,.1), mgp=pmpg, tcl=-.25, xaxs="i")
bp <- barplot(wtxt, 
              xlab="AAS sites per window",ylab=NA,
              axes=FALSE, axisnames = FALSE)
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
axis(2)
mtext("# of windows",2,1.6)
legend("topright",
       paste(sum(wcnt>1,na.rm=TRUE), "windows with >1 AAS"), bty="n")
dev.off()

plotdev(file.path(pfig.path,paste0("windows_AAS_log")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.15,.1), mgp=pmpg, tcl=-.25, xaxs="i")
bp <- barplot(wtxt, 
              xlab="AAS sites per window",ylab=NA,
              axes=FALSE, axisnames = FALSE, log="y", ylim=c(1, max(wtxt)))
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
options(scipen=999)
axis(2, at=10^(1:10), las=2)
axis(2, at=c(1:10,(1:10)*10,(1:10)*100, (1:10)*1000), tcl=-.125,
     labels=FALSE)
options(scipen=0)
mtext("# of windows",2,1.6)
legend("topright",
       paste(sum(wcnt>1,na.rm=TRUE), "windows with >1 AAS"), bty="n")
dev.off()

## split by protein and window
wins <- lapply(split(wsite$pos, wsite$windowID), range)

## only realy windows (length>1)
wlen <- unlist(lapply(wins, diff))+1
wins <- wins[wlen>1]
wlen <- wlen[wlen>1]

wins <- as.data.frame(do.call(rbind, wins))
wins$ensembl <- unlist(lapply(strsplit(rownames(wins),"-"), function(x) x[1]))

write.table(window.file,
            x=wins[,c(3,1,2)],
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## WINDOW LENGTH DISTRIBUTION,

plotdev(file.path(pfig.path,paste0("windows_length")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
hist(nchar(unique(bdat$BP)), breaks=0:max(wlen), border=2, freq=FALSE,
     axes=FALSE, xlab=NA, ylab=NA, main=NA)
legend("topright","BP length", col=2, lty=1, bty="n")
par(new=TRUE)
hst <- hist(wlen, breaks=0:max(wlen), axes=TRUE, main=NA,
            xlab="window length", ylab="count")
dev.off()

plotdev(file.path(pfig.path,paste0("windows_length_zoom")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
hst <- hist(wlen, breaks=0:max(wlen), axes=TRUE, main=NA,
            xlab="window length", ylab="count", xlim=c(0,100))
dev.off()


plotdev(file.path(pfig.path,paste0("windows_length_log")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
hst <- hist(log10(wlen), axes=FALSE, main=NA,
            xlab="window lengths", ylab="count")
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125,
     labels=FALSE)
axis(2)
dev.off()


## ALL WINDOWS - WITH RAAS filter

osite <-site[order(site$ensembl, site$pos),]
osite <- osite[osite$RAAS.median > -1,] ## FILTER

wsite <- split(osite, osite$ensembl) # split by protein

## use function in saap_utils to get windows
wsite <- lapply(wsite, get.windows, maxdst=maxdst)
wsite <- do.call(rbind, wsite)
rownames(wsite) <- rownames(osite)

wsite$windowID <- paste(wsite$ensembl, wsite$window, sep="-")

wcnt <- table(wsite$windowID)
wtxt <- wcnt
wtxt[wtxt>max.saap] <- paste0(">", max.saap)
wtxt <- table(wtxt)[c(as.character(1:max.saap),paste0(">", max.saap))]

plotdev(file.path(pfig.path,paste0("windows_highRAAS_AAS")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.15,.1), mgp=pmpg, tcl=-.25, xaxs="i")
bp <- barplot(wtxt, 
              xlab="AAS sites per window",ylab=NA,
              axes=FALSE, axisnames = FALSE)
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
axis(2)
mtext("# of windows",2,1.6)
legend("topright",
       paste(sum(wcnt>1,na.rm=TRUE), "windows with >1 AAS"), bty="n")
dev.off()

plotdev(file.path(pfig.path,paste0("windows_AAS_highRAAS_log")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.25,.1), mgp=pmpg, tcl=-.25, xaxs="i")
bp <- barplot(wtxt, 
              xlab="AAS sites per window",ylab=NA,
              axes=FALSE, axisnames = FALSE, log="y")
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
options(scipen=999)
axis(2, at=10^(1:10), las=2)
axis(2, at=c(1:10,(1:10)*10,(1:10)*100, (1:10)*1000), tcl=-.125,
     labels=FALSE)
options(scipen=0)
mtext("# of windows",2,1.6)
legend("topright",
       paste(sum(wcnt>1,na.rm=TRUE), "windows with >1 AAS"), bty="n")
figlabel("RAAS>0.1", pos="topright", region="plot")
dev.off()

## split by protein and window
wins <- lapply(split(wsite$pos, wsite$windowID), range)

## only realy windows (length>1)
wlen <- unlist(lapply(wins, diff))+1
wins <- wins[wlen>1]
wlen <- wlen[wlen>1]

wins <- as.data.frame(do.call(rbind, wins))
wins$ensembl <- unlist(lapply(strsplit(rownames(wins),"-"), function(x) x[1]))

##write.table(window.file,
##            x=wins[,c(3,1,2)],
## sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## WINDOW LENGTH DISTRIBUTION,

bplen <- nchar(unique(bdat$BP)) 
wmax <- max(bplen, max(wlen))

plotdev(file.path(pfig.path,paste0("windows_highRAAS_length")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
hist(bplen, breaks=0:wmax, border=2, freq=FALSE,
     axes=FALSE, xlab=NA, ylab=NA, main=NA)
legend("topright","BP length", col=2, lty=1, bty="n")
par(new=TRUE)
hst <- hist(wlen, breaks=0:wmax, axes=TRUE, main=NA,
            xlab="window length", ylab="count")
dev.off()

plotdev(file.path(pfig.path,paste0("windows_highRAAS_length_zoom")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
hst <- hist(wlen, breaks=0:wmax, axes=TRUE, main=NA,
            xlab="window length", ylab="count", xlim=c(0,100))
dev.off()


plotdev(file.path(pfig.path,paste0("windows_highRAAS_length_log")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
hst <- hist(log10(wlen), axes=FALSE, main=NA,
            xlab="window lengths", ylab="count")
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125,
     labels=FALSE)
axis(2)
dev.off()

