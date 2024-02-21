
source("~/work/mistrans/scripts/saap_utils.R")

library(Biostrings) # for blosum62
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

addPoints <- function(ovl, value="median") {
    z <- t(apply(ovl[[value]], 2, rev))
    cx <- c(z)
    mx <- max(cx,na.rm=TRUE)
    mn <- min(cx,na.rm=TRUE)
    cx <- 3*(cx-mn)/(mx-mn)
    points(x = rep(1:ncol(z), nrow(z)),
           y = rep(nrow(z):1, each = ncol(z)), cex = cx)
}

## TODO:
## * which AA are missing from mapped file and why? likely
## because BP didnt match, generate from SAAP/BP,
## * global or local background distribution?

####!!!
## TODO 20240215:
## * rows: plot only largest contributor! Q->G, T->V,
##   TODO: find largest RAAS effect size.
## * RAAS distributions: show all and scale/color by pval
## * formal approach to missing data: map measured peptides to the
##   full peptide space (main peptides or all possible).
## * plot by protein class, e.g only Albumins.
## * CHECK HANDLING OF DUPLICATE SAAP here and in tmt retrieval


## FROM->TO BY AA PROPERTY CLASSES
## DONT FILTER FOR CODONS, use hdat

## AA PROPERTY CLASSES
## https://en.wikipedia.org/wiki/Amino_acid#/media/File:ProteinogenicAminoAcids.svg
AAPROP <- list(charged=c("R","H","K","D","E"),
               polar=c("S","T","N","Q"),
               special=c("C","U","G","P"),
               hydrophobic=c("A","V","I","L","M","F","Y","W"))
tmp <- unlist(AAPROP)
AAPROP <- sub("[0-9]+","", names(tmp))
names(AAPROP) <- tmp
aaprop <- sub("hydrophobic", "hphobic", AAPROP)

###  CODONS
aa <- unique(GENETIC_CODE)
CODONS <- rep("", length(aa))
for ( i in seq_along(aa) )
    CODONS[i] <- paste(names(which(GENETIC_CODE==aa[i])), collapse=";")
names(CODONS) <- aa
ACODONS <- paste0(names(CODONS),": ", CODONS)
names(ACODONS) <- aa

## 3<->1 letter code
AAC <- seqinr::aaa( seqinr::a())
names(AAC) <-  seqinr::a()
AAC <- AAC[AAC!="Stp"]

## plot colors
docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))
gcols <- grey.colors(n=100, start=.9, end=0)

#### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")

in.file <- file.path(out.path,"saap_mapped.tsv")
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")

### PARAMETERES

p.min <- 1e-10
p.txt <- 1e-5

exclude.albumin <- TRUE # FALSE #  
exclude.frequent <- FALSE # TRUE # 
frequent <- c("Q","W","T","S")
only.unique <- FALSE

LAB <- "all"
fig.path <- file.path(proj.path,"figures","saap_raasprofiles")
if (  exclude.albumin ) {
    fig.path <- paste0(fig.path,"_noalb")
    LAB <- "no Alb./Hemog."
}
if ( exclude.frequent ) {
    tmp <- paste(frequent,collapse=",")
    fig.path <- paste0(fig.path,"_", gsub(",","",tmp))
    LAB <- paste("no", tmp)
}
if ( only.unique ) {
    fig.path <- paste0(fig.path,"_uniqueSAAP")
    LAB <- paste("unique SAAP")
}
dir.create(fig.path, showWarnings=FALSE)

### START

## PARSE & FILTER DATA
dat <- read.delim(in.file)
tmtf <- read.delim(tmt.file)

## convert to logical
tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)
dat$Hemoglobin.Albumin <- as.logical(dat$Hemoglobin.Albumin)

#### TODO: 


## fuse all excluded tags for each unique SAAP
alls <- rbind(dat[,c("Keep.SAAP","SAAP")],
              tmtf[,c("Keep.SAAP","SAAP")])
alls <- split(alls$Keep.SAAP, alls$SAAP)

## inconsistents
##alls[which(unlist(lapply(alls, function(x) length(unique(x))))>1)]

## UNIFY FILTER COLUMNS
## remove if tagged so for any dataset
keep <- unlist(lapply(alls, all))
dat$keep <- keep[dat$SAAP]
tmtf$keep <- keep[tmtf$SAAP]

## remove excluded
hdat <- dat[which(dat$keep),]

## get raw RAAS data TMT level
## remove excluded
cat(paste("removing", sum(!tmtf$Keep.SAAP),
          "tagged as false positive on TMT level\n"))
tmtf <- tmtf[tmtf$keep,]
## exclude NA or Inf
rm <- is.na(tmtf$RAAS) | is.infinite(tmtf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from TMT level\n"))
tmtf <- tmtf[!rm,]

## CLEAN PROTEIN LEVEL
tmtsaap <- unique(tmtf$SAAP)
rm <- !hdat$SAAP%in%tmtsaap

cat(paste("removing", sum(rm), "SAAP without values at TMT level\n"))
hdat <- hdat[!rm,]


## first find mutated AA pairs
## (column AASin input mixes I/L)
## and split mutated AA pairs into from/to
## TODO: get this from the columns in the mapped file, once clear
## why missing
saaps <- strsplit(hdat$SAAP,"")
bases <- strsplit(hdat$BP, "")
fromtol <- lapply(1:length(saaps), function(i) {
    pos <- which(saaps[[i]]!=bases[[i]])
    c(from=bases[[i]][pos], to=saaps[[i]][pos])
})
fromto <- do.call(rbind, fromtol)

## replace from/to columns
hdat$from <- fromto[,1]
hdat$to <- fromto[,2]

## just useful for plots instead of raw codons
hdat$aacodon <- paste(hdat$from, hdat$codon, sep="-")

## fromto column
hdat$fromto <- apply(fromto, 1, paste0, collapse=":")

## generate columns by AA property
hdat$pfrom <- aaprop[fromto[,1]]
hdat$pto <- aaprop[fromto[,2]]
hdat$pfromto <- paste0(hdat$pfrom,":",hdat$pto)

hdat$frompto <- paste0(hdat$from, ":", hdat$pto)

## STRUCTURE: add bins for numeric values
## iupred3
iupred3.bins <- cut(hdat$iupred3, breaks=seq(0,1,.2))
rsrt <- levels(iupred3.bins)
hdat$iupred3.bins <- as.character(iupred3.bins)
hdat$iupred3.bins[is.na(hdat$iupred3.bins)] <- "na"
## anchor2
anchor2.bins <- cut(hdat$anchor2, breaks=seq(0,1,.2))
rsrt <- levels(anchor2.bins)
hdat$anchor2.bins <- as.character(anchor2.bins)
hdat$anchor2.bins[is.na(hdat$anchor2.bins)] <- "na"
## secondary structure by names
ssrt <- c(E="sheet", H="helix", C="coil")
hdat$sstruc <- ssrt[hdat$s4pred]
hdat$sstruc[hdat$sstruc==""] <- "na"

## MAP protein level info to TMT Data
idx <- match(tmtf$SAAP, hdat$SAAP)
tmtf$from <- hdat$from[idx]
tmtf$to <- hdat$to[idx]
tmtf$fromto <- hdat$fromto[idx]
tmtf$pfrom <- hdat$pfrom[idx]
tmtf$pto <- hdat$pto[idx]
tmtf$pfromto <- hdat$pfromto[idx]
tmtf$frompto <- hdat$frompto[idx]
tmtf$codon <- hdat$codon[idx]
tmtf$aacodon <- hdat$aacodon[idx]
tmtf$iupred3.bins <- hdat$iupred3.bins[idx]
tmtf$anchor2.bins <- hdat$anchor2.bins[idx]
tmtf$sstruc <- hdat$sstruc[idx]

tmtf$albumin <- hdat$Hemoglobin.Albumin[idx]

### FILTER
if (  exclude.albumin ) {
    rmsaap <- tmtf$SAAP[tmtf$albumin]
    tmtf <- tmtf[!tmtf$SAAP%in%rmsaap,]
}
if ( exclude.frequent ) {
    rmsaap <- tmtf$SAAP[tmtf$from%in%frequent]
    tmtf <- tmtf[!tmtf$SAAP%in%rmsaap,]
}
if ( only.unique ) {
    ## FOR TESTING
    ## just take the first of each SAAP per dataset
    tmtf <- tmtf[!duplicated(paste(tmtf$SAAP,tmtf$Dataset)),]
}

### TODO: merge unique SAAP here to test effects of multiple
### measurements.

## from this point on we only work with tmtf, which
## now should have all information

## sorting and labels
uds <- sort(unique(tmtf$Dataset))
uds <- c("Healthy",uds[uds!="Healthy"])

srt <- c("charged","polar","hphobic","special")
srt <- paste(rep(srt,each=4), rep(srt,4), sep=":")

axex <- ftlabels(srt) # axis labels with arrows

dpath <- file.path(fig.path,"dists")
dir.create(dpath, showWarnings=FALSE)
fname <- file.path(dpath,paste0("AAprop_wtests_"))
source("~/work/mistrans/scripts/saap_utils.R")
ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
                    delog=TRUE, rows="pfromto", cols="Dataset",
                    bg=TRUE, row.srt=srt, col.srt=uds,
                    use.test=t.test, do.plots=TRUE, xlab="TMT level RAAS",
                    fname=fname, verb=1)


## common histogram with color by sign and
## first, calculate densities
dns <- ovw$ALL
for ( i in 1:nrow(ovw$p.value) ) {
    for ( j in 1:ncol(ovw$p.value) ) {
        rw <- rownames(ovw$p.value)[i]
        cl <- colnames(ovw$p.value)[j]
        vls <- ovw$ALL[[rw]][[cl]]
        if (length(vls)>1 ) {
            dns[[rw]][[cl]] <- density(vls)
        }
    }
}
## get max
mxs <- rep(NA, length(dns))
for ( i in seq_along(dns) )
    mxs[i] <- max(unlist(lapply(dns[[i]],
                                function(z)
                                    ifelse("y"%in%names(z),
                                           max(unlist(z["y"])),-10))))

png(file.path(fig.path,paste0("AAprop_densities.png")),
    res=300, width=4, height=2, units="in")
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(tmtf$RAAS, col="#77777755", border=NA, 
     xlim=c(-6,4),
     xlab="TMT level RAAS", ylab=NA, main=NA, axes=FALSE)
dns <- ovw$ALL
par(new=TRUE)
plot(1, col=NA, xlim=c(-6,4), ylim=c(0,.75),#max(mxs)),
     xlab=NA, ylab=NA, axes=FALSE)
for ( i in 1:nrow(ovw$p.value) ) {
    for ( j in 1:ncol(ovw$p.value) ) {
        rw <- rownames(ovw$p.value)[i]
        cl <- colnames(ovw$p.value)[j]
        vls <- ovw$ALL[[rw]][[cl]]
        if (length(vls)>1 ) {
            dns[[rw]][[cl]] <- density(vls)
            if ( ovw$p.value[i,j] <1e-3 )
                lines(dns[[rw]][[cl]],
                      col=ifelse(ovw$sign[i,j]==-1,"blue","red"),
                      lwd=-log10(ovw$p.value[i,j])/-log10(p.min))
            ## TODO: annotate
            
        }# else abline(v=vls)
    }
}
axis(2, labels=NA)
mtext("density",2,.5)
axis(1)
dev.off()



png(file.path(fig.path,paste0("AAprop_wtests_legend.png")),
    res=300, width=2, height=2, units="in")
par(mai=c(.6,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlapsLegend(p.min=p.min, p.txt=p.txt, type=2, col=ttcols)
dev.off()

png(file.path(fig.path,paste0("AAprop_wtests.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1, ylab=NA, xlab="", col=ttcols, show.total=TRUE)
axis(2, length(axex):1, labels=axex, las=2)
##addPoints(ovw)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

## local raas colors
png(file.path(fig.path,"raas_profile_colors.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(c(ovw$mean),
                          q=0.01, colf=viridis::viridis,
                          n=10,
                          xlim=c(-4,1),
                          plot=TRUE,
                          mai=c(.5,.5,.1,.1), xlab="All TMT level RAAS")
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()


png(file.path(fig.path,paste0("AAprop_raas_mean.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="mean", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1, ylab=NA, xlab="", text.cex=.8)
axis(2, length(axex):1, labels=axex, las=2)
mtext("mean TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()


png(file.path(fig.path,paste0("AAprop_raas_median.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="median", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1, ylab=NA, xlab="", text.cex=.8)
axis(2, length(axex):1, labels=axex, las=2)
mtext("median TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

### median vs. p-value as a measure of effect size
plot(ovw$count, ovw$median, xlab="count", ylab="RAAS median")
plot(ovw$count, -log(ovw$p.value), xlab="count", ylab=expression(-log(p)))
plot(ovw$median, -log(ovw$p.value))

png(file.path(fig.path,paste0("AAprop_volcano.png")),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.75,.1,.75), mgp=c(1.3,.3,0), tcl=-.25)
volcano(ovw, cut=100, p.txt=5, v.txt=c(-Inf,-1), density=TRUE,
        xlab="mean TMT RAAS", value="mean")
abline(v=mean(ovw$mean,na.rm=TRUE))
dev.off()

png(file.path(fig.path,paste0("AAprop_volcano_mean.png")),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.75,.1,.75), mgp=c(1.3,.3,0), tcl=-.25)
volcano(ovw, cut=100, p.txt=5, v.txt=c(-Inf,-1), density=TRUE,
        xlab="mean TMT RAAS", value="mean")
abline(v=mean(ovw$mean,na.rm=TRUE))
dev.off()

png(file.path(fig.path,paste0("AAprop_count_unique.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
cnt <- ovw$unique
cnt[cnt==0] <- NA
txt <- ovw$unique
txt[txt==0] <- ""
txt.col <- ifelse(ovw$unique>quantile(c(ovw$unique),.95),"white","black")
image_matrix(cnt, col=gcols, axis=1:2,
             text=txt, text.col=txt.col, ylab="",
             xlab="", text.cex=.8)
##axis(4, nrow(ovw$p.value):1, labels=ACODONS[rownames(ovw$p.value)], las=2)
mtext("mistranslated", 2, 3.5, adj=-1.7, cex=1.2)
mtext("# unique SAAP", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()


ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="from", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)
png(file.path(fig.path,paste0("fromAA_wtests.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="",
             col=ttcols, show.total=TRUE)
mtext("from", 2, 1.3)
##addPoints(ovw)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="to", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)
png(file.path(fig.path,paste0("toAA_wtests.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="",
             col=ttcols, show.total=TRUE)
mtext("to", 2, 1.3)
##addPoints(ovw)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()


## plot 16 plots for pfrom-to combos, for each to all AA
for ( ptype in unique(tmtf$pfromto) ) {

    dtmt <- tmtf[tmtf$pfromto==ptype,]
    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="fromto", cols="Dataset",
                       bg=TRUE, value="RAAS", col.srt=uds,
                       use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                       verb=1)
    ## calculate optimal figure height: result fields + figure margins (mai)
    nh <- nrow(ovw$p.value) *.2 + 1.4
    nw <- ncol(ovw$p.value) *.4 + 1.15
    
    if ( nrow(ovw$p.value) < 5 ) {
        nh <- nh*2; nw <- nw*2
    }
    
    ## plot figure
      png(file.path(fig.path,paste0("type_",ptype, "_AA_wtests.png")),
          height=nh, width=nw, res=300, units="in")
    par(mai=c(1,.75,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:2, ylab=NA, xlab="",
                 col=ttcols, show.total=TRUE)
    mtext("from:to", 2, 2)
    figlabel(LAB, pos="bottomleft", cex=1)
    figlabel(ptype, pos="bottomright", cex=1.5)
    dev.off()
}

## plot 20 from plots for pfrom-to combos, for each to all AA
for ( ptype in unique(tmtf$to) ) {

    dtmt <- tmtf[tmtf$to==ptype,]
    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="fromto", cols="Dataset",
                       bg=TRUE, value="RAAS", col.srt=uds,
                       use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                       verb=0)
    ## calculate optimal figure height: result fields + figure margins (mai)
    nh <- nrow(ovw$p.value) *.2 + 1.4
    nw <- ncol(ovw$p.value) *.4 + 1.15
    
    if ( nrow(ovw$p.value) < 5 ) {
        nh <- nh*2; nw <- nw*2
    }
    
    ## plot figure
      png(file.path(fig.path,paste0("toAA_",ptype, "_AA_wtests.png")),
          height=nh, width=nw, res=300, units="in")
    par(mai=c(1,.75,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:2, ylab=NA, xlab="",
                 col=ttcols, show.total=TRUE)
    mtext("from:to", 2, 2)
    figlabel(LAB, pos="bottomleft", cex=1)
    figlabel(paste("to",ptype), pos="bottomright", cex=1.5)
    dev.off()
}

## plot 20 from plots for pfrom-to combos, for each to all AA
for ( ptype in unique(tmtf$from) ) {

    dtmt <- tmtf[tmtf$from==ptype,]
    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="fromto", cols="Dataset",
                       bg=TRUE, value="RAAS", col.srt=uds,
                       use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                       verb=0)
    ## calculate optimal figure height: result fields + figure margins (mai)
    nh <- nrow(ovw$p.value) *.2 + 1.4
    nw <- ncol(ovw$p.value) *.4 + 1.15
    
    if ( nrow(ovw$p.value) < 5 ) {
        nh <- nh*2; nw <- nw*2
    }
    
    ## plot figure
    png(file.path(fig.path,paste0("fromAA_",ptype, "_AA_wtests.png")),
        height=nh, width=nw, res=300, units="in")
    par(mai=c(1,.75,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:2, ylab=NA, xlab="",
                 col=ttcols, show.total=TRUE)
    mtext("from:to", 2, 2)
    figlabel(LAB, pos="bottomleft", cex=1)
    figlabel(paste("from",ptype), pos="bottomright", cex=1.5)
    dev.off()
}

## plot 20 from plots for AA->AAprop
for ( ptype in unique(tmtf$from) ) {

    dtmt <- tmtf[tmtf$from==ptype,]
    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="frompto", cols="Dataset",
                       bg=TRUE, value="RAAS", col.srt=uds,
                       use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                       verb=0)
    ## calculate optimal figure height: result fields + figure margins (mai)
    nh <- nrow(ovw$p.value) *.2 + 1.4
    nw <- ncol(ovw$p.value) *.4 + 1.4
    
    if ( nrow(ovw$p.value) < 2 ) {
        nh <- nh*2; nw <- nw*2
    }
    
    ## plot figure
      png(file.path(fig.path,paste0("pfromAA_",ptype, "_AA_wtests.png")),
          height=nh, width=nw, res=300, units="in")
    par(mai=c(1,1,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:2, ylab=NA, xlab="",
                 col=ttcols, show.total=TRUE)
        figlabel(LAB, pos="bottomleft", cex=1)
    figlabel(ptype, pos="bottomright", cex=1.5)
    dev.off()
}

## by AA->AA
for ( ds in uds ) {

    dtmt <- tmtf[tmtf$Dataset==ds,]

    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="to", cols="from",
                       bg=FALSE, value="RAAS", 
                       use.test=t.test, do.plots=FALSE, 
                       verb=0)

    png(file.path(fig.path,paste0(ds,"_AA_wtests.png")),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
                 show.total=TRUE)
    figlabel(LAB, pos="bottomleft", cex=1)
    figlabel(pos="bottomright", text=ds, cex=2, font=2)
    dev.off()
    
    png(file.path(fig.path,paste0(ds,"_AA_volcano.png")),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.75,.1,.75), mgp=c(1.3,.3,0), tcl=-.25)
    volcano(ovw, cut=100, p.txt=3, v.txt=c(-Inf,-1), density=TRUE,
            xlab="mean TMT RAAS", value="mean")
    dev.off()

    png(file.path(fig.path,paste0(ds,"_raas_profile_colors.png")),
                  res=300, width=3, height=3, units="in")
    par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
    lraas.col <- selectColors(c(ovw$mean),
                              q=.05, colf=viridis::viridis,

                              n=100,
                              xlim=c(-4,1),
                              plot=TRUE,
                              mai=c(.5,.5,.1,.1), xlab="All TMT level RAAS")
    axis(1, at=seq(-4,4,.5), labels=FALSE)
    figlabel(LAB, pos="bottomleft", cex=1)
    dev.off()
    
    png(file.path(fig.path,paste0(ds,"_AA_raas_mean.png")),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotovl(ovw, value="mean", text="count", txt.cut=-2,
            col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
            axis=1:2, ylab=NA, xlab="", text.cex=.8)
    mtext("mean TMT-level RAAS", 3, .5, cex=1.2)
    figlabel(LAB, pos="bottomleft", cex=1)
    figlabel(ds, pos="bottomright", cex=1.5)
    dev.off()
    
    png(file.path(fig.path,paste0(ds,"_AA_count_unique.png")),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    cnt <- ovw$unique
    cnt[cnt==0] <- NA
    txt <- ovw$unique
    txt[txt==0] <- ""
    txt.col <- ifelse(ovw$unique>quantile(c(ovw$unique),.95),"white","black")
    image_matrix(cnt, col=gcols, axis=1:2,
                 text=txt, text.col=txt.col, ylab="",
                 xlab="", text.cex=.8)
    ##axis(4, nrow(ovw$p.value):1, labels=ACODONS[rownames(ovw$p.value)], las=2)
    mtext("mistranslated", 2, 3.5, adj=-1.7, cex=1.2)
    mtext("# unique SAAP", 3, .5, cex=1.2)
    figlabel(LAB, pos="bottomleft", cex=1)
    dev.off()
}


## by codon 
ctmt <- tmtf[tmtf$codon!="",]
for ( ds in uds ) {

    dtmt <- tmtf[ctmt$Dataset==ds,]

    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="to", cols="aacodon",
                       bg=FALSE, use.test=t.test, do.plots=FALSE, 
                       verb=0)

    png(file.path(fig.path,paste0(ds,"_codon_wtests.png")),
        res=300, width=14, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.7,.5,.7,2.7))
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:3, ylab=NA, xlab="", col=ttcols,
                 show.total=FALSE)
    axis(4, nrow(ovw$p.value):1,
         labels=ACODONS[rownames(ovw$p.value)], las=2)
    figlabel(LAB, pos="topright", cex=1.5)
    figlabel(pos="bottomright", text=ds, cex=2, font=2)
    dev.off()
    
    png(file.path(fig.path,paste0(ds,"_codon_volcano.png")),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.75,.1,.75), mgp=c(1.3,.3,0), tcl=-.25)
    volcano(ovw, cut=100, p.txt=3, v.txt=c(-Inf,-1), density=TRUE,
            xlab="mean TMT RAAS", value="mean")
    figlabel(LAB, pos="topright", cex=1)
    dev.off()

    png(file.path(fig.path,paste0(ds,"_codon_count_unique.png")),
        res=300, width=14, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.7,.5,.7,2.7))
    cnt <- ovw$unique
    cnt[cnt==0] <- NA
    txt <- ovw$unique
    txt[txt==0] <- ""
    txt.col <- ifelse(ovw$unique>quantile(c(ovw$unique),.99),"white","black")
    image_matrix(cnt, col=gcols, axis=1:2,
                 text=txt, text.col=txt.col, ylab="",
                 xlab="", text.cex=.8)
    ##axis(4, nrow(ovw$p.value):1, labels=ACODONS[rownames(ovw$p.value)], las=2)
    mtext("mistranslated", 2, 3.5, adj=-1.7, cex=1.2)
    mtext("# unique SAAP", 3, .5, cex=1.2)
    figlabel(LAB, pos="topright", cex=1.5)
    figlabel(pos="bottomright", text=ds, cex=2, font=2)
    dev.off()
}


### FIND MOST FREQUENT
## TODO: instead find strongest
## TODO: remove most frequent per Dataset, eg. N:M vs. T:V in PDAC

ptype <- unique(tmtf$pfromto)
pmost <- sapply(ptype, function(x) {
    names(which.max(table(tmtf$fromto[tmtf$pfromto==x])))
})[srt]

tdat <- tmtf[tmtf$fromto%in%pmost,]

ovw <- raasProfile(x=tdat, id="SAAP", 
                   rows="fromto", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=pmost, col.srt=uds,
                   use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)

png(file.path(fig.path,paste0("AAprop_wtests_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
mtext("most frequent", 2, 3, cex=1.5)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

png(file.path(fig.path,paste0("AAprop_raas_mean_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="mean", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1:2, ylab=NA, xlab="", text.cex=.8)
mtext("most frequent", 2, 3, cex=1.5)
mtext("mean TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

png(file.path(fig.path,paste0("AAprop_raas_median_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="median", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1:2, ylab=NA, xlab="", text.cex=.8)
mtext("most frequent", 2, 3, cex=1.5)
mtext("median TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

## remove most frequent
tdat <- tmtf[!tmtf$fromto%in%pmost,]

ovw <- raasProfile(x=tdat, id="SAAP", 
                   rows="pfromto", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=srt, col.srt=uds,
                   use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)

png(file.path(fig.path,paste0("AAprop_wtests_wo_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
mtext("w/o freq.", 1, 2, cex=1.5, adj=-1)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

png(file.path(fig.path,paste0("AAprop_raas_mean_wo_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="mean", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1:2, ylab=NA, xlab="", text.cex=.8)
mtext("w/o freq.", 1, 2, cex=1.5, adj=-1)
mtext("mean TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

png(file.path(fig.path,paste0("AAprop_raas_median_wo_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="median", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1:2, ylab=NA, xlab="", text.cex=.8)
mtext("w/o freq.", 1, 2, cex=1.5, adj=-1)
mtext("median TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()


### BY STRUCTURAL FEATURES

## IUPRED3
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="iupred3.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0, fname=file.path(fig.path,"iupred3_"))

png(file.path(fig.path,paste0("structure_iupred3.png")),
    res=300, width=4, height=2, units="in")
par(mai=c(.7,1,.01,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
mtext("IUPRED3", 2, 3.5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=.7)
dev.off()

## ANCHOR2
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="anchor2.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0, fname=file.path(fig.path,"anchor2_"))

png(file.path(fig.path,paste0("structure_anchor2.png")),
    res=300, width=4, height=2, units="in")
par(mai=c(.7,1,.01,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
mtext("ANCHOR2", 2, 3.5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=.7)
dev.off()

## S4 PRED SECONDARY STRUCTURE
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="sstruc", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(ssrt), col.srt=uds,
                   use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0, fname=file.path(fig.path,"s4pred_"))

png(file.path(fig.path,paste0("structure_s4pred.png")),
    res=300, width=4, height=1.5, units="in")
par(mai=c(.7,1,.01,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
mtext("S4PRED", 2, 3.5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=.7)
dev.off()


