
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
dir.create(fig.path, showWarnings=FALSE)

### START

## PARSE & FILTER DATA
dat <- read.delim(in.file)
tmtf <- read.delim(tmt.file)


dat$Hemoglobin.Albumin <- as.logical(dat$Hemoglobin.Albumin)

#### TODO: missing values

## NOTE: inconsistent remove tags
if ( FALSE ) {
    saap <- "LGMFNIQHGK" # flagged Keep=FALSE in tmt but not in hdat?
    saap <- "WYNLAIGSTGPWLK"
    saap <- "HIADLAGNSEVILVVPAFNVINGGSHAGNK"
    tmtf[which(tmtf$SAAP==saap),c("Keep.SAAP","Dataset")]
    dat[which(dat$SAAP==saap),c("remove","Keep.SAAP","Dataset")]
}

## remove excluded
hdat <- dat[which(dat$Keep.SAAP),]


## get raw RAAS data TMT level
tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)
## remove excluded
cat(paste("removing", sum(!tmtf$Keep.SAAP),
          "tagged as false positive on TMT level\n"))
tmtf <- tmtf[tmtf$Keep.SAAP,]
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

## fromto column
hdat$fromto <- apply(fromto, 1, paste0, collapse=":")

## generate columns by AA property
hdat$pfrom <- aaprop[fromto[,1]]
hdat$pto <- aaprop[fromto[,2]]
hdat$pfromto <- paste0(hdat$pfrom,":",hdat$pto)


#### FILTER

if (  exclude.albumin ) {
    ## TODO: also exclude from TMF level file
    hdat <- hdat[!hdat$Hemoglobin.Albumin,]
}

if ( exclude.frequent ) {

    ## TODO: MAKE SURE TO REMOVE THE SAME TYPES!!
    ## this is probably not required since we
    ## only look up SAAP that appear in hdat!

    rmsaap <- hdat$SAAP[hdat$from%in%frequent]
    tmtf <- tmtf[!tmtf$SAAP%in%rmsaap,]

    hdat <- hdat[!hdat$from%in%frequent,]
}

## find which are missing
##length(grep("NA",hdat$pfromto))
##adat <- hdat[grep("NA",hdat$pfromto, invert=TRUE),]

## sorting and labels
uds <- sort(unique(hdat$Dataset))
uds <- c("Healthy",uds[uds!="Healthy"])

srt <- c("charged","polar","hphobic","special")
srt <- paste(rep(srt,each=4), rep(srt,4), sep=":")

axex <- ftlabels(srt) # axis labels with arrows

dpath <- file.path(fig.path,"dists")
dir.create(dpath, showWarnings=FALSE)
fname <- file.path(dpath,paste0("AAprop_wtests_"))
source("~/work/mistrans/scripts/saap_utils.R")
ovw <- raasProfile(x=hdat, id="SAAP", values=tmtf, 
                   rows="pfromto", cols="Dataset",
                   bg=TRUE, vid="RAAS", row.srt=srt, col.srt=uds,
                   use.test=t.test, do.plots=TRUE, xlab="TMT level RAAS",
                   fname=fname, verb=0)
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
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

## local raas colors
png(file.path(fig.path,"raas_profile_colors.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(c(ovw$median),
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
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()


png(file.path(fig.path,paste0("AAprop_raas_median.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="median", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1, ylab=NA, xlab="", text.cex=.8)
axis(2, length(axex):1, labels=axex, las=2)
mtext("median TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

### median vs. p-value as a measure of effect size
plot(ovw$count, ovw$median, xlab="count", ylab="RAAS median")
plot(ovw$count, -log(ovw$p.value), xlab="count", ylab=expression(-log(p)))
plot(ovw$median, -log(ovw$p.value))

png(file.path(fig.path,paste0("AAprop_volcano.png")),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.75,.1,.75), mgp=c(1.3,.3,0), tcl=-.25)
volcano(ovw, cut=100, p.txt=5, v.txt=c(-2.5,-1), density=TRUE,
        xlab="median TMT RAAS", value="median")
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
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()


ovw <- raasProfile(x=hdat, id="SAAP", values=tmtf, 
                   rows="from", cols="Dataset",
                   bg=TRUE, vid="RAAS", col.srt=uds,
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
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

ovw <- raasProfile(x=hdat, id="SAAP", values=tmtf, 
                   rows="to", cols="Dataset",
                   bg=TRUE, vid="RAAS", col.srt=uds,
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
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()


## plot 16 plots for pfrom-to combos, for each to all AA
for ( ptype in unique(hdat$pfromto) ) {

    ddat <- hdat[hdat$pfromto==ptype,]
    ovw <- raasProfile(x=ddat, id="SAAP", values=tmtf, 
                       rows="fromto", cols="Dataset",
                       bg=TRUE, vid="RAAS", col.srt=uds,
                       use.test=t.test, do.plots=FALSE, xlab="TMT level RAAS",
                       verb=0)
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
    figlabel(LAB, pos="bottomleft", cex=1.5)
    figlabel(ptype, pos="bottomright", cex=1.5)
    dev.off()
}

## plot 20 from plots for pfrom-to combos, for each to all AA
for ( ptype in unique(hdat$to) ) {

    ddat <- hdat[hdat$to==ptype,]
    ovw <- raasProfile(x=ddat, id="SAAP", values=tmtf, 
                       rows="fromto", cols="Dataset",
                       bg=TRUE, vid="RAAS", col.srt=uds,
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
    figlabel(LAB, pos="bottomleft", cex=1.5)
    figlabel(ptype, pos="bottomright", cex=1.5)
    dev.off()
}

## plot 20 from plots for pfrom-to combos, for each to all AA
for ( ptype in unique(hdat$from) ) {

    ddat <- hdat[hdat$from==ptype,]
    ovw <- raasProfile(x=ddat, id="SAAP", values=tmtf, 
                       rows="fromto", cols="Dataset",
                       bg=TRUE, vid="RAAS", col.srt=uds,
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
    figlabel(LAB, pos="bottomleft", cex=1.5)
    figlabel(ptype, pos="bottomright", cex=1.5)
    dev.off()
}

## plot 20 from plots for AA->AAprop
hdat$frompto <- paste0(hdat$from, ":", hdat$pto)
for ( ptype in unique(hdat$from) ) {

    ddat <- hdat[hdat$from==ptype,]
    ovw <- raasProfile(x=ddat, id="SAAP", values=tmtf, 
                       rows="frompto", cols="Dataset",
                       bg=TRUE, vid="RAAS", col.srt=uds,
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
        figlabel(LAB, pos="bottomleft", cex=1.5)
    figlabel(ptype, pos="bottomright", cex=1.5)
    dev.off()
}

## by AA->AA
for ( ds in uds ) {
    ddat <- hdat[hdat$Dataset==ds,]
    dtmt <- tmtf[tmtf$Dataset==ds,]
    dtmt <- split(dtmt$RAAS,  dtmt$SAAP)

    ovw <- raasProfile(x=ddat, id="SAAP", values=dtmt, 
                       rows="to", cols="from",
                       bg=FALSE, vid="RAAS", 
                       use.test=w.test, do.plots=FALSE, 
                       verb=0)

    png(file.path(fig.path,paste0(ds,"_AA_wtests.png")),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
                 show.total=TRUE)
    figlabel(LAB, pos="bottomleft", cex=1.5)
    figlabel(pos="bottomright", text=ds, cex=2, font=2)
    dev.off()
    
    png(file.path(fig.path,paste0(ds,"_raas_profile_colors.png")),
                  res=300, width=3, height=3, units="in")
    par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
    lraas.col <- selectColors(c(ovw$median),
                              q=.05, colf=viridis::viridis,

                              n=100,
                              xlim=c(-4,1),
                              plot=TRUE,
                              mai=c(.5,.5,.1,.1), xlab="All TMT level RAAS")
    axis(1, at=seq(-4,4,.5), labels=FALSE)
    figlabel(LAB, pos="bottomleft", cex=1)
    dev.off()
    
    png(file.path(fig.path,paste0(ds,"_AA_raas_median.png")),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotovl(ovw, value="median", text="count", txt.cut=-2,
            col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
            axis=1:2, ylab=NA, xlab="", text.cex=.8)
    mtext("median TMT-level RAAS", 3, .5, cex=1.2)
    figlabel(LAB, pos="bottomleft", cex=1.5)
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
    figlabel(LAB, pos="bottomleft", cex=1.5)
    dev.off()
}


## by codon 
cdat <- hdat[hdat$codon!="",]
## add useful aa/codon tag
cdat$aacodon <- paste(cdat$from,cdat$codon, sep="-")


for ( ds in uds ) {
    ddat <- cdat[cdat$Dataset==ds,]
    dtmt <- tmtf[tmtf$Dataset==ds,]
    dtmt <- split(dtmt$RAAS,  dtmt$SAAP)

    ovw <- raasProfile(x=ddat, id="SAAP", values=dtmt, 
                       rows="to", cols="aacodon",
                       bg=FALSE, use.test=w.test, do.plots=FALSE, 
                       verb=0)

    png(file.path(fig.path,paste0(ds,"_codon_wtests.png")),
        res=300, width=14, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.7,.5,.7,2.7))
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:3, ylab=NA, xlab="", col=ttcols,
                 show.total=FALSE)
    axis(4, nrow(ovw$p.value):1, labels=ACODONS[rownames(ovw$p.value)],
         las=2)
    figlabel(LAB, pos="topright", cex=1.5)
    figlabel(pos="bottomright", text=ds, cex=2, font=2)
    dev.off()
 
}

### FIND MOST FREQUENT
## TODO: instead find strongest
## TODO: remove most frequent per Dataset, eg. N:M vs. T:V in PDAC

ptype <- unique(hdat$pfromto)
pmost <- sapply(ptype, function(x) {
    names(which.max(table(hdat$fromto[hdat$pfromto==x])))
})[srt]

tdat <- hdat[hdat$fromto%in%pmost,]

ovw <- raasProfile(x=tdat, id="SAAP", values=tmtf, 
                   rows="fromto", cols="Dataset",
                   bg=TRUE, vid="RAAS", row.srt=pmost, col.srt=uds,
                   use.test=w.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)

png(file.path(fig.path,paste0("AAprop_wtests_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
mtext("most frequent", 2, 3, cex=1.5)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

png(file.path(fig.path,paste0("AAprop_raas_mean_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="mean", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1:2, ylab=NA, xlab="", text.cex=.8)
mtext("most frequent", 2, 3, cex=1.5)
mtext("mean TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

png(file.path(fig.path,paste0("AAprop_raas_median_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="median", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1:2, ylab=NA, xlab="", text.cex=.8)
mtext("most frequent", 2, 3, cex=1.5)
mtext("median TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

## remove most frequent
tdat <- hdat[!hdat$fromto%in%pmost,]

ovw <- raasProfile(x=tdat, id="SAAP", values=tmtf, 
                   rows="pfromto", cols="Dataset",
                   bg=TRUE, vid="RAAS", row.srt=srt, col.srt=uds,
                   use.test=w.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)

png(file.path(fig.path,paste0("AAprop_wtests_wo_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
mtext("w/o freq.", 1, 2, cex=1.5, adj=-1)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

png(file.path(fig.path,paste0("AAprop_raas_mean_wo_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="mean", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1:2, ylab=NA, xlab="", text.cex=.8)
mtext("w/o freq.", 1, 2, cex=1.5, adj=-1)
mtext("mean TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

png(file.path(fig.path,paste0("AAprop_raas_median_wo_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotovl(ovw, value="median", text="count", txt.cut=-2,
        col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
        axis=1:2, ylab=NA, xlab="", text.cex=.8)
mtext("w/o freq.", 1, 2, cex=1.5, adj=-1)
mtext("median TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()


### TEST IUPRED3 by AA CLASS

filt <- hdat$Dataset=="PDAC" #hdat$Dataset!="Healthy"
boxplot(hdat$iupred3[filt] ~ hdat$fromto[filt], las=2)
boxplot(hdat$anchor2[filt] ~ hdat$fromto[filt], las=2)

boxplot(hdat$iupred3[filt] ~ hdat$pfromto[filt], las=2)
boxplot(hdat$anchor2[filt] ~ hdat$pfromto[filt], las=2)


### TEST HPHOBIC->SPECIAL HIGH RAAS IN CCRCC

filt <- hdat$Dataset=="CCRCC" & hdat$pfromto=="hphobic:special"
cat(paste(nrow(hdat[filt,]), "unique hits\n"))
table(hdat$fromto[filt])
##A:G A:P L:P Y:G 
##  8   2   1   2 
ridx <- (which(tmtf$SAAP%in%hdat$SAAP[filt] & tmtf$Dataset=="CCRCC"))
tmtf$RAAS[ridx]

## most A:G at mean RAAS, extreme RAAS
## comes from SAAP with three RAAS: L:P
ridx <- (which(tmtf$SAAP%in%hdat$SAAP[filt] & tmtf$Dataset=="CCRCC"))
tmtf$RAAS[ridx]
filt <- hdat$Dataset=="CCRCC" & hdat$fromto=="L:P"
ridx <- (which(tmtf$SAAP%in%hdat$SAAP[filt] & tmtf$Dataset=="CCRCC"))
tmtf$RAAS[ridx]

hdat[filt,] ## PROTEIN NOT MATCHED, missing from codons etc.
## YTLPPGVDPTKVSSSLSPEGTLTVEAPMPK
## BLAST: 
##  	estrogen receptor-related protein [Homo sapiens] 	Homo sapiens 	91.8 	91.8 	100% 	2e-19 	96.67% 	78 	AAB20722.1
## Select seq pdb|3Q9P|A 	Chain A, Heat shock protein beta-1 [Homo sapiens] 	Homo sapiens 	91.8 	91.8 	100% 	2e-19 	96.67% 	85 	3Q9P_A
## Select seq emb|CAA34498.1| 	p24k-1 (AA 1-91) [Homo sapiens] 	Homo sapiens 	91.8 	91.8 	100% 	3e-19 	96.67% 	91 	CAA34498.1

ridx <- (which(tmtf$SAAP%in%hdat$SAAP[filt] & tmtf$Dataset=="CCRCC"))
tmtf$RAAS[ridx]

## get max
dtmt <- tmtf[ tmtf$Dataset=="CCRCC",]
tmt <- split(dtmt$RAAS, dtmt$SAAP) 

filt <- hdat$Dataset=="CCRCC" & hdat$pfromto=="hphobic:special"

sid <- names(which.max(unlist(lapply(tmt[hdat$SAAP[filt]], function (x) max(x)))))
##sid <- "IIHLGGK"
hdat[hdat$SAAP==sid,]
## Q96IR2 · ZN845_HUMAN <- max RAAS in set of 17 hphobic->special in CCRCC
## Protein Zinc finger protein 845

## NOTE: maximal value not in L:P

## NOTE: only 7 RAAS at 8 hits for A:G
filt <- hdat$Dataset=="CCRCC" & hdat$fromto=="A:G"
table(hdat$SAAP[filt])
which(tmtf$SAAP%in%hdat$SAAP[filt] & tmtf$Dataset=="CCRCC")
hdat$SAAP[which(!hdat$SAAP[filt]%in%tmtf$SAAP[tmtf$Dataset=="CCRCC"])]


## SIMPLE BOXPLOTS

boxplot(df$RAAS ~ df$pfrom)

#### MIXED LINEAR MODEL:
## TODO: does this help to answer the question,
## whether amino acid properties or individual AA best explain
## Dataset-specific high RAAS
library(lme4)
library(lmerTest)

## TODO: full TMT matric
idx <- match(tmtf$SAAP, hdat$SAAP)
tmtf$from <- hdat$from[idx]
tmtf$to <- hdat$to[idx]
tmtf$pfrom <- hdat$pfrom[idx]
tmtf$pto <- hdat$pto[idx]

## TODO: can we use this?

df <- data.frame(RAAS=tmtf$RAAS,
                 Dataset=tmtf$Dataset,
                 from=hdat$from[idx],
                 to=hdat$to[idx],
                 pfrom=hdat$pfrom[idx],
                 pto=hdat$pto[idx],
                 stringsAsFactors = TRUE)

fprp <- RAAS ~ pfrom*pto + (1|Dataset)
mprp <- lmer(fprp, data = df, REML=TRUE,
             control = lmerControl(optimizer="bobyqa",
                                   optCtr = list(maxfun = 1e9)))
summary(mprp)
image_matrix(vcov(mprp), col=viridis::viridis(100), axis=1:2)

faa <- RAAS ~ from*to + (1|Dataset)
mfaa <- lmer(faa, data = df, REML=TRUE,
             control = lmerControl(optimizer="bobyqa",
                                   optCtr = list(maxfun = 1e9)))
summary(mfaa)
X <- model.matrix(fprp, df)

##summary(frn_diff_parsimoniuous)
      
## raas_value ~ globaleeigenschaft1*globaleeigenschaft2*globaleeigenschaft3 + (1+globaleeigenschaft1*globaleeigenschaft2*globaleeigenschaft3|Aminosäure)

##ich benutz lme4 für lmer() und lmerTest damit da p-werte dabei stehen
##alles ohne klammern und ohne | sind fixed effects
##in der klammer dann die random effects

##parsimonious_m_formula <- 
##value ~ 1 + delay + valence + pe_unsigned_zero + delay:valence + 
##    delay:pe_unsigned_zero + valence:pe_unsigned_zero + delay:valence:pe_unsigned_zero + 
##    (1 + delay + valence + delay:valence + pe_unsigned_zero + 
##        delay:pe_unsigned_zero | id) + (1 | electrode)

##frn_diff_parsimoniuous <- lmer(parsimonious_m_formula,
##                               data = data, REML=T,
##                               control = lmerControl(optimizer="bobyqa", optCtr = list(maxfun = 1e9)))
##summary(frn_diff_parsimoniuous)
