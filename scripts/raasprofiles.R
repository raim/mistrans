
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

## TODO 20240222
## * plotovl: plot all function,
## * fall back on median,
## * re-blast saap, re-annotate,
## * Healthy: untangle tissues and collapse cancer.

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

use.test <- t.test

healthy <- FALSE

exclude.albumin <- TRUE # FALSE # 
only.unique <- FALSE # TRUE # 

exclude.frequent <- FALSE # TRUE # 
frequent <- c("Q","W","T","S")

LAB <- "all"
fig.path <- file.path(proj.path,"figures","raasprofiles")
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
    fig.path <- paste0(fig.path,"_unique")
    LAB <- paste0(LAB, ", unique SAAP")
}

SETID <- ifelse(healthy,"tissues","cancer")


## figure output paths
dir.create(fig.path, showWarnings=FALSE)
## folder for detailed distributions
dpath <- file.path(fig.path,"dists")
dir.create(dpath, showWarnings=FALSE)

### START

## PARSE & FILTER DATA
dat <- read.delim(in.file)
tmtf <- read.delim(tmt.file)

## convert to logical
tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)
dat$Hemoglobin.Albumin <- as.logical(dat$Hemoglobin.Albumin)

#### TODO: find out which/how many TMT level SAAP are missing
## from the protein level file, and why.


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
cat(paste("removing", sum(!dat$keep),
          "tagged as false positive on protein level\n"))
hdat <- dat[which(dat$keep),]

## get raw RAAS data TMT level
## remove excluded
cat(paste("removing", sum(!tmtf$keep),
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
tmtf$iupred3 <- hdat$iupred3[idx]
tmtf$anchor2 <- hdat$anchor2[idx]
tmtf$iupred3.bins <- hdat$iupred3.bins[idx]
tmtf$anchor2.bins <- hdat$anchor2.bins[idx]
tmtf$sstruc <- hdat$sstruc[idx]

tmtf$albumin <- hdat$Hemoglobin.Albumin[idx]

### CHOOSE DATASETS
if ( healthy ) {
    ds <- tmtf$Dataset
    tmtf$Dataset[ds=="Healthy"] <- tmtf$TMT.Tissue[ds=="Healthy"]
    tmtf$Dataset[ds!="Healthy"] <- "Cancer"
}

### MEAN AND MEDIAN RAAS
## per data set and unique SAAP
usaap <- paste0(tmtf$SAAP,"/",tmtf$BP,"/",tmtf$Dataset)
araas <- split(tmtf$RAAS, usaap)

raas.median <- unlist(lapply(araas, median))
raas.emean <- unlist(lapply(araas, mean))
raas.mean <- unlist(lapply(araas, function(x) {log10(mean(10^x))}))

tmtf$unique <- usaap
tmtf$RAAS.mean <- raas.mean[usaap]
tmtf$RAAS.median <- raas.median[usaap]

raas.bins <- cut(tmtf$RAAS.median, breaks=seq(-6,4,.5))
raas.srt <- levels(raas.bins)
tmtf$raas.bins <- as.character(raas.bins)
tmtf$raas.bins[is.na(tmtf$raas.bins)] <- "na"

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
    tmtf$RAAS.orig <- tmtf$RAAS
    tmtf$RAAS <- tmtf$RAAS.mean
}

### TODO: merge unique SAAP here to test effects of multiple
### measurements.

## from this point on we only work with tmtf, which
## now should have all information

## sorting and labels
uds <- sort(unique(tmtf$Dataset))
uds <- c("Healthy",uds[uds!="Healthy"])
uds <- c("Cancer",uds[uds!="Cancer"])
uds <- uds[uds%in%tmtf$Dataset]

srt <- c("charged","polar","hphobic","special")
srt <- paste(rep(srt,each=4), rep(srt,4), sep=":")

axex <- ftlabels(srt) # axis labels with arrows


## by AA properties
fname <- file.path(dpath,paste0(SETID,"_AAprop_wtests_"))
ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
                    delog=TRUE, rows="pfromto", cols="Dataset",
                    bg=TRUE, row.srt=srt, col.srt=uds,
                    use.test=use.test, do.plots=TRUE, xlab="TMT level RAAS",
                    fname=fname, verb=1)

## local raas colors
png(file.path(fig.path,paste0(SETID,"raas_profile_colors.png")),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(c(ovw$median),
                          q=0.01, colf=viridis::viridis,
                          n=10, xlim=c(-4,1), plot=TRUE,
                          mai=c(.5,.5,.1,.1), xlab="All TMT level RAAS")
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

## legend for all wtests
png(file.path(fig.path,paste0("wtests_legend.png")),
    res=300, width=2, height=2, units="in")
par(mai=c(.6,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlapsLegend(p.min=p.min, p.txt=p.txt, type=2, col=ttcols)
dev.off()

## plot all
plotProfiles(ovw, fname=file.path(fig.path,paste0("AAprop_",SETID)),
             mai=c(.8,1.5,.5,.5),
             ttcols=ttcols, value="median",
             rlab=LAB, llab="",
             vcols=lraas.col$col, vbrks=lraas.col$breaks,
             gcols=gcols)

## TODO: median vs. p-value as a measure of effect size
if ( interactive() ) {
    plot(ovw$count, ovw$median, xlab="count", ylab="RAAS median")
    plot(ovw$count, -log(ovw$p.value), xlab="count", ylab=expression(-log(p)))
    plot(ovw$median, -log(ovw$p.value))
}

## by "from" amino acid
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="from", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)
plotProfiles(ovw, fname=file.path(fig.path,paste0("fromAA_",SETID)),
             mtxt="encoded AA",
             mai=c(.8,.5,.5,.5),
             ttcols=ttcols, value="median",
             rlab=LAB, llab="",
             vcols=lraas.col$col, vbrks=lraas.col$breaks,
             gcols=gcols)

## by "to" amino acid
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="to", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)
plotProfiles(ovw, fname=file.path(fig.path,paste0("toAA_",SETID)),
             mai=c(.8,.5,.5,.5),
             mtxt="substituted AA", ttcols=ttcols, value="median",
             rlab=LAB, llab="", vcols=lraas.col$col, vbrks=lraas.col$breaks,
             gcols=gcols)


## plot 16 plots for pfrom-to combos, for each to all AA
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="fromto", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=1)

## plot only signficant
ovwp <- sortOverlaps(ovw, axis=2, p.min=p.txt, cut=TRUE)
if ( nrow(ovwp$p.value)>0 )
    plotProfiles(ovwp,
                 fname=file.path(fig.path,paste0("AA_",SETID,"_cut")),
                 mai=c(.8,.6,.5,.5), ttcols=ttcols, value="median",
                 rlab=LAB, llab=ptype,
                 vcols=lraas.col$col, vbrks=lraas.col$breaks,
                 gcols=gcols)
for ( ptype in unique(tmtf$pfromto) ) {

    ## sort by to
    qsrt <- sort(unique(tmtf$fromto[tmtf$pfromto==ptype]))
    ft <- unlist(lapply(strsplit(qsrt,":"), function(x) x[2]))
    qsrt <- qsrt[order(ft)]
                 
    ovwp <- sortOverlaps(ovw, srt=qsrt)

    plotProfiles(ovwp, fname=file.path(fig.path,paste0(ptype,"_",SETID)),
                 mai=c(.8,.6,.5,.5), ttcols=ttcols, value="median",
                 rlab=LAB, llab=ptype,
                 vcols=lraas.col$col, vbrks=lraas.col$breaks,
                 gcols=gcols)
    ## cut
    ovwp <- sortOverlaps(ovwp, axis=2, p.min=1e-5, cut=TRUE)
    if ( nrow(ovwp$p.value)>0 )
        plotProfiles(ovwp,
                     fname=file.path(fig.path,paste0(ptype,"_",SETID,"_cut")),
                     mai=c(.8,.6,.5,.5), ttcols=ttcols, value="median",
                     rlab=LAB, llab=ptype,
                     vcols=lraas.col$col, vbrks=lraas.col$breaks,
                     gcols=gcols)
}



## by codon->AA
ctmt <- tmtf[tmtf$codon!="",]
for ( ds in uds ) {


    ## TODO: understand background

    dtmt <- ctmt[ctmt$Dataset==ds,]

    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="to", cols="aacodon",
                       bg=FALSE, use.test=use.test, do.plots=FALSE, 
                       verb=0)

    ## TODO: re-create previous codon plots
    plotProfiles(ovw, fname=file.path(fig.path,paste0("codon_",SETID,"_",ds)),
                 fw=.2, mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
                 rlab=LAB, llab=ds,
                 vcols=lraas.col$col, vbrks=lraas.col$breaks,
                 gcols=gcols)
}

## by AA->AA
for ( ds in uds ) {

    ## TODO: understand background

    dtmt <- tmtf[tmtf$Dataset==ds,]

    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="to", cols="from",
                       bg=FALSE, value="RAAS", 
                       use.test=use.test, do.plots=FALSE, 
                       verb=1)

    plotProfiles(ovw, fname=file.path(fig.path,paste0("AA_",SETID,"_",ds)),
                 mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
                 rlab=LAB, llab=ds,
                 vcols=lraas.col$col, vbrks=lraas.col$breaks,
                 gcols=gcols)

}



### BY STRUCTURAL FEATURES

for ( ds in uds ) {

    tmtd <- tmtf
    if ( ds!="all" )
        tmtd <- tmtf[tmtf$Dataset==ds,]
    ovl <- clusterCluster(paste0(tmtd$raas.bins), tmtd$iupred3.bins,
                          cl1.srt=rev(raas.srt))
    par(mai=c(1,1,.5,.5), mgp=c(3,.3,0), tcl=-.25)
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt,
             xlab="IUPRED3", ylab="median RAAS")
    
    plotCor(tmtd$RAAS.median, tmtd$iupred3, xlab="median RAAS", ylab="IUPRED3")

    ovl <- clusterCluster(paste0(tmtd$raas.bins), tmtd$anchor2.bins,
                          cl1.srt=rev(raas.srt))
    par(mai=c(1,1,.5,.5), mgp=c(3,.3,0), tcl=-.25)
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt,
             xlab="ANCHOR2", ylab="median RAAS")
    
    plotCor(tmtd$RAAS.median, tmtd$anchor2, xlab="median RAAS", ylab="ANCHOR2")

}

## IUPRED3 vs. RAAS bins
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="iupred3.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)


## IUPRED3
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="iupred3.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=TRUE, xlab="TMT level RAAS",
                   verb=0, 
                   fname=file.path(dpath,paste0("iupred3_",SETID,"_")))

plotProfiles(ovw, fname=file.path(fig.path,paste0("structure_iupred3_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             rlab=LAB, mtxt="IUPRED3", mtxt.line=3.3,
             vcols=lraas.col$col, vbrks=lraas.col$breaks,
             gcols=gcols)

## ANCHOR2
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="anchor2.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0, fname=file.path(fig.path,"anchor2_"))

plotProfiles(ovw, fname=file.path(fig.path,paste0("structure_anchor2_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             rlab=LAB, mtxt="ANCHOR2", mtxt.line=3.3,
             vcols=lraas.col$col, vbrks=lraas.col$breaks,
             gcols=gcols)

## S4 PRED SECONDARY STRUCTURE
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="sstruc", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(ssrt), col.srt=uds,
                   use.test=use.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0, fname=file.path(fig.path,"s4pred_"))

plotProfiles(ovw, fname=file.path(fig.path,paste0("structure_s4pred_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             rlab=LAB, mtxt="SPRED4", mtxt.line=3.3,
             vcols=lraas.col$col, vbrks=lraas.col$breaks,
             gcols=gcols)


