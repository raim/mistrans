
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

## SORT CODONS BY AA PROPERTY
CODL <- strsplit(CODONS, ";")[names(aaprop)]
CODL <- CODL[unlist(lapply(CODL, function(x) !is.null(x)))]
CODL <- lapply(CODL, sort)
COD.SRT <- paste(GENETIC_CODE[unlist(CODL)], unlist(CODL), sep="-")

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

gen.path <- "~/data/mammary/"
codon.file <- file.path(gen.path,"processedData","coding_codons.tsv")

in.file <- file.path(out.path,"saap_mapped3.tsv")
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")

### PARAMETERES

RAAS.MIN <- -4
RAAS.MAX <-  0

p.min <- 1e-10
p.txt <- 1e-5

p.dot <- p.txt
dot.sze <- c(.1,1.5)

use.test <- t.test

healthy <- FALSE # TRUE # 

## TODO: extracellular is mostly album/globin - analyze
exclude.extracellular <- FALSE # TRUE # 
exclude.albumin <- FALSE # TRUE # 
only.unique <- FALSE # TRUE # 
include.kr <- FALSE # TRUE # 

exclude.frequent <- FALSE # TRUE # 
frequent <- c("Q","W","T","S")

LAB <- "all"
fig.path <- file.path(proj.path,"figures","raasprofiles3")
if (  exclude.albumin ) {
    fig.path <- paste0(fig.path,"_noalb")
    LAB <- "-Alb./Hemog."
}
if ( exclude.frequent ) {
    tmp <- paste(frequent,collapse=",")
    fig.path <- paste0(fig.path,"_", gsub(",","",tmp))
    LAB <- paste0("-", tmp)
}
if ( only.unique ) {
    fig.path <- paste0(fig.path,"_unique")
    LAB <- paste0(LAB, ", unique SAAP")
}
if ( exclude.extracellular ) {
    fig.path <- paste0(fig.path,"_no_extracellular")
    LAB <- paste0(LAB, " -extracell.")
}
if ( include.kr ) {
    fig.path <- paste0(fig.path,"_wKR")
    LAB <- paste0(LAB, " +K/R")
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

if ( include.kr ) {
    tmtf$Keep.SAAP[grep("R|K", tmtf$AAS)] <- TRUE
}

## tag albumin/hemoglobin
## TODO: compare and fuse with protein level annotation
##       as.logical(dat$Hemoglobin.Albumin)
dat$Hemoglobin.Albumin <- dat$albumin|dat$globin 
dat$Keep.SAAP <- !dat$IG

#### TODO: find out which/how many TMT level SAAP are missing
## from the protein level file, and why.

### UNIFY FILTER COLUMNS

## fuse all excluded tags for each unique SAAP
alls <- rbind(dat[,c("Keep.SAAP","SAAP")],
              tmtf[,c("Keep.SAAP","SAAP")])
alls <- split(alls$Keep.SAAP, alls$SAAP)

## inconsistents
##alls[which(unlist(lapply(alls, function(x) length(unique(x))))>1)]

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

ina <- which(is.na(idx))
if ( length(ina)>0 ) {
    cat(paste("TODO:", length(ina), "missing from unique saap file.\n"))
    tmtf <- tmtf[-ina,]
    idx <- idx[-ina]
    
}


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

tmtf$gene <- hdat$gene[idx]
tmtf$transcript <- hdat$transcript[idx]
tmtf$protein <- hdat$protein[idx]

tmtf$albumin <- hdat$Hemoglobin.Albumin[idx]
tmtf$extracellular <- hdat$extracellular[idx]

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

## TODO: albumin vs. globin
if ( FALSE ) {
    ovl <- clusterCluster(cl1=tmtf$extracellular, cl2=tmtf$albumin)
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt)
    dev.off()
}

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

    ## just take the first of each SAAP/BP per Dataset
    tmtf <- tmtf[!duplicated(tmtf$unique),]
    tmtf$RAAS.orig <- tmtf$RAAS
    tmtf$RAAS <- tmtf$RAAS.median #TODO: use delogged mean?
}
if (  exclude.extracellular ) {
    rmsaap <- tmtf$SAAP[tmtf$extracellular]
    tmtf <- tmtf[!tmtf$SAAP%in%rmsaap,]
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

### START ANALYSIS

## color legends

## RAAS COLORS
png(file.path(fig.path,paste0("legend_raas.png")),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(tmtf$RAAS,
                          mn=RAAS.MIN, mx=RAAS.MAX,colf=viridis::viridis,
                          n=50, plot=TRUE,
                          mai=c(.5,.5,.1,.1),
                          xlab=expression(TMT~level~log[10]*RAAS))
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

## globally used RAAS colors!!
vcols <- lraas.col$col
vbrks <- lraas.col$breaks

## legend for all two-sided statistics
png(file.path(fig.path,paste0("legend_wtests.png")),
    res=300, width=2, height=2, units="in")
par(mai=c(.6,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlapsLegend(p.min=p.min, p.txt=p.txt, type=2, col=ttcols)
dev.off()

## legend for dot plot
pp <- seq(0, -log10(p.dot), length.out=6)
rs <- seq(RAAS.MIN,RAAS.MAX+1, length.out=6)

pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))

colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- rs
    
ovw <- list()
ovw$p.value <- t(10^-pm)
ovw$median <- t(rm)
source("~/work/mistrans/scripts/saap_utils.R")

mai <- c(.5,.5,.1,.1)
fh <- fw <- .2
nh <- nrow(ovw$p.value) *fh + mai[1] + mai[3]
nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]
segmenTools::plotdev(file.path(fig.path,"legend_dotplot"),
                     height=nh, width=nw, res=300)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovw, value="median",
           vbrks=vbrks,
           vcols=vcols, 
           dot.sze=dot.sze, p.dot=p.dot, axis=1:2,
           ylab=expression(log[10]~p),
           xlab=expression(log[10]~RAAS))
dev.off()

           
## global distribution by cancer type
           
ylm <- range(tmtf$RAAS)
par(mfcol=c(1,2))
boxplot(tmtf$RAAS ~ factor(tmtf$Dataset, levels=uds), ylim=ylm)
tmtu <- tmtf[!duplicated(tmtf$unique),]
boxplot(tmtu$RAAS.median ~ factor(tmtu$Dataset, levels=uds), ylim=ylm)
## TODO: pairs of same SAAP
tmtl <- split(tmtf$Dataset, tmtf$SAAP)

## RAAS profiles by AA properties
fname <- file.path(dpath,paste0("AAprop_",SETID,"_wtests_"))
ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
                    delog=TRUE, rows="pfromto", cols="Dataset",
                    bg=TRUE, row.srt=srt, col.srt=uds,
                    use.test=use.test, do.plots=TRUE,
                    xlab=expression(TMT~level~log[10]*RAAS),
                    fname=fname, verb=0)
ovwp <- sortOverlaps(ovw, axis=2, p.min=p.txt, cut=TRUE)



## plot all
source("~/work/mistrans/scripts/saap_utils.R")
plotProfiles(ovw, fname=file.path(fig.path,paste0("AAprop_",SETID)),
             mai=c(.8,1.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="",
             vcols=vcols, vbrks=vbrks,
             gcols=gcols, verb=1)
if ( nrow(ovwp$p.value)>0 )
    plotProfiles(ovwp,
                 fname=file.path(fig.path,paste0("AAprop_",SETID,"_cut")),
                 mai=c(.8,1.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="",
                 vcols=vcols, vbrks=vbrks,
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
                   use.test=use.test, do.plots=FALSE,
                   xlab=expression(TMT~level~log[10]*RAAS),
                   verb=0)
plotProfiles(ovw, fname=file.path(fig.path,paste0("fromAA_",SETID)),
             mtxt="encoded AA",
             mai=c(.8,.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="",
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## by "to" amino acid
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="to", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=expression(TMT~level~log[10]*RAAS),
                   verb=0)
plotProfiles(ovw, fname=file.path(fig.path,paste0("toAA_",SETID)),
             mai=c(.8,.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             mtxt="substituted AA", ttcols=ttcols, value="median",
             rlab=LAB, llab="", vcols=vcols, vbrks=vbrks,
             gcols=gcols)


## plot 16 plots for pfrom-to combos, for each to all AA
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="fromto", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=TRUE,
                   fname=file.path(dpath,paste0("AA_",SETID,"_")),
                   xlab=expression(TMT~level~log[10]*RAAS),
                   verb=0)

## plot only signficant
ovwp <- sortOverlaps(ovw, axis=2, p.min=p.txt, cut=TRUE)
if ( nrow(ovwp$p.value)>0 )
    plotProfiles(ovwp,
                 fname=file.path(fig.path,paste0("AA_",SETID,"_cut")),
                 mai=c(.8,.6,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="",
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols)
for ( ptype in unique(tmtf$pfromto) ) {

    ## sort by to
    qsrt <- sort(unique(tmtf$fromto[tmtf$pfromto==ptype]))
    ft <- unlist(lapply(strsplit(qsrt,":"), function(x) x[2]))
    qsrt <- qsrt[order(ft)]
                 
    ovwp <- sortOverlaps(ovw, srt=qsrt)

    plotProfiles(ovwp, fname=file.path(fig.path,paste0(ptype,"_",SETID)),
                 mai=c(.8,.6,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=ptype,
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols)
    ## cut
    ovwp <- sortOverlaps(ovwp, axis=2, p.min=p.txt, cut=TRUE)
    if ( nrow(ovwp$p.value)>1 ) #TODO: fix error with nrow=1 in 
        plotProfiles(ovwp,
                     fname=file.path(fig.path,paste0(ptype,"_",SETID,"_cut")),
                     mai=c(.8,.6,.5,.5), ttcols=ttcols, value="median",
                     p.min=p.min, p.txt=p.txt,
                     dot.sze=dot.sze, p.dot=p.dot,
                     rlab=LAB, llab=ptype,
                     vcols=vcols, vbrks=vbrks,
                     gcols=gcols)
    
}



## by codon->AA
ctmt <- tmtf[tmtf$codon!="",]
codons <- read.delim(codon.file, row.names=1)
for ( ds in c(uds,"all") ) {

    ## TODO: add codon frequency two-sided bar plot

    dtmt <- ctmt
    if ( ds!="all" )
        dtmt <- ctmt[ctmt$Dataset==ds,]


    ## NOTE: bg=FALSE, no column-wise background !
    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="to", cols="aacodon",
                       bg=FALSE, use.test=use.test, do.plots=FALSE, 
                       verb=0)

    ## TODO: expand ovw to ALL codons and AA to align plots!

    
    ## sort by amino acid property!
    aasrt <- names(aaprop)[names(aaprop)%in%rownames(ovw$p.value)]
    ovw <- sortOverlaps(ovw, axis=2, srt=aasrt, cut=TRUE)


    ## CUSTOM CODON FREQUENCY PLOTS

    ## CODON COUNTS AT THE AAS OF UNIQUE SAAP
    lcodt <- apply(ovw$unique, 2, sum)
    names(lcodt) <- sub(".*-","",names(lcodt))
    ## AAS codon frequencies per AA
    lcodl <- split(lcodt, GENETIC_CODE[names(lcodt)])
    
    lcodl <- lapply(lcodl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
    lcodl <- lcodl[names(aaprop)] ## ##SORT BY AA PROP

    lcods <- unlist(lapply(lcodl, function(x) x/sum(x)))

    ## CODON COUNTS IN ALL UNIQUE TRANSCRIPTS 
    cod  <- codons[unique(dtmt$transcript),]
    codt <- apply(cod,2,sum)
    ## transcript codon frequencies per AA
    codl <- split(codt, GENETIC_CODE[sub(".*\\.","",names(codt))])
    codl <- codl
    cods <- unlist(lapply(codl, function(x) x/sum(x)))

    ## sort (AA property ->codon frequency)
    cdsrt <- sub("\\.","-",names(unlist(lcodl)))
    ovw <- sortOverlaps(ovw, axis=1, srt=cdsrt, cut=FALSE)

    ## codon class lines
    cdn.cnt <- table(sub("-.*","",colnames(ovw$p.value)))[names(aaprop)]
    cdn.cnt <- cumsum(cdn.cnt[!is.na(cdn.cnt)])
    
    ## SORT CODON FREQUENCIES as in main plots
    names(cods) <- sub("\\.","-", names(cods))
    names(lcods) <- sub("\\.","-", names(lcods))
    lcods <- lcods[cdsrt]
    cods <- cods[cdsrt]

    ## TODO: re-create previous codon plots
    plotProfiles(ovw, fname=file.path(fig.path,paste0("codon_",SETID,"_",ds)),
                 fw=.2, mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=ds, mtxt="substituted AA",
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols, col.lines=cdn.cnt)

    ## TODO: better AA (y-axis) pre-sorting by property
    ovwp <- sortOverlaps(ovw, axis=1, p.min=p.txt, cut=TRUE)
    if ( nrow(ovwp$p.value)>1 & ncol(ovwp$p.value)>0 )
        plotProfiles(ovwp,
                     fname=file.path(fig.path,paste0("codon_",
                                                     SETID,"_",ds,"_cut")),
                     fw=.2, mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
                     p.min=p.min, p.txt=p.txt,
                     dot.sze=dot.sze, p.dot=p.dot,
                     rlab=LAB, llab=ds,
                     vcols=vcols, vbrks=vbrks,
                     gcols=gcols)


    mai <- c(.05,.5,.1,.5)
    fw <- .2
    nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]
    y <- log2(lcods/cods[names(lcods)]) # log2 ratio of codon frequencies
    ylm <- max(abs(y),na.rm=TRUE)

    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_ratio_old")),
            type="png", res=300, width=nw,height=.75)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plot(1:length(lcods), y, type="h", axes=FALSE, xlab=NA, lwd=3,
         ylab=expression(lg[2]~r),
         xlim=c(0.5,length(lcods)+.5), col=2, ylim=c(-ylm,ylm))
    abline(h=0)
    for ( ax in c(2,4) ) {
        axis(ax, at=seq(-5,5,.5), labels=FALSE)
        axis(ax, at=c(-.5,.5), cex.axis=1, las=2)
    }
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    dev.off()
    mai.bar <- mai
    mai.bar[c(2,4)] <- mai.bar[c(2,4)] +.05
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,"_codons_ratio")),
            type="png", res=300, width=nw,height=.75)
    par(mai=mai.bar, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    bp <- barplot(y, axes=FALSE, xlab=NA, beside=TRUE,
                  ylab=expression(lg[2]~r), las=2, xaxt='n', 
                  width=.5, space=.5, ylim=c(-ylm,ylm))
    ##abline(v=bp[2,cdn.cn]+  (bp[1,cdn.cn+1]-bp[2,cdn.cn])/2)
    cdn.cn <- head(cdn.cnt, length(cdn.cnt)-1)
    abline(v=bp[cdn.cn,]+unique(diff(bp[,1]))/2)
    ##abline(v=par("usr")[1]+unique(diff(bp)[,1])*cdn.cnt, lwd=1, xpd=TRUE)
    ##x2 <- par("usr")[2]-par("usr")[1]
    ##rel.cnt <- cdn.cnt/max(cdn.cnt)
    ##rel.cnt <- head(rel.cnt, length(rel.cnt)-1)
    ##abline(v=par("usr")[1]+x2*rel.cnt, lwd=1, xpd=TRUE)
    abline(h=0)
    for ( ax in c(2,4) ) {
        axis(ax, at=seq(-5,5,.1), labels=FALSE, tcl=par("tcl")/2, line=.2)
        axis(ax, at=seq(-5,5,.5), labels=FALSE, line=.2)
        axis(ax, at=c(-1,-.5,0,.5,1),
             label=c("-1","0.5","0","0.5","1"), cex.axis=1, las=2, line=.2)
    }
    dev.off()
    
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,"_codons_barplot")),
            type="png", res=300, width=nw, height=.75)
    par(mai=mai.bar, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    bp  <- barplot(rbind(lcods,cods[names(lcods)]),
                   axes=FALSE, xlab=NA, beside=TRUE,
                   ylab="", las=2, xaxt='n')
    cdn.cn <- head(cdn.cnt, length(cdn.cnt)-1)
    abline(v=bp[2,cdn.cn]+  (bp[1,cdn.cn+1]-bp[2,cdn.cn])/2)
    ##x2 <- par("usr")[2]
    ##rel.cnt <- cdn.cnt/max(cdn.cnt)
    ##rel.cnt <- head(rel.cnt, length(rel.cnt)-1)
    ##abline(v=.5+x2*rel.cnt, lwd=1, xpd=TRUE)
    axis(2, las=2, line=.2)
    axis(4, las=2, line=.2)
    dev.off()
    

    ## hypergeo tests
    aaa <- unique( GENETIC_CODE)
    aam <- matrix(1, ncol=length(lcodt), nrow=1)
    colnames(aam) <- names(lcodt)
    aac <- aap <- aam
    aac[] <- 0
    
    for ( aa in unique(aaa) ) {
        alc <- names(GENETIC_CODE)[GENETIC_CODE==aa]
        
        if ( length(alc)<1 ) next
        for ( j in  seq_along(alc) ) {
            
            wcd <- alc[ j] # white balls
            bcd <- alc[-j] # black balls
            
            if ( !wcd%in%names(lcodt) ) next
            
            m <- codt[wcd]; # number of white balls
            n <- sum(codt[bcd], na.rm=TRUE); # number of black balls
            q <- lcodt[wcd] # white balls drawn - CODON OF INTEREST
            k <- sum(lcodt[alc], na.rm=TRUE) # number of balls drawn
            
            aac[,wcd] <- q
            aam[,wcd] <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)
            aap[,wcd] <- phyper(q=q, m=m, n=n, k=k, lower.tail=TRUE)
            ##cat(paste("testing", aa, wcd, signif(aam[,wcd]), "\n"))
        }
    }
    
    colnames(aam) <- colnames(aap) <- colnames(aac) <-
        paste0(GENETIC_CODE[colnames(aam)],"-", colnames(aam))

    ## sort
    aam <- aam[,cdsrt,drop=FALSE]
    aap <- aap[,cdsrt,drop=FALSE]
    aac <- aac[,cdsrt,drop=FALSE]

    
    ovl <- list()
    ovl$p.value <- aam
    ovl$count <- aac
    
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_enriched")),
            type="png", res=300, width=nw,height=.25)
    par(mai=c(.05,.5,.05,.5), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotOverlaps(ovl, p.min=1e-5, p.txt=1e-2, axis=3, xlab=NA, ylab=NA,
                 col=upcols, text.cex=.7)
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    box()
    dev.off()
    
    ovl$p.value <- aap
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_deprived")),
            type="png", res=300, width=nw,height=.25)
    par(mai=c(.05,.5,.05,.5), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotOverlaps(ovl, p.min=1e-5, p.txt=1e-2, axis=3, xlab=NA, ylab=NA,
                 col=docols, text.cex=.7)
    box()
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    dev.off()


    ovl$p.value <- rbind(aam,aap)
    rownames(ovl$p.value) <- c("more","less")
    ovl$count <- rbind(aac,aac)
    ovl$count[] <- 0
    ovl$count[ovl$p.value<1e-3] <- 1 
    ovl$sign <- rbind(rep( 1,length(aam)),
                      rep(-1,length(aam)))
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_hypergeo")),
            type="png", res=300, width=nw, height=.5+.45)
    par(mai=c(.5,.5,.05,.5), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, axis=2, xlab=NA, ylab=NA,
                 text.cex=.7, col=ttcols)
    axis(1, at=1:length(aac), labels=aac, las=2)
    box()
    abline(v=.5+cdn.cnt, lwd=1, xpd=FALSE)
    dev.off()
    

}



## by AA->AA
for ( ds in uds ) {

    dtmt <- tmtf[tmtf$Dataset==ds,]

    ## NOTE: bg=FALSE, no column-wise background!
    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="to", cols="from",
                       bg=FALSE, value="RAAS", 
                       use.test=use.test, do.plots=FALSE, 
                       verb=0)

    plotProfiles(ovw, fname=file.path(fig.path,paste0("AA_",SETID,"_",ds)),
                 mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=ds,
                 vcols=vcols, vbrks=vbrks,
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
             xlab="IUPRED3", ylab=expression(median~log[10]*RAAS))
    
    plotCor(tmtd$RAAS.median, tmtd$iupred3, xlab="median RAAS", ylab="IUPRED3")

    ovl <- clusterCluster(paste0(tmtd$raas.bins), tmtd$anchor2.bins,
                          cl1.srt=rev(raas.srt))
    par(mai=c(1,1,.5,.5), mgp=c(3,.3,0), tcl=-.25)
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt,
             xlab="ANCHOR2", ylab=expression(median~log[10]*RAAS))
    
    plotCor(tmtd$RAAS.median, tmtd$anchor2,
            xlab=expression(median~log[10]*RAAS), ylab="ANCHOR2")

}

## TODO: IUPRED3 vs. RAAS bins
##ovw <- raasProfile(x=tmtf, id="SAAP", 
##                   rows="iupred3.bins", cols="Dataset",
##                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
##                   col.srt=uds,
##                   use.test=use.test, do.plots=FALSE, xlab="TMT level RAAS",
##                   verb=0)
##

## IUPRED3
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="iupred3.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=expression(TMT~level~log[10]*RAAS),
                   verb=0, 
                   fname=file.path(dpath,paste0("iupred3_",SETID,"_")))

plotProfiles(ovw, fname=file.path(fig.path,paste0("structure_iupred3_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, mtxt="IUPRED3", mtxt.line=3.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## ANCHOR2
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="anchor2.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=expression(median~log[10]*RAAS),
                   verb=0, fname=file.path(fig.path,"anchor2_"))

plotProfiles(ovw, fname=file.path(fig.path,paste0("structure_anchor2_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, mtxt="ANCHOR2", mtxt.line=3.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## S4 PRED SECONDARY STRUCTURE
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="sstruc", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(ssrt), col.srt=uds,
                   use.test=use.test, do.plots=TRUE,
                   xlab=expression(TMT~level~log[10]*RAAS),
                   verb=0, fname=file.path(dpath,"s4pred_"))

plotProfiles(ovw, fname=file.path(fig.path,paste0("structure_s4pred_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, mtxt="SPRED4", mtxt.line=3.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)


dev.off()

if ( FALSE ) {

    ## developing dot plots for combined effect/p plot
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt, txt.col=NA,
                 text.cex=.8, axis=1:2, ylab=NA, xlab=NA,
                 col=ttcols, show.total=TRUE)
    mai=c(.8,.9,.5,.5)
    value <- "median"
    vcols <- vcols
    vbrks <- vbrks
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
    navals <- ovw[[value]]
    navals[] <- NA
    image_matrix(navals, breaks=vbrks,
                 col=vcols, axis=1, xlab=NA, ylab=NA)
    p <- -log10(ovw$p.value)
    p[p>-log10(p.min)] <- -log10(p.min)
    z <- p/-log10(p.min)
    points(x = rep(1:ncol(z), nrow(z)),
           y = rep(nrow(z):1, each= ncol(z)),
           cex=1.5*c(t(z)), pch=19,
           col=num2col(t(ovw[[value]]),limits=range(vbrks),
                       colf=viridis::viridis, n=length(vcols)))
    axis(2, length(axex):1, labels=axex, las=2)
    axis(4, at=nrow(ovw$num.query):1, labels=ovw$num.query[,1],las=2)
    axis(3, at=1:ncol(ovw$num.target), labels=ovw$num.target[1,],las=2)
    if ( !missing(mtxt) ) mtext(mtxt, 2, mtxt.line)
    if ( !missing(llab) ) figlabel(llab, pos="bottomleft", cex=1.2)
    if ( !missing(rlab) ) figlabel(rlab, pos="bottomright", cex=.8)

    
}
