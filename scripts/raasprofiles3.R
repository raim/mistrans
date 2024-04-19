
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

## AA colors; most divergent from
## https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
n <- length(CODL)
set.seed(1)
aa.cols <- sample(color, n)
names(aa.cols) <- names(CODL)
## clustal
aa.cols[] <- "gray"
aa.cols[c("G","P","S","T")] <- "orange"
aa.cols[c("H","K","R")] <- "red"
aa.cols[c("F","W","Y")] <- "blue"
aa.cols[c("I","L","M","V")] <- "green"

aa.pchs <- rep(5, length(aa.cols))
names(aa.pchs) <- names(aa.cols)
aa.pchs[c("G","P","S","T")] <- 1:4
aa.pchs[c("H","K","R")] <- 1:3
aa.pchs[c("F","W","Y")] <- 1:3
aa.pchs[c("I","L","M","V")] <- 1:4

## shapely/rasmol
## via 
aa.cols[c("D","E")]   <- "#E60A0A"  # acidid
aa.cols[c("K","R")]   <- "#145AFF"  # basic
aa.cols[c("S","T")]   <- "#FA9600"  # S/T - P-targets
aa.cols[c("C","M")]   <- "#E6E600"  # sulfur-containing
aa.cols[c("P")]     <- "#DC9682"    # stiff backbone 
aa.cols[c("F","Y")]   <- "#3232AA"  # aromatic
aa.cols[c("W")]     <- "#B45AB4"    # large aromatic
aa.cols[c("H")]     <- "#8282D2"    # aromatic, nitrogen
aa.cols[c("N","Q")]   <- "#00DCDC"  # nitrogen-containing
aa.cols[c("G")]     <- "#EBEBEB"    # just backbone
aa.cols[c("L","V","I")] <- "#0F820F"# branched chaine
aa.cols[c("A")]     <- "#C8C8C8"    # minimal residue

aa.pchs[c("D","E")]   <- c(19,4)
aa.pchs[c("C","M")]   <- c(4,19)
aa.pchs[c("K","R")]   <- c(19,4)
aa.pchs[c("S","T")]   <- c(19,4)
aa.pchs[c("F","Y")]   <- c(19,4)
aa.pchs[c("N","Q")]   <- c(19,4)
aa.pchs[c("G")]     <- 19 #light grey
aa.pchs[c("L","V","I")] <- c(19,2,4)
aa.pchs[c("A")]     <- 19
aa.pchs[c("W")]     <- 4
aa.pchs[c("H")]     <- 19
aa.pchs[c("P")]     <- 4

aa.srt <- rev(names(sort(aa.cols)))

## color by codon frequency
if ( FALSE ) {
    aa.cols <-
        unlist(lapply(CODL, length))
    aa.pchs[] <- 1
}


## plot colors
docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))
gcols <- grey.colors(n=100, start=.9, end=0)

#### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")

mam.path <- "~/data/mammary/"
codon.file <- file.path(mam.path,"processedData","coding_codons.tsv")
dana14.file <- file.path(mam.path,"originalData","dana14_codons.csv")

in.file <- file.path(out.path,"saap_mapped3.tsv")
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")

## protein complexes
humap.file <- file.path(mam.path,"originalData","humap2_complexes_20200809.txt")
corum.file <- file.path(mam.path,"originalData","humanComplexes.txt")

## gene name mappings

## uniprot/ensembl/name mappings
uni2ens.file <- file.path(mam.path,"originalData","uniprot_ensembl.dat")
uni2nam.file <- file.path(mam.path,"originalData","uniprot_name.dat")

## genome feature file
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

### PARAMETERS

## TODO: use for raasProfile and aaProfile calls
p.adjust <- "none" ## multiple hypothesis testing

RAAS.MIN <- -4
RAAS.MAX <-  1
colors <- "arno" # "inferno" #"rocket" # "viridis" # 

COLF <- get(colors)

p.min <- 1e-10
p.txt <- 1e-5

p.dot <- p.min # p.txt
dot.sze <- c(.3,2)

use.test <- t.test # w.test # 

## CANCERS OR HEALTHY TISSUES
healthy <- FALSE # TRUE #  


## DATA FILTERS

## TODO: extracellular is mostly album/globin - analyze
exclude.nterm <- FALSE # TRUE # 
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
    LAB <- paste0(LAB, "+K/R")
}
if (  exclude.nterm ) {
    fig.path <- paste0(fig.path,"_nterm")
    LAB <- paste0(LAB, "-N3")
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

## TODO: INCONSISTENT tags between my protein mapping and Shiri's TMT file
##alls[which(unlist(lapply(alls, function(x) length(unique(x))))>1)]

## remove if tagged so for any dataset
keep <- unlist(lapply(alls, all)) # all keep tags must be tree
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


tmtf$pos <- hdat$pos[idx]
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

tmtf$name <- hdat$name[idx]
tmtf$gene <- hdat$gene[idx]
tmtf$transcript <- hdat$transcript[idx]
tmtf$protein <- hdat$protein[idx]
tmtf$ensembl <- hdat$ensembl[idx]
tmtf$MANE.protein <- hdat$MANE.protein[idx]

tnm <- tmtf$name
tnm[tnm==""] <- tmtf$ensembl[tnm==""]
tmtf$name <- tnm

tmtf$albumin <- hdat$Hemoglobin.Albumin[idx]
tmtf$extracellular <- hdat$extracellular[idx]

tmtf$site <- hdat$site[idx]
tmtf$rsite <- hdat$site[idx]/nchar(hdat$BP)[idx]
tmtf$rsite.bins <- cut(tmtf$rsite, seq(0,1,.2))

### CHOOSE DATASETS
if ( healthy ) {
    ds <- tmtf$Dataset
    tmtf$Dataset[ds=="Healthy"] <- tmtf$TMT.Tissue[ds=="Healthy"]
    tmtf$Dataset[ds!="Healthy"] <- "Cancer"
}

### MEAN AND MEDIAN RAAS

## STATS PER DATA SET and UNIQUE SAAP
usaap <- paste0(tmtf$SAAP,"/",tmtf$BP,"/",tmtf$Dataset)
araas <- split(tmtf$RAAS, usaap)

raas.emedian <- unlist(lapply(araas, median))
raas.median <- unlist(lapply(araas, function(x) {log10(median(10^x))}))
raas.emean <- unlist(lapply(araas, mean))
raas.mean <- unlist(lapply(araas, function(x) {log10(mean(10^x))}))

tmtf$unique <- usaap
tmtf$RAAS.mean <- raas.emean[usaap]
tmtf$RAAS.median <- raas.emedian[usaap]

raas.bins <- cut(tmtf$RAAS.median, breaks=c(-6,-4,-2,-1,0,3))
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
if (  exclude.nterm ) {
    rmsaap <- tmtf$site <4
    tmtf <- tmtf[!rmsaap,]
}


## from this point on we only work with tmtf, which
## now should have all information

## sorting and labels
uds <- sort(unique(tmtf$Dataset))
uds <- c("Healthy",uds[uds!="Healthy"])
uds <- c("Cancer",uds[uds!="Cancer"])
uds <- uds[uds%in%tmtf$Dataset]

## additionally all and cancer-only
auds <- c(uds,"all")
if ( SETID=="cancer" )
    auds <- c(auds, "cancer")

srt <- c("charged","polar","hphobic","special")
srt <- paste(rep(srt,each=4), rep(srt,4), sep=":")

axex <- ftlabels(srt) # axis labels with arrows


## RAAS COLORS
png(file.path(fig.path,paste0("legend_raas.png")),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(tmtf$RAAS,
                          mn=RAAS.MIN, mx=RAAS.MAX,colf=COLF,
                          n=50, plot=TRUE,
                          mai=c(.5,.5,.1,.1),
                          xlab=expression(TMT~level~log[10]*RAAS))
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(colors, pos="bottomleft", cex=1)
figlabel(LAB, pos="bottomright", cex=1)
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
rs <- seq(RAAS.MIN,RAAS.MAX, length.out=6)

pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))

colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- round(rs,1)
    
ovw <- list()
ovw$p.value <- t(10^-pm)
ovw$median <- t(rm)

plab <- expression(log[10]~p)
if ( p.adjust=="q" )
    plab <- expression(log[10]~q)

mai <- c(.5,.5,.1,.1)
fh <- fw <- .2
nh <- nrow(ovw$p.value) *fh + mai[1] + mai[3]
nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]

for ( colstyle in c("viridis","rocket","inferno","arno") ) {

    scols <- get(colstyle, mode="function")(length(vcols))
    
    segmenTools::plotdev(file.path(fig.path,paste0("legend_dotplot_",colstyle)),
                         height=nh, width=nw, res=300)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
    dotprofile(ovw, value="median",
               vbrks=vbrks,
               vcols=scols, 
               dot.sze=dot.sze, p.dot=p.dot, axis=1:2,
               ylab=plab,
               xlab=expression(log[10]~RAAS))
    figlabel(colstyle, pos="bottomleft", cex=1)
    dev.off()
}
segmenTools::plotdev(file.path(fig.path,paste0("legend_dotplot")),
                     height=nh, width=nw, res=300)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovw, value="median",
           vbrks=vbrks,
           vcols=vcols, 
           dot.sze=dot.sze, p.dot=p.dot, axis=1:2,
           ylab=plab,
           xlab=expression(log[10]~RAAS))
figlabel(colors, pos="bottomleft", cex=1)
dev.off()

## legend for dot plot
pp <- seq(0, -log10(p.dot), length.out=3)
rs <- c(-4,-2,0,1) #seq(RAAS.MIN,RAAS.MAX, length.out=3)
pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))
colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- round(rs,1)
ovw <- list(p.value=t(10^-pm),
            median=t(rm))

mai <- c(.4,.5,.05,.06)
fh <- fw <- .2
nh <- nrow(ovw$p.value) *fh + mai[1] + mai[3]
nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]
segmenTools::plotdev(file.path(fig.path,paste0("legend_dotplot_tight")),
                     height=nh, width=nw, res=300)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovw, value="median",
           vbrks=vbrks,
           vcols=vcols, 
           dot.sze=dot.sze, p.dot=p.dot, axis=1:2,
           ylab=plab,
           xlab=NA)
mtext(expression(log[10]~RAAS), 1, 1.1, adj=-.1)
dev.off()

## UNIPROT <-> ENSEMBL MAPPING
## NOTE: ~90 duplicated ensembl IDs
uni2ens <- read.delim(uni2ens.file, header=FALSE)
uni2ens[,2] <- sub("\\..*", "", uni2ens[,2]) # remove ensembl version tag
uni2ens[,1] <- sub("-.*", "", uni2ens[,1]) # remove uniprot version tag
## remove duplicate ensembl IDs
## TODO: is the list sorted? best uniprot hit?
uni2ens <- uni2ens[!duplicated(uni2ens[,2]),]

## UNIPROT <-> NAME MAPPING
## TODO: add MANE column
uni2nam <- read.delim(uni2nam.file, header=FALSE)
uni2nam[,1] <- sub("-.*", "", uni2nam[,1]) # remove uniprot version tag

## ENSEMBL <-> NAME MAPPING
## TODO: add MANE column
genes <- read.delim(feature.file)
genes <- genes[genes$proteins!="" & !is.na(genes$proteins),]
genes <- genes[genes$name!="" & !is.na(genes$name),]
ptl <- strsplit(genes$proteins, ";")
ens2nam <- rep(genes$name, lengths(ptl))
names(ens2nam) <- unlist(ptl)
pnms <- ens2nam


### START ANALYSIS

## PROTEINS

## unique BP/SAAP specific RAAS
source("~/work/mistrans/scripts/saap_utils.R")

## median RAAS per unique protein position
posl <- split(tmtf$RAAS, paste(tmtf$ensembl, tmtf$pos))
posl <- as.data.frame(listProfile(posl, y=tmtf$RAAS))
colnames(posl) <- paste0("RAAS.", colnames(posl))
posl$ensembl <- sub(" .*", "", rownames(posl))

## MEDIAN RAAS FOR EACH UNIQUE PROTEIN POSITION
posl <- split(posl$RAAS.median, posl$ensembl)

## protein median of median RAAS
pbstat <- listProfile(posl, y=tmtf$RAAS)


## volcano
covl <- list()
covl$median <- as.matrix(pbstat$median)
covl$p.value <- pbstat[,"p",drop=FALSE]
## replace by gene names
##rownames(covl$p.value) <-
##    pnms[rownames(covl$p.value)]

dense2d(pbstat$median, log10(pbstat$n))

plotdev(file.path(fig.path,paste0("volcano_proteins_unique")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.1,.5), mgp=c(1.2,.3,0), tcl=-.25)
res <- volcano(covl, value="median",
               p.txt=-log10(0.01), v.txt=c(-2,-2), cut=5, mid=-2,
               ids=pnms, xlab=expression(log[10](bar(RAAS))))
dev.off()


## ALL RAAS per protein
ptl <- split(tmtf$RAAS, tmtf$ensembl)
ptstat <- listProfile(ptl, y=tmtf$RAAS)

## TODO: nice plots, try to use raasProfiler!
covl <- list()
covl$median <- as.matrix(ptstat$median)
covl$p.value <- ptstat[,"p",drop=FALSE]

plotdev(file.path(fig.path,paste0("volcano_proteins")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.1,.5), mgp=c(1.2,.3,0), tcl=-.25)
res <- volcano(covl, value="median",
               p.txt=12, v.txt=c(-2,-1), cut=50, mid=0,
               ids=pnms, xlab=expression(log[10](bar(RAAS))))
abline(v=0)
dev.off()

## PROTEIN WINDOWS - 50 AA
plen <- dat$len[!duplicated(dat$ensembl)]
names(plen) <- dat$ensembl[!duplicated(dat$ensembl)]
plen <- plen[!is.na(plen)]

## generate protein windows of length 50
## NOTE: last window to end is <=50
pwins <- sapply(plen, function(x) unique(c(seq(0,x,50),x)))
pwins <- lapply(pwins, function(x) cbind(start=1+x[1:(length(x)-1)],
                                         end=x[2:length(x)]))
pwins <- data.frame(ID=rep(names(pwins), lengths(pwins)/2),
                    do.call(rbind, pwins))
## TODO: collect unique and site-specific median RAAS for each window
## and do tests.


## PROTEIN COMPLEXES

## TODO: first boil down RAAS to median per site to avoid
## over-estimation by single sites, eg. Q->G in PSMA1.

## rename names by ensembl in position-wise RAAS medians

## distribution in protein complexes
humap <- read.csv(humap.file)
hulst <- strsplit(humap[,3]," ")
huids <- humap[,1]

humap <- read.delim(corum.file)
hulst <- strsplit(humap$subunits.UniProt.IDs.,";")
huids <- sub(" complex$","",humap[,2])

names(hulst)<- huids

## replace complex uniprot genes by ALL ensembl proteins that map to it
huens <- lapply(hulst, function(x) {
    unique(uni2ens[uni2ens[,1]%in%x,2])
})


## COMPLEXES - MEDIAN SITE RAAS
## replace complex ensembl IDs with rows in RAAS table
huraas <- lapply(huens, function(x) c(na.omit(pbstat[x,"median"])))

## RAAS statistics
hustat <- listProfile(huraas, y=tmtf$RAAS)
hustat <- as.data.frame(hustat)

## TODO: nice plots, try to use raasProfiler!
covl <- list()
covl$median <- as.matrix(hustat$median)
covl$p.value <- hustat[,"p",drop=FALSE]

plotdev(file.path(fig.path,paste0("volcano_corum_unique")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.1,.5), mgp=c(1.2,.3,0), tcl=-.25)
res <- volcano(covl, value="median",
               p.txt=3, v.txt=c(-Inf,-1.5), cut=50, mid=-1,
               xlab=expression(log[10](bar(RAAS))))
abline(v=0)
dev.off()


## COMPLEXES - UNIQUE RAAS
## replace complex ensembl IDs with rows in RAAS table
huraas <- lapply(huens, function(x) tmtf$RAAS[tmtf$ensembl%in%x])

## RAAS statistics
hustat <- listProfile(huraas, y=tmtf$RAAS)
hustat <- as.data.frame(hustat)

## TODO: nice plots, try to use raasProfiler!
covl <- list()
covl$median <- as.matrix(hustat$median)
covl$p.value <- hustat[,"p",drop=FALSE]

plotdev(file.path(fig.path,paste0("volcano_corum")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.1,.5), mgp=c(1.2,.3,0), tcl=-.25)
res <- volcano(covl, value="median",
               p.txt=20, v.txt=c(-Inf,-1), cut=50, mid=0,
               xlab=expression(log[10](bar(RAAS))))
abline(v=0)
dev.off()

## NOTE: EZR, MSN, 

## investigate
humap[grep("Spliceosome, E", humap$ComplexName),c("ComplexName",
                                               "Complex.comment",
                                               "subunits.Gene.name.")]
humap[grep("Drosha", humap$ComplexName),c("ComplexName",
                                               "Complex.comment",
                                               "subunits.Gene.name.")]

## GLOBAL DISTRIBUTION BY CANCER TYPE
           
ylm <- range(tmtf$RAAS)
plotdev(file.path(fig.path,paste0("RAAS_distribution")),
        type="png", res=300, width=3,height=2)
par(mai=c(.7,.5,.1,.1), mgp=c(1.2,.3,0), tcl=-.25, xaxs="i")
boxplot(tmtf$RAAS ~ factor(tmtf$Dataset, levels=uds), ylim=ylm, las=2,
        xlab=NA, ylab=expression(log[10]*RAAS[TMT]), cex=.5, pch=19,
        pars=list(outcol="#00000055"), axes=FALSE)
for ( ax in c(2,4) ) axis(ax)
axis(1, at=1:length(uds), labels=uds, las=2)
figlabel(LAB, pos="bottomright", cex=.9)
dev.off()

plotdev(file.path(fig.path,paste0("RAAS_distribution_density")),
        type="png", res=300, width=3,height=2)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
plot(density(tmtf$RAAS), ylim=c(0,.5), col=NA,
     xlab=expression(log[10]*RAAS[TMT]), main=NA, xlim=ylm)
for ( i in seq_along(uds) )
    lines(density(tmtf$RAAS[tmtf$Dataset==uds[i]]), col=i, lwd=2)
legend("topright", uds, col=seq_along(uds), lty=1, seg.len=.5, lwd=2, bty="n",
       ncol=1, cex=.8, y.intersp=.9)
figlabel(LAB, pos="bottomleft", cex=.8)
dev.off()

tmtu <- tmtf[!duplicated(tmtf$unique),]
plotdev(file.path(fig.path,paste0("RAAS_distribution_unique")),
        type="png", res=300, width=3,height=2)
par(mai=c(.7,.5,.1,.1), mgp=c(1.2,.3,0), tcl=-.25, xaxs="i")
boxplot(tmtu$RAAS.median ~ factor(tmtu$Dataset, levels=uds), ylim=ylm, las=2,
        xlab=NA, ylab=expression(log[10]*bar(RAAS[unique])),
        cex=.5, pch=19, pars=list(outcol="#00000055"), axes=FALSE)
for ( ax in c(2,4) ) axis(ax)
axis(1, at=1:length(uds), labels=uds, las=2)
figlabel(LAB, pos="bottomright", cex=.9)
dev.off()

plotdev(file.path(fig.path,paste0("RAAS_distribution_unique_density")),
        type="png", res=300, width=3,height=2)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
plot(density(tmtu$RAAS.median), ylim=c(0,.5), col=NA,
     xlab=NA, main=NA, xlim=ylm)
mtext(expression(log[10]*bar(RAAS[unique])), 1, 1.5)
for ( i in seq_along(uds) )
    lines(density(tmtu$RAAS.median[tmtu$Dataset==uds[i]]), col=i, lwd=2)
legend("topright", uds, col=seq_along(uds), lty=1, seg.len=.5, lwd=2, bty="n",
       ncol=1, cex=.8, y.intersp=.9)
figlabel(LAB, pos="bottomleft", cex=.8)
dev.off()

## TODO: pairs of same SAAP
##tmtl <- split(tmtf$Dataset, tmtf$SAAP)

## RAAS profiles by AAS site in peptides
## TODO: median RAAS of unique BP/SAAP vs. rel.position,
##       similar to iupred3
for ( ds in auds ) {

    tmtd <- tmtf
    if ( ds!="all" ) {
        if ( ds!="cancer" )
            tmtd <- tmtf[tmtf$Dataset==ds,]
        else
            tmtd <- tmtf[tmtf$Dataset!="Healthy",]
    }
    plotdev(file.path(fig.path,paste0("AASsite_",SETID,"_cor_",ds)),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tmtd$RAAS.median, tmtd$rsite,
            xlab=NA,
            ylab="relative AAS position in peptide")
    mtext(expression(log[10]*bar(RAAS[unique])), 1, 1.6)
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    dev.off()
}

fname <- file.path(dpath,paste0("AASsite_",SETID,"_wtests_"))
ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
                    delog=TRUE, rows="rsite.bins", cols="Dataset",
                    bg=TRUE, row.srt=rev(rsrt), col.srt=uds,
                    use.test=use.test, do.plots=FALSE,
                    xlab=expression(TMT~level~log[10]*RAAS),
                    fname=fname, verb=0)
plotProfiles(ovw, fname=file.path(fig.path,paste0("AASsite_",SETID)),
             mai=c(.8,1,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             ttcols=ttcols, value="median",
             mtxt="relative AAS position", mtxt.line=3.5,
             rlab=LAB, llab="",
             vcols=vcols, vbrks=vbrks,
             gcols=gcols, verb=1)
  
## RAAS profiles by AA properties
fname <- file.path(dpath,paste0("AAprop_",SETID,"_wtests_"))
ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
                    delog=TRUE, rows="pfromto", cols="Dataset",
                    bg=TRUE, row.srt=srt, col.srt=uds,
                    use.test=use.test, do.plots=TRUE,
                    xlab=expression(TMT~level~log[10]*RAAS),
                    fname=fname, verb=0)
ovwp <- sortOverlaps(ovw, axis=2, p.min=p.min, cut=TRUE)
ovwp <- sortOverlaps(ovwp, axis=2, srt=sort(rownames(ovwp$p.value)))



## plot all
plotProfiles(ovw, fname=file.path(fig.path,paste0("AAprop_",SETID)),
             mai=c(.7,1.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="",
             vcols=vcols, vbrks=vbrks,
             gcols=gcols, verb=1)
if ( nrow(ovwp$p.value)>0 )
    plotProfiles(ovwp,
                 fname=file.path(fig.path,paste0("AAprop_",SETID,"_cut")),
                 mai=c(.7,1.5,.5,.5), ttcols=ttcols, value="median",
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

## by AA property transition
for ( ptype in unique(tmtf$pfromto) ) {

    ## sort by to
    qsrt <- sort(unique(tmtf$fromto[tmtf$pfromto==ptype]))
    ft <- unlist(lapply(strsplit(qsrt,":"), function(x) x[2]))
    qsrt <- qsrt[order(ft)]
                 
    ovwp <- sortOverlaps(ovw, srt=qsrt)

    plotProfiles(ovwp, fname=file.path(fig.path,paste0(ptype,"_",SETID)),
                 mai=c(.8,1.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="", mtxt=ptype, mtxt.line=3,
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols)
    ## cut
    ovwp <- sortOverlaps(ovwp, axis=2, p.min=p.txt, cut=TRUE)
    if ( nrow(ovwp$p.value)>1 ) #TODO: fix error with nrow=1 in 
        plotProfiles(ovwp,
                     fname=file.path(fig.path,paste0(ptype,"_",SETID,"_cut")),
                     mai=c(.8,1.5,.5,.5), ttcols=ttcols, value="median",
                     p.min=p.min, p.txt=p.txt,
                     dot.sze=dot.sze, p.dot=p.dot,
                     rlab=LAB, llab="",mtxt=ptype, mtxt.line=3,
                     vcols=vcols, vbrks=vbrks,
                     gcols=gcols)
    
}


## plot only significantly high RAAS
ovwp <-ovw
##ovwp$p.value[ovw$sign<0] <- 1
source("~/work/mistrans/scripts/saap_utils.R")
## only significant
ovwp <- sortOverlaps2(ovwp, axis=2, p.min=p.min, cut=TRUE)
## re-sort by AA properties
newsrt <- sort(rownames(ovwp$p.value))
ovwp <- sortOverlaps2(ovwp, axis=2, srt=newsrt)
par(family="sans") 
if ( nrow(ovwp$p.value)>0 )
    plotProfiles(ovwp,
                 fname=file.path(fig.path,paste0("AA_",SETID,"_cut")),
                 mai=c(.7,.5,.5,.5), fw=.2, fh=.2,
                 ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="", ffam="monospace",
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols)

## TODO: plot on p-value correction
plot(ovw$p.value, qvalue::qvalue(c(ovw$p.value))$qvalues)


## by AA->AA
for ( ds in auds ) {

    tmtd <- tmtf
    if ( ds!="all" ) {
        if ( ds!="cancer" )
            tmtd <- tmtf[tmtf$Dataset==ds,]
        else
            tmtd <- tmtf[tmtf$Dataset!="Healthy",]
    }

    ## NOTE: bg=FALSE, no column-wise background!
    ovw <- raasProfile(x=tmtd, id="SAAP", 
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





## BY CODON->AA
## TODO:
## * expand ovw to ALL codons and AA to align plots,
## * plot log2 odds ratio vs. 
ctmt <- tmtf[tmtf$codon!="",]

## load global gene-wise codon counts
codons <- read.delim(codon.file, row.names=1)

## Dana and Tuller 2014/2015
## decoding time/sec -> 1/decoding rate
if ( file.exists(dana14.file) )
    decode <- 1/read.csv(dana14.file,
                         row.names=1)[,"H..sapiens5.HEK293",drop=FALSE]

## global analysis: codon frequency vs. RAAS
## CODON COUNTS IN ALL UNIQUE TRANSCRIPTS 
cod  <- codons[unique(ctmt$transcript),]
codt <- apply(cod,2,sum)
## transcript codon frequencies per AA
codl <- split(codt, GENETIC_CODE[sub(".*\\.","",names(codt))])
    codl <- lapply(codl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
    codl <- codl[names(aaprop)[names(aaprop)%in%names(codl)]] ## SORT BY AA PROP
codf <- lapply(codl, function(x) x/sum(x)) # codon frequency
cods <- unlist(codf)
names(cods) <- sub("\\.","-",names(cods))
## GLOBAL CODON SORTING by background frequencies
## (AA property ->codon frequency)
CODSRT <- sub("\\.","-",names(unlist(codl)))
if ( !include.kr )
    CODSRT <-
        CODSRT[grep("[KR]", CODSRT, invert=TRUE)]

## codon rank - TODO: remove or improve this
codr <- lapply(codl, function(x) round(6*rank(x)/length(x))) # codon rank
codr$W[1] <- 0
codr$M[1] <- 0
## codon rank
codr <- unlist(codr)
names(codr) <- sub("\\.", "-", names(codr))
ctmt$codon.rank <- as.character(codr[ctmt$aacodon])
    
## codon frequency bins
ctmt$codon.bins <- cut(cods[ctmt$aacodon], seq(0,1,.1))

## TODO: do by 2|3, 4 or 6 codons
cod.bins <- cut(cods[ctmt$aacodon], seq(0,1,.1))
ncod <- unlist(lapply(CODL, length))


if ( interactive() & exists("decode", mode="numeric") ) {
    ## TODO: analyze per AA - positive trend
    ## globally: negative trend
    plotCor(cods, decode[sub(".*-","", names(cods)),1], density=FALSE,
            col=NA, xlab="codon frequency, all AAS transcripts",
            ylab=expression(decoding~rate(codons/s)))
    points(cods, decode[sub(".*-","", names(cods)),1], lwd=1, cex=1,
           col=aa.cols[sub("-.*","",names(cods))],
           pch=aa.pchs[sub("-.*","",names(cods))])
    legend("top",unique(sub("-.*","",names(cods))),
           col=aa.cols[unique(sub("-.*","",names(cods)))],
           pch=aa.pchs[unique(sub("-.*","",names(cods)))],
           pt.cex=1, ncol=2, cex=.8, y.intersp=.75,
           bty="n")
}


## remove single codon!
ltmt <- ctmt#[!ctmt$from%in%c("W","M","H","D","E","N","Q","F","Y"),]

## NOTE: bg=FALSE, no column-wise background !
ovw <- raasProfile(x=ltmt, id="SAAP", 
                   rows="codon.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   row.srt=rev(levels(ctmt$codon.bins)),
                   use.test=use.test, do.plots=FALSE,
                   fname=file.path(dpath,paste0("codon_bins_",SETID,"_")),
                   xlab=expression(TMT~level~log[10]*RAAS),
                   verb=1)
plotProfiles(ovw, fname=file.path(fig.path,paste0("codon_bins_",SETID)),
             mai=c(.8,1,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             mtxt="codon frequency", mtxt.line=3.5,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="", vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## CODON FREQUENCY RANK
ovw <- raasProfile(x=ltmt, id="SAAP", 
                   rows="codon.rank", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   ##row.srt=rev(levels(ctmt$codon.bins)),
                   use.test=use.test, do.plots=FALSE,
                   fname=file.path(dpath,paste0("codon_rank_",SETID,"_")),
                   xlab=expression(TMT~level~log[10]*RAAS),
                   verb=1)
plotProfiles(ovw, fname=file.path(fig.path,paste0("codon_rank_",SETID)),
             mai=c(.8,.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             mtxt="codon rank",# mtxt.line=3.5,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="", vcols=vcols, vbrks=vbrks,
             gcols=gcols)



## per codon profiles
for ( ds in auds ) {

    ## TODO: add codon frequency two-sided bar plot

    dtmt <- ctmt
    if ( ds!="all" ) {
        if ( ds!="cancer" )
            dtmt <- ctmt[ctmt$Dataset==ds,]
        else
            dtmt <- ctmt[ctmt$Dataset!="Healthy",]
    }

    ## NOTE: bg=FALSE, no column-wise background !
    dtmt$all <- "all"
    ova <- raasProfile(x=dtmt, id="SAAP", 
                       rows="all", cols="aacodon",
                       col.srt=CODSRT, filter=FALSE,
                       bg=FALSE, use.test=use.test, do.plots=FALSE, 
                       verb=0)

    ## NOTE: bg=FALSE, no column-wise background !
    ## TODO: move this outside loop and use rows as bg
    ovd <- raasProfile(x=dtmt, id="SAAP", 
                       rows="Dataset", cols="aacodon",
##                       row.srt=uds,
                       col.srt=CODSRT, filter=FALSE,
                       bg=TRUE, bg.dir="row",
                       use.test=use.test, do.plots=FALSE, 
                       verb=0)
    nr <- uds[uds%in%rownames(ovd$p.value)]
    if ( length(nr)>1 )
        ovd <- sortOverlaps(ovd, axis=2, srt=nr)



    ## NOTE: bg=FALSE, no column-wise background !
    ovw <- raasProfile(x=dtmt, id="SAAP", 
                       rows="to", cols="aacodon",
                       col.srt=CODSRT, filter=FALSE,
                       bg=FALSE, use.test=use.test, do.plots=FALSE, 
                       verb=0)


    
    ## sort by amino acid property via shapely colors!
    ##aasrt <- names(aaprop)[names(aaprop)%in%rownames(ovw$p.value)]
    aasrt <- aa.srt[aa.srt%in%rownames(ovw$p.value)]
    ovw <- sortOverlaps(ovw, axis=2, srt=aasrt, cut=TRUE)


    ## CUSTOM CODON FREQUENCY PLOTS

    ## CODON COUNTS AT THE AAS OF UNIQUE SAAP
    lcodt <- apply(ovw$unique, 2, sum)
    names(lcodt) <- sub(".*-","",names(lcodt))
    ## AAS codon frequencies per AA
    lcodl <- split(lcodt, GENETIC_CODE[names(lcodt)])
    
    lcodl <- lapply(lcodl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
    lcodl <- lcodl[names(aaprop)] ## ##SORT BY AA PROP

    ## AA-specific codon frequency
    lcods <- unlist(lapply(lcodl, function(x) x/sum(x)))

    ## CODON COUNTS IN ALL UNIQUE TRANSCRIPTS 
    cod  <- codons[unique(dtmt$transcript),]
    codt <- apply(cod,2,sum)
    ## transcript codon frequencies per AA
    codl <- split(codt, GENETIC_CODE[sub(".*\\.","",names(codt))])
    codl <- lapply(codl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
    codl <- codl[names(aaprop)[names(aaprop)%in%names(codl)]] ## SORT BY AA PROP
    cods <- unlist(lapply(codl, function(x) x/sum(x)))

    ## sort by background (AA property ->codon frequency)
    cdsrt <- sub("\\.","-",names(unlist(codl)))

    ## filter available
    cdsrt <- cdsrt[cdsrt%in%colnames(ovw$p.value)]

    ## sort overlaps
    ##ovw <- sortOverlaps(ovw, axis=1, srt=cdsrt, cut=FALSE)
    ##ovd <- sortOverlaps(ovd, axis=1, srt=cdsrt, cut=FALSE)
    ##ova <- sortOverlaps(ova, axis=1, srt=cdsrt, cut=FALSE)
    
    ## codon frequency bins
    codon.bins <- cut(lcods[dtmt$aacodon], seq(0,1,.2))

    ## codon class lines
    cdn.cnt <- table(sub("-.*","",colnames(ovw$p.value)))[names(aaprop)]
    cdn.cnt <- cumsum(cdn.cnt[!is.na(cdn.cnt)])
    
    ## SORT CODON FREQUENCIES as in main plots
    ## TODO: constant sorting of  ALL codons 
    names(cods) <- sub("\\.","-", names(cods))
    names(lcods) <- sub("\\.","-", names(lcods))
    lcods <- lcods[cdsrt]
    cods <- cods[cdsrt]

    ## log2 ratio of frequencies
    lratio <- (lcods/cods[names(lcods)]) # log2 ratio of codon frequencies
    rlm <- max(abs(lratio),na.rm=TRUE)

 
    
    plotProfiles(ova, fname=file.path(fig.path,paste0("codon_",
                                                      SETID,"_",ds,
                                                      "_all")),
                 fw=.2, mai=c(.05,.5,.05,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=ds, mtxt="",
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols, col.lines=cdn.cnt)


    plotProfiles(ovd, fname=file.path(fig.path,paste0("codon_",
                                                      SETID,"_",ds,"_Dataset")),
                 fw=.2, mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=ds, mtxt="",
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols, col.lines=cdn.cnt)


    ## TODO: re-create previous codon plots
    plotProfiles(ovw, fname=file.path(fig.path,paste0("codon_",SETID,"_",ds)),
                 fw=.2, mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=ds, mtxt="substituted AA",
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols, col.lines=cdn.cnt)

### CODON FREQUENCY ANALYSIS

    mai <- c(.05,.5,.1,.5)
    fw <- .2
    nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]

    ## sort by amino acid property via shapely colors!
    ##aasrt <- names(aaprop)[names(aaprop)%in%rownames(ovw$p.value)]
    aasrt <-
        aa.srt[aa.srt%in%sub("-.*","",names(cods))]
 
 
    ## relation codon base frequency and ratio AAS
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_ratio")),
           type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(cods, lratio, density=FALSE,
            ylab=expression(f[AAS]/f[bg]),
            xlab=expression(f[bg]), xlim=c(0,1),
            col=NA,pch=19, lwd=2, cex=.6)
    points(cods, lratio, lwd=2, cex=1,
           col=aa.cols[sub("-.*","",names(cods))],
           pch=aa.pchs[sub("-.*","",names(cods))])
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    ## higher frequency tends to be above diagonal
    ## TODO: better AA colors
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies")),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(cods, lcods, density=FALSE,
            ylab=expression(codon~frequency~f[AAS]),
            xlab=expression(codon~frequency~f[bg]), xlim=c(0,1), ylim=c(0,1),
            col=NA,pch=19, lwd=2, cex=.6)
    points(cods, lcods, lwd=2, cex=1,
           col=aa.cols[sub("-.*","",names(cods))],
           pch=aa.pchs[sub("-.*","",names(cods))])
    ##abline(a=0,b=1, col=2)
    legend("bottomright",aasrt,
           col=aa.cols[aasrt],
           pch=aa.pchs[aasrt],
           pt.cex=1, lwd=2, lty=NA, seg.len=0, ncol=2, cex=.8, y.intersp=.75,
           bty="n")
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_diff")),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(cods, lcods-cods, density=FALSE,
         ylab=expression(f[AAS]-f[bg]),
         xlab=expression(f[bg]), xlim=c(0,1),
         col=NA,pch=19, lwd=2, cex=.6)
    points(cods, lcods-cods, lwd=2, cex=1,
           col=aa.cols[sub("-.*","",names(cods))],
           pch=aa.pchs[sub("-.*","",names(cods))])
    if ( FALSE )
        legend("bottomright",aasrt,
               col=aa.cols[aasrt],
               pch=aa.pchs[aasrt],
               pt.cex=1, lwd=2, lty=NA, ncol=2, cex=.6, y.intersp=.75,
               bty="n")
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()
    
    ## TODO: get median RAAS for each codon and plot against codon frequency

    ram <- sapply(cdsrt, function(cl) median(dtmt$RAAS[dtmt$aacodon==cl]))
    ram <- ram[names(lcods)]
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_raas")),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(ram, lcods, density=FALSE, xlab=expression(log[10]*bar(RAAS)),
            ylab=expression(codon~frequency~f[AAS]), col=NA)
    points(ram, lcods, lwd=2, cex=1,
           col=aa.cols[sub("-.*","",names(cods))],
           pch=aa.pchs[sub("-.*","",names(cods))])
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_raas_nocol")),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(ram, lcods, density=FALSE, xlab=expression(log[10]*bar(RAAS)),
            ylab=expression(codon~frequency~f[AAS]), col="#00000099", pch=19,
            lwd=0,cex=1)
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_raas_slopes")),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(ram, lcods, density=FALSE, xlab=expression(log[10]*bar(RAAS)),
         ylab=expression(codon~frequency~f[AAS]), col=NA)
    slps <- rep(NA, length(CODL))
    for ( k in 1:length(CODL) ) {
        aa <- names(CODL)[k]
        cds <- paste0(aa,"-",CODL[[k]])
        if ( sum(cds%in%names(lcods))<2 ) next
        rm <- ram[cds]
        cd <- lcods[cds]

        cd <- cd[!is.na(rm)]
        rm <- rm[!is.na(rm)]
        if ( length(rm)<1 ) next
  
        points(rm, cd, lwd=2, cex=1,
               col=aa.cols[aa], pch=aa.pchs[aa])
        if ( length(rm)<2 ) next
        cdf <- lm(cd ~rm)
        slps[k] <- coef(cdf)[2]
        lines(rm, predict(cdf), col=aa.cols[aa])
    }
    dev.off()
    
    ## plot slopes
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_raas_slopes_bynum")),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(unlist(lapply(CODL,length)), slps, col=aa.cols[names(CODL)],
         pch=aa.pchs[names(CODL)], lwd=2,
         xlab="number of codons per/AA")
    abline(h=0)
    dev.off()
    plotdev(file.path(fig.path,
                      paste0("codon_",SETID,"_",ds,
                             "_codons_frequencies_raas_slopes_bynum_zoom")),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(unlist(lapply(CODL,length)), slps, col=aa.cols[names(CODL)],
         pch=aa.pchs[names(CODL)], lwd=2, ylim=c(-1.5,0),
         xlab="number of codons per/AA")
    abline(h=0)
    dev.off()
    
    ## TODO: mean has no correlation to RAAS, median does!
    if ( interactive() ) {
        plotCor(ova$median[,names(lcods)], lcods, density=FALSE, col=NA)
        points(ova$median[,names(lcods)], lcods, lwd=1, cex=1,
               col=aa.cols[sub("-.*","",names(cods))],
               pch=aa.pchs[sub("-.*","",names(cods))])
        plotCor(ova$mean[,names(lcods)], lcods, density=FALSE, col=NA)
        points(ova$mean[,names(lcods)], lcods, lwd=1, cex=1,
               col=aa.cols[sub("-.*","",names(cods))],
               pch=aa.pchs[sub("-.*","",names(cods))])
        plotCor(ova$mean[,names(lcods)], ova$median[,names(lcods)],
                density=FALSE, col=NA)
        points(ova$mean[,names(lcods)], ova$median[,names(lcods)], lwd=1, cex=1,
               col=aa.cols[sub("-.*","",names(cods))],
               pch=aa.pchs[sub("-.*","",names(cods))])
    }

    if ( exists("decode", mode="numeric") ) {
        plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                          "_codons_frequencies_raas_dana14")),
                type="png", res=300, width=3,height=3)
        par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
        plotCor(ram, decode[sub(".*-","",names(ram)),1],
                density=FALSE, xlab=expression(log[10]*bar(RAAS)),
                ylab=expression(decoding~rate/(codons/s)), col=NA)
        points(ram, decode[sub(".*-","",names(ram)),1], lwd=2, cex=1,
               col=aa.cols[sub("-.*","",names(cods))],
               pch=aa.pchs[sub("-.*","",names(cods))])
        figlabel(ds, pos="bottomleft", font=2, cex=1.2)
        figlabel(LAB, pos="bottomright", cex=.7)
        for ( ax in 3:4 )  axis(ax, labels=FALSE)
        dev.off()
    }

    ramm <- matrix(ram,nrow=1)
    colnames(ramm) <- names(ram)
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_raas")),
            type="png", res=300, width=nw, height=.25+.1)
    par(mai=c(.05,.5,.05,.5), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    image_matrix(ramm, col=vcols, breaks=vbrks,
                 axis=NA, las=2, ylab=NA)
    mtext(expression(bar(RAAS)), 2, .15, las=2, cex=1)
    box()
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    abline(v=.5+cdn.cnt, lwd=1, xpd=FALSE, col="white")
    dev.off()
    
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequency_raas")),
            type="png", res=300, width=nw,height=3)
    par(mai=c(.8,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    boxplot(dtmt$RAAS ~ factor(dtmt$aacodon, levels=cdsrt),
            las=2, xlab="", ylab=expression(TMT~level~log[10]*RAAS),
            cex=.5, pch=19, pars=list(outcol="#00000055"))
    axis(3, at=1:length(cdsrt), labels=round(lcods[cdsrt],2), las=2)
    axis(4)
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_ratio_old")),
            type="png", res=300, width=nw,height=.75)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plot(1:length(lcods), lratio, type="h", axes=FALSE, xlab=NA, lwd=3,
         ylab=expression(lg[2]~r),
         xlim=c(0.5,length(lcods)+.5), col=2, ylim=c(-rlm,rlm))
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
    bp <- barplot(lratio, axes=FALSE, xlab=NA, beside=TRUE,
                  ylab=expression(lg[2]~r), las=2, xaxt='n', 
                  width=.5, space=.5, ylim=c(-rlm,rlm))
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
    aam <- aam[,CODSRT,drop=FALSE]
    aap <- aap[,CODSRT,drop=FALSE]
    aac <- aac[,CODSRT,drop=FALSE]

    ## construct overlap class
    ovl <- list()
    pvl <- aam
    pvl[aap<aam] <- aap[aap<aam]
    ovl$p.value <- pvl
    ovl$count <- aac
    ovl$sign <- ifelse(aam<aap, 1, -1)
    
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_hypergeo")),
            type="png", res=300, width=nw, height=.25+.1)
    par(mai=c(.05,.5,.05,.5), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, axis=NA, xlab=NA, ylab=NA,
                 text.cex=.7, col=ttcols)
    ##axis(1, at=1:length(aac), labels=aac, las=2)
    mtext(expression(f), 2, .25, las=2)
    box()
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    dev.off()

    ## two roes
    ovl$p.value <- rbind(aam,aap)
    rownames(ovl$p.value) <- c("more","less")
    ovl$count <- rbind(aac,aac)
    ovl$count[1,aam>=aap] <- 0
    ovl$count[2,aam< aap] <- 0 
    ovl$sign <- rbind(rep( 1,length(aam)),
                      rep(-1,length(aam)))
    plotdev(file.path(fig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_hypergeo2")),
            type="png", res=300, width=nw, height=.35+.1)
    par(mai=c(.05,.5,.05,.5), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, axis=2, xlab=NA, ylab=NA,
                 text.cex=.7, col=ttcols)
    ##axis(1, at=1:length(aac), labels=aac, las=2)
    ##axis(1, at=1:length(aac), labels=aac, las=2)
    ##mtext(expression(f), 2, .25, las=2)
    box()
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    dev.off()
    
}




### BY STRUCTURAL FEATURES
for ( ds in auds ) {

    tmtd <- tmtf
    if ( ds!="all" ) {
        if ( ds!="cancer" )
            tmtd <- tmtf[tmtf$Dataset==ds,]
        else
            tmtd <- tmtf[tmtf$Dataset!="Healthy",]
    }
    ovl <- clusterCluster(tmtd$iupred3.bins, paste0(tmtd$raas.bins), 
                          cl1.srt=c(rev(levels(iupred3.bins)),"na"),
                          cl2.srt=raas.srt)
    plotdev(file.path(fig.path,paste0("structure_iupred3_",SETID,"_ovl_",ds)),
            type="png", res=300, width=3,height=3)
    par(mai=c(.9,.9,.5,.5), mgp=c(3.3,.3,0), tcl=-.25)
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt,
                 xlab=NA, ylab=NA,
                 show.total=TRUE)
    mtext("disordered score", 2, 3.5)
    mtext(expression(log[10]*bar(RAAS[unique])), 1, 3.5)
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    dev.off()
    
    plotdev(file.path(fig.path,paste0("structure_iupred3_",SETID,"_cor_",ds)),
            type="png", res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tmtd$RAAS.median, tmtd$iupred3,
            xlab=NA, ylab="disordered score")
    mtext(expression(log[10]*bar(RAAS[unique])), 1, 1.6)
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    dev.off()
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
             rlab=LAB, mtxt="disordered score", mtxt.line=3.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## ANCHOR2
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="anchor2.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=expression(log[10]*bar(RAAS)),
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
                   verb=0, fname=file.path(dpath,"structure_s4pred_"))

plotProfiles(ovw, fname=file.path(fig.path,paste0("structure_s4pred_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, mtxt="SPRED4", mtxt.line=3.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)



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
