
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

## axis labels
xl.hlfm <- expression(protein~"half-life"/h)
xl.hlf <- expression(protein~"half-life"/h)

## overrule specific y-axis label
xl.prota <- xl.raas

### ADDITIONAL DATA
## protein half-lives - @Mathieson2018
math18.file <- file.path(mam.path,"originalData",
                         "41467_2018_3106_MOESM5_ESM.xlsx")

## 20S targets - @Pepelnjak2024
pepe24.file <- file.path(mam.path,"originalData",
                         "44320_2024_15_moesm1_esm.xlsx")

##  @Watson2023 - T/osmo
## NOTE/TODO: unused since this is data from mouse cells;
## perhaps useful when mapped to human orthologs.
w23prot.file <- file.path(mam.path,"originalData",
                          "41586_2023_6626_MOESM4_ESM.xlsx")
w23pprot.file <- file.path(mam.path,"originalData",
                           "41586_2023_6626_MOESM5_ESM.xlsx")

## @Yang2022: thermal stability prediction
## https://structure-next.med.lu.se/ProTstab2/
protstab.file <- file.path(mam.path,"originalData",
                           "ProTstab2_human.csv")

## @Savitski2014: Tracking cancer drugs in living cells by thermal
## profiling of the proteome
thermo.file <- file.path(mam.path, "originalData",
                         "savitski14_tableS11.xlsx")
thatp.file <- file.path(mam.path, "originalData",
                         "savitski14_tableS3.xlsx")


### START ANALYSIS


#### PROTEINS and COMPLEXES

## TODO: align use of raasProfile vs. listProfile

## MEDIAN SITE AND PROTEIN RAAS

### median raas per unique mane protein site

sitl <- split(tmtf$RAAS, paste(tmtf$ensembl, tmtf$pos))
site <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## mod. column names and add protein and site info
colnames(site) <- paste0("RAAS.", colnames(site))
site$ensembl <- sub(" .*", "", rownames(site))
site$pos <- as.numeric(sub(".* ", "", rownames(site)))

## add gene names
site$name <- ens2nam [site$ensembl]

## add uniprot id
site$uniprot <- unlist(lapply(ens2u[site$ensembl], paste, collapse=";"))

## RAAS COLOR
site$RAAS.color <- num2col(site$RAAS.median,
                           limits=c(RAAS.MIN, RAAS.MAX), colf=arno)

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



if ( interactive() ) {
    ## test alternative measures of consistent RAAS
    ## CV:
    dense2d(site$RAAS.cv, -log10(site$RAAS.p.value),
            xlab="CV", ylab=expression(-log10(p)))
    dense2d(site$RAAS.cv, site$RAAS.sd, xlab="CV", ylab="SD")
    hist(site$RAAS.sd)
}

## ORDER SITES

## order proteins by RAAS
ptstat$rank <- rank(ptstat$median)

## make sure its ordered
## TODO: ORDER PROTEINS BY SIZE OR NUMBER OF AAS
site <-site[order(site$ensembl, site$pos),]

## order by protein RAAS rank (for hotspot plot)
site$rank <- ptstat[site$ensembl,"rank"]
site <-site[order(site$rank, site$pos),]

##
nsites <- as.numeric(sub(".*\\.","",tagDuplicates(site$ensembl)))
nsites[is.na(nsites)] <- 1
site$n <- nsites

### COMPARE dataset sizes

fmat <- matrix(NA, ncol=2, nrow=nrow(tmtf))
##fmat[,3] <- tmtf$Dataset
fmat[,2] <- tmtf$Dataset
fmat[,1] <- tmtf$Dataset
##fmat[duplicated(tmtf$unique),2] <- "na"
fmat[duplicated(paste(tmtf$BP,tmtf$SAAP)),1] <- "na"

cnts <- apply(fmat, 2, table)
cnts <- do.call(rbind,lapply(cnts, function(x) x[uds]))
rownames(cnts) <- c("BP/SAAP", "#RAAS")
if ( interactive() ) {
    barplot(cnts, beside=TRUE, legend=TRUE)
    
    hist(bdat$n, breaks=100)
}


### COLLECT PROTEIN DATA

#### TODO: why do we have n=747 at exon but n=713 for
#### the log intensity plot on exon/intron?

## PROTEIN ABUNDANCES via Razor.protein.precursor.intensity

## remove brackets
alst <- sub("\\[","", sub("\\]","", tmtf$Razor.protein.precursor.intensity))
## split multiple values and convert to numeric
alst <- lapply(alst, function(x) as.numeric(unlist(strsplit(x, ","))))
## remove zeros and take unique value
alst <- lapply(alst, function(x) unique(x[x!=0]))
## there are still many entries with multiple values
table(lengths(alst))
alst <- unlist(lapply(alst, mean))

## MEDIAN of PROTEINS
alst <- split(alst, tmtf$ensembl)
alst <- lapply(alst, function(x) x[!is.na(x)])

alst <- lapply(alst, as.numeric)

barplot(table(lengths(alst)))

## TODO: analyze stats before just taking the median
alst <- lapply(alst, median, na.rm=TRUE)
alst <- unlist(alst)

ptstat$intensity <- alst[rownames(ptstat)]


## RELATIVE POSITON

## PROTEIN LENGTH from saap_mapped4.tcv
plen <- split(hdat$len, hdat$ensembl)
plen <- lapply(plen, unique)
table(lengths(plen)) # check: all should be the same
plen <- unlist(lapply(plen, function(x) x[1]))

ptstat$length <- plen[rownames(ptstat)]

## THERMOSTABILITY @Savitski2014
therm <- as.data.frame(read_xlsx(thermo.file, sheet=2))
mlt <- therm$meltP_Jurkat
names(mlt) <- therm[,2]

ptstat$Tmelt <- mlt[match(pnms[rownames(ptstat)], names(mlt))]

## PROTEIN HALF-LIFE @Mathieson2018
hlvd <- readxl::read_xlsx(math18.file)

## mean half-live over all replicates and cell types
## TODO: consider halflife distributions instead of just taking mean
## NOTE: THIS includes mouse data, but they correlate.
cidx <- grep("half_life", colnames(hlvd), value=TRUE)
hlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
names(hlv) <- unlist(hlvd[,1])

ptstat$halflife <- hlv[match(pnms[rownames(ptstat)], names(hlv))]

## @Savitski2014 - difference of melting T w/o ATP
## TODO: analyze this better
therm <- as.data.frame(read_xlsx(thatp.file, sheet=2))
mlt <- as.numeric(therm$diff_meltP_Exp1)
names(mlt) <- therm[,2]

ptstat$DeltaTmelt <- mlt[match(pnms[rownames(ptstat)], names(mlt))]

## ProtStab2 Prediction
protstab <- read.csv(protstab.file)
## use refseq ID w/o version number as row names
rownames(protstab) <- sub("\\.[0-9]*","", protstab$id)

## get local reduced refseq mapping
## TODO: use full n<->n mapping and maximize yield
ps2ens <- refseq2ens[refseq2ens[,1] %in%rownames(ptstat),]
ps2ens <- ps2ens[ps2ens[,2] %in%rownames(protstab),]
protstab$ensembl <- ps2ens[match(rownames(protstab),ps2ens[,2]),1]


ptstat$ProtStab2 <- protstab$Human_predict_Tm[match(rownames(ptstat),
                                                    protstab$ensembl)]

## GET WHOLE PROTEIN MEAN IUPRED3 SCORE (from saap_peptides4.tsv
iu3 <- split(hdat$iupred3.protein, hdat$ensembl)
## QC: all protein level
table(unlist(lapply(iu3, function(x) length(unique))))
iu3 <- unlist(lapply(iu3, unique))

ptstat$iupred3 <- iu3[rownames(ptstat)]

## @Pepelnjak2024: in vitro 20S targets
## TODO: simplify mapping
p20 <- data.frame(readxl::read_xlsx(pepe24.file, sheet=2))

xl.20s <- expression(median~log[2]("20S"/control))
xl.20p <- expression("20S target:"~sign~"x"~log[10](p))

## MEDIAN LG2FC PER PROTEIN
p20s <- p20[!is.na(p20$Log2.FC.),]
p20l <- split(p20s[,c("Log2.FC.")],
              p20s[,"UniprotID"])
p20p <- unlist(lapply(p20l, median, na.rm=TRUE))

## median -log10(p)*sign
p20s$pscale <- sign(p20s$Log2.FC.) * -log10(p20s$p.value)
p20l <- split(p20s[,c("pscale")],
              p20s[,"UniprotID"])
p20pv <- unlist(lapply(p20l, median, na.rm=TRUE))


## get uniprot mapping
## TODO: akugb with uni2e and ens2u  in init script.
p2ens <- uni2ens
p2ens <- p2ens[p2ens[,2]%in%tmtf$ensembl,]
p2ens <- p2ens[p2ens[,1]%in%names(p20p),]

cat(paste(sum(duplicated(p2ens[,1])),"uniprot with multiple ensembl\n"))
cat(paste(sum(duplicated(p2ens[,2])),"uniprot with multiple ensembl\n"))
p2el <- unlist(split(p2ens[,2], p2ens[,1]))

e2u <- names(p2el)
names(e2u) <- p2el

if ( length(unique(lengths(p2el)))>1 )
    stop("non-unique uniprot 2 ensembl mapping")


lg2fc <- p20p
names(lg2fc) <- p2el[names(lg2fc)]
ptstat$p20_lgfc <- lg2fc[rownames(ptstat)]

pval <- p20pv
names(pval) <- p2el[names(pval)]
ptstat$p20_pval <- pval[rownames(ptstat)]


## write-out collected protein results
write.table(cbind(protein=rownames(ptstat), ptstat),
            file=file.path(pfig.path, "proteins_raas.tsv"), sep="\t",
            na="", row.names=FALSE, quote=FALSE)

### HOTSPOTS

## TODO:
## * RAAS profiles by sliding window

pns <- ptstat$n
pns[pns>15] <- ">15"
paas <- table(pns)[c(as.character(1:15),">15")]
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
pns[pns>15] <- ">15"
paas <- table(pns)[c(as.character(1:15),">15")]
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
pns[pns>15] <- ">15"
paas <- table(pns)[c(as.character(1:15),">15")]
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
pns[pns>15] <- ">15"
paas <- table(pns)[c(as.character(1:15),">15")]

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
bp <- barplot(paas, xlab="SAAP per BP",
              ylab=NA, axes=FALSE, axisnames = FALSE)
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
axis(2)
mtext("# of peptides", 2, 1.6)
legend("topright",
       paste(sum(bpstat$nsaap>1,nr.rm=TRUE), "peptides with >1 SAAP"), bty="n")
dev.off()

plotdev(file.path(pfig.path,paste0("hotspots_SAAP_per_peptide_log")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp <- barplot(paas, xlab="SAAP per BP",ylab=NA,
              axes=FALSE, axisnames = FALSE, log="y")
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
options(scipen=999)
axis(2, at=10^(1:10), las=2)
axis(2, at=c(1:10,(1:10)*10,(1:10)*100, (1:10)*1000), tcl=-.125,
     labels=FALSE)
options(scipen=0)
mtext("# of peptides",2,1.6)
legend("topright",
       paste(sum(bpstat$nsaap>1,nr.rm=TRUE), "peptides with >1 SAAP"), bty="n")
dev.off()


pns <- bpstat$n
pns[pns>15] <- ">15"
paas <- table(pns)[c(as.character(1:15),">15")]

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

### POSITION ON PROTEIN

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
    axis(1, at=log10(rep(1:10, 7) * 10^rep(-4:2, each=10)), tcl=-.125, labels=FALSE)
    dev.off()
}

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

filters <- list(all=1:nrow(bdat),
                high=which(bdat$RAAS> -1))
fnms <- c(all="all", high="RAAS> 0.1")

for ( i in length(filters):1 ) {

    nm <- names(filters)[i]
    rdat <- bdat[filters[[i]],]
    
    dsts <- matrix(NA, nrow=nrow(rdat), ncol=nrow(rdat))
    same <- dsts ## same basepeptide filter
    for ( i in 1:nrow(rdat) ) {

        dst <- abs(rdat$pos[i] - rdat$pos)
        dst[rdat$ensembl!=rdat$ensembl[i]] <- Inf
        dsts[i,] <- dst
        same[i,] <- rdat$BP[i] == rdat$BP
        
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
    plot(hst$mids+.5, hst$counts, type="l", xlim=c(0,bpmx*2),
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
    plot(hst$mids+.5, hst$counts, type="l", xlim=c(0,bpmx),
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
    abline(v=10)
    dev.off()
    

}


### AUTOCORRELATIONS 

## TODO: do we need to clean data for this?

## * concatenated proteins
## * concatenated BP.

N <- 150

## autocorrelation of of site median RAAS with concatenated protein sequences

## * initiate vector for concatenated proteins with median RAAS,
## * set actual RAAS values,
## * calculate auto-correlation.

osite <-site[order(site$ensembl, site$pos),]
ositel <- split(osite, osite$ensembl)
##ositel <- ositel[unlist(lapply(ositel, nrow))>5]


## RAAS

## expand each protein to full sequence vector
## initialized to NA (not: global median RAAS)
##globalraas <- median(tmtf$RAAS, na.rm=TRUE)
araas <- lapply(ositel, function(x) {
    aras <- rep(NA, plen[unique(x$ensembl)])
    aras[x$pos] <- x$RAAS.median # TODO: max(abs(median-RAAS))
    aras
})
araas <- unlist(araas)

## calculate autocorrelation
acor <- rep(NA, N)
for ( k in 1:N ) {
    rng <- 1:(length(araas)-N)
    acor[k] <- cor(araas[rng], araas[rng+k], use="pairwise.complete")
}



if ( interactive() ) {

    ## TODO: different handling of NA in R's acf
    bcor <- acf(araas, lag.max=N, na.action = na.pass)
    ##abline(v=10, col=2)
    
    plot(bcor$acf[1:N+1], type="l")
    lines(acor, col=2)
    plotCor(acor, bcor$acf[1:N+1])
    
    bfft <- fft(bcor$acf[1:N+1])
    dc <- bfft[1]/N # TODO: compare with mean cor
    plot(abs(bfft)/N, type="h")

    afft <- fft(acor)
    dc <- afft[1]/N # TODO: compare with mean cor
    plot(abs(afft)/N, type="h")
}


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

## TODO: fourier with correct period
if ( interactive() ) {
    nfft <- fft(ncor)[1:(N/2)]
    dc <- nfft[1]/N # TODO: compare with mean cor
    tosc <- N/1:(N/2-1)
    plot(tosc, abs(nfft[2:length(nfft)])/N, type="h", xlim=c(2,6))
}

### NOTE: R's acf is equivalent to manual pearson correlation
### but behaves differently with NAs
if ( interactive() ) {

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



## high RAAS measurement frequency
rsite <- osite[osite$RAAS.median> .1,]
rsitel <- split(rsite, rsite$ensembl)
rraas <- lapply(rsitel, function(x) {
    aras <- rep(0, plen[unique(x$ensembl)])
    aras[x$pos] <- x$RAAS.n
    aras
})
rraas <- unlist(rraas)

## calculate autocorrelation
rcor <- rep(NA, N)
for ( k in 1:N ) {
    rng <- 1:(length(rraas)-N)
    rcor[k] <- cor(rraas[rng], rraas[rng+k])
}

## plot auto-correlations

plotdev(file.path(pfig.path,paste0("autocorrelation_protein_RAAS")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, acor, type="l", xlab="lag k", ylab=expression(C[RAAS]))
points(1:N, acor, cex=.25)
abline(v=10)
dev.off()

plotdev(file.path(pfig.path,paste0("autocorrelation_protein_frequency")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, ncor, type="l", xlab="lag k", ylab=expression(C[frequency]))
points(1:N, ncor, cex=.25)
dev.off()

plotdev(file.path(pfig.path,
                  paste0("autocorrelation_protein_frequency_highRAAS")),
        type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, rcor, type="l", xlab="lag k", ylab=expression(C[frequency]))
points(1:N, rcor, cex=.25)
figlabel(expression(RAAS > 0.1), pos="topright", region="plot")
dev.off()

## AUTO-CORRELATION FOR CONCATENATED BASE PEPTIDES

## split by BP + site in BP
bpl <- split(tmtf$RAAS, paste(tmtf$BP, tmtf$site))
bpe <- listProfile(bpl, y=tmtf$RAAS, use.test=use.test, min=3)

## mod. column names and add protein and site info
colnames(bpe) <- paste0("RAAS.", colnames(bpe))
bpe$BP <- sub(" .*", "", rownames(bpe))
bpe$site <- as.numeric(sub(".* ", "", rownames(bpe)))
bpe$len <- nchar(bpe$BP)

## ORDER by BP and position and split by BP
obpe <- bpe[order(bpe$BP, bpe$site),]
obpl <- split(obpe, obpe$BP)
##ositel <- ositel[unlist(lapply(ositel, nrow))>5]


## RAAS

## expand each protein to full sequence vector
## initialized to NA 
araas <- lapply(obpl, function(x) {
    aras <- rep(NA, unique(x$len))
    aras[x$site] <- x$RAAS.median
    aras
})
araas <- unlist(araas)

## calculate autocorrelation
acor <- rep(NA, N)
for ( k in 1:N ) {
    rng <- 1:(length(araas)-N)
    acor[k] <- cor(araas[rng], araas[rng+k], use="pairwise.complete")
}

## MEASUREMENT FREQUENCY

## expand each protein to full sequence vector
## initialized to 0 and add RAAS counts
nraas <- lapply(obpl, function(x) {
    aras <- rep(0, unique(x$len))
    aras[x$site] <- x$RAAS.n
    aras
})
nraas <- unlist(nraas)

## calculate autocorrelation
ncor <- rep(NA, N)
for ( k in 1:N ) {
    rng <- 1:(length(nraas)-N)
    ncor[k] <- cor(nraas[rng], nraas[rng+k])
}

plotdev(file.path(pfig.path,paste0("autocorrelation_BP_RAAS")),
            type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, acor, type="l", xlab="lag k", ylab=expression(C[RAAS]))
points(1:N, acor, cex=.25)
abline(v=10)
dev.off()

plotdev(file.path(pfig.path,paste0("autocorrelation_BP_frequency")),
            type=ftyp, res=300, width=4, height=2)
par(mai=c(.5,.5,.1,.1), mgp=pmpg, tcl=-.25)
plot(1:N, ncor, type="l", xlab="lag k", ylab=expression(C[frequency]))
points(1:N, ncor, cex=.25)
abline(v=10)
dev.off()

## TODO: fourier of auto-correlation function


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
wtxt[wtxt>15] <- ">15"
wtxt <- table(wtxt)[c(as.character(1:15),">15")]

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
wtxt[wtxt>15] <- ">15"
wtxt <- table(wtxt)[c(as.character(1:15),">15")]

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

