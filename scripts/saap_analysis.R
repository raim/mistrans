
## TODO
## 2012: Jonathan Weissman, Ignolia, Cell paper: MOUSE embryonic stem cells,
## ribosome density at mutated position.
##  https://doi.org/10.1016/j.cell.2011.10.002
## Bacterial mistranslation v ribosome density.

## ribosome density: ribosome density per codon divided
## by the mean ribosome density of transcript.


library(viridis)
library(segmenTools)
library(Biostrings) # for blosum62
data(BLOSUM62)
data(PAM250)
options(stringsAsFactors=FALSE)

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"

gen.path <- file.path(mam.path, "originalData")
dat.path <- file.path(proj.path,"originalData")
fig.path <- file.path(proj.path,"figures","saap_analysis")
out.path <- file.path(proj.path,"processedData")

feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

## SAAP mapped to protein and codons
in.file <- file.path(out.path,"saap_mapped.tsv")
dat <- read.delim(in.file)

### HOTSPOTS
## get ALL proteins - from project mammary
genes <- read.delim(feature.file)
##genes <- genes[genes$type=="protein_coding",]

## look for proteins with more than one mutation
## in close vicinity
hotspots <- split(dat$pos, f=dat$protein)
hotspots <- lapply(hotspots, sort)

mutn <- unlist(lapply(hotspots, length))

idx <-  which.max(mutn)
#idx <- which(names(hotspots)=="ENSP00000506126")
id <- sub("_.*","",names(hotspots)[idx])
name <- genes$name[grep(id, genes$protein)]
if ( length(name)==0 ) name <- id
hpos <- table(hotspots[[idx]])

png(file.path(fig.path,"hotspots_example.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(x=as.numeric(names(hpos)), y=hpos, type="h",
     main=name, xlab="position in protein",
     ylab="# SAAP")
axis(2)
dev.off()

#### PLOTS

raas.col <- "Mean_precursor_RAAS"
hdat <- dat#[dat[,raas.col]> -1,]

###  CODONS
aa <- unique(GENETIC_CODE)
CODONS <- rep("", length(aa))
for ( i in seq_along(aa) )
    CODONS[i] <- paste(names(which(GENETIC_CODE==aa[i])), collapse=";")
names(CODONS) <- aa

## TODO:
## * store codon context,-1,+1,
## * check consistency of dat$from and dat$codon,
##   -> requires to account for genomic mutations?
## * load codon usage table, AA usage table,


## codon positions
cpos <- strsplit(hdat$codon,"")
## https://pubmed.ncbi.nlm.nih.gov/11164038/ A vs. U in 2nd

c1 <- unlist(lapply(cpos, function(x) x[1]))
c2 <- unlist(lapply(cpos, function(x) x[2]))
c3 <- unlist(lapply(cpos, function(x) x[3]))


boxplot(hdat[,raas.col] ~ c1, ylab=raas.col)
boxplot(hdat[,raas.col] ~ c2, ylab=raas.col)
boxplot(hdat[,raas.col] ~ c3, ylab=raas.col)

cfrq <- rbind("1"=table(c1),
              "2"=table(c2),
              "3"=table(c3))
cfrq <- cfrq[,c("A","T","G","C")]
png(file.path(fig.path,"codons_type.png"),
    res=300, width=3, height=3, units="in")
par(mai=c(.5,.25,.2,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(cfrq,beside=TRUE, legend.text=rownames(cfrq),
        args.legend=list(x="top",ncol=3, inset=c(-.02,-.1), bty="n",
                         title="codon position"))
dev.off()
png(file.path(fig.path,"codons_pos.png"),
    res=300, width=3, height=3, units="in")
par(mai=c(.5,.25,.2,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(t(cfrq),beside=TRUE, legend.text=colnames(cfrq),
        args.legend=list(x="top",ncol=4, inset=c(-.02,-.1), bty="n"),
        xlab="codon position")
dev.off()



fta <- table(hdat$to, hdat$from)

image_matrix(fta, axis=1:2, ylab="mistranslated", xlab="encoded")

png(file.path(fig.path,"codons_all.png"),
    res=300, width=15, height=5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
par(mai=c(.7,.5,.7,3))
ftc <- table(hdat$to,paste(hdat$from,hdat$codon, sep="-"))
txt <- ftc
txt[txt==0] <- ""
ftc[ftc==0] <- NA
image_matrix(ftc, axis=1:3, text=txt, ylab="mistranslated",
             text.cex=.7, xlab="")
axis(4, nrow(ftc):1, labels=CODONS[rownames(ftc)], las=2)
dev.off()


### POSITION v LENGTH

## size cut off
size.cutoff <- c(2.5e2,1e3)
small <- hdat$len <  size.cutoff[1]
large <- hdat$len >= size.cutoff[1] & hdat$len <= size.cutoff[2]
huge  <- hdat$len >  size.cutoff[2]



hist(hdat$pos, breaks=seq(0,4e4,100), xlim=c(0,3000), xlab="absolute position",
     main=NA)
hist(hdat$len, breaks=seq(0,4e4,100), border=2, add=TRUE)
legend("topright", c("BP position","protein length"), col=c(1,2), lty=1)


png(file.path(fig.path,"position_length_total.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(log10(hdat$len),
     xlab=expression(protein~length/log[10]), breaks=50)
abline(v=log10(size.cutoff))
dev.off()


png(file.path(fig.path,"position_length_relative.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(hdat$rpos, log10(hdat$len), colf=viridis::viridis,
        ##len, log="y",
        xlab="relative position of AAS in protein",
        ylab=expression(protein~length/log[10]))
abline(h=log10(size.cutoff))
dev.off()

png(file.path(fig.path,"position_length_absolute.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i",
    xpd=TRUE)
dense2d(hdat$len, hdat$pos, colf=viridis::viridis,
        xlim=c(0,1500), ylim=c(0,1500),
        ylab=expression(SAAP~position),
        xlab=expression(protein~length), nbin=512)
par(xpd=FALSE)
##abline(a=0, b=1)
abline(a=0, b=.5)
abline(v=size.cutoff)
dev.off()

png(file.path(fig.path,"position_hist_absolute.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i",
    xpd=TRUE)
hist(hdat$pos)
dev.off()


png(file.path(fig.path,"position_ecdf_relative_short.png"),
    res=300, width=2.5, height=2.5, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
cdf <- ecdf(hdat$rpos[small]) 
plot(seq(0,1,.01), cdf(seq(0,1,.01)), xlim=c(0,1), ylim=c(0,1),
     xlab="relative position in protein",
     ylab=expression(ecdf(x)), col=2, type="l", main=NA)
abline(a=0,b=1, col=1)
dev.off()

png(file.path(fig.path,"position_ecdf_relative_mid.png"),
    res=300, width=2.5, height=2.5, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
cdf <- ecdf(hdat$rpos[large])
plot(seq(0,1,.01), cdf(seq(0,1,.01)), xlim=c(0,1), ylim=c(0,1),
     xlab="relative position in protein",
     ylab=expression(ecdf(x)), col=2, type="l", main=NA)
abline(a=0,b=1, col=1)
dev.off()

png(file.path(fig.path,"position_ecdf_relative_long.png"),
    res=300, width=2.5, height=2.5, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
cdf <- ecdf(hdat$rpos[huge])
plot(seq(0,1,.01), cdf(seq(0,1,.01)), xlim=c(0,1), ylim=c(0,1),
     xlab="relative position in protein",
     ylab=expression(ecdf(x)), col=2, type="l", main=NA)
abline(a=0,b=1, col=1)
dev.off()

png(file.path(fig.path,"position_ecdf_absolute.png"),
    res=300, width=2.5, height=2.5, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
cdf <- ecdf(hdat$pos)
x <- seq(0,max(hdat$pos),1)
plot(x, cdf(x), xlim=c(0,1000), xlab="absolute position in protein",
     ylab=expression(ecdf(x)), type="l", col=2)
abline(v=8.5, col=1)
dev.off()

png(file.path(fig.path,"position_ecdf_absolute_zoom.png"),
    res=300, width=2.5, height=2.5, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
cdf <- ecdf(hdat$pos)
plot(seq(0,100,.01), cdf(seq(0,100,.01)),xlim=c(0,30), ylim=c(0,.05),
     xlab="absolute position in protein",
     ylab=expression(ecdf(x)), type="l", col=2)
abline(v=10, col=1)
dev.off()

png(file.path(fig.path,"position_length_absolute_log.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
dense2d(log10(hdat$len), log10(hdat$pos),
        colf=viridis::viridis, ylim=c(1,4), xlim=c(1,4),
        ylab=expression(SAAP~position/log[10]),
        xlab=expression(protein~length/log[10]))
abline(a=0,b=1)
dev.off()

png(file.path(fig.path,"position_hist_relative.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(hdat$rpos, breaks=seq(0,1,0.1),
     main=paste(sum(!is.na(hdat$rpos)), "unique SAAP"),
     xlab="relative position of AAS in protein")
loc.sze <- table(hdat$rpos>.5)
dev.off()

## BACKGROUND: main peptides
load(file.path(out.path, "mapped_peptides.rda"))
png(file.path(fig.path,"position_hist_relative_main_peptides.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
hist(unlist(mpos), breaks=seq(0,1,0.1),
     main=paste(sum(!is.na(unlist(mpos))), "main peptides"),
     xlab="relative position of AAS in protein")
par(new=TRUE)
cdf <- ecdf(unlist(mpos))
plot(seq(0,1,.01), cdf(seq(0,1,.01)), type="l", axes=FALSE, col=2,
     xlab=NA, ylab=NA)
abline(a=0,b=1)
mtext(expression(ecdf(x)), 4, 1.3)
axis(4)
dev.off()

## SMALL v LARGE PROTEINS
png(file.path(fig.path,"position_hist_relative_short.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
hist(hdat$rpos[small], breaks=seq(0,1,0.1),
     main=paste0(sum(small,na.rm=TRUE),
                 " short proteins <",size.cutoff[1]," aa"),
     xlab="relative position of AAS in protein")
par(new=TRUE)
cdf <- ecdf(hdat$rpos[small])
plot(seq(0,1,.01), cdf(seq(0,1,.01)), type="l", axes=FALSE, col=2,
     xlab=NA, ylab=NA)
abline(a=0,b=1)
mtext(expression(ecdf(x)), 4, 1.3)
axis(4)
dev.off()
png(file.path(fig.path,"position_hist_relative_long.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
hist(hdat$rpos[huge], breaks=seq(0,1,0.1),
     main=paste0(sum(huge,na.rm=TRUE),
                 " long proteins >",size.cutoff[2]," aa"),
     xlab="relative position of AAS in protein")
par(new=TRUE)
cdf <- ecdf(hdat$rpos[huge])
plot(seq(0,1,.01), cdf(seq(0,1,.01)), type="l", axes=FALSE, col=2,
     xlab=NA, ylab=NA)
abline(a=0,b=1)
mtext(expression(ecdf(x)), 4, 1.3)
axis(4)
dev.off()
png(file.path(fig.path,"position_hist_relative_mid.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
hist(hdat$rpos[large], breaks=seq(0,1,0.1),
     main=paste0(sum(large,na.rm=TRUE),
                 " proteins >",size.cutoff[1]," aa, <",size.cutoff[2]," aa"),
     xlab="relative position of AAS in protein")
par(new=TRUE)
cdf <- ecdf(hdat$rpos[large])
plot(seq(0,1,.01), cdf(seq(0,1,.01)), type="l", axes=FALSE, col=2,
     xlab=NA, ylab=NA)
abline(a=0,b=1)
mtext(expression(ecdf(x)), 4, 1.3)
axis(4)
dev.off()


raas.col <- "Mean_precursor_RAAS"
png(file.path(fig.path,"position_RAAS_absolute.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(hdat$pos, hdat[,raas.col], colf=viridis::viridis,#xlim=c(0,1e4),
        ylim=range(hdat[,raas.col], na.rm=TRUE),
        xlab="absolute position of AAS in protein",
        ylab=raas.col)
dev.off()

png(file.path(fig.path,"position_RAAS_relative.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(hdat$rpos, hdat[,raas.col], colf=viridis::viridis,
        xlim=range(hdat$rpos, na.rm=TRUE),
        ylim=range(hdat[,raas.col], na.rm=TRUE),
        xlab="relative position of AAS in protein",
        ylab=raas.col)
dev.off()

png(file.path(fig.path,"protein_length_RAAS.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hdat[,raas.col], hdat$len, colf=viridis::viridis, na.rm=TRUE,
        ylim=c(0,1e4),
        xlim=range(hdat[,raas.col], na.rm=TRUE),
        ylab="protein length",
        xlab=raas.col)
dev.off()

png(file.path(fig.path,"protein_length_RAAS_log.png"),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hdat[,raas.col], log10(hdat$len), colf=viridis::viridis,na.rm=TRUE,
        xlim=range(hdat[,raas.col], na.rm=TRUE),
        ylab=expression(log[10](protein~length)),
        xlab=raas.col)
dev.off()




#### REPLACED AA SIMILARITIES

## first find mutated AA pairs
## (column AASin input mixes I/L)
## and split mutated AA pairs into from/to

## TODO: get this from the columns in the mapped file!
saaps <- strsplit(dat$SAAP,"")
bases <- strsplit(dat$BP, "")
fromto <- lapply(1:length(saaps), function(i) {
    pos <- which(saaps[[i]]!=bases[[i]])
    c(from=bases[[i]][pos], to=saaps[[i]][pos])
})


## analyze AA similarity by different matrices
## TODO: remove extra columns
simmats <- c("BLOSUM62", "PAM250")

## reporter vs. precursor RAAS

png(file.path(fig.path,paste0("raas_precursor_reporter.png")),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(dat$Mean_reporter_RAAS, dat$Mean_precursor_RAAS, na.rm=TRUE,
        xlab="mean reporter RAAS", ylab="mean precursor RAAS",
        colf=viridis::viridis)
abline(a=0, b=1, col=1, lty=2, lwd=2)
dev.off()

for ( i in seq_along(simmats) ) {

    mid <- simmats[i]
    MAT <- get(mid)[AA_STANDARD,AA_STANDARD]

    ylab <- "similarity between replaced AA" #"BLOSUM62 similarity"

    ## similarity of replaced AA
    sim <- unlist(lapply(fromto, function(x) MAT[x[1],x[2]]))

    rid <- "Mean_precursor_RAAS"
    raas <- unlist(dat[,rid])
    
    brks <- c(-5,-3,-2,-1,0,1,4)
    raas_bins <- cut(raas, brks)
    raas_dual <- rep("RAAS>0", length(raas))
    raas_dual[raas < 0] <- "RAAS<0"
    
    df <- data.frame(raas=raas,
                     bins=raas_bins,
                     dual=raas_dual,
                     sim=sim)
    df <- df[!is.na(raas),]# & !is.infinite(raas),]
    dff <- df[!is.infinite(df$raas),]
    
    
    
    png(file.path(fig.path,paste0("AAS_",mid,"_raas_dense.png")),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(dff$sim, dff$raas, ylab=ylab, xlab=rid, colf=viridis::viridis)
    dev.off()

    png(file.path(fig.path,paste0("AAS_",mid,"_raas_boxplot.png")),
        res=300, width=5, height=3.5, units="in")
    par(mai=c(.75,.5,.5,.1), mgp=c(1.3,.3,0), tcl=-.25)
    boxplot(sim ~ bins, data=df, ylab=ylab,
            xlab=NA, las=2, at=seq_along(levels(df$bins)))
    axis(3, at=seq_along(levels(df$bins)), labels=table(df$bins),las=2)
    mtext(rid, 1, 2.5)
    dev.off()
    

    png(file.path(fig.path,paste0("AAS_",mid,"_raas_boxplot_dual.png")),
        res=300, width=3.75, height=3.75, units="in")
    b62 <- cbind.data.frame(MAT[lower.tri(MAT)],
                            mid)
    ddff <- data.frame(sim=c(dff[,"sim"], b62[,1]),
                       dual=c(dff[,"dual"],b62[,2]))
    ddff[,2] <- as.factor(ddff[,2])
    ylim <- range(ddff[,1])
    ylim[2]  <- ylim[2]*1.1
    par(mai=c(.5,.5,.25,0), mgp=c(1.3,.3,0), tcl=-.25)
    boxplot(sim ~ dual, data=ddff, xlab=NA,
            ylab=ylab, axes=FALSE, at=1:3, box=FALSE, ylim=ylim)
    axis(2)
    axis(1, at=1:3, labels=paste0(levels(ddff[,2]), "\n", table(ddff[,2])),
         mgp=c(0,1.5,0))
    ## TODO: significance bars
    tt <- wilcox.test(ddff$sim[ddff$dual=="RAAS<0"],
                      ddff$sim[ddff$dual=="RAAS>0"])
    text(2.5, par("usr")[4]*.9, paste("p =",signif(tt$p.value,1)),
         xpd=TRUE, pos=3)
    arrows(x0=2, y0=par("usr")[4]*.9, x1=3, angle=90,
           code=3, length=.05, xpd=TRUE)
    tt <- wilcox.test(ddff$sim[ddff$dual==mid],
                      ddff$sim[ddff$dual=="RAAS<0"])
    text(1.5, par("usr")[4], paste("p =",signif(tt$p.value,1)), xpd=TRUE, pos=3)
    arrows(x0=1, y0=par("usr")[4], x1=2, angle=90, code=3, length=.05, xpd=TRUE)
    dev.off()
    
    png(file.path(fig.path,paste0("AAS_",mid,"_similarities_hist.png")),
        res=300, width=5, height=3.5, units="in")
    hist(MAT[lower.tri(MAT)])
    dev.off()
}
