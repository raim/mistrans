
library(readxl)
library(segmenTools)
options(stringsAsFactors=FALSE)
options(scipen = 0)

data.path <- "~/data/"
mam.path <- file.path(data.path, "mammary")
yeast.path <- file.path(data.path, "yeast")

fig.path <- file.path(mam.path,"figures","halflives")
dir.create(fig.path)
ftyp <- "png"
corW <- corH <- 2.5
pmai <- c(.4,.4,.2,.2)
pmpg <- c(1,.2,0)

## yeast protein halflives - @Christiano2014
yfeature.file <- file.path(yeast.path,"feature_R64-1-1_20110208_allData.csv")
yhlf.file <- file.path(yeast.path,"processedData","christiano14.csv")
## human protein halflives - @Mathieson2018
hhlf.file <- file.path(mam.path,"originalData",
                       "41467_2018_3106_MOESM5_ESM.xlsx")
hlen.file <- file.path(mam.path,"processedData",
                       "protein_length.tsv")
## orthologous genes via feature file
hfeature.file <- file.path(mam.path,"features_GRCh38.110.tsv")


## @Pepelnjak2024: 20S targets
pepe24.file <- file.path(mam.path,"originalData",
                         "44320_2024_15_moesm1_esm.xlsx")

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

## uniprot/refseq/ensembl/name mappings
uni2ens.file <- file.path(mam.path,"originalData","uniprot_ensembl.dat")
uni2nam.file <- file.path(mam.path,"originalData","uniprot_name.dat")
refseq.file <- file.path(mam.path,"originalData",
                         "ensembl_refseq_20240528.tsv.gz")

## human genes - restrict to protein coding and MANE collection
hgenes <- read.delim(hfeature.file)
filter <- hgenes$type=="protein_coding" & hgenes$MANE!=""
hgenes <- hgenes[filter,]

## map human:yeast genes ortholog column!
ygns <- strsplit(hgenes$yeast, ";")
names(ygns) <- hgenes$ID

olen <- unlist(lapply(ygns, length))
table(olen)

## map yeast cohorts to human genes
## TAKE ONLY FIRST OF YEAST ORTHOLOGS
ymap <- unlist(lapply(strsplit(hgenes$yeast, ";"), function(x) x[1]))
names(ymap) <- hgenes$name
y2h <- hgenes$name[!is.na(ymap)]
names(y2h) <- ymap[!is.na(ymap)]

hhld <-  readxl::read_xlsx(hhlf.file)
## mean half-live over all replicates and cell types
## TODO: consider distribution
hhld <- hhld[,-grep("Mouse", colnames(hhld))]
cidx <- grep("half_life", colnames(hhld), value=TRUE)

hhlf <- apply(hhld[,cidx], 1, mean, na.rm=TRUE)
names(hhlf) <- unlist(hhld[,1])

##hhlf[hhlf>1000] <- 1000
## yeast half-lives
yhld <- read.delim(yhlf.file)
yhld$hgene <-
    y2h[yhld[,1]]
yhld <- yhld[!is.na(yhld$hgene),]
yhld <- yhld[!duplicated(yhld$hgene),]
rownames(yhld) <- yhld$hgene

## UNIPROT <-> ENSEMBL MAPPING
## NOTE: ~90 duplicated ensembl IDs
uni2ens <- read.delim(uni2ens.file, header=FALSE)
uni2ens[,2] <- sub("\\..*", "", uni2ens[,2]) # remove ensembl version tag
uni2ens[,1] <- sub("-.*", "", uni2ens[,1]) # remove uniprot version tag
## 1-many lists
uni2e <- split(uni2ens[,2], uni2ens[,1])
ens2u <- split(uni2ens[,1], uni2ens[,2])
## REFSEQ <-> ENSEMBL MAPPING
refseq2ens <- read.delim(refseq.file, header=TRUE)


### AXIS LABELS
xl.hlf <- expression(protein~"half-life"/h)
xl.hhlf <- expression(human~protein~"half-life"/h)
xl.yhlf <- expression(yeast~protein~"half-life"/h)
xl.hdeg <- expression(human~degradation/h^-1)
xl.ydeg <- expression(yeast~degradation/h^-1)

xl.hlgd <- expression(human~degradation/log(d^-1))
xl.ylgd <- expression(yeast~degradation/log(d^-1))

## get half-lives
yhlf <- yhld[names(hhlf),"t1.2..hours."]
yhlf[yhlf==">= 100"] <- 120
yhlf <- as.numeric(yhlf)

# deg rate in days
hdeg <- log(2)/(hhlf)
ydeg <- log(2)/(yhlf)

hlgd <- log(log(2)/(hhlf/24))
ylgd <- log(log(2)/(yhlf/24))


plotdev(file.path(fig.path,"protein_halflives_human_yeast_zoom"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(hhlf, yhlf, xlim=c(0,200), ylim=c(0,20),
        xlab=xl.hhlf,
        ylab=xl.yhlf, title=TRUE, cor.legend=FALSE)#, legpos="topright")
##abline(v=48)
figlabel("zoom", pos="bottomleft", font=2, cex=1.2)
dev.off()

plotdev(file.path(fig.path,"protein_halflives_human_yeast"),
        type=ftyp, height=3, width=4, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(hhlf, yhlf, 
        xlab=xl.hhlf,
        ylab=xl.yhlf, title=TRUE, cor.legend=FALSE)#, legpos="topright")
##abline(v=48)
dev.off()

plotdev(file.path(fig.path,"protein_degradation_human_yeast"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(hdeg, ydeg,
        xlab=xl.hdeg,
        ylab=xl.ydeg, title=TRUE, cor.legend=FALSE)#, legpos="topright")
showl <- (hlgd < -4 & ylgd < -1)
showh <- (ylgd > 1 & hlgd > 0) | ylgd> 3
showb <- (hlgd < -4) | (hlgd > -1 & ylgd < 0)
showby <- ylgd < -1 & hlgd < -1
show <- TRUE

dev.off()

plotdev(file.path(fig.path,"protein_degradation_human_yeast_zoom"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(hdeg, ydeg,
        xlab=xl.hdeg,
        ylab=xl.ydeg, ylim=c(0,.3), title=TRUE, cor.legend=FALSE)#, legpos="topright")
showl <- (hlgd < -4 & ylgd < -1)
showh <- (ylgd > 1 & hlgd > 0) | ylgd> 3
showb <- (hlgd < -4) | (hlgd > -1 & ylgd < 0)
showby <- ylgd < -1 & hlgd < -1
show <- TRUE

dev.off()

## log as in multi-species paper 
plotdev(file.path(fig.path,"protein_degradation_human_yeast_log"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(hlgd, ylgd,
        xlab=xl.hlgd,
        ylab=xl.ylgd, title=TRUE, cor.legend=FALSE)#, legpos="topright")
showl <- (hlgd < -4 & ylgd < -1)
showh <- (ylgd > 1 & hlgd > 0) | ylgd> 3
showb <- (hlgd < -4) | (hlgd > -1 & ylgd < 0)
showby <- ylgd < -1 & hlgd < -1
show <- TRUE
if ( show ) {
    shadowtext(hlgd[showl], ylgd[showl], labels=names(hlgd)[showl],
               xpd=TRUE, pos=4, cex=.7, col=1)
    shadowtext(hlgd[showh], ylgd[showh], labels=names(hlgd)[showh],
               xpd=TRUE, pos=4, cex=.7, col=1)
    shadowtext(hlgd[showb], ylgd[showb], labels=names(hlgd)[showb],
               xpd=TRUE, pos=4, cex=.7, col=1)
    shadowtext(hlgd[showby], ylgd[showby], labels=names(hlgd)[showby],
               xpd=TRUE, pos=3, cex=.7, col=1, srt=90)
}
##abline(v=48)
dev.off()

if ( interactive() ) {
    hist(hhlf,breaks=100)
    hist(log(2)/hhlf, breaks=100)
    hist(log(log(2)/(hhlf/24)))
}

### HALF-LIFE and LENGTH

## YEAST
ygenes <- read.delim(yfeature.file)
ygenes <- ygenes[ygenes$type=="gene",]
yhld <- read.delim(yhlf.file)

yhlf <- yhld[match(ygenes$ID, yhld$ENSG), "t1.2..hours."]
yhlf[which(yhlf==">= 100")] <- 120
yhlf <- as.numeric(yhlf)



## HUMAN
hgenes <- read.delim(hfeature.file)
filter <- hgenes$type=="protein_coding" & hgenes$MANE!=""
hgenes <- hgenes[filter,]

## thermostability data
therm <- as.data.frame(readxl::read_xlsx(thermo.file, sheet=2))
mlt <- therm$meltP_Jurkat
names(mlt) <- therm[,2]

## half life data
hhld <-  as.data.frame(readxl::read_xlsx(hhlf.file))
## mean half-live over all replicates and cell types, incl Mouse neurons
## TODO: consider distribution
##hhld <- hhld[,-grep("Mouse", colnames(hhld))]
##hhld <- hhld[,grep("Hepatocytes", colnames(hhld))]
cidx <- grep("half_life", colnames(hhld), value=TRUE)
hhlf <- apply(hhld[,cidx], 1, mean, na.rm=TRUE)
names(hhlf) <-   hhld[,1]

## protein lengths
hlen <- read.delim(hlen.file, header=FALSE, row.names=1)
rownames(hlen) <- sub("\\.[0-9]*", "", rownames(hlen))

## TODO: predicted stability
protstab <- read.csv(protstab.file)
## use refseq ID w/o version number as row names
rownames(protstab) <- sub("\\.[0-9]*","", protstab$id)

## get local reduced refseq mapping
## TODO: use full n<->n mapping and maximize yield
ps2ens <- refseq2ens[refseq2ens[,1] %in%  hgenes$MANE.protein,]
ps2ens <- ps2ens[ps2ens[,2] %in%rownames(protstab),]
protstab$ensembl <- ps2ens[match(rownames(protstab),ps2ens[,2]),1]
## match to RAAS table
idx <- match(hgenes$MANE.protein, protstab$ensembl)
pstab <- protstab$Human_predict_Tm[idx]
names(pstab) <- hgenes$MANE.protein

## expand proteins and gene names
hmap <- cbind(hgenes$name, hgenes$MANE.protein)
hmap <- cbind.data.frame(hmap,
                         len=hlen[hmap[,2],1],
                         hlf=hhlf[hmap[,1]],
                         mlt=mlt[hmap[,1]],
                         mlt_predicted=pstab[hmap[,2]])

## remove duplicated gene names AFTER SORTING FOR LENGTH
## NOTE: correlation much stronger when taking the longest proteoform
## for each gene. MANE PROTEIN INSTEAD
#hmap <- hmap[order(hmap$len, decreasing=TRUE),]
#hmap <- hmap[!duplicated(hmap[,1]),]

### TODO: add MANE.protein column to human feature file
## or generate protein file

## plot half-lives vs lengths

plotdev(file.path(fig.path,"protein_halflive_length_yeast"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ygenes$p.length), log10(yhlf), round=3,
        ylim=log10(c(.1,2e3)),  xlim=log10(c(50,3e4)),
        xlab="protein length", ylab="mean protein half-life/h", axes=FALSE,
        title=TRUE, cor.legend=FALSE)
for ( ax in 1:2 ) {
    axis(ax, at=0:9, labels=10^(0:9))
    axis(ax, at=log10(rep(1:10, 6) * 10^rep(-1:4, each=10)),
         tcl=-.125, labels=FALSE)
}
figlabel("yeast", pos="bottomleft", cex=1.2, font=2)
dev.off()

plotdev(file.path(fig.path,"protein_halflive_length_human"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(hmap$len), log10(hmap$hlf), round=3,
        ylim=log10(c(1,3e4)), xlim=log10(c(50,3e4)),
        xlab="protein length", ylab="mean protein half-life/h", axes=FALSE,
        title=TRUE, cor.legend=FALSE)
for ( ax in 1:2 ) {
    axis(ax, at=0:9, labels=10^(0:9))
    axis(ax, at=log10(rep(1:10, 6) * 10^rep(-1:4, each=10)),
         tcl=-.125, labels=FALSE)
}
figlabel("human", pos="bottomleft", cex=1.2, font=2)
dev.off()

plotdev(file.path(fig.path,"protein_halflive_melting_human"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(hmap$mlt, log10(hmap$hlf), round=3,
        ##ylim=log10(c(1,3e4)), xlim=log10(c(50,3e4)),
        xlab="melting T/°C", ylab="mean protein half-life/h", axes=FALSE,
        title=TRUE, cor.legend=FALSE)
for ( ax in 2 ) {
    axis(ax, at=0:9, labels=10^(0:9))
    axis(ax, at=log10(rep(1:10, 6) * 10^rep(-1:4, each=10)),
         tcl=-.125, labels=FALSE)
}
axis(1)
figlabel("human", pos="bottomleft", cex=1.2, font=2)
dev.off()

plotdev(file.path(fig.path,"protein_length_melting_human"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(hmap$mlt, log10(hmap$len), round=3,
        ##ylim=log10(c(1,3e4)), xlim=log10(c(50,3e4)),
        xlab="melting T/°C", ylab="protein length", axes=FALSE,
        title=TRUE, cor.legend=FALSE)
for ( ax in 2 ) {
    axis(ax, at=0:9, labels=10^(0:9))
    axis(ax, at=log10(rep(1:10, 6) * 10^rep(-1:4, each=10)),
         tcl=-.125, labels=FALSE)
}
axis(1)
figlabel("human", pos="bottomleft", cex=1.2, font=2)
dev.off()

plotdev(file.path(fig.path,"protein_melting_length_human"),
        type=ftyp, height=corH, width=corW, res=200)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(hmap$len), hmap$mlt, round=3,
        ##ylim=log10(c(1,3e4)), xlim=log10(c(50,3e4)),
        ylab="melting T/°C", xlab="protein length", axes=FALSE,
        title=TRUE, cor.legend=FALSE)
for ( ax in 1 ) {
    axis(ax, at=0:9, labels=10^(0:9))
    axis(ax, at=log10(rep(1:10, 6) * 10^rep(-1:4, each=10)),
         tcl=-.125, labels=FALSE)
}
axis(2)
figlabel("human", pos="bottomleft", cex=1.2, font=2)
dev.off()


## NOTE: low correlation of melting to predicted
## but both correlate well with length
plotCor(hmap$mlt, hmap$mlt_predicted)
plotCor(log(hmap$len),hmap$mlt_predicted)
plotCor(log(hmap$len),hmap$mlt)

hdat <- hmap[,3:ncol(hmap)]
##pca <- prcomp(hdat, scale=TRUE)


