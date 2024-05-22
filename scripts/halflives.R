
library(segmenTools)
options(stringsAsFactors=FALSE)

data.path <- "~/data/"
mam.path <- file.path(data.path, "mammary")
yeast.path <- file.path(data.path, "yeast")

fig.path <- file.path(mam.path,"figures","halflives")
dir.create(fig.path)
ftyp <- "png"

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

hhld <-  readxl::read_xlsx(math18.file)
## mean half-live over all replicates and cell types
## TODO: consider distribution
hhld <- hhld[,-grep("Mouse", colnames(hhld))]
hhld <- hhld[,grep("Hepatocytes", colnames(hhld))]
cidx <- grep("half_life", colnames(hhld), value=TRUE)
hhlf <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
names(hhlf) <- unlist(hlvd[,1])

##hhlf[hhlf>1000] <- 1000

yhld <- read.delim(yhlf.file)
yhld$hgene <-
    y2h[yhld[,1]]
yhld <- yhld[!is.na(yhld$hgene),]
yhld <- yhld[!duplicated(yhld$hgene),]
rownames(yhld) <- yhld$hgene

xl.hlf <- expression(protein~"half-life"/h)
xl.hhlf <- expression(human~protein~"half-life"/h)
xl.yhlf <- expression(yeast~protein~"half-life"/h)
xl.hdeg <- expression(human~degradation/h^-1)
xl.ydeg <- expression(yeast~degradation/h^-1)

xl.hlgd <- expression(human~degradation/log(d^-1))
xl.ylgd <- expression(yeast~degradation/log(d^-1))


yhlf <- yhld[names(hhlf),"t1.2..hours."]
yhlf[yhlf==">= 100"] <- 120
yhlf <- as.numeric(yhlf)

# deg rate in days
hdeg <- log(2)/(hhlf)
ydeg <- log(2)/(yhlf)

hlgd <- log(log(2)/(hhlf/24))
ylgd <- log(log(2)/(yhlf/24))


plotdev(file.path(fig.path,"protein_halflives_human_yeast_zoom"),
        type=ftyp, height=3, width=4, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hhlf, yhlf, xlim=c(0,200), ylim=c(0,20),
        xlab=xl.hhlf,
        ylab=xl.yhlf)#, legpos="topright")
##abline(v=48)
dev.off()

plotdev(file.path(fig.path,"protein_halflives_human_yeast"),
        type=ftyp, height=3, width=4, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hhlf, yhlf, 
        xlab=xl.hhlf,
        ylab=xl.yhlf)#, legpos="topright")
##abline(v=48)
dev.off()

plotdev(file.path(fig.path,"protein_degradation_human_yeast"),
        type=ftyp, height=4, width=4, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hdeg, ydeg,
        xlab=xl.hdeg,
        ylab=xl.ydeg)#, legpos="topright")
showl <- (hlgd < -4 & ylgd < -1)
showh <- (ylgd > 1 & hlgd > 0) | ylgd> 3
showb <- (hlgd < -4) | (hlgd > -1 & ylgd < 0)
showby <- ylgd < -1 & hlgd < -1
show <- TRUE

dev.off()

plotdev(file.path(fig.path,"protein_degradation_human_yeast_zoom"),
        type=ftyp, height=4, width=4, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hdeg, ydeg,
        xlab=xl.hdeg,
        ylab=xl.ydeg, ylim=c(0,.3))#, legpos="topright")
showl <- (hlgd < -4 & ylgd < -1)
showh <- (ylgd > 1 & hlgd > 0) | ylgd> 3
showb <- (hlgd < -4) | (hlgd > -1 & ylgd < 0)
showby <- ylgd < -1 & hlgd < -1
show <- TRUE

dev.off()

## log as in multi-species paper 
plotdev(file.path(fig.path,"protein_degradation_human_yeast_log"),
        type=ftyp, height=4, width=4, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hlgd, ylgd,
        xlab=xl.hlgd,
        ylab=xl.ylgd)#, legpos="topright")
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

plotCor(log10(ygenes$p.length), log10(yhlf), round=3)


## HUMAN
hgenes <- read.delim(hfeature.file)
filter <- hgenes$type=="protein_coding" & hgenes$MANE!=""
hgenes <- hgenes[filter,]

hhld <-  readxl::read_xlsx(math18.file)
## mean half-live over all replicates and cell types
## TODO: consider distribution
hhld <- hhld[,-grep("Mouse", colnames(hhld))]
hhld <- hhld[,grep("Hepatocytes", colnames(hhld))]
cidx <- grep("half_life", colnames(hhld), value=TRUE)
hhlf <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
names(hhlf) <- unlist(hlvd[,1])

## protein lengths
hlen <- read.delim(hlen.file, header=FALSE, row.names=1)
rownames(hlen) <- sub("\\.[0-9]*", "", rownames(hlen))

## expand proteins and gene names
hmap <- strsplit(hgenes$proteins,";")
hmap <- cbind(rep(hgenes$name, lengths(hmap)),
              unlist(hmap))
hmap <- cbind.data.frame(hmap,
                         len=hlen[hmap[,2],1],
                         hlf=hhlf[hmap[,1]])

## remove duplicated gene names AFTER SORTING FOR LENGTH
## NOTE: correlation much stronger when taking the longest proteoform
## for each gene.
## NOTE/TODO: select MANE PROTEIN INSTEAD
hmap <- hmap[order(hmap$len, decreasing=TRUE),]
hmap <- hmap[!duplicated(hmap[,1]),]

plotCor(log10(hmap$len), log10(hmap$hlf), round=3,
        xlab="protein length", ylab="protein half-life/h")

### TODO: add MANE.protein column to human feature file
## or generate protein file
