
library(segmenTools)
options(stringsAsFactors=FALSE)

data.path <- "~/data/"
mam.path <- file.path(data.path, "mammary")
yeast.path <- file.path(data.path, "yeast")

fig.path <- file.path(mam.path,"figures","halflives")
dir.create(fig.path)
ftyp <- "png"

## yeast protein halflives - @Christiano2014
yhlf.file <- file.path(yeast.path,"processedData","christiano14.csv")
## human protein halflives - @Mathieson2018
hhlf.file <- file.path(mam.path,"originalData",
                       "41467_2018_3106_MOESM5_ESM.xlsx")
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

yhlf <- yhld[names(hhlf),"t1.2..hours."]
yhlf[yhlf==">= 100"] <- 100
yhlf <- as.numeric(yhlf)

plotdev(file.path(fig.path,"protein_halflives_human_yeast_zoom"),
        type=ftyp, height=3, width=4, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hhlf, yhlf, xlim=c(0,200), ylim=c(0,20),
        xlab="protein halflive/h, human",
        ylab="protein halflive/h, yeast")#, legpos="topright")
##abline(v=48)
dev.off()

plotdev(file.path(fig.path,"protein_halflives_human_yeast"),
        type=ftyp, height=3, width=4, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(hhlf, yhlf, 
        xlab="protein halflive/h, human",
        ylab="protein halflive/h, yeast")#, legpos="topright")
##abline(v=48)
dev.off()

plotdev(file.path(fig.path,"protein_degradation_human_yeast"),
        type=ftyp, height=3, width=4, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(log(log(2)/(hhlf/24)), log(log(2)/(yhlf/24)),
        xlab="log deg. rate/d, human",
        ylab="log deg. rate/d, yeast")#, legpos="topright")
##abline(v=48)
dev.off()

hist(hhlf,breaks=100)
hist(log(2)/hhlf, breaks=100)
hist(log(log(2)/(hhlf/24)))
