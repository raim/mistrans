library(segmenTools)

##library(RiboProfiling)
##library(ribosomeProfilingQC)
##library(riboSeqR)
##library(RiboDIPA,ORFik)
## see https://github.com/zhpn1024/ribotish

## attempt to get the correct path for input data;
## these should be downloaded or generated as outlined
## in README.md
proj.path <- file.path(Sys.getenv("DECODE")) 
if ( proj.path=="" ) # author's local path
    proj.path <- "/home/raim/data/mistrans"
SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")

add.path <- file.path(proj.path,"additionalData")
mam.path <- file.path("/home/raim/data/mammary")

rseq.path <- file.path(proj.path, "figures", "riboseq")
dir.create(rseq.path)

ftyp <- "png"


ribo.file <- file.path(proj.path, "processedData",
                       "Iwasaki19_All.RiboProElong.bed")

ribo.file <- file.path(mam.path, "processedData",
                       "gwipsvizRiboseq.bed")

source("~/programs/segmenTools/R/parsers.R")
ribo <- bed2coor(ribo.file, header = c("chr", "start", "end", "score"))

## strand info by reverted start/end? would be illegal for bed, i think
any(ribo$end<ribo$start)

## fraction of multiple position values
sum(ribo$end==ribo$start)/nrow(ribo) #93% are single position values
sum(ribo$end>ribo$start)/nrow(ribo) # 6% cover several positions


## EXPAND RANGES TO FULL POSITIONS
lens <- ribo$end-ribo$start+1

## exponential decrease of covered positions, mostly 2
barplot(table(ribo$end-ribo$start+1))

## expand positions
riba <- data.frame(chr=rep(ribo$chr, lens),
                   start=rep(ribo$start, lens),
                   score=rep(ribo$score, lens))

head(riba)

## TODO: loop over diff(riba$start==0) and increase each by 1, until none are
## left
while( any(diff(riba$start)==0) ) {

    idx <- which(diff(riba$start)==0)+1
    riba$start[idx] <- riba$start[idx]+1
}
    
## number of non-consecutive runs
sum(diff(riba$start)!=1)

## investigate a range with many reads
plot(riba$start, log10(riba$score), xlim=c(6e5, 8e5), type="h")
## zoom in 
plot(riba$start, log10(riba$score), xlim=c(628e3, 635e3), type="h")

plotdev(file.path(rseq.path,paste0("riboseq_example")),
        res=300, type=ftyp, width=4, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
plot(riba$start, log10(riba$score), xlim=c(6298e2, 6299e2), type="h",
     ylab=expression(log[10](count)), xlab="chromosome position")
axis(1, at=1:1e6, labels=FALSE)
dev.off()

## chromosome length index
chr.file <- file.path(add.path,"sequenceIndex.csv.gz")
chrMap <- read.delim(chr.file)
chrS <- c(0,cumsum(as.numeric(chrMap[,3])))
chrL <- chrMap$length
chrIdx <- chrMap[,1]
names(chrIdx) <- chrMap[,2]

## convert AAS chromosome to index
##riba$chr <- chrIdx[riba$chr]
##riba$chr.orig <- riba$chr
##riba <- coor2index(riba, chrS=chrS)
##riba$chr <- riba$chr.orig

## codon period

chr <- 3

rm(all);
gc()

all <- rep(0, chrL[chr])
all[riba$start[riba$chr==chr]] <- log10(riba$score[riba$chr==chr])

racf <- acf(all, lag.max=150, plot=FALSE)


plotdev(file.path(rseq.path,paste0("riboseq_autocor")),
        res=300, type=ftyp, width=4, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
plot(racf$lag, racf$acf, type="h", axes=FALSE, ylim=c(0,1),
     ylab=expression(ACF), xlab="lag")
axis(2)
axis(1)#, at=seq(0,297,3), las=2)
dev.off()

## TODO:
## * load transcript coordinates and get total count,
## * get relative count for each codon,
## * plot relative count by RAAS.
