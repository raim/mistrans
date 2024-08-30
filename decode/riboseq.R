library(segmenTools)

## ANALYZE RIBOSEQ DATA

SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source(file.path(SRC.PATH, "raas_init.R"))

## additionally required data files
mam.path <- file.path(Sys.getenv("MAMDATA")) 
if ( mam.path=="" ) # author's local path
    mam.path <- "/home/raim/data/mammary"

## data to be added to additionalData folder
cds.file <- file.path(mam.path,"originalData", "transcript_coordinates.tsv")

## data from additionalData folder
chr.file <- file.path(add.path,"sequenceIndex.csv.gz")


## bed files generated from downloaded bigwig files 
ribo.file <- file.path(proj.path, "processedData",
                       "Iwasaki19_All.RiboProElong.bed")
ribo.file <- file.path(mam.path, "processedData",
                       "gwipsvizRiboseq.bed.gz")

## figure output path
rseq.path <- file.path(proj.path, "figures", "riboseq")
dir.create(rseq.path)

##source("~/programs/segmenTools/R/parsers.R")
cat(paste("loading riboseq file", ribo.file, "\n"))
ribo <- bed2coor(ribo.file, header = c("chr", "start", "end", "score"))

## some interactive QC and exploration
if ( interactive() ) { 
    ## strand info by reverted start/end? would be illegal for bed, i think
    any(ribo$end<ribo$start)

    ## fraction of multiple position values
    sum(ribo$end==ribo$start)/nrow(ribo) #93% are single position values
    sum(ribo$end>ribo$start)/nrow(ribo) # 6% cover several positions

    ## exponential decrease of covered positions, mostly 2
    barplot(table(ribo$end-ribo$start+1))
}

## EXPAND BED RANGES TO FULL POSITIONS
lens <- ribo$end-ribo$start+1

## expand positions
riba <- data.frame(chr=rep(ribo$chr, lens),
                   start=rep(ribo$start, lens),
                   score=rep(ribo$score, lens))
## loop over diff(riba$start==0) and increase each by 1, until none are left
while( any(diff(riba$start)==0) ) {
    idx <- which(diff(riba$start)==0)+1
    riba$start[idx] <- riba$start[idx]+1
}
    

## some interactive QC and exploration
if ( interactive() ) {
    ## number of non-consecutive runs
    sum(diff(riba$start)!=1)

    ## investigate a range with many reads
    plot(riba$start, log10(riba$score), xlim=c(6e5, 8e5), type="h")
    ## zoom in 
    plot(riba$start, log10(riba$score), xlim=c(628e3, 635e3), type="h")
}

## plot an example region
plotdev(file.path(rseq.path,paste0("riboseq_example")),
        res=300, type=ftyp, width=4, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
plot(riba$start, log10(riba$score), xlim=c(6298e2, 6299e2), type="h",
     ylab=expression(log[10](count)), xlab="chromosome position")
axis(1, at=1:1e6, labels=FALSE)
dev.off()

## CHROMOSOME LENGTH INDEX
## required for segmenTool's indexing scheme
chrMap <- read.delim(chr.file)
chrS <- c(0,cumsum(as.numeric(chrMap[,3])))
chrL <- chrMap$length
chrIdx <- chrMap[,1]
names(chrIdx) <- chrMap[,2]


## TEST CODON PERIOD
chr <- 3

if ( exists("all") ) rm(all);
gc()

## expand to full vector
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
mtext(paste("chromosome",chr), 3, 0)
dev.off()


### SUMMARIZE DATA FOR TRANSCRIPTS

## transcript coordinates
cds <- read.delim(cds.file)
cds <- cds[cds$transcript%in%genes$MANE,]

## NOTE: bed file had no strand info
## use transcripts w/o strand info
cdsi <- cds[,c("transcript","chr","start","end")]#, "strand")] 
cdsi$chr <- chrIdx[cdsi$chr]
cdsi <- coor2index(cdsi, chrS=chrS)


## convert chromosome to index
riba$chr <- chrIdx[riba$chr]
riba$chr.orig <- riba$chr
riba <- coor2index(riba, chrS=chrS)
riba$chr <- riba$chr.orig

## large vector for full chromosome of riboseq values
## NOTE: high memory usage!
rall <- rep(0, max(chrS))
rall[riba$start] <- riba$score

## get riboseq counts for CDS
## TODO: convert to GRanges and use this?
cdsl <- split(cdsi[,c("start","end")], cdsi$transcript)

## FOR EACH TRANSCRIPT:
## get length and sum all riboseq values
rvals <- lapply(cdsl, function(x) {
    coors <- unlist(apply(x,1, function(y) y[1]:y[2]))

    ncds <- length(coors)
    navail <- sum(rall[coors]>0)
    ## sum riboseq values if present
    vals <- sum(rall[coors])
    ## return
    c(len=ncds, n=navail, val=vals)
})
rvals <- as.data.frame(do.call(rbind, rvals) )

## per CDS nucleotide
rvals$valn <- rvals$val/rvals$len

## GET AAS VALUES

FILTER <- csite$fromto=="Q:G"

aas <- csite[FILTER,c("chr","coor")] # ignore strand,since no strand info in input
aas$chr <- chrIdx[aas$chr]
aai <- coor2index(aas, chrS=chrS)

## per nucleotide over codon
## NOTE: introduces errors for splice sites in codon

range <- -1:1
aar <- rep(0, nrow(aai))
for ( r in range )
    aar <- aar + rall[aai$coor + r]
aar <- aar/length(range)


## relative value
aan <- aar/rvals[csite$transcript[FILTER], "valn"]

plotdev(file.path(rseq.path,paste0("riboseq_raas_norm")),
        res=300, type=ftyp, width=2.5, height=2.5)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,.3,0), tcl=-.25)
plotCor(log2(aan), csite$RAAS[FILTER],
        xlab=expression(log[2](count[AAS]/count[transcript])), ylab=xl.raas,
        cor.legend=FALSE, title=TRUE)
dev.off()

plotdev(file.path(rseq.path,paste0("riboseq_raas_raw")),
        res=300, type=ftyp, width=2.5, height=2.5)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,.3,0), tcl=-.25)
plotCor(log10(aar), csite$RAAS[FILTER],
        xlab=expression(log[10](count[AAS])), ylab=xl.raas,
        cor.legend=FALSE, title=TRUE)
dev.off()

if ( FALSE ) {
    ## slow but less memory: 
    rvals <- lapply(cdsl, function(x) {
        
        ## get transcript CDS coordinates
        coors <- unlist(apply(x,1, function(y) y[1]:y[2]))
        if ( any(duplicated(coors)) ) stop("overlapping cds")
        
        ## count coordinates that are present
        idx <- which(riba$start%in%coors)
        ncds <- length(coors)
        navail <- length(idx)
        ## sum riboseq values if present
        vals <- ifelse(navail>0, sum(riba[idx,"score"]), 0)
        ## return
        c(len=ncds, n=navail, val=vals)
    })
}

## TODO:
## * load transcript coordinates and get total count,
## * get relative count for each codon,
## * plot relative count by RAAS.
