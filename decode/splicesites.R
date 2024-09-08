library(segmenTools)
library(stringr)

### ANALYZE SPLICE SITES
### TODO: extract sequences and clustering

dat.path <- "/home/raim/data"
mam.path <- file.path(dat.path, "mammary")
mis.path <- file.path(dat.path, "mistrans","processedData")
fig.path <- file.path(dat.path, "mistrans","figures","splicing")

dir.create(fig.path)

ftyp <- "png"

feature.file <- file.path(mam.path,"features_GRCh38.110.tsv.gz")
cdsmap <- file.path(mam.path,"processedData","protein_cds_structure.dat")
cdspos <- file.path(mam.path,"processedData","protein_cds_coordinates.tsv")
pfasta <- file.path(mis.path,"all_proteins.fa")
tfasta <- file.path(mam.path,"processedData","coding.fa")
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")


## list of all MANE proteins
genes <- read.delim(feature.file)
genes <- genes[genes$proteins!="" & !is.na(genes$proteins),]
##genes <- genes[genes$name!="" & !is.na(genes$name),]
noname <- genes$name=="" | is.na(genes$name)
genes$name[noname] <- genes$ID[noname]
ptl <- strsplit(genes$proteins, ";")
ens2nam <- rep(genes$name, lengths(ptl))
names(ens2nam) <- unlist(ptl)
pnms <- ens2nam

MANES <- genes$MANE.protein
MANES <- MANES[MANES!=""]

## protein sequences

fas <- readFASTA(pfasta, grepID=TRUE)
fas <- fas[names(fas)%in%MANES]
seql <- lapply(fas, function(x) x$seq)



## get transcript exon map
cds <- read.delim(cdsmap,header=FALSE, row.names=1)
cdl <- strsplit(cds[,1],";")
names(cdl) <- rownames(cds)
cdl <- lapply(cdl, as.numeric)
## genomic position of CDS
cpos <- read.delim(cdspos, header=FALSE, sep=" ")
colnames(cpos) <- c("ID","chr","start","end","strand")
posl <- split(cpos[,2:5], cpos[,1])
## NOTE: reverse order for minus strand
posl <- lapply(posl, function(x)
    x[order(x$start, decreasing=x$strand=="-"),])
## splice sites in protein coordinates
PRS <- lapply(cdl, function(x) x/3)

AA <- unique(Biostrings::GENETIC_CODE)
AA <- AA[AA!="*"]
## "[KR]..Q" strongest of the motifs
patterns <- unique(c(AA, "H","G","AQ","[KR]A","[KR][KR]","A", "[KR]", "Q", "[KR]A[A-Z]?Q", "[KR]..Q"))

pids <- fs::path_sanitize(patterns)
names(pids) <- patterns

for ( pat in patterns ) {

    pid <- pids[pat]


    ## extract actual PATTERN 
    kraqs <- stringr::str_extract_all(seql, pat)

    ## get positions of PATTERN
    kpos <- gregexpr(pat, seql)
    names(kpos) <- names(seql)
    
    ## FIND DISTANCES OF THE MOTIF (START) TO CLOSEST
    ## SPLICE SITE
    
    ## splice sites in AA coordinates
    pdl <- PRS[names(seql)]
    
    ## remove all sequences without PATTERN
    rmk <- unlist(lapply(kpos, function(x) x[1]==-1))
    kpos <- kpos[!rmk]
    pdl <- pdl[!rmk]
    kraqs <- kraqs[!rmk]
    ## loop through all proteins with PATTERN
    ssd <- exl <- NULL
    for ( i in seq_along(kpos) ) {
        hits <- kpos[[i]]
        cds <- pdl[[i]]
        len <- diff(c(0,cds)) #TODO: correct exon length?
        ssd[[i]] <- exl[[i]] <- rep(NA, length(hits))
        for ( j in seq_along(hits) ) {
            ssdst <- hits[j] - cds[unique(1:(length(cds)-1))]
            idx <- which.min(abs(ssdst))
            ssd[[i]][j] <- ssdst[idx] # minimal distance to splice site
            exon <- which(cds>=hits[j])[1] # in which exon is the start?
            exl[[i]][j] <- len[exon]
        }
    }
    names(ssd) <- names(exl) <- names(kpos)
    
    ## match lengths
    klen <- lapply(kpos, function(x) attributes(x)$match.length)
    
    ## check positions
    kn <- unlist(lapply(ssd, function(x) any(x==0)))
    kid <- names(which(kn)[1])
    
    cdl[[kid]]
    pdl[[kid]]
    ssd[[kid]]
    kpos[[kid]]
    
    ## NOTE: table on fractions doesn't work, needs rounding!
    sst <- table(round(unlist(ssd),2))
    dst <- as.numeric(names(sst))
    class(sst) <- "numeric"
    
    plotdev(file.path(fig.path,paste0(pid, "_absolutepos_bar")),
            height=3.5, width=5, res=300, type=ftyp)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(dst, sst, type="h", xlim=c(-40,40),
         xlab="distance to closest splice site", ylab="Frequency")
    legend('topright', pat)
    ##points(dst,sst)
    axis(3, labels=FALSE)
    dev.off()
    plotdev(file.path(fig.path,paste0(pid, "_absolutepos_bar_zoom")),
            height=3.5, width=5, res=300, type=ftyp)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(dst, sst, type="h", xlim=c(-10,10),
         xlab="distance to closest splice site", ylab="Frequency", axes=FALSE)
    ##points(dst,sst)
    axis(1, at=-10:10)
    axis(2)
    axis(3, labels=FALSE)
    axis(3, at=-10:10, labels=FALSE)
    legend('topright', pat)
    dev.off()
    plotdev(file.path(fig.path,paste0(pid, "_absolutepos_smooth")),
            height=3.5, width=5, res=300, type=ftyp)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(dst, ma(sst,10), type="l", xlim=c(-40,40),
         xlab="distance to closest splice site",
         ylab="Frequency (moving average)")
    axis(3, labels=FALSE)
    legend('topright', pat)
    dev.off()


    sst <- table(round(unlist(ssd),0))
    dst <- as.numeric(names(sst))
    class(sst) <- "numeric"
    plotdev(file.path(fig.path,paste0(pid, "_absolutepos_bar_round")),
            height=3.5, width=5, res=300, type=ftyp)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(dst, sst, type="h", xlim=c(-40,40),
         xlab="distance to closest splice site", ylab="Frequency")
    axis(1, at=-100:100, tcl=-.125, labels=FALSE)
    axis(3, labels=FALSE)
    legend('topright', pat)
    dev.off()
    
    
    rpos <- unlist(ssd)/unlist(exl)
    plotdev(file.path(fig.path,paste0(pid, "_relativepos_hist")),
            height=3.5, width=5, res=300, type=ftyp)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    hist(rpos, breaks=100,
         xlab="distance to closest splice site / exon length", main=NA)
    axis(3, labels=FALSE)
    legend('topright', pat)
    dev.off()
    
}

## extract splice sites

tfas <- readFASTA(tfasta, grepID=TRUE)
tfas <- tfas[names(tfas)%in%genes$MANE]
tseq <- lapply(tfas, function(x) x$seq)

## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
## reverse map transcript-protein
pamrt <-rownames(trmap)
names(pamrt) <- trmap[,1]

## rename transcripts
names(tseq) <- pamrt[names(tseq)]

## get all splice site sequences

rng <- 18
ass <- list()
asp <- list()
for ( i in seq_along(tseq) ) {
    pid <- names(tseq)[i]
    n <- length(cdl[[pid]])-1
    if ( n==0 ) next
    ass[[pid]] <- rep("", n)
    asp[[pid]] <- rep(NA, n)
    for ( j in 1:n ) {
        pos <- cdl[[pid]][j]
        ass[[pid]][j] <- substr(tseq[[pid]], pos-rng+1, pos+rng)
        asp[[pid]][j] <- pos%%3 # codon frame
    }
}

## add separator for unlist-derived counter
names(ass) <- names(asp) <- paste0(names(ass), "_")

asp <- unlist(asp)
ass <- unlist(ass)

## TODO: translate using frame info
## and plot sequence logos for each frame separately

## remove all non-full length splice sites (close to ends)
rms <- nchar(ass) != rng*2
asp <- asp[!rms]
ass <- ass[!rms]

asm <- do.call(rbind, strsplit(ass, ""))
colnames(asm) <- c(-rng:-1,1:rng)

ssf <- rbind(A=apply(asm, 2, function(x) sum(x=="A")),
             C=apply(asm, 2, function(x) sum(x=="C")),
             G=apply(asm, 2, function(x) sum(x=="G")),
             T=apply(asm, 2, function(x) sum(x=="T")))

plotdev(file.path(fig.path,paste0("splicesite_nucleotides")),
        height=3.5, width=7.5, res=300, type=ftyp)
par(mai=c(.5,.75,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
bp <- barplot(ssf, las=2, beside=TRUE, col=1:4, legend=TRUE)
mtext("distance from splice site",1, 1.3)
dev.off()

## NOTE:
## consensus by phase
##  0: [AC]AG - AAG: K, CAG: Q,
## +1: AG[AG]: R,
## +2: G[GA][TA] - GG[TA]: G, GAT: E, GAA: D. 
