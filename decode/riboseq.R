library(segmenTools)

## ANALYZE RIBOSEQ DATA

SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source(file.path(SRC.PATH, "raas_init.R"))

do.autocor <- !interactive()

## figure output path
rseq.path <- file.path(fig.path, "riboseq") #,RID)
dir.create(rseq.path)

## additionally required data files
mam.path <- file.path(Sys.getenv("MAMDATA")) 
if ( mam.path=="" ) # author's local path
    mam.path <- "/home/raim/data/mammary"

## data to be added to additionalData folder
cds.file <- file.path(mam.path,"originalData", "transcript_coordinates.tsv")

## data from additionalData folder
chr.file <- file.path(add.path,"sequenceIndex.csv.gz")

motif.file <- file.path(proj.path, "processedata", "saap_motifs.tsv")



## CHROMOSOME LENGTH INDEX
## required for segmenTool's indexing scheme
chrMap <- read.delim(chr.file)
chrS <- c(0,cumsum(as.numeric(chrMap[,3])))
chrL <- chrMap$length
chrIdx <- chrMap[,1]
names(chrIdx) <- chrMap[,2]

## TRANSCRIPT COORDINATES (MANE only)
cds <- read.delim(cds.file)
cds <- cds[cds$transcript%in%genes$MANE,]

## NOTE: bed file had no strand info
## use transcripts w/o strand info
cdsi <- cds[,c("transcript","chr","start","end")]#, "strand")] 
cdsi$chr <- chrIdx[cdsi$chr]
cdsi <- coor2index(cdsi, chrS=chrS)

## exons as list per transcript
cdsl <- split(cdsi[,c("start","end")], cdsi$transcript)
        

## RIBO-SEQ DATA SETS TO ANALYZE
## bed files generated from downloaded bigwig files
dataSets <- c("Iwasaki19_All.RiboProElong",
              "human_eIF3b_bound_40S.RiboProElong",
              "gwipsvizRiboseq",
              "Iwasaki16_All.RiboCov",
              "Chen20_All.RiboCov"
              )

ds <- dataSets[1]
cors <- NULL
for ( ds in dataSets ) {

    ribo.file <- file.path(proj.path, "processedData",
                           paste0(ds, ".bed.gz"))

    RID <- sub("\\.bed$", "", sub("\\.gz$","",basename(ribo.file)))
    
    cat(paste("loading riboseq data set", ds, "\n"))
    ribo <- bed2coor(ribo.file, header = c("chr", "start", "end", "score"))
    
    ## some interactive QC and exploration
    if ( FALSE ) { 
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
    if ( FALSE ) {
        ## number of non-consecutive runs
        sum(diff(riba$start)!=1)
        
        ## investigate a range with many reads
        plot(riba$start, log10(riba$score), xlim=c(6e5, 8e5), type="h")
        ## zoom in 
        plot(riba$start, log10(riba$score), xlim=c(628e3, 635e3), type="h")
    }
    
    ## plot an example region
    plotdev(file.path(rseq.path,paste0(RID,"_example")),
        res=300, type=ftyp, width=4, height=3)
    par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
    plot(riba$start, log10(riba$score), xlim=c(6298e2, 6299e2), type="h",
         ylab=expression(log[10](count)), xlab="chromosome position")
    axis(1, at=1:1e6, labels=FALSE)
    dev.off()

    
    ## TEST CODON PERIOD
    if ( do.autocor ) {
        chr <- 3
        
        if ( exists("all") ) rm(all);
        gc()
        
        ## expand to full vector
        all <- rep(0, chrL[chr])
        all[riba$start[riba$chr==chr]] <- log10(riba$score[riba$chr==chr])
        
        racf <- acf(all, lag.max=150, plot=FALSE)
        
        plotdev(file.path(rseq.path,paste0(RID,"_autocor")),
                res=300, type=ftyp, width=4, height=3)
        par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
        plot(racf$lag, racf$acf, type="h", axes=FALSE, ylim=c(0,1),
             ylab=expression(ACF), xlab="lag")
        axis(2)
        axis(1)#, at=seq(0,297,3), las=2)
        mtext(paste("chromosome",chr), 3, 0)
        dev.off()
    }
    
### SUMMARIZE DATA FOR TRANSCRIPTS
    ## TODO: convert to GRanges and use this?
    

    ## convert chromosome to index
    riba$chr <- chrIdx[riba$chr]
    riba$chr.orig <- riba$chr
    riba <- coor2index(riba, chrS=chrS)
    riba$chr <- riba$chr.orig
    
    ## large vector for full chromosome of riboseq values
    ## NOTE: high memory usage!
    rall <- rep(0, max(chrS))
    rall[riba$start] <- riba$score
    

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
    
    usite <- site
    
    aas <- usite[,c("chr","coor")] # ignore strand,since no strand info in input
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
    aan <- aar/rvals[usite$transcript, "valn"]
    
    ## clean up to avoid memory problems
    rm(rall)
    rm(ribo)
    rm(riba)
    gc()
    
    ## PLOT
    
    ## interactive QC and data exploration
    if ( FALSE ) {
        ##boxplot(log2(aan) ~ usite$fromto)
        plotCor(as.numeric(usite$iupred3), log(aan))
        plotCor(log10(as.numeric(usite$protein.intensity)), log(aan))
        plotCor(as.numeric(usite$MMSeq2), log(aan)) # <- SIGNIFICANT
        
    }
    
    
    ## subset to certain types?
    FILTER <- rep(TRUE, nrow(usite)) # usite$fromto=="Q:G" # 

    ## MAIN: correlation RAAS to relative read count at AAS
    plotdev(file.path(rseq.path,paste0(RID,"_raas_norm")),
            res=300, type=ftyp, width=2.5, height=2.5)
    par(mai=c(.5,.5,.25,.15), mgp=c(1.4,.3,0), tcl=-.25)
    cr <- plotCor(log2(aan[FILTER]), usite$RAAS.median[FILTER],
                  xlab=expression(log[2](C[AAS]/C[transcript])), ylab=xl.raas,
                  cor.legend=FALSE, title=TRUE)
    figlabel(RID, pos="bottomleft", cex=.7)
    dev.off()

    ## STORE
    cors[[RID]] <- cr
    
    plotdev(file.path(rseq.path,paste0(RID,"_raas_raw")),
            res=300, type=ftyp, width=2.5, height=2.5)
    par(mai=c(.5,.5,.25,.15), mgp=c(1.4,.3,0), tcl=-.25)
    plotCor(log10(aar[FILTER]), usite$RAAS.median[FILTER],
            xlab=expression(log[10](count[AAS])), ylab=xl.raas,
            cor.legend=FALSE, title=TRUE)
    dev.off()
    
    ## transcripts
    ## only for non-0 AAS
    
    faas <- FILTER & aar!=0
    
    plotdev(file.path(rseq.path,paste0(RID,"_raas_transcript")),
            res=300, type=ftyp, width=2.5, height=2.5)
    par(mai=c(.5,.5,.25,.15), mgp=c(1.4,.3,0), tcl=-.25)
    plotCor(log10(rvals[usite$transcript[faas], "valn"]),
            usite$RAAS.median[faas],
            xlab=expression(log[10](count[transcript])), ylab=xl.raas,
            cor.legend=FALSE, title=TRUE)
    dev.off()
    
    
    plotdev(file.path(rseq.path,paste0(RID,"_aas_transcript")),
            res=300, type=ftyp, width=2.5, height=2.5)
    par(mai=c(.5,.5,.25,.15), mgp=c(1.4,.3,0), tcl=-.25)
    plotCor(log10(rvals[usite$transcript[FILTER], "valn"]),
            log10(aar[FILTER]),
            xlab=expression(log[10](count[transcript])),
            ylab=expression(log[10](count[AAS])),
            cor.legend=FALSE, title=TRUE)
    dev.off()

    
    ## TODO:
    ## * load motif classification and analyze specifically by motif
    ## * dotplot, riboseq count bins vs. Datasets
    
    ## expand to tmtu
    
    tmtu <- tmtf
    
    idx <- match(paste(tmtu$BP, tmtu$SAAP), paste(usite$BP, usite$SAAP))
    
    ## TODO: why missing?
    ina <- which(is.na(idx))
    
    ## REMOVE MISSING BP/SAAP
    if ( length(ina)>0 ) {
        
        cat(paste("TODO:", length(ina), "BP/SAAP missing from BP/SAAP file.\n"))
        tmtu <- tmtu[-ina,]
        idx <- idx[-ina]    
    }
    riboseq <- rep(NA, nrow(tmtu))
    riboseq <- log2(aan[idx])
    riboseq[is.infinite(riboseq)] <- NA
    rmx <- 3
    riboseq[riboseq>  rmx] <-  rmx
    riboseq[riboseq< -rmx] <- -rmx
    riboseq.bins <- cut(riboseq,
                        breaks=seq(min(riboseq,na.rm=TRUE),
                                   max(riboseq,na.rm=TRUE), length.out=6))
    levels(riboseq.bins) <- c(levels(riboseq.bins), "NA")
    riboseq.bins[is.na(riboseq.bins)] <- "NA"
    tmtu$riboseq.bins <-riboseq.bins
    
    ovls <- raasProfile(x=tmtu, id="SAAP", 
                        rows="riboseq.bins", cols="Dataset",
                        row.srt=rev(levels(riboseq.bins)),
                        col.srt=uds,
                        bg=TRUE, value="RAAS",
                        use.test=use.test, do.plots=FALSE,
                        xlab=xl.raas,
                        verb=0)
    
    omai <- c(.75,1.25,.5,.5)
    nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
    nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

    plotdev(file.path(rseq.path,paste0(RID,"_dotplot")),
            height=nh, width=nw, res=300, type=ftyp, bg="white")
    par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
    dotprofile(ovls, value="median",
               p.dot=p.dot, dot.sze=dot.sze, vbrks=abrks, vcols=acols, 
               xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8, axis=1:2)
    ##mtext(expression(log[2](count~ratio)),2, 3.3)
    dev.off()

}

## save results
save(cors, file=file.path(out.path, "riboseq_correlations.rda"))
