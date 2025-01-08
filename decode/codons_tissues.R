
## CODON FREQUENCY ANALYSIS of AMINO ACID SUBSTITUTION SITES

## TODO 20241220:

## * load tissue-wise transcript counts, table EV2 by Eraslan et
##   al. 2019 10.15252/msb.20188513,
## * calculate scaled codon frequencies: transcript*count for each tissue,
## * correlate with tissue-specific RAAS values.

SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source(file.path(SRC.PATH, "raas_init.R"))

    

## local output path
ctfig.path <- file.path(fig.path,"codons", "tissues")
dir.create(ctfig.path, showWarnings=FALSE)

### START ANALYSIS


### CODON FREQUENCIES vs. RAAS

## ONLY USE RAAS WITH ASSIGNED CODONS
ctmt <- tmtf[tmtf$codon!="",]


## add columns for codon positions
cpos <- strsplit(ctmt$codon,"")
ctmt$pos1 <- unlist(lapply(cpos, function(x) x[1]))
ctmt$pos2 <- unlist(lapply(cpos, function(x) x[2]))
ctmt$pos3 <- unlist(lapply(cpos, function(x) x[3]))

## LOAD GLOBAL GENE-WISE CODON COUNTS
codons <- read.delim(codon.file, row.names=1)

## Wu et al. 2019: Codon Stability Coefficient
csc <- read.csv(wu19.file, row.names=1)

## Dana and Tuller 2014/2015
## decoding time/sec -> 1/decoding rate
decod <- 1/read.csv(dana14.file,
                    row.names=1)[,"H..sapiens5.HEK293",drop=FALSE]

## relative rate per AA
decodl <-
    split(decod, GENETIC_CODE[rownames(decod)])
decodl <- lapply(decodl, function(x) (x-min(x[,1]))/(max(x[,1]-min(x[,1]))))
decodl <- do.call(rbind, decodl)
decodl[is.na(decodl)] <- 1
rownames(decodl) <- sub("M", "M.ATG", rownames(decodl))
rownames(decodl) <- sub("W", "W.TGG", rownames(decodl))
rownames(decodl) <- sub(".*\\.", "", rownames(decodl))


## CODON COUNT and FREQUENCY IN ALL UNIQUE TRANSCRIPTS 

if ( any(!ctmt$transcript%in%rownames(codons)) )
    stop("some codon frequencies not found, rerun calculation")

## codon counts in all mapped transcripts
## NOTE: this is used as our background frequency
cod  <- codons[unique(ctmt$transcript),]
codt <- apply(cod,2,sum) # total count

### CODON SORTING

## PER AA
## transcript codon frequencies per AA and SORTING by frequency and AA prop
codl <- split(codt, GENETIC_CODE[names(codt)])
codl <- lapply(codl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT

if ( FALSE ) # 20241225: to keep stop codon, TODO: sort Fbg below!
    codl <- codl[aap.srt[aap.srt%in%names(codl)]] ## SORT AA BY PROP

## LOCAL CODON SORTING by background frequencies
## (AA property ->codon frequency)
codon.srt <- sub("\\.","-",names(unlist(codl)))
## remove K|R codons
codon.srt <- codon.srt[grep("[KR]", codon.srt, invert=TRUE)]

## codon frequencies per AA
Fbg <- lapply(codl, function(x) x/sum(x)) # codon frequency
Fbg <- unlist(Fbg)
names(Fbg) <- sub("\\.","-",names(Fbg))



### TISSUE-SPECIFIC CODON FREQUENCIES (202412)
## for each tissue multiple codon per transcript count by
## tissue-specific transcript count, and calculate a tissue-specific
## codon frequency.

## tissue-specific transcript counts
eraslan19 <- read.delim(eraslan19.file, row.names=3)

## TODO: 500 transcripts from tissue-specific data sets are missing from the
## the codon count table! why?
## gene IDs seem available, but transcript IDs aren't, map from gene names
## TODO: load full gene transcript mapping

## remap via gene ID to transcript MANE IDs
eraslan19$transcript <- genes$MANE[match(eraslan19$EnsemblGeneID, genes$ID)]

## TODO: messages of missing numbers 
sum(is.na(eraslan19$transcript)) # only 37 missing
sum(eraslan19$transcript!=rownames(eraslan19), na.rm=TRUE) # 5250 different transcripts!
sum(duplicated(eraslan19$transcript)) # 139 duplicated transcripts


## TODO: use MANE transcript match only where transcript is missing
## from codon table!
##rownames(eraslan19) <- eraslan19$transcript

tids <- data.frame(TID=rownames(eraslan19),
                   GID=eraslan19$EnsemblGeneID,
                   name=eraslan19$GeneName,
                   TMANE=eraslan19$transcript)

tcnts <- eraslan19[,grep("exonic",colnames(eraslan19))]
colnames(tcnts) <- sub("_exonicMRNA", "",  colnames(tcnts))

## mean counts over replicates
## TODO: inspect standard deviations
tiss <- unique(sub("_.*", "", colnames(tcnts)))
tcnt <- sapply(tiss, function(x)
    apply(tcnts[,grep(x,colnames(tcnts)),drop=FALSE],1,mean))

## find missing matches and fill up via gene ID transcripts
## TODO: record that cs. 500 IDs were not found and were replace
## by MANE transcript
tids$match <- tids$TID
missing <- !tids$match%in%rownames(codons)
tids$match[missing] <- tids$TMANE[missing]

missing <- !tids$match%in%rownames(codons)

## remove 29 still missing

tids <- tids[!missing,]
tcnt <- tcnt[!missing,]
rownames(tcnt) <- tids$match

## list of tissue-specific codons counts

tcodons <- codons[rownames(tcnt),]
tcodcnts <- lapply(tiss, function(x) tcodons*tcnt[,x])
names(tcodcnts) <- tiss

source(file.path(SRC.PATH, "raas_utils.R"))
Fbg2 <- codonFrequencies(cod, sort=TRUE, unlist=TRUE)
if ( any(Fbg!=Fbg2[names(Fbg)]) )
    stop("new codon frequency function doesn't work")

## AA-specific codon frequencies
tcodfreq <- lapply(tcodcnts, codonFrequencies, sort=TRUE, unlist=TRUE)


## sort and combine tissue-specific codon frequencies
tmpnms <- names(tcodfreq[[1]])
tcodfreq <- lapply(tcodfreq, function(x) x[tmpnms])
tcodfreq <- do.call(rbind, tcodfreq)

if ( interactive() ) {
    hist(tcodfreq[,1]) ## NOTE: large differences in stop codon usage
    hist(tcodfreq[,4])
}

plotdev(file.path(ctfig.path,paste0("tissues_codons_frequencies")),
        type=ftyp, res=300, width=12,height=6)
par(mai=c(.65,1.5,.1,.1), mgp=c(1.3, .3, 0), tcl=-.25, family='monospace')
image_matrix(rbind(tcodfreq,
                   raw=Fbg[match(colnames(tcodfreq), names(Fbg))]),
             axis=1:2, col=viridis(100),
             breaks=seq(0,1,length.out=101), xlab=NA, ylab=NA)
abline(v=.5+cumsum(table(sub("-.*","", colnames(tcodfreq)))),
       col='black', lwd=.8, xpd=TRUE)
abline(v=.5+cumsum(table(sub("-.*","", colnames(tcodfreq)))),
       col='white', lwd=1)
dev.off()

## TODO: huge heatmap of per transcript codon frequencies.
if ( FALSE ) {
    rcodons <- codons/apply(codons,1,sum)
    plotdev(file.path(ctfig.path,paste0("all_codons_frequencies")),
            type=ftyp, res=300, width=12,height=100)
    par(mai=c(.65,1.5,.1,.1), mgp=c(1.3, .3, 0), tcl=-.25, family='monospace')
    image_matrix(rcodons,
             axis=1:2, col=viridis(100),
             breaks=seq(0,1,length.out=101), xlab=NA, ylab=NA)
    abline(v=.5+cumsum(table(sub("-.*","", colnames(tcodfreq)))),
           col='white', lwd=1)
    dev.off()
}

## REPORT AVAILABLE AND MISSING TISSIES
havefreq <- tolower(rownames(tcodfreq))
haveraas <- unique(ctmt$TMT.Tissue[ctmt$Dataset=="Healthy"])

cat(paste("codon frequencies available but no RAAS:",
          paste(havefreq[!havefreq%in%haveraas], collapse="; "), "\n"))
cat(paste("RAAS available but no codon frequencies:",
          paste(haveraas[!haveraas%in%havefreq], collapse="; "), "\n"))

## map different tissue names between data sets
tmap <- setNames(havefreq, havefreq)
tmap["salivarygland"] <- "saliva-secretinggland"
tmap["brain"] <- "cerebralcortex"
tmap["fat"] <- "adiposetissue"
tmap["thyroid"] <- "thyroidgland"
tmap["appendices"] <- "vermiformappendix"
tmap["prostate"] <- "prostategland"

havefreq <- tmap[tolower(rownames(tcodfreq))]
haveraas <- unique(ctmt$TMT.Tissue[ctmt$Dataset=="Healthy"])

cat(paste("codon frequencies available but no RAAS:",
          paste(havefreq[!havefreq%in%haveraas], collapse="; "), "\n"))
cat(paste("RAAS available but no codon frequencies:",
          paste(haveraas[!haveraas%in%havefreq], collapse="; "), "\n"))

## CORRELATE TO RAAS

## global scaled codon frequency
acnt <- apply(tcnt, 1, sum)
acodcnts <- tcodons*acnt
afrq <- codonFrequencies(acodcnts)

## global median RAAS over codons
araas <- sapply(names(afrq), function(cl)
    median(10^ctmt$RAAS[which(ctmt$aacodon==cl)]))

## NOTE: SCALED VALUES REPRODUCE NON-SCALED VERSION
plotdev(file.path(ctfig.path,paste0("codons_raas_freq_all")),
            type=ftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(afrq, log10(araas[names(afrq)]), 
        title=TRUE, cor.legend=FALSE, density=FALSE, pch=1,
        xlab="codon frequency", ylab=xl.raas)
figlabel("all tissues", pos='bottomleft')
dev.off()

for ( i in 1:nrow(tcodfreq) ) {

    tid <- tmap[tolower(rownames(tcodfreq)[i])]

    ## get codon frequency for this tissue
    tfrq <- tcodfreq[i,]

    ## get codon-specific RAAS for this tissue
    dtmt <- ctmt[ctmt$TMT.Tissue==tid,]
    if ( nrow(dtmt)==0 ) {
        cat(paste(tid, "NOT IN DATA\n"))
        next
    }
    traas <- sapply(names(tfrq), function(cl)
        median(10^dtmt$RAAS[dtmt$aacodon==cl]))

    plotdev(file.path(ctfig.path,paste0("codons_raas_freq_", tid)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tfrq, log10(traas), 
            title=TRUE, cor.legend=FALSE, density=FALSE, pch=1,
            xlab="codon frequency", ylab=xl.raas)
    figlabel(tid, pos='bottomleft')
    dev.off()
    plotdev(file.path(ctfig.path,paste0("codons_raas_freq_diff_", tid)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tfrq - afrq, log10(traas/araas),
            title=TRUE, cor.legend=FALSE, density=FALSE, pch=1,
            xlab=expression(f[tissue]-f[all]),
            ylab=expression(log[10](RAAS[tissue]/RAAS[all])))
    figlabel(tid, pos='bottomleft')
    dev.off()

    plotdev(file.path(ctfig.path,paste0("codons_freq_", tid)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tfrq, afrq, title=TRUE, cor.legend=FALSE, density=FALSE, pch=1,
            xlab=bquote(f[tissue]),
            ylab=expression(f[all]))
    figlabel(tid, pos='bottomleft')
    dev.off()
    plotdev(file.path(ctfig.path,paste0("codons_raas_", tid)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(log10(traas), log10(araas), 
            title=TRUE, cor.legend=FALSE, density=FALSE, pch=1,
            xlab=bquote(log[10](RAAS[tissue])),
            ylab=expression(log[10](RAAS[all])))
    figlabel(tid, pos='bottomleft')
    dev.off()
}
