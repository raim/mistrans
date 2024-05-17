
## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
source("~/work/mistrans/scripts/raasprofiles3_init.R")


cfig.path <- file.path(fig.path,"codons")
dir.create(cfig.path, showWarnings=FALSE)

### START ANALYSIS


### CODON FREQUENCIES vs. RAAS

## BY CODON->AA
## TODO:
## * CODON POSITION ANALYSIS.

## only use RAAS with assigned codons
ctmt <- tmtf[tmtf$codon!="",]

## add columns for codon positions
cpos <- strsplit(ctmt$codon,"")
ctmt$pos1 <- unlist(lapply(cpos, function(x) x[1]))
ctmt$pos2 <- unlist(lapply(cpos, function(x) x[2]))
ctmt$pos3 <- unlist(lapply(cpos, function(x) x[3]))

## LOAD GLOBAL GENE-WISE CODON COUNTS
codons <- read.delim(codon.file, row.names=1)

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


## Wu et al. 2019: Codon Stability Coefficient
csc <- read.csv(wu19.file, row.names=1)

## Gingold et al. 2014: tRNAs in proliferation vs. differentiation vs. other
## via trna file
trnas <- read.delim(trna.file)
## remove unclassified
## TODO: clean this mapping in the tRNA annotation script,
## and specifically also handle Met(e) vs. Met(i)
trnas <- trnas[trnas$CL.gingold14!="",]

## expand codons
tcodons <- strsplit(trnas$codons, ",")
tclasses <- trnas$CL.gingold14
tclasses[tclasses==""] <- "na."
tclasses <- rep(tclasses, lengths(tcodons))
tcodons <- unlist(tcodons)

## TODO: assign each codon by majority vote
codtab <- table(tcodons, tclasses)
ccls <- colnames(codtab)[apply(codtab,1,which.max)]
names(ccls) <- rownames(codtab)

ctmt$ccls <- ccls[ctmt$codon]

## global analysis: codon frequency vs. RAAS

## CODON COUNT and FREQUENCY IN ALL UNIQUE TRANSCRIPTS 

## codon counts in all mapped transcripts
## NOTE: this is used as our background frequency
cod  <- codons[unique(ctmt$transcript),]
codt <- apply(cod,2,sum) # total count

if ( any(!ctmt$transcript%in%rownames(codons)) )
    stop("some codon frequencies not found, rerun calculation")

### CODON SORTING

## PER AA
## transcript codon frequencies per AA and SORTING by frequency and AA prop
codl <- split(codt, GENETIC_CODE[sub(".*\\.","",names(codt))])
codl <- lapply(codl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
codl <- codl[aap.srt[aap.srt%in%names(codl)]] ## SORT BY AA PROP

## LOCAL CODON SORTING by background frequencies
## (AA property ->codon frequency)
codon.srt <- sub("\\.","-",names(unlist(codl)))
if ( !include.kr )
    codon.srt <- codon.srt[grep("[KR]", codon.srt, invert=TRUE)]

## codon frequencies per AA
Fbg <- lapply(codl, function(x) x/sum(x)) # codon frequency
Fbg <- unlist(Fbg)
names(Fbg) <- sub("\\.","-",names(Fbg))



## median codon RAAS:
## NOTE: median over all measurements
Craas <- sapply(COD.SRT, function(cl)
    log10(median(10^ctmt$RAAS[ctmt$aacodon==cl])))

## calculate unique SAAP/BP medians before
## taking codon median RAAS
site <- split(ctmt$RAAS, ctmt$unique.site)
site <- listProfile(site, y=ctmt$RAAS, use.test=use.test, min=3)
site.codons <- lapply(split(ctmt$aacodon,  ctmt$unique.site), unique)
if ( length(table(lengths(site.codons)))!=1 )
    stop("multiple codons per protein site")
site <- cbind(site, codon=unlist(site.codons))

## median of site medians
## NOTE: correlation to codon frequency is lost
## when using the median of medians.
Craas.site <- sapply(COD.SRT, function(cl)
    log10(median(10^site$median[site$codon==cl])))

## add columns for codon positions
cpos <- strsplit(sub("[A-Z]-","",site$codon),"")
site$pos1 <- unlist(lapply(cpos, function(x) x[1]))
site$pos2 <- unlist(lapply(cpos, function(x) x[2]))
site$pos3 <- unlist(lapply(cpos, function(x) x[3]))


## CODON FREQUENCIES at AAS
lcodt <- table(site$codon)
names(lcodt) <- sub(".*-","",names(lcodt))
## AAS codon frequencies per AA
lcodl <- split(lcodt, GENETIC_CODE[names(lcodt)])

lcodl <- lapply(lcodl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
lcodl <- lcodl[names(aaprop)] ## ##SORT BY AA PROP

## store for sanity check:
## why more codons via raasProfile unique counter than here?
lcodt.global <- lcodt
Fbg.lst <- codl
Faas.lst <- lcodl

## AA-specific codon frequency
Faas <- unlist(lapply(lcodl, function(x) x/sum(x)))
names(Faas) <- sub("\\.","-", names(Faas))

## inspect
if ( interactive() ) { # properly plotted below loop
    plotCor(Craas.site[names(Craas)], Craas)
    plotCor(Craas.site[names(Faas)], Faas)
    plotCor(Craas[names(Faas)], Faas)
    plotCor(Craas[names(Fbg)], Fbg)
    plotCor(Fbg[names(Faas)], Faas, xlim=c(0,1), ylim=c(0,1))
}

## codon rank - TODO: remove or improve this
codr <- lapply(codl, function(x) round(6*rank(x)/length(x))) # codon rank
codr$W[1] <- 0
codr$M[1] <- 0
## codon rank
codr <- unlist(codr)
names(codr) <- sub("\\.", "-", names(codr))
ctmt$codon.rank <- as.character(codr[ctmt$aacodon])
    
## codon frequency bins
ctmt$codon.bins <- cut(Fbg[ctmt$aacodon], seq(0,1,.1))

## TODO: do by 2|3, 4 or 6 codons
cod.bins <- cut(Fbg[ctmt$aacodon], seq(0,1,.1))
ncod <- unlist(lapply(CODL, length))



## Dana and Tuller 2014: decoding times
## TODO: analyze per AA - positive trend
## globally: negative trend
plotdev(file.path(cfig.path,paste0("codon_frequencies_decoding_dana14")),
        type=ftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, decod[sub(".*-","", names(Fbg)),1], density=FALSE,
        col=NA, xlab="codon frequency, all AAS transcripts",
        ylab=expression(decoding~rate(codons/s)))
points(Fbg, decod[sub(".*-","", names(Fbg)),1], lwd=1, cex=1,
       col=aa.cols[sub("-.*","",names(Fbg))],
       pch=aa.pchs[sub("-.*","",names(Fbg))])
dev.off()

plotdev(file.path(cfig.path,paste0("codon_frequencies_reldecoding_dana14")),
        type=ftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, decodl[sub(".*-","", names(Fbg)),1], density=FALSE,
        col=NA, xlab="codon frequency, all AAS transcripts",
        ylab=expression(relative~decoding~rate))
points(Fbg, decodl[sub(".*-","", names(Fbg)),1], lwd=1, cex=1,
       col=aa.cols[sub("-.*","",names(Fbg))],
       pch=aa.pchs[sub("-.*","",names(Fbg))])
dev.off()

## Wu et al. 2019: Codon Stability Coefficient
plotdev(file.path(cfig.path,paste0("codon_frequencies_CSC_wu19")),
        type=ftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, csc[sub(".*-","", names(Fbg)),"X293T_endo"], density=FALSE,
        col=NA, xlab="codon frequency, all AAS transcripts",
        ylab="Codon Stability Coefficient",
        legpos="bottomright")
points(Fbg, csc[sub(".*-","", names(Fbg)),"X293T_endo"], lwd=1, cex=1,
       col=aa.cols[sub("-.*","",names(Fbg))],
       pch=aa.pchs[sub("-.*","",names(Fbg))])
dev.off()





## PER CODON RAAS PROFILES per Dataset

## matrix to store codon frequencies in bg and at AAS
allcodons <- sub("\\.","-",names(unlist(codl)))
codon.results <- matrix(NA, nrow=length(COD.SRT), ncol=length(auds))
colnames(codon.results) <- auds
rownames(codon.results) <- COD.SRT
codon.Fbg <- codon.Faas <- codon.raas <- codon.results

    source("~/work/mistrans/scripts/saap_utils.R")
    source("~/programs/segmenTools/R/plotUtils.R")
    ##ftyp="pdf"

for ( ds in auds ) {

    ## TODO: add codon frequency two-sided bar plot

    dtmt <- ctmt
    if ( ds!="all" ) {
        if ( ds!="cancer" )
            dtmt <- ctmt[ctmt$Dataset==ds,]
        else
            dtmt <- ctmt[ctmt$Dataset!="Healthy",]
    }
    dsLAB <- ds
    if ( ds=="all" ) dsLAB <- ""

    ## NOTE: bg=FALSE, no column-wise background !
    dtmt$all <- "all"
    ova <- raasProfile(x=dtmt, id="unique.site", 
                       rows="all", cols="aacodon",
                       col.srt=codon.srt, filter=FALSE,
                       bg=FALSE, use.test=use.test, do.plots=FALSE, 
                       verb=0)

    ## NOTE: bg=TRUE :column-wise background
    ## TODO: move this outside loop and use rows as bg
    ovd <- raasProfile(x=dtmt, id="unique.site", 
                       rows="Dataset", cols="aacodon",
                       col.srt=codon.srt, filter=FALSE,
                       bg=TRUE, bg.dir="row",
                       use.test=use.test, do.plots=FALSE, 
                       verb=0)
    nr <- uds[uds%in%rownames(ovd$p.value)]
    if ( length(nr)>1 )
        ovd <- sortOverlaps(ovd, axis=2, srt=nr)



    ## CODON FREQUENCY ANALYSIS - by "to" amino acid

    ## NOTE: bg=FALSE, no column-wise background !
    ovw <- raasProfile(x=dtmt, id="unique.site", 
                       rows="to", cols="aacodon",
                       col.srt=codon.srt, filter=FALSE,
                       bg=FALSE, use.test=use.test, do.plots=FALSE, 
                       verb=0)

    ## TODO: why do the frequencies derived from ovw$unique differ from
    ## those calculated above?

    
    ## sort by amino acid property via shapely colors!
    ##aasrt <- names(aaprop)[names(aaprop)%in%rownames(ovw$p.value)]
    aasrt <- aap.srt[aap.srt%in%rownames(ovw$p.value)]
    ovw <- sortOverlaps(ovw, axis=2, srt=aasrt, cut=TRUE)


    ## CUSTOM CODON FREQUENCY PLOTS

    ## CODON COUNTS AT THE AAS OF UNIQUE BP
    ## NOTE: this was wrong, since the same codon
    ## with different substitutions is counted several times
    ## lcodt <- apply(ovw$unique, 2, sum)
    ## names(lcodt) <- sub(".*-","",names(lcodt))


    dsite <- split(dtmt$RAAS, dtmt$unique.site)
    dsite <- listProfile(dsite, y=dtmt$RAAS, use.test=use.test, min=3)
    dsite.codons <- lapply(split(dtmt$aacodon,  dtmt$unique.site), unique)
    if ( length(table(lengths(dsite.codons)))!=1 )
        stop("multiple codons per protein site")
    dsite <- cbind(dsite, codon=unlist(dsite.codons))
    
    ## codon frequencies at AAS
    lcodt <- table(dsite$codon)
    names(lcodt) <- sub(".*-","",names(lcodt))
  
    ## AAS codon frequencies per AA
    lcodl <- split(lcodt, GENETIC_CODE[names(lcodt)])
    
    lcodl <- lapply(lcodl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
    lcodl <- lcodl[aap.srt] ## ##SORT BY AA PROP

    ## AA-specific codon frequency
    Faas.ds <- unlist(lapply(lcodl, function(x) x/sum(x)))
    names(Faas.ds) <- sub("\\.","-", names(Faas.ds))

    ## CODON COUNTS IN ALL UNIQUE TRANSCRIPTS
    ## DATASET-SPECIFIC
    cod  <- codons[unique(dtmt$transcript),]
    codt <- apply(cod, 2, sum)

    ## transcript codon frequencies per AA
    codl <- split(codt, GENETIC_CODE[sub(".*\\.","",names(codt))])
    codl <- lapply(codl, sort, decreasing=TRUE) # SORT BY MOST FREQUENT
    codl <- codl[aap.srt[aap.srt%in%names(codl)]] # SORT BY AA PROP
    Fbg.ds <- unlist(lapply(codl, function(x) x/sum(x)))
    names(Fbg.ds) <- sub("\\.","-", names(Fbg.ds))

    ## why is global codon count at AAS lower than that via raasProfile?
    if ( ds=="all" ) {
        Fbg.ds.lst <- codl
        Faas.ds.lst <- lcodl
        tmp <- cbind(unlist(lcodt.global[names(lcodt)]), lcodt)
    }

    ## MEDIAN CODON RAAS
    ## get median RAAS for each codon and plot against codon frequency
    Craas.ds <- sapply(codon.srt, function(cl)
        log10(median(10^dtmt$RAAS[dtmt$aacodon==cl])))

    codon.Fbg[names(Fbg.ds),ds] <- Fbg.ds
    codon.Faas[names(Faas.ds),ds] <- Faas.ds
    codon.raas[names(Craas.ds), ds] <- Craas.ds

    
    ## CODON SORTING and PLOT SETTINGS

    ## TODO: simplify sorting and do globally, same for all
    
    ## sort by background (AA property ->codon frequency)
    ## i.e. as codl
    cdsrt <- sub("\\.","-",names(unlist(codl)))

    ## FILTER AVAILABLE
    cdsrt <- cdsrt[cdsrt%in%colnames(ovw$p.value)]
    
    ## codon frequency bins
    codon.bins <- cut(Faas.ds[dtmt$aacodon], seq(0,1,.2))

    ## codon class lines in plots
    cdn.cnt <- table(sub("-.*","",colnames(ovw$p.value)))[aap.srt]
    cdn.cnt <- cumsum(cdn.cnt[!is.na(cdn.cnt)])
    
    ## SORT CODON FREQUENCIES as in main plots
    ## TODO: constant sorting of  ALL codons 
    Faas.ds <- Faas.ds[cdsrt]
    Fbg.ds <- Fbg.ds[cdsrt]
    Craas.ds <- Craas.ds[cdsrt]

    ## log2 ratio of frequencies
    Fratio <- (Faas.ds/Fbg.ds) # log2 ratio of codon frequencies

    ## axis limit for ratio
    rlm <- max(abs(Fratio),na.rm=TRUE)

    
    plotProfiles(ova, fname=file.path(cfig.path,paste0("codon_",
                                                      SETID,"_",ds,
                                                      "_all")),
                 fw=.2, mai=c(.05,.75,.05,.6), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=dsLAB, ftyp=ftyp, mtxt="",
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols, col.lines=cdn.cnt, ffam="monospace")


    plotProfiles(ovd, fname=file.path(cfig.path,paste0("codon_",
                                                      SETID,"_",ds,"_Dataset")),
                 fw=.2, mai=c(.05,.75,.05,.6), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=dsLAB, ftyp=ftyp, mtxt="",
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols, col.lines=cdn.cnt, ffam="monospace")
    maimain <- c(.6,.75,.6,.6) # re-used for codon axis
    plotProfiles(ovw, fname=file.path(cfig.path,paste0("codon_",SETID,"_",ds)),
                 fw=.2, mai=maimain, ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=dsLAB, ftyp=ftyp,
                 mtxt="Incorporated AA", mtxt.cex=1.5,
                 vcols=vcols, vbrks=vbrks,
                 axis1.col=aap.cols[sub("-.*","",colnames(ovw$p.value))],
                 axis2.col=aap.cols[rownames(ovw$p.value)],
                 gcols=gcols, col.lines=cdn.cnt, ffam="monospace", verb=1)
        
### CODON FREQUENCY ANALYSIS

   ## sort by amino acid property via shapely colors!
    ##aasrt <- names(aaprop)[names(aaprop)%in%rownames(ovw$p.value)]
    aasrt <-
        aap.srt[aap.srt%in%sub("-.*","",names(Fbg.ds))]


    ## CODON FREQUENCY CORRELATIONS
    ## relation codon base frequency and ratio AAS
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_ratio")),
           type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(Fbg.ds, Fratio, density=FALSE,
            ylab=expression(f[AAS]/f[bg]),
            xlab=expression(f[bg]), xlim=c(0,1),
            col=NA,pch=19, lwd=2, cex=.6, line.methods="ols")
    points(Fbg.ds, Fratio, lwd=2, cex=1,
           col=aa.cols[sub("-.*","",names(Fbg.ds))],
           pch=aa.pchs[sub("-.*","",names(Fbg.ds))])
    figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    ## higher frequency tends to be above diagonal
    ## TODO: better AA colors
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(Fbg.ds, Faas.ds, density=FALSE,
            ylab=expression(codon~frequency~f[AAS]),
            xlab=expression(codon~frequency~f[bg]), xlim=c(0,1), ylim=c(0,1),
            col=NA,pch=19, lwd=2, cex=.6, line.methods="ols")
    points(Fbg.ds, Faas.ds, lwd=2, cex=1,
           col=aa.cols[sub("-.*","",names(Fbg.ds))],
           pch=aa.pchs[sub("-.*","",names(Fbg.ds))])
    ##abline(a=0,b=1, col=2)
    legend("bottomright",aasrt,
           col=aa.cols[aasrt],
           pch=aa.pchs[aasrt],
           pt.cex=1, lwd=2, lty=NA, seg.len=0, ncol=2, cex=.8, y.intersp=.75,
           bty="n")
    figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_diff")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(Fbg.ds, Faas.ds-Fbg.ds, density=FALSE,
         ylab=expression(f[AAS]-f[bg]),
         xlab=expression(f[bg]), xlim=c(0,1),
         col=NA,pch=19, lwd=2, cex=.6, line.methods="ols")
    points(Fbg.ds, Faas.ds-Fbg.ds, lwd=2, cex=1,
           col=aa.cols[sub("-.*","",names(Fbg.ds))],
           pch=aa.pchs[sub("-.*","",names(Fbg.ds))])
    if ( FALSE )
        legend("bottomright",aasrt,
               col=aa.cols[aasrt],
               pch=aa.pchs[aasrt],
               pt.cex=1, lwd=2, lty=NA, ncol=2, cex=.6, y.intersp=.75,
               bty="n")
    figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_raas")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(Craas.ds, Faas.ds, ylim=c(0,1),
            density=FALSE, xlab=xl.raas,
            ylab=expression(codon~frequency~f[AAS]), col=NA, line.methods="ols")
    points(Craas.ds, Faas.ds, lwd=2, cex=1, 
           col=aa.cols[sub("-.*","",names(Fbg.ds))],
           pch=aa.pchs[sub("-.*","",names(Fbg.ds))])
    figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()
    
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_raas_nocol")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(Craas.ds, Faas.ds, ylim=c(0,1),
            density=FALSE, xlab=xl.raas,
            ylab=expression(codon~frequency~f[AAS]), col="#00000099", pch=19,
            lwd=0,cex=1) #, line.methods="")
    ##figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
    ##figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_raas_slopes")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(Craas.ds, Faas.ds, xlab=xl.raas, ylim=c(0,1),
         ylab=expression(codon~frequency~f[AAS]), col=NA)
    slps <- rep(NA, length(CODL))
    for ( k in 1:length(CODL) ) {
        aa <- names(CODL)[k]
        cds <- paste0(aa,"-",CODL[[k]])
        if ( sum(cds%in%names(Faas.ds))<2 ) next
        rm <- Craas.ds[cds]
        cd <- Faas.ds[cds]

        cd <- cd[!is.na(rm)]
        rm <- rm[!is.na(rm)]
        if ( length(rm)<1 ) next
  
        points(rm, cd, lwd=2, cex=1,
               col=aa.cols[aa], pch=aa.pchs[aa])
        if ( length(rm)<2 ) next
        cdf <- lm(cd ~rm)
        slps[k] <- coef(cdf)[2]
        lines(rm, predict(cdf), col=aa.cols[aa])
    }
    dev.off()
    
    ## plot slopes
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequencies_raas_slopes_bynum")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(unlist(lapply(CODL,length)), slps, col=aa.cols[names(CODL)],
         pch=aa.pchs[names(CODL)], lwd=2,
         xlab="number of codons per/AA")
    abline(h=0)
    dev.off()
    plotdev(file.path(cfig.path,
                      paste0("codon_",SETID,"_",ds,
                             "_codons_frequencies_raas_slopes_bynum_zoom")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(unlist(lapply(CODL,length)), slps, col=aa.cols[names(CODL)],
         pch=aa.pchs[names(CODL)], lwd=2, ylim=c(-1.5,0),
         xlab="number of codons per/AA")
    abline(h=0)
    dev.off()
    
    ## TODO: why doesn't mean correlate to RAAS, while median does!
    if ( interactive() ) {
        fds <- Faas.ds[names(Faas.ds)%in%colnames(ova$median)]
        plotCor(ova$median[,names(fds)], fds, density=FALSE, col=NA)
        points(ova$median[,names(fds)], fds, lwd=1, cex=1,
               col=aa.cols[sub("-.*","",names(fds))],
               pch=aa.pchs[sub("-.*","",names(fds))])
        plotCor(ova$mean[,names(fds)], fds, density=FALSE, col=NA)
        points(ova$mean[,names(fds)], fds, lwd=1, cex=1,
               col=aa.cols[sub("-.*","",names(fds))],
               pch=aa.pchs[sub("-.*","",names(fds))])
        plotCor(ova$mean[,names(fds)], ova$median[,names(fds)],
                density=FALSE, col=NA)
        points(ova$mean[,names(fds)],
               ova$median[,names(fds)], lwd=1, cex=1,
               col=aa.cols[sub("-.*","",names(fds))],
               pch=aa.pchs[sub("-.*","",names(fds))])
    }

    ## COMPARE TO EXTERNAL MEASURES

    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                       "_codons_frequencies_raas_wu19")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(Craas.ds, csc[sub(".*-","",names(Craas.ds)),"X293T_endo"],
            density=FALSE, xlab=xl.raas,
            ylab="Codon Stability Coefficient", col=NA)
    points(Craas.ds, csc[sub(".*-","",names(Craas.ds)),"X293T_endo"],
           lwd=2, cex=1, col=aa.cols[sub("-.*","",names(Craas.ds))],
           pch=aa.pchs[sub("-.*","",names(Craas.ds))])
    figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
        figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()

    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                       "_codons_frequencies_raas_dana14")),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(Craas.ds, decod[sub(".*-","",names(Craas.ds)),1],
            density=FALSE, xlab=xl.raas,
            ylab=expression(decoding~rate/(codons/s)), col=NA)
    points(Craas.ds, decod[sub(".*-","",names(Craas.ds)),1], lwd=2, cex=1,
           col=aa.cols[sub("-.*","",names(Craas.ds))],
           pch=aa.pchs[sub("-.*","",names(Craas.ds))])
    figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
        figlabel(LAB, pos="bottomright", cex=.7)
    for ( ax in 3:4 )  axis(ax, labels=FALSE)
    dev.off()


    ## full codon dotplots!
    mai <- c(.05,.75,.1,.6)
    fw <- .2
    nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]

    source("~/work/mistrans/scripts/saap_utils.R")
    maix <- maimain
    maix[1] <- maix[3] <- 0.5
    nhm <- nrow(ovw$p.value) *fw + maix[1] + maix[3]
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_axis")),
            type=ftyp, res=300, width=1, height=nhm)
    par(mai=maimain, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotCodonAxis(CODP[rownames(ovw$p.value)], aa.cols=aap.cols, cex=2)
    dev.off()
    
     ## codon RAAS aligned with codon dotplot 
    ramm <- matrix(Craas.ds,nrow=1)
    colnames(ramm) <- names(Craas.ds)
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_raas")),
            type=ftyp, res=300, width=nw, height=.25+.1)
    par(mai=c(.05,.75,.05,.6), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    image_matrix(ramm, col=vcols, breaks=vbrks,
                 axis=NA, las=2, ylab=NA)
    mtext(expression(bar(RAAS)), 2, .15, las=2, cex=1)
    box()
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    abline(v=.5+cdn.cnt, lwd=1, xpd=FALSE, col="white")
    dev.off()
    
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_frequency_raas")),
            type=ftyp, res=300, width=nw,height=2.9)
    par(mai=c(1,.75,.05,.6), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i",
        family="monospace")
    tmp <- boxplot(dtmt$RAAS ~ factor(dtmt$aacodon, levels=cdsrt),
                   las=2, xlab="", ylab=xl.raas, axes=FALSE,
                   cex=.5, pch=19, pars=list(outcol="#00000055"))
    ##axis(3, at=1:length(cdsrt), labels=round(Faas.ds[cdsrt],2), las=2)
    axis(4); axis(2)
    Map(axis, 1, at=1:length(tmp$names), labels=tmp$names, las=2,
        col.axis=aap.cols[sub("-.*","",tmp$names)], col=NA, font=2)
    abline(v=.5+cdn.cnt, lwd=1, xpd=FALSE)
    axis(1, at=c(0,.5+cdn.cnt), labels=FALSE, tcl=-.5)
    figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    text(c(3.25,13,25,42),rep(-9.5, 5), cex=1.2, font=2,
         labels=aaprop.srt, col=aaprop.cols[aaprop.srt], xpd=TRUE)
    mtext("mRNA codons at amino acid substitution sites", 1, 4, cex=1.5)
    dev.off()

    
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_ratio_old")),
            type=ftyp, res=300, width=nw,height=.75)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plot(1:length(Faas.ds), Fratio, type="h", axes=FALSE, xlab=NA, lwd=3,
         ylab=expression(lg[2]~r),
         xlim=c(0.5,length(Faas.ds)+.5), col=2, ylim=c(-rlm,rlm))
    abline(h=0)
    for ( ax in c(2,4) ) {
        axis(ax, at=seq(-5,5,.5), labels=FALSE)
        axis(ax, at=c(-.5,.5), cex.axis=1, las=2)
    }
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    dev.off()
    
    mai.bar <- mai
    mai.bar[c(2,4)] <- mai.bar[c(2,4)] +.05
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,"_codons_ratio")),
            type=ftyp, res=300, width=nw,height=.75)
    par(mai=mai.bar, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    bp <- barplot(Fratio, axes=FALSE, xlab=NA, beside=TRUE,
                  ylab=expression(lg[2]~r), las=2, xaxt='n', 
                  width=.5, space=.5, ylim=c(-rlm,rlm))
    ##abline(v=bp[2,cdn.cn]+  (bp[1,cdn.cn+1]-bp[2,cdn.cn])/2)
    cdn.cn <- head(cdn.cnt, length(cdn.cnt)-1)
    abline(v=bp[cdn.cn,]+unique(diff(bp[,1]))/2)
    ##abline(v=par("usr")[1]+unique(diff(bp)[,1])*cdn.cnt, lwd=1, xpd=TRUE)
    ##x2 <- par("usr")[2]-par("usr")[1]
    ##rel.cnt <- cdn.cnt/max(cdn.cnt)
    ##rel.cnt <- head(rel.cnt, length(rel.cnt)-1)
    ##abline(v=par("usr")[1]+x2*rel.cnt, lwd=1, xpd=TRUE)
    abline(h=0)
    for ( ax in c(2,4) ) {
        axis(ax, at=seq(-5,5,.1), labels=FALSE, tcl=par("tcl")/2, line=.2)
        axis(ax, at=seq(-5,5,.5), labels=FALSE, line=.2)
        axis(ax, at=c(-1,-.5,0,.5,1),
             label=c("-1","0.5","0","0.5","1"), cex.axis=1, las=2, line=.2)
    }
    dev.off()
    
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",
                                       ds,"_codons_barplot")),
            type=ftyp, res=300, width=nw, height=.75)
    par(mai=mai.bar, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    bp  <- barplot(rbind(Faas.ds,Fbg.ds[names(Faas.ds)]),
                   axes=FALSE, xlab=NA, beside=TRUE,
                   ylab="", las=2, xaxt='n')
    cdn.cn <- head(cdn.cnt, length(cdn.cnt)-1)
    abline(v=bp[2,cdn.cn]+  (bp[1,cdn.cn+1]-bp[2,cdn.cn])/2)
    ##x2 <- par("usr")[2]
    ##rel.cnt <- cdn.cnt/max(cdn.cnt)
    ##rel.cnt <- head(rel.cnt, length(rel.cnt)-1)
    ##abline(v=.5+x2*rel.cnt, lwd=1, xpd=TRUE)
    axis(2, las=2, line=.2)
    axis(4, las=2, line=.2)
    dev.off()


    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",
                                       ds,"_codons_barplot_onlyfbg")),
            type=ftyp, res=300, width=nw, height=1)
    par(mai=mai.bar, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    bp  <- barplot(rbind(Fbg.ds[names(Faas.ds)]),
                   axes=FALSE, xlab=NA, beside=TRUE,
                   ylab="", las=2, xaxt='n')
    cdn.cn <- head(cdn.cnt, length(cdn.cnt)-1)
    abline(v=bp[,cdn.cn]+unique(diff(bp[1,]))/2)
   
    ##x2 <- par("usr")[2]
    ##rel.cnt <- cdn.cnt/max(cdn.cnt)
    ##rel.cnt <- head(rel.cnt, length(rel.cnt)-1)
    ##abline(v=.5+x2*rel.cnt, lwd=1, xpd=TRUE)
    for ( ax in c(2,4) ) {
        axis(ax, at=c(0,.5,1), las=2, line=.2)
        axis(ax, at=seq(0,1,.1), labels=FALSE, line=.2, tcl=par("tcl")/2)
    }
    mtext("codon\nfrequency", 2, 2)
    dev.off()
    

    ## hypergeo tests
    aaa <- unique( GENETIC_CODE)
    aam <- matrix(1, ncol=length(lcodt), nrow=1)
    colnames(aam) <- names(lcodt)
    aac <- aap <- aam
    aac[] <- 0
    
    for ( aa in unique(aaa) ) {
        alc <- names(GENETIC_CODE)[GENETIC_CODE==aa]
        
        if ( length(alc)<1 ) next
        for ( j in  seq_along(alc) ) {
            
            wcd <- alc[ j] # white balls
            bcd <- alc[-j] # black balls
            
            if ( !wcd%in%names(lcodt) ) next
            
            m <- codt[wcd]; # number of white balls
            n <- sum(codt[bcd], na.rm=TRUE); # number of black balls
            q <- lcodt[wcd] # white balls drawn - CODON OF INTEREST
            k <- sum(lcodt[alc], na.rm=TRUE) # number of balls drawn
            
            aac[,wcd] <- q
            aam[,wcd] <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)
            aap[,wcd] <- phyper(q=q, m=m, n=n, k=k, lower.tail=TRUE)
            ##cat(paste("testing", aa, wcd, signif(aam[,wcd]), "\n"))
        }
    }
    
    colnames(aam) <- colnames(aap) <- colnames(aac) <-
        paste0(GENETIC_CODE[colnames(aam)],"-", colnames(aam))

    ## sort
    csrt <- codon.srt[codon.srt%in%colnames(aam)]
    aam <- aam[,csrt,drop=FALSE]
    aap <- aap[,csrt,drop=FALSE]
    aac <- aac[,csrt,drop=FALSE]

    ## construct overlap class
    ovl <- list()
    pvl <- aam
    pvl[aap<aam] <- aap[aap<aam]
    ovl$p.value <- pvl
    ovl$count <- aac
    ovl$sign <- ifelse(aam<aap, 1, -1)
    
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_hypergeo")),
            type=ftyp, res=300, width=nw, height=.25+.1)
    par(mai=c(.05,.75,.05,.6), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, axis=NA, xlab=NA, ylab=NA,
                 text.cex=.7, col=ttcols)
    ##axis(1, at=1:length(aac), labels=aac, las=2)
    mtext(expression(f), 2, .25, las=2)
    box()
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    dev.off()

    ## two roes
    ovl$p.value <- rbind(aam,aap)
    rownames(ovl$p.value) <- c("more","less")
    ovl$count <- rbind(aac,aac)
    ovl$count[1,aam>=aap] <- 0
    ovl$count[2,aam< aap] <- 0 
    ovl$sign <- rbind(rep( 1,length(aam)),
                      rep(-1,length(aam)))
    plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                      "_codons_hypergeo2")),
            type=ftyp, res=300, width=nw, height=.35+.1)
    par(mai=c(.05,.75,.05,.6), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, axis=2, xlab=NA, ylab=NA,
                 text.cex=.7, col=ttcols)
    ##axis(1, at=1:length(aac), labels=aac, las=2)
    ##axis(1, at=1:length(aac), labels=aac, las=2)
    ##mtext(expression(f), 2, .25, las=2)
    box()
    abline(v=.5+cdn.cnt, lwd=1, xpd=TRUE)
    dev.off()
    
}

## CODON FREQUENCY ANALYSIS

## sanity check
## data set 'all' values should correspond to
## the values calculated above the loop
if ( interactive() ) {
    plotCor(codon.raas[names(Craas),"all"], Craas)
    plotCor(codon.Fbg[names(Fbg),"all"], Fbg)
    plotCor(codon.Faas[names(Faas),"all"], Faas)
}

## RESULT TABLE
codon.freq <- data.frame(codon=sub(".*-","",names(Fbg)),
                         AA=sub("-.*","",names(Fbg)),
                         frequency_all=Fbg,
                         frequency_AAS=Faas[names(Fbg)],
                         RAAS=Craas[names(Fbg)])
write.table(codon.freq, file=file.path(cfig.path, "codon_frequencies.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE, na="")

fftyp <- ftyp # "pdf" # 

## exclude the single codon AA (M and W) for correlation analysis
fbg <- Fbg
fbg[fbg==1] <- NA

plotdev(file.path(cfig.path,paste0("codons_faas_both_correct")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(codon.Faas[names(Faas),"all"], Faas, round=4,
        xlab=expression(f[AAS]~"NOT via raasProfile function"),
        ylab=expression(f[AAS]~"from unique protein AAS sites"),
        density=FALSE, pch=19, col="#00000099")
dev.off()


## TODO:
## get p-values for TLS - https://stats.stackexchange.com/questions/166483/p-values-for-regression-coefficients-in-total-least-squares-regression

plotdev(file.path(cfig.path,paste0("codons_faas_fbg")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
plotCor(Fbg, Faas[names(Fbg)], ##outliers=which(Fbg==1),
        xlab=expression("Relative codon frequency,"~f[bg]),
        xlim=c(0,1), ylim=c(0,1),
        ylab=expression("Relative codon frequency,"~f[AAS]),
        line.methods=c("ols"), density=FALSE, pch=19, col="#00000099")
#abline(a=0, b=1, col="gray")
dev.off()

plotdev(file.path(cfig.path,paste0("codons_faas_fbg_ratio_woMW")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, Faas[names(Fbg)]/Fbg, outliers=which(Fbg==1),
        legpos="bottomright",
        xlab=expression("Relative codon frequency"~f[bg]),
        xlim=c(0,1), ##ylim=c(0,1),
        ylab=expression(f[AAS]/f[bg]),
        line.methods=c("ols","tls"), density=FALSE, pch=19, col="#00000099")
dev.off()
plotdev(file.path(cfig.path,paste0("codons_faas_fbg_ratio_withMW")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, Faas[names(Fbg)]/Fbg, ##outliers=which(Fbg==1),
        legpos="bottomright",
        xlab=expression("Relative codon frequency"~f[bg]),
        xlim=c(0,1), ##ylim=c(0,1),
        ylab=expression(f[AAS]/f[bg]),
        line.methods=c("ols","tls"), density=FALSE, pch=19, col="#00000099")
dev.off()

plotdev(file.path(cfig.path,paste0("codons_raas_fbg_woMW")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, Craas[names(Fbg)], xlim=c(0,1), outliers=which(Fbg==1),
        line.methods=c("tls","ols"), density=FALSE, pch=19,
        col="#00000099",
        xlab=expression("Relative codon frequency"~f[bg]), ylab=xl.raas)
dev.off()

plotdev(file.path(cfig.path,paste0("codons_raas_fbg_withMW")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, Craas[names(Fbg)], xlim=c(0,1), ##outliers=which(Fbg==1),
        line.methods=c("tls","ols"), density=FALSE, pch=19, col="#00000099",
        xlab=expression("Relative codon frequency"~f[bg]), ylab=xl.raas)
dev.off()

## NOTE:
## * do not exclude M and W,
## * show fbg instead of fAAS vs. RAAS,
## * ignore fAAS > fbg.

## official publication figures
## no lines, to tls; colors for supplement

plotdev(file.path(cfig.path,paste0("codons_raas_fbg")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
plotCor(Fbg, Craas[names(Fbg)], xlim=c(0,1), ##outliers=which(Fbg==1),
        line.methods=c("ols"), density=FALSE, pch=19,
        col="#00000099",
        xlab=expression("Relative codon frequency"~f[bg]), ylab=NA)
mtext(xl.raas, 2, 1.2)
dev.off()

## Wu et al. 2019: Codon Stability Coefficient
plotdev(file.path(cfig.path,paste0("codons_raas_wu19")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
plotCor(csc[sub(".*-","",names(Craas)),"X293T_endo"], Craas, 
        line.methods=c("ols"), density=FALSE, ylab=NA,
        xlab=NA, col=NA)
mtext("Codon Stability Coefficient", 1, 1.2)
points(csc[sub(".*-","",names(Craas)),"X293T_endo"], Craas, 
       lwd=NA, cex=1, pch=19, col="#00000099")
mtext(xl.raas, 2, 1.2)
dev.off()

## Dana and Tuller 2014: Decoding Rate
plotdev(file.path(cfig.path,paste0("codons_raas_dana14")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
pc <- plotCor(decod[sub(".*-","",names(Craas)),1], Craas, 
        line.methods=c("ols"), density=FALSE, ylab=NA,
        xlab=expression(decoding~rate/(codons/s)), col=NA,
        legpos="bottomright")
mtext(xl.raas, 2, 1.2)
points(decod[sub(".*-","",names(Craas)),1], Craas, 
       lwd=NA, cex=1, pch=19, col="#00000099")
dev.off()
## Dana and Tuller 2014: Relative Decoding Rate
plotdev(file.path(cfig.path,paste0("codons_raas_dana14_rel")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
pc <- plotCor(decodl[sub(".*-","",names(Craas)),1], Craas, 
        line.methods=c("ols"), density=FALSE, ylab=NA,
        xlab=expression(scaled~decoding~rate), col=NA,
        legpos="bottomright")
mtext(xl.raas, 2, 1.2)
points(decodl[sub(".*-","",names(Craas)),1], Craas, 
       lwd=NA, cex=1, pch=19, col="#00000099")
dev.off()

## encoded amino acid for coloring
aatmp <- sub("-.*", "", names(Fbg))
aasort <- names(sort(aa.cols))
aasort <- aasort[aasort%in%aatmp]

plotdev(file.path(cfig.path,paste0("codons_raas_fbg_byAA")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, Craas[names(Fbg)], xlim=c(0,1),
        line.methods=c("tls","ols"), density=FALSE, pch=19,
        xlab=expression("Relative codon frequency"~f[bg]),
        ylab=xl.raas, col=NA)
points(Fbg, Craas[names(Fbg)], col=aa.cols[aatmp], pch=aa.pchs[aatmp], lwd=2)
dev.off()


plotdev(file.path(cfig.path,paste0("codons_faas_fbg_ratio_byAA")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, Faas[names(Fbg)]/Fbg,
        legpos="bottomright",
        xlab=expression("Relative codon frequency, MTP"~f[bg]),
        xlim=c(0,1), col=NA,
        ylab=expression(f[AAS]/f[bg]),
        line.methods=c("ols","tls"), density=FALSE, pch=19)
points(Fbg, Faas[names(Fbg)]/Fbg, col=aa.cols[aatmp], pch=aa.pchs[aatmp])
dev.off()


plotdev(file.path(cfig.path,paste0("codons_faas_fbg_byAA")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, Faas[names(Fbg)], 
        xlab=expression("Relative codon frequency, MTP"~f[bg]),
        xlim=c(0,1), ylim=c(0,1), col=NA,
        ylab=expression("Relative codon frequency, AAS"~f[AAS]),
        line.methods=c("ols","tls"), density=FALSE, pch=19)
##abline(a=0, b=1, col="gray")
points(Fbg, Faas[names(Fbg)], col=aa.cols[aatmp], pch=aa.pchs[aatmp])
legend("bottomright",aasort,
       col=aa.cols[aasort],
       pch=aa.pchs[aasort],
       pt.cex=1, lwd=1, lty=NA, seg.len=0.1, ncol=2, cex=.8, y.intersp=.8,
       bty="n")
dev.off()


plotdev(file.path(cfig.path,paste0("codons_raas_wu19_byAA")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(csc[sub(".*-","",names(Craas)),"X293T_endo"], Craas, 
        line.methods=c("ols","tls"), density=FALSE, ylab=xl.raas,
        xlab="Codon Stability Coefficient", col=NA)
points(csc[sub(".*-","",names(Craas)),"X293T_endo"], Craas, 
       lwd=2, cex=1, col=aa.cols[sub("-.*","",names(Craas))],
       pch=aa.pchs[sub("-.*","",names(Craas))])
dev.off()

plotdev(file.path(cfig.path,paste0("codons_raas_dana14_byAA")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(decod[sub(".*-","",names(Craas)),1], Craas, 
        density=FALSE, ylab=xl.raas,
        xlab=expression(decoding~rate/(codons/s)), col=NA,
        legpos="bottomright")
points(decod[sub(".*-","",names(Craas)),1], Craas, 
       lwd=2, cex=1,
       col=aa.cols[sub("-.*","",names(Craas))],
       pch=aa.pchs[sub("-.*","",names(Craas))])
dev.off()

plotdev(file.path(cfig.path,paste0("codons_raas_dana14_rel_byAA")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(decodl[sub(".*-","",names(Craas)),1], Craas, 
        density=FALSE, ylab=xl.raas,
        xlab=expression(scaled~decoding~rate), col=NA,
        legpos="bottomright")
points(decodl[sub(".*-","",names(Craas)),1], Craas, 
       lwd=2, cex=1,
       col=aa.cols[sub("-.*","",names(Craas))],
       pch=aa.pchs[sub("-.*","",names(Craas))])
dev.off()


### NOT USED: remove or archive?

## CODON BINS by Dataset
ovw <- raasProfile(x=ctmt, id="SAAP", 
                   rows="codon.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   row.srt=rev(levels(ctmt$codon.bins)),
                   use.test=use.test, do.plots=FALSE,
                   fname=file.path(dpath,paste0("codon_bins_",SETID,"_")),
                   xlab=xl.raas,
                   verb=1)
plotProfiles(ovw, fname=file.path(cfig.path,paste0("codon_bins_",SETID)),
             mai=c(.8,1,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             mtxt="codon frequency", mtxt.line=3.5,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="", ftyp=ftyp,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## CODON FREQUENCY RANK by Dataset
ovw <- raasProfile(x=ctmt, id="SAAP", 
                   rows="codon.rank", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   ##row.srt=rev(levels(ctmt$codon.bins)),
                   use.test=use.test, do.plots=FALSE,
                   fname=file.path(dpath,paste0("codon_rank_",SETID,"_")),
                   xlab=xl.raas,
                   verb=1)
plotProfiles(ovw, fname=file.path(cfig.path,paste0("codon_rank_",SETID)),
             mai=c(.8,.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             mtxt="codon rank",# mtxt.line=3.5,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="", ftyp=ftyp,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

### CODON POSITION
ovd <- raasProfile(x=ctmt, id="unique.site", 
                   rows="pos2", cols="Dataset",
                   col.srt=uds, filter=FALSE,
                   bg=TRUE, #bg.dir="row",
                   use.test=use.test, do.plots=FALSE, 
                   verb=0)
par(mai=c(1,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovd, value="median", vcols=vcols, vbrks=vbrks, axis=1:2,
           p.dot=1e-10, dot.sze=c(.3,2), show.total=TRUE)

boxplot(ctmt$RAAS ~ ctmt$pos1)
boxplot(ctmt$RAAS ~ ctmt$pos2)
boxplot(ctmt$RAAS ~ ctmt$pos3)

## per unique site
boxplot(site$median ~ site$pos1)
boxplot(site$median ~ site$pos2)
boxplot(site$median ~ site$pos3)

### CODON CLASS
ovd <- raasProfile(x=ctmt, id="unique.site", 
                   rows="ccls", cols="Dataset",
                   col.srt=uds, filter=FALSE,
                   bg=TRUE, #bg.dir="row",
                   use.test=use.test, do.plots=FALSE, 
                   verb=0)
par(mai=c(1,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovd, value="median", vcols=acols, vbrks=abrks, axis=1:2,
           p.dot=1e-10, dot.sze=c(.3,2), show.total=TRUE)

ctmt$CCLS <-
    toupper(sub("differentiation","diff.",
                sub("proliferation", "proli.", ctmt$ccls)))

for ( ds in auds ) {
    for ( type in c("all",unique(ctmt$from)) ) {
        dtmt <- ctmt
        dsl <- ""
        if ( ds!="all" ) {
            dsl <- ds
            if ( ds!="cancer" )
                dtmt <- ctmt[ctmt$Dataset==ds,]
            else
                dtmt <- ctmt[ctmt$Dataset!="Healthy",]
        }

        
        if ( type!="all" )
            dtmt <- dtmt[dtmt$from==type,]
        if ( nrow(dtmt)<2 ) next

        
        ## comparison types
        cnt <- sort(table(dtmt$CCLS))
        ts <- names(cnt)
        if ( length(ts)==3 ) {
            d1 <- "DIFF."
            d2 <- "PROLI."
        } else {
            d1 <- ts[1]
            d2 <- ts[2]
        }
        
        tt <- NULL
        if ( sum(dtmt$CCLS==d1, na.rm=TRUE)>1 &
             sum(dtmt$CCLS==d2, na.rm=TRUE)>1 ) {
            tt <- t.test(dtmt$RAAS[dtmt$CCLS==d1],
                         dtmt$RAAS[dtmt$CCLS==d2])
        } else next
        if ( tt$p.value >0.05 ) next

        cat(paste(ds, "plotting", type, "\n"))

        plotdev(file.path(cfig.path,
                          paste0("codons_classes_gingold14_",ds,"_",type)),
                type=ftyp, res=300, width=2.5,height=3.5)
        par(mai=c(.6,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25, xpd=TRUE)
        bp <- boxplot(dtmt$RAAS ~ dtmt$CCLS, xlab=NA,
                      ylab=expression(log[10](RAAS)), axes=FALSE,
                      ylim=c(-5.7,3.3))
        if ( min(table(dtmt$CCLS)) < 50 )
            stripchart(RAAS ~ CCLS, vertical = TRUE, data = dtmt, 
                       method = "jitter", add = TRUE, pch = 20,
                       col="#00000077", cex=1,
                       axes=FALSE)
        axis(2)
        axis(1, at=1:length(bp$names), labels=bp$names, las=2, col=NA,
             mgp=c(0,0,0))
        axis(3,  at=1:length(bp$names), labels=cnt[bp$names], las=1, col=NA,
             mgp=c(0,-1,0), cex.axis=.8)
        ##mtext("codon classes\nGingold et al. 2014", 1, 3)
        if ( !is.null(tt) ) {
            x0 <- which(bp$names==d1)
            x1 <- which(bp$names==d2)
            arrows(x0=x0, x1=x1,
                   y0=par("usr")[4], code=3, angle=90, length=.05, xpd=TRUE)
            text((x0+x1)/2,
                 par("usr")[4], paste0("p=",signif(tt$p.value,1)),
                 pos=3, xpd=TRUE)
        }
        figlabel(dsl, pos="bottomleft", font=2, cex=1.2) 
        figlabel(type, pos="bottomright", font=2, cex=1.2) 
        dev.off()
    }
}
