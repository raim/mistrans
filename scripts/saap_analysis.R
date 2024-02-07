
## TODO
## 2012: Jonathan Weissman, Ignolia, Cell paper: MOUSE embryonic stem cells,
## ribosome density at mutated position.
##  https://doi.org/10.1016/j.cell.2011.10.002
## Bacterial mistranslation v ribosome density.

## ribosome density: ribosome density per codon divided
## by the mean ribosome density of transcript.

## TODO: RAAS vs. codon->AA profiles
## * combine patient and TMT levels RAAS,
## * indicate points in volcano,

## * CODONS: unique SAAP


## * get more AA similarity measures, eg. Grentham's
## https://en.wikipedia.org/wiki/Amino_acid_replacement#Physicochemical_distances

source("~/work/mistrans/scripts/saap_utils.R")

library(readxl)
library(viridis)
library(segmenTools)
library(Biostrings) # for blosum62
data(BLOSUM62)
data(PAM250)
library(seqinr) # just for 1/3 letter code
library(grantham) # grantham distance

options(stringsAsFactors=FALSE)

###  CODONS
aa <- unique(GENETIC_CODE)
CODONS <- rep("", length(aa))
for ( i in seq_along(aa) )
    CODONS[i] <- paste(names(which(GENETIC_CODE==aa[i])), collapse=";")
names(CODONS) <- aa
ACODONS <- paste0(names(CODONS),": ", CODONS)
names(ACODONS) <- aa

## 3<->1 letter code
AAC <- seqinr::aaa( seqinr::a())
names(AAC) <-  seqinr::a()
AAC <- AAC[AAC!="Stp"]

## generate Grantham distance
## TODO: store several matrices in project folder and just load them instead
## of generating or using a package, e.g.
## https://en.wikipedia.org/wiki/Amino_acid_replacement#Physicochemical_distances

gmat <- matrix(NA, ncol=length(AAC), nrow=length(AAC))
colnames(gmat) <- rownames(gmat) <- names(AAC)
for ( i in 1:nrow(gmat) )
    for ( j in 1:ncol(gmat) )
        gmat[i,j] <- unlist(grantham_distance(x = AAC[names(AAC)[i]],
                                              y = AAC[names(AAC)[j]],
                                              method="exact")[,3])
GRANTHAM <- gmat

NTS <- c("A","T","G","C")
## note: same as genomebrowser
ntcols <- c(A=rgb(85/255,107/255,47/255), ## green
            T=rgb(1,0,0), ## red
            G=rgb(1,.65,0),## orange
            C=rgb(0,0,1)) ## blue

cd.col <- c("darkgray", rgb(t(col2rgb(2))/255), rgb(t(col2rgb(3))/255))
names(cd.col) <- 1:3


#### START

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"

gen.path <- file.path(mam.path, "originalData")
dat.path <- file.path(proj.path,"originalData")
fig.path <- file.path(proj.path,"figures","saap_analysis")
out.path <- file.path(proj.path,"processedData")

in.file <- file.path(out.path,"saap_mapped.tsv")
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

##tmt.file <- file.path(proj.path,"originalData",
##                      "All_filtered_SAAP_TMTlevel_quant_df.xlsx")
##pat.file <- file.path(proj.path,"originalData",
##                      "All_filtered_SAAP_patient_level_quant_df.xlsx")
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")
pat.file <- file.path(proj.path,"originalData",
                      "All_SAAP_patient_level_quant_df.txt")

## analysis parameters
use.test <- stats::t.test # w.test

## plot colors
docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))

gcols <- grey.colors(n=100, start=.9, end=0)

## SAAP mapped to protein and codons
dat <- read.delim(in.file)
dat$Hemoglobin.Albumin <-
    as.logical(dat$Hemoglobin.Albumin)
PRAAS <- "Mean.precursor.RAAS" # precursor ion/experiment level
MRAAS <- "RAAS.median"


## TODO: calculate mean/median RAAS here from raw data
##       at the required level

## get ALL proteins - from project mammary
genes <- read.delim(feature.file)
##genes <- genes[genes$type=="protein_coding",]



### HOTSPOTS

## look for proteins with more than one mutation
## in close vicinity
hotspots <- split(dat$pos, f=dat$ensembl)
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

## QC: mean over proteins!


#### PLOTS

## get raw RAAS data patient level
patf <- read.delim(pat.file)
patf$Keep.SAAP <-  as.logical(patf$Keep.SAAP)
## remove excluded
cat(paste("removing", sum(!patf$Keep.SAAP),
          "tagged as false positive on patient level\n"))
patf <- patf[patf$Keep.SAAP,]
## exclude NA or Inf
rm <- is.na(patf$RAAS) | is.infinite(patf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from patient level\n"))
patf <- patf[!rm,]
pat <- split(patf$RAAS, patf$SAAP)

## FILTERED TABLES

## add useful aa/codon tag
dat$aacodon <- paste(dat$from,dat$codon, sep="-")



hdat <- dat[!(dat$remove | dat$Hemoglobin.Albumin),]
cdat <- hdat[hdat$codon!="",] # NOTE: this filters out duplicated SAAP



## get raw RAAS data TMT level
tmtf <- read.delim(tmt.file)
tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)
## remove excluded
cat(paste("removing", sum(!tmtf$Keep.SAAP),
          "tagged as false positive on TMT level\n"))
tmtf <- tmtf[tmtf$Keep.SAAP,]
## exclude NA or Inf
rm <- is.na(tmtf$RAAS) | is.infinite(tmtf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from TMT level\n"))
tmtf <- tmtf[!rm,]
tmt <- split(tmtf$RAAS, tmtf$SAAP)

## get RAAS coloring scheme
## TODO: use tmt level RAAS
png(file.path(fig.path,"raas_colors.png"),
    res=300, width=3, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
raas.col <- selectColors(unlist(tmt),
                         q=0, mx=1, mn=-4,
                         colf=viridis::viridis, plot=TRUE,
                         mai=c(.5,.5,.1,.1),
                         xlab="All TMT level RAAS",
                         n=40)
dev.off()


### CODONS vs. RAAS
## LOOP OVER DATASETS
for ( set in c(unique(cdat$Dataset),"all") ) {

    set.path <- file.path(fig.path, set)
    dir.create(set.path, showWarnings=FALSE)

    cat(paste(set,"\n"))

    ## EXPERIMENT FILTER
    if ( set=="all" ) {
        s.cdat <- cdat
        s.tmtf <- tmtf
    } else {
        s.cdat <- cdat[cdat$Dataset==set,]
        s.tmtf <- tmtf[tmtf$Dataset==set,]
    }
    s.tmt <- split(s.tmtf$RAAS, s.tmtf$SAAP)

    png(file.path(set.path,"raas_distribution.png"),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
    hist(unlist(pat), add=FALSE, border=2, main=NA)
    hist(unlist(s.tmt), add=TRUE, border=3)
    hist(s.cdat$RAAS.median, add=TRUE)
    legend("topright", c("patient","TMT","unique SAAP w codon"),
           col=c(2,3,1), lty=1)
    dev.off()

    ##
    to.lst <- list()
    for ( tol in sort(unique(s.cdat$to)) )
        to.lst[[tol]] <- unlist(tmt[s.cdat$SAAP[s.cdat$to == tol]])
    png(file.path(set.path,"raas_to.png"),
        res=300, width=6, height=3, units="in")
    par(mai=c(.5,.5,.5,.15), mgp=c(1.4,.3,0), tcl=-.25)
    boxplot(to.lst,ylab="all TMT RAAS")
    axis(3, at=1:length(to.lst), label=unlist(lapply(to.lst,length)),las=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    ## TTESTS

    ova <- raasProfile(x=s.cdat, id="SAAP", values=s.tmt,
                   rows="to", cols="aacodon", use.test=t.test, do.plots=FALSE)

    source("~/programs/segmenTools/R/clusterTools.R")
    ##TODO: find bug in plotOverlaps for LSCC
    png(file.path(set.path,"codons_all_ttests.png"),
        res=300, width=14, height=5, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.7,.5,.7,2.7))
    tmp <- plotOverlaps(ova, p.min=1e-10, p.txt=1e-5,
                 text.cex=.8, axis=1:3, ##type="unique",
                 ylab="mistranslated AA",
                 xlab="", col=ttcols, cut=FALSE)#, rmz =FALSE, short=FALSE)
    axis(4, nrow(ova$p.value):1, labels=ACODONS[rownames(ova$p.value)], las=2)
    figlabel(pos="bottomright", text="p value/t-test", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()

    
    png(file.path(set.path,"codons_all_ttests_legend.png"),
        res=300, width=2, height=2, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlapsLegend(p.min=1e-10, p.txt=1e-5, type=2, col=ttcols)
    dev.off()
    
    if ( !interactive() ) {
        
        ova <- raasProfile(x=s.cdat, id="SAAP", values=s.tmt,
                           rows="to", cols="aacodon", use.test=w.test,
                           do.plots=FALSE)
        
        png(file.path(set.path,"codons_all_wtests.png"),
            res=300, width=14, height=5, units="in")
        par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
        par(mai=c(.7,.5,.7,2.7))
        plotOverlaps(ova, p.min=1e-10, p.txt=1e-5,
                     text.cex=.7, axis=1:3, ##type="unique",
             ylab="mistranslated AA",
             xlab="", col=ttcols)#, rmz =FALSE, short=FALSE)
        axis(4, nrow(ova$p.value):1, labels=ACODONS[rownames(ova$p.value)],
             las=2)
        figlabel(pos="bottomright", text="p value/wilcox", cex=1.5, font=2)
        figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
        dev.off()
        
    }
    
    ## median RAAS
    png(file.path(set.path,"codons_all_raas_median.png"),
        res=300, width=14, height=5, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.7,.5,.7,2.7))
    txt <- ova$count
    txt[txt==0] <- ""
    txt.col <- ifelse(ova$median< -2, "white","black")
    image_matrix(ova$median, col=raas.col$col, cut=TRUE, breaks=raas.col$breaks,
                 axis=1:3, text=txt, text.col=txt.col, ylab="mistranslated",
                 xlab="", text.cex=.8)
    axis(4, nrow(ova$p.value):1, labels=ACODONS[rownames(ova$p.value)], las=2)
    figlabel(pos="bottomright", text="median RAAS", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    ## mean RAAS
    png(file.path(set.path,"codons_all_raas_mean.png"),
        res=300, width=14, height=5, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.7,.5,.7,2.7))
    txt <- ova$count
    txt[txt==0] <- ""
    txt.col <- ifelse(ova$mean< -2, "white","black")
    image_matrix(ova$mean, col=raas.col$col, cut=TRUE, breaks=raas.col$breaks,
                 axis=1:3, text=txt, text.col=txt.col, ylab="mistranslated",
             xlab="", text.cex=.8)
    axis(4, nrow(ova$p.value):1, labels=ACODONS[rownames(ova$p.value)], las=2)
    figlabel(pos="bottomright", text="mean RAAS", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    ## total count
    png(file.path(set.path,"codons_all_count.png"),
        res=300, width=14, height=5, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.7,.5,.7,2.7))
    cnt <- ova$count
    cnt[cnt==0] <- NA
    txt <- ova$count
    txt[txt==0] <- ""
    txt.col <- ifelse(ova$count>300,"white","black")
    image_matrix(cnt, col=gcols, axis=1:3,
             text=txt, text.col=txt.col, ylab="mistranslated",
             xlab="", text.cex=.8)
    axis(4, nrow(ova$p.value):1, labels=ACODONS[rownames(ova$p.value)], las=2)
    figlabel(pos="bottomright", text="total count", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    png(file.path(set.path,"codons_all_count_unique.png"),
        res=300, width=14, height=5, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.7,.5,.7,2.7))
    cnt <- ova$unique
    cnt[cnt==0] <- NA
    txt <- ova$unique
    txt[txt==0] <- ""
    txt.col <- ifelse(ova$unique>100,"white","black")
    image_matrix(cnt, col=gcols, axis=1:3,
                 text=txt, text.col=txt.col, ylab="mistranslated",
                 xlab="", text.cex=.8)
    axis(4, nrow(ova$p.value):1, labels=ACODONS[rownames(ova$p.value)], las=2)
    figlabel(pos="bottomright", text="unique SAAP count", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    df <- data.frame(freq=c(ova$count), raas=c(ova$median))
    png(file.path(set.path,"codons_all_raas_frequency.png"),
        res=300, width=3, height=3, units="in")
    par(mai=c(.5,.5,.2,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(df$raas, log10(df$freq), ylab="codon->AA count",
            xlab="median RAAS", axes=FALSE)
    axis(1)
    axis(2, at=log10(c(1:10,(1:10)*10,1:10*100)), labels=FALSE, tcl=par("tcl")/2)
    axis(2, at=log10(10^(0:5)), labels=10^(0:5))
    dev.off()


    ## VOLCANO PLOT

    ##source("~/work/mistrans/scripts/saap_utils.R")
    png(file.path(set.path,"codons_all_ttests_volcano.png"),
        res=300, width=5, height=3, units="in")
    par(mai=c(.5,.5,.1,.5), mgp=c(1.3,.3,0), tcl=-.25)
    volcano(ova, xlab="median RAAS", value="median",
            p.txt=5, v.txt=c(-4,-1),
            mid=mean(unlist(s.tmt)))
    abline(v=mean(unlist(s.tmt)))
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    

    
    ## t-tests over AA
    
    ova <- raasProfile(x=s.cdat, id="SAAP", values=s.tmt,
                       rows="to", cols="from", use.test=t.test, do.plots=FALSE)
    
    png(file.path(set.path,"AA_all_ttests.png"),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ova, p.min=1e-10, p.txt=1e-5,
                 text.cex=.7, axis=1:2,
                 ylab="mistranslated AA",
                 xlab="encoded AA", col=ttcols,
                 show.total=TRUE)#, rmz =FALSE, short=FALSE)
    figlabel(pos="bottomright", text="p value/t-test", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    png(file.path(set.path,"AA_all_ttests_legend.png"),
        res=300, width=2, height=2, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlapsLegend(p.min=1e-10, p.txt=1e-5, type=2, col=ttcols)
    dev.off()
    
    if ( !interactive() ) {
        
        ova <- raasProfile(x=s.cdat, id="SAAP", values=s.tmt,
                           rows="to", cols="from", use.test=w.test,
                           do.plots=FALSE)
        
        png(file.path(set.path,"AA_all_wtests.png"),
            res=300, width=5, height=5, units="in")
        par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        plotOverlaps(ova, p.min=1e-10, p.txt=1e-5,
                     text.cex=.7, axis=1:2,
                     ylab="mistranslated AA",
                     xlab="encoded AA", col=ttcols,
                     show.total=TRUE)#, rmz =FALSE, short=FALSE)
        figlabel(pos="bottomright", text="p value/wilcox", cex=1.5, font=2)
        figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
        dev.off()
        
        png(file.path(set.path,"AA_all_wtests_raas.png"),
            res=300, width=5, height=5, units="in")
        par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        ovr <- ova
        ovr$count <- round(ovr$median,1)
        ovr$count[is.na(ovr$count)] <- ""
        plotOverlaps(ovr, p.min=1e-10, p.txt=1e-5,
             text.cex=.7, axis=1:2,
             ylab="mistranslated AA",
             xlab="encoded AA", col=ttcols,
             show.total=TRUE, rmz =FALSE, short=FALSE)
        figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
        figlabel(pos="bottomright", text="p value/wilcox", cex=1.5, font=2)
        dev.off()
        
    }
    
    
    ## VOLCANO PLOT
    png(file.path(set.path,"AA_all_ttests_volcano.png"),
        res=300, width=5, height=3, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    volcano(ova, xlab="median RAAS", value="median",
            p.txt=5, v.txt=c(-4,-1),
            mid=mean(unlist(s.tmt)))
    abline(v=median(unlist(s.tmt)))
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    

    ## median RAAS
    png(file.path(set.path,"AA_all_raas_median.png"),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
    txt <- ova$count
    txt[txt==0] <- ""
    txt.col <- ifelse(ova$median< -2, "white","black")
    image_matrix(ova$median, col=raas.col$col, cut=TRUE, breaks=raas.col$breaks,
                 axis=1:2, text=txt, text.col=txt.col, ylab="mistranslated",
                 xlab="encoded AA", , text.cex=.8)
    figlabel(pos="bottomright", text="median RAAS", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    ## mean RAAS
    png(file.path(set.path,"AA_all_raas_mean.png"),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
    txt <- ova$count
    txt[txt==0] <- ""
    txt.col <- ifelse(ova$mean< -2, "white","black")
    image_matrix(ova$mean, col=raas.col$col, cut=TRUE, breaks=raas.col$breaks,
                 axis=1:2, text=txt, text.col=txt.col, ylab="mistranslated",
                 xlab="encoded AA", text.cex=.8)
    figlabel(pos="bottomright", text="mean RAAS", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    ## total count
    png(file.path(set.path,"AA_all_count.png"),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
    cnt <- ova$count
    cnt[cnt==0] <- NA
    txt <- ova$count
    txt[txt==0] <- ""
    txt.col <- ifelse(ova$count>300,"white","black")
    image_matrix(cnt, col=gcols, axis=1:2,
                 text=txt, text.col=txt.col, ylab="mistranslated",
                 xlab="encoded AA", text.cex=.8)
    figlabel(pos="bottomright", text="total count", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    png(file.path(set.path,"AA_all_count_unique.png"),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
    cnt <- ova$unique
    cnt[cnt==0] <- NA
    txt <- ova$unique
    txt[txt==0] <- ""
    txt.col <- ifelse(ova$unique>100,"white","black")
    image_matrix(cnt, col=gcols, axis=1:2,
                 text=txt, text.col=txt.col, ylab="mistranslated",
                 xlab="encoded AA", text.cex=.8)
    figlabel(pos="bottomright", text="unique SAAP count", cex=1.5, font=2)
    figlabel(pos="bottomleft", text=set, cex=1.5, font=2)
    dev.off()
    
    df <- data.frame(freq=c(ova$count), raas=c(ova$median))
    png(file.path(set.path,"AA_all_raas_frequency.png"),
        res=300, width=3, height=3, units="in")
    par(mai=c(.5,.5,.2,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(df$raas, log10(df$freq), ylab="AA->AA count",
            xlab="median RAAS", axes=FALSE)
    axis(1)
    axis(2, at=log10(c(1:10,(1:10)*10,1:10*100)), labels=FALSE, tcl=par("tcl")/2)
    axis(2, at=log10(10^(0:5)), labels=10^(0:5))
    dev.off()
    
    
    

    
###  CODONS
    
    ## use DiffLogo as for Behle et al.! 
    library(seqLogo)
    require(ggplot2)
    require(ggseqlogo)
    require(DiffLogo)
    ## with pwm1: AAS enriched codon,
    ##      pwm2: global human, local proteins, local peptides.
    ## dlogo <- createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
    ## dlogo <- enrichDiffLogoObjectWithPvalues(dlogo,n1=n1, n2=n2)
    ## see 20201204_RM_topA_promoters.R for plot
    

    ## TODO:
    ## * store codon context,-1,+1,
    ## * check consistency of dat$from and dat$codon,
    ##   -> requires to account for genomic mutations?
    ## * load codon usage table, AA usage table,
    
    
    
    
    
    ## codon positions
    cpos <- strsplit(s.cdat$codon,"")
    ## https://pubmed.ncbi.nlm.nih.gov/11164038/ A vs. U in 2nd
    c1 <- factor(unlist(lapply(cpos, function(x) x[1])), levels=NTS)
    c2 <- factor(unlist(lapply(cpos, function(x) x[2])), levels=NTS)
    c3 <- factor(unlist(lapply(cpos, function(x) x[3])), levels=NTS)
    cm <- cbind(c1=as.character(c1), c2=as.character(c2), c3=as.character(c3))


    cfrq <- rbind("1"=table(c1),
                  "2"=table(c2),
              "3"=table(c3))
    cfrq <- cfrq[,NTS]
    
    pwm <- t(cfrq)
    pwm <- t(t(pwm)/apply(pwm,2,sum))
    
    png(file.path(set.path,"codons_logo.png"),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.05,.05), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    seqLogo::seqLogo(makePWM(pwm), ic.scale=FALSE)
    dev.off()
    ##DiffLogo::seqLogo(makePWM(pwm), sparse=TRUE, ylim=c(0,.3))
    
    ## TODO: difflogo for each codon, AAS vs. rest of dataset,
    ## sequence logos of AAs and nucleotides surrounding AAS.
    
      

    ## GET CODON-LEVEL RAAS MEANS from tmt level
    cp.lst <- list()
    for ( pos in 1:3 ) {
        cp.lst[[paste0("c",pos)]] <- list()
        for ( nt in names(ntcols) ) 
        cp.lst[[paste0("c",pos)]][[nt]] <-
            unlist(s.tmt[s.cdat$SAAP[cm[,pos]==nt]])
    }
    
    png(file.path(set.path,"codons_pos.png"),
        res=300, width=3, height=1.5, units="in")
    par(mai=c(.25,.5,.2,.25), mgp=c(1.3,.3,0), tcl=-.25)
    barplot(t(cfrq),beside=TRUE, legend.text=colnames(cfrq),
            args.legend=list(x="top",ncol=4, inset=c(-.02,-.1), bty="n"),
            xlab="codon position", col=ntcols)
    dev.off()
    
    png(file.path(set.path,"codons_pos_raas.png"),
        res=300, width=3, height=3, units="in")
    par(mai=c(.5,.5,.2,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    par(mai=c(.5,.5,.15,.25))
    boxplot(cp.lst[[1]], ylab="all TMT RAAS", boxwex=.75, at=1:4,
            names=NA, xlab=NA, axes=FALSE, xlim=c(0.25,14), 
            col=ntcols, cex=.5)
    axis(2)
    boxplot(cp.lst[[2]], ylab=NA, add=TRUE, at=1:4 + 4.5, boxwex=.75,
            names=NA, axes=FALSE, col=ntcols, cex=.5)
    boxplot(cp.lst[[3]], ylab=NA, add=TRUE, at=1:4 + 9, boxwex=.75,
            names=NA, axes=FALSE, col=ntcols, cex=.5)
    axis(1, at=c(2.5,7,11.5), labels=1:3, col=NA)
    mtext("codon position", 1, 1.3)
    axis(4)
    dev.off()
    
    png(file.path(set.path,"codons_type.png"),
        res=300, width=3, height=1.5, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    par(mai=c(0.1,.5,.2,.1))
    barplot(cfrq,beside=TRUE, legend.text=rownames(cfrq), col=cd.col,
            args.legend=list(x="top",ncol=3, inset=c(-.02,-.1), bty="n",
                             title="codon position"))
    dev.off()
    png(file.path(set.path,"codons_type_raas.png"),
        res=300, width=3, height=3, units="in")
    par(mai=c(.5,.5,.2,.1), mgp=c(1.3,.3,0), tcl=-.25)
    par(mai=c(.5,.5,.15,.1))
    boxplot(cp.lst[[1]], ylab="all TMT RAAS", at=1:4 -.5, boxwex=.18,
            names=rep(1,4), xlab=NA, xlim=c(.4,4.05), axes=FALSE, col=cd.col[1])
    axis(2)
    boxplot(cp.lst[[2]], ylab=NA, add=TRUE, at=1:4 -.3, boxwex=.18,
            names=rep(2,4), axes=FALSE, col=cd.col[2])
    boxplot(cp.lst[[3]], ylab=NA, add=TRUE, at=1:4 -.1, boxwex=.18,
            names=rep(3,4), axes=FALSE, col=cd.col[3])
    axis(1, at=1:4 -.3, mgp=c(10,1.3,0), tcl=0, labels=levels(c1))
    dev.off()
    
    ## RAAS by all codons
    aa.lst <- list()
    for ( aac in sort(unique(s.cdat$aacodon)) )
        aa.lst[[aac]] <- unlist(s.tmt[s.cdat$SAAP[s.cdat$aacodon==aac]])
    
    
    png(file.path(set.path,"codons_raas.png"),
        res=300, width=12.1, height=5, units="in")
    par(mai=c(.7,.5,.5,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    boxplot(aa.lst, las=2,
            xlab=NA, ylab="RAAS")
    tb <- table(s.cdat$aacodon)
    tt <- unlist(lapply(aa.lst, length))
    ##axis(3, at=1:length(tt), labels=tb, las=2)
    dev.off()
    
#### REPLACED AA SIMILARITIES

    ## TODO:
    ## * plot for several selections, strong difference for unfiltered 22k,
    ## * do on protein-level?
    ## * SWITCH TO RAAS DISTRIBUTION.


    ## first find mutated AA pairs
    ## (column AASin input mixes I/L)
    ## and split mutated AA pairs into from/to
    ## TODO: get this from the columns in the mapped file!
    saaps <- strsplit(s.cdat$SAAP,"")
    bases <- strsplit(s.cdat$BP, "")
    fromto <- lapply(1:length(saaps), function(i) {
        pos <- which(saaps[[i]]!=bases[[i]])
        c(from=bases[[i]][pos], to=saaps[[i]][pos])
    })
    


    ## analyze AA similarity by different matrices
    ## TODO: remove extra columns
    simmats <- c("BLOSUM62", "PAM250","GRANTHAM")
    for ( i in seq_along(simmats) ) {

        mid <- simmats[i]
        MAT <- get(mid)[AA_STANDARD,AA_STANDARD]
        
        ylab <- "similarity between replaced AA" #"BLOSUM62 similarity"
        
        ## similarity of replaced AA
        sim <- unlist(lapply(fromto, function(x) MAT[x[1],x[2]]))
        
        rid <- "Mean.precursor.RAAS"
        raas <- unlist(s.cdat[,rid])
        
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
        
        
        
        png(file.path(set.path,paste0("AAS_",mid,"_raas_dense.png")),
        res=300, width=3.5, height=3.5, units="in")
        par(mai=c(.5,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25)
        plotCor(dff$sim, dff$raas, ylab=ylab, xlab=rid, colf=viridis::viridis)
        dev.off()
        
        png(file.path(set.path,paste0("AAS_",mid,"_raas_boxplot.png")),
            res=300, width=5, height=3.5, units="in")
        par(mai=c(.75,.5,.5,.1), mgp=c(1.3,.3,0), tcl=-.25)
        boxplot(sim ~ bins, data=df, ylab=ylab,
                xlab=NA, las=2, at=seq_along(levels(df$bins)))
        axis(3, at=seq_along(levels(df$bins)), labels=table(df$bins),las=2)
        mtext(rid, 1, 2.5)
        dev.off()
        
        
        png(file.path(set.path,paste0("AAS_",mid,"_raas_boxplot_dual.png")),
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
        text(1.5, par("usr")[4], paste("p =",signif(tt$p.value,1)),
             xpd=TRUE, pos=3)
        arrows(x0=1, y0=par("usr")[4], x1=2, angle=90, code=3,
               length=.05, xpd=TRUE)
        dev.off()
        
        png(file.path(set.path,paste0("AAS_",mid,"_similarities_hist.png")),
            res=300, width=5, height=3.5, units="in")
        hist(MAT[lower.tri(MAT)])
        dev.off()
    }
}

for ( set in c(unique(cdat$Dataset),"all") ) {

    set.path <- file.path(fig.path, set)
    dir.create(set.path, showWarnings=FALSE)

    cat(paste(set,"\n"))

    ## EXPERIMENT FILTER
    if ( set=="all" ) {
        s.hdat <- hdat
        s.tmtf <- tmtf
    } else {
        s.hdat <- hdat[cdat$Dataset==set,]
        s.tmtf <- tmtf[tmtf$Dataset==set,]
    }
    s.tmt <- split(s.tmtf$RAAS, s.tmtf$SAAP)


### SECONDARY STRUCTURE

    ## TODO: load background

    ## explore s4pred results
    ss.aas <- table(s.hdat$s4pred)
    ss.aas <- ss.aas/sum(ss.aas)
    
    ## TODO: find out what's missing
    ss.aas <- ss.aas[names(ss.aas)!=""]
    
    ss.tot <- apply(s.hdat[,paste0(names(ss.aas),".protein")],2, sum,na.rm=TRUE)
    ss.tot <- ss.tot/sum(ss.tot)
    
    ss.tab <- rbind(AAS=ss.aas, total=ss.tot)
    
    ss.nms <- c(C=expression(coil),
                E=expression(beta),
                H=expression(alpha))
    
    ## slight enrichment of beta/alpha
    png(file.path(set.path,"structure_s4pred.png"),
        res=300, width=3, height=3, units="in")
    par(mai=c(.25,.5,.15,.25), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    barplot(ss.tab, beside=TRUE, legend=TRUE, ylab="fraction",
            names.arg=ss.nms[colnames(ss.tab)])
    axis(4)
    dev.off()
    
    
    ss.lst <- list()
    for ( ssl in c("C","E","H") )
        ss.lst[[ssl]] <- unlist(s.tmt[s.hdat$SAAP[cdat$s4pred==ssl]])
    
    png(file.path(set.path,"structure_s4pred_raas.png"),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.25), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    boxplot(ss.lst, ylab="all TMT RAAS", xlab=NA, las=2, axes=FALSE)
    axis(2)
    axis(4)
    axis(1, at=1:length(ss.lst), labels=ss.nms[names(ss.lst)], las=1)
    figlabel("AAS s4pred", pos="bottomleft")
    dev.off()
    
    
    ## GET STRUCTURE SCORES
    
    iup <- s.hdat$iupred3
    ## TODO: why negative IUPRED3??
    if ( min(iup,na.rm=TRUE)<0 )
        iup <- iup - min(iup,na.rm=TRUE)
    iubg <- s.hdat$iupred3.protein
    
    anc <- s.hdat$anchor2
    anbg <- s.hdat$anchor2.protein
    
    png(file.path(set.path,"structure_iupred3_bg.png"),
        res=300, width=3, height=3, units="in")
    par(mai=c(.5,.5,.1,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotCor(iup, iubg, xlab="iupred3 score  at AAS site",
            ylab="mean iupred3 score of protein") 
    dev.off()
    

    ## bin score to compare to full list of RAAS
    iup.bins <- cut(iup, breaks=seq(0,1,.2))
    levels(iup.bins)
    iup.lst <- list()
    for ( iupl in levels(iup.bins) )
        iup.lst[[iupl]] <-
            unlist(s.tmt[s.hdat$SAAP[as.character(iup.bins)==iupl]])

    png(file.path(set.path,"structure_iupred3_raas.png"),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.75,.5,.5,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    boxplot(iup.lst, ylab="all TMT RAAS", xlab=NA,
            las=2)
    axis(3, at=1:length(iup.lst),
         labels=unlist(lapply(iup.lst, length)), las=2)
    figlabel("AAS iupred3", pos="bottomleft")
    dev.off()
    
    ## slight positive trend of unstructured/anchor vs RAAS
    png(file.path(set.path,"structure_iupred3.png"),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    hist(iubg, border=NA, col=NA, breaks=seq(0,1,.05), xlab="iupred3 score",
         main=NA)
    hist(iubg, col=1, border=1, add=TRUE, breaks=seq(0,1,.05))
    hist(iup, border=2, col=paste0(rgb(t(col2rgb(2))/255),77),
         add=TRUE, breaks=seq(0,1,.05))
    legend("topright", c("AAS", "mean of protein",
                         paste0("p=",signif(wilcox.test(iup, iubg)$p.value,2))),
           col=c(2,1,NA), lty=1, bty="n")
    dev.off()       
    
    
    ## anchor2
    
    png(file.path(set.path,"structure_anchor2_bg.png"),
        res=300, width=3, height=3, units="in")
    par(mai=c(.5,.5,.1,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotCor(anc, anbg, xlab="anchor2 score  at AAS site",
            ylab="mean anchor2 score of protein") 
    dev.off()
    
    png(file.path(set.path,"structure_anchor2.png"),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    hist(anc, border=NA, col=NA, breaks=seq(0,1,.05), main=NA,
         xlab="anchor2 score")
    hist(anbg, col=1, border=1, add=TRUE, breaks=seq(0,1,.05))
    hist(anc, border=2, col=paste0(rgb(t(col2rgb(2))/255),77),
         add=TRUE, breaks=seq(0,1,.05), main=NA)
    legend("topright", c("AAS", "mean of protein",
                         paste0("p=",signif(wilcox.test(anc, anbg)$p.value,2))),
           col=c(2,1,NA), lty=1, bty="n")
    dev.off()
    
    ## bin score to compare to full list of RAAS
    anc.bins <- cut(anc, breaks=seq(0,1,.2))
    levels(anc.bins)
    anc.lst <- list()
    for ( ancl in levels(anc.bins) )
        anc.lst[[ancl]] <-
            unlist(s.tmt[s.hdat$SAAP[as.character(anc.bins)==ancl]])
    
    png(file.path(set.path,"structure_anchor2_raas.png"),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.75,.5,.5,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    boxplot(anc.lst, ylab="all TMT RAAS", xlab=NA, las=2)
    axis(3, at=1:length(anc.lst),
         labels=unlist(lapply(anc.lst, length)), las=2)
    figlabel("AAS anchor2", pos="bottomleft")
    dev.off()
    
    
    ## normalized AAS/protein
    png(file.path(set.path,"structure_anchor2_raas_norm.png"),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotCor(s.hdat$RAAS.median, anc/anbg, xlab="median RAAS",
            ylab="normalized anchor2 score, AAS/protein")
    dev.off()
    
    png(file.path(set.path,"structure_iupred3_raas_norm.png"),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plotCor(s.hdat$RAAS.median, iup/iubg, xlab="median RAAS",
            ylab="normalized iupred3 score, AAS/protein", ylim=c(0,4))
    dev.off()
    
    
### POSITION v LENGTH
    
    ## size cut off
    size.cutoff <- c(2.5e2,1e3)
    small <- s.hdat$len <  size.cutoff[1]
    large <- s.hdat$len >= size.cutoff[1] & s.hdat$len <= size.cutoff[2]
    huge  <- s.hdat$len >  size.cutoff[2]
    
    hist(s.hdat$pos, breaks=seq(0,4e4,100), xlim=c(0,3000),
         xlab="absolute position",
         main=NA)
    hist(s.hdat$len, breaks=seq(0,4e4,100), border=2, add=TRUE)
    legend("topright", c("BP position","protein length"), col=c(1,2), lty=1)
    
    
    png(file.path(set.path,"position_length_total.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    hist(log10(s.hdat$len),
         xlab=expression(protein~length/log[10]), breaks=50)
    abline(v=log10(size.cutoff))
    dev.off()
    
    
    png(file.path(set.path,"position_length_relative.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(s.hdat$rpos, log10(s.hdat$len), colf=viridis::viridis,
            ##len, log="y",
            xlab="relative position of AAS in protein",
            ylab=expression(protein~length/log[10]))
    abline(h=log10(size.cutoff))
    dev.off()
    
    png(file.path(set.path,"position_length_absolute.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i",
        xpd=TRUE)
    dense2d(s.hdat$len, s.hdat$pos, colf=viridis::viridis,
            xlim=c(0,1500), ylim=c(0,1500),
            ylab=expression(SAAP~position),
            xlab=expression(protein~length), nbin=512)
    par(xpd=FALSE)
    ##abline(a=0, b=1)
    abline(a=0, b=.5)
    abline(v=size.cutoff)
    dev.off()
    
    png(file.path(set.path,"position_hist_absolute.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i",
        xpd=TRUE)
    hist(s.hdat$pos)
    dev.off()
    
    
    png(file.path(set.path,"position_ecdf_relative_short.png"),
        res=300, width=2.5, height=2.5, units="in")
    par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
    cdf <- ecdf(s.hdat$rpos[small]) 
    plot(seq(0,1,.01), cdf(seq(0,1,.01)), xlim=c(0,1), ylim=c(0,1),
         xlab="relative position in protein",
         ylab=expression(ecdf(x)), col=2, type="l", main=NA)
    abline(a=0,b=1, col=1)
    dev.off()
    
    png(file.path(set.path,"position_ecdf_relative_mid.png"),
        res=300, width=2.5, height=2.5, units="in")
    par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
    cdf <- ecdf(s.hdat$rpos[large])
    plot(seq(0,1,.01), cdf(seq(0,1,.01)), xlim=c(0,1), ylim=c(0,1),
         xlab="relative position in protein",
         ylab=expression(ecdf(x)), col=2, type="l", main=NA)
    abline(a=0,b=1, col=1)
    dev.off()
    
    png(file.path(set.path,"position_ecdf_relative_long.png"),
        res=300, width=2.5, height=2.5, units="in")
    par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
    cdf <- ecdf(s.hdat$rpos[huge])
    plot(seq(0,1,.01), cdf(seq(0,1,.01)), xlim=c(0,1), ylim=c(0,1),
         xlab="relative position in protein",
         ylab=expression(ecdf(x)), col=2, type="l", main=NA)
    abline(a=0,b=1, col=1)
    dev.off()
    
    png(file.path(set.path,"position_ecdf_absolute.png"),
        res=300, width=2.5, height=2.5, units="in")
    par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
    cdf <- ecdf(s.hdat$pos)
    x <- seq(0,max(s.hdat$pos, na.rm=TRUE),1)
    plot(x, cdf(x), xlim=c(0,1000), xlab="absolute position in protein",
         ylab=expression(ecdf(x)), type="l", col=2)
    abline(v=8.5, col=1)
    dev.off()
    
    png(file.path(set.path,"position_ecdf_absolute_zoom.png"),
        res=300, width=2.5, height=2.5, units="in")
    par(mai=c(.5,.5,.15,.15), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
    cdf <- ecdf(s.hdat$pos)
    plot(seq(0,100,.01), cdf(seq(0,100,.01)),xlim=c(0,30), ylim=c(0,.05),
         xlab="absolute position in protein",
         ylab=expression(ecdf(x)), type="l", col=2)
    abline(v=10, col=1)
    dev.off()
    
    png(file.path(set.path,"position_length_absolute_log.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i")
    dense2d(log10(s.hdat$len), log10(s.hdat$pos),
            colf=viridis::viridis, ylim=c(1,4), xlim=c(1,4),
            ylab=expression(SAAP~position/log[10]),
            xlab=expression(protein~length/log[10]))
    abline(a=0,b=1)
    dev.off()
    
    png(file.path(set.path,"position_hist_relative.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    hist(s.hdat$rpos, breaks=seq(0,1,0.1),
         main=paste(sum(!is.na(s.hdat$rpos)), "unique SAAP"),
         xlab="relative position of AAS in protein")
    loc.sze <- table(s.hdat$rpos>.5)
    dev.off()
    


    ## BACKGROUND: main peptides
    load(file.path(out.path, "mapped_peptides.rda"))
    png(file.path(set.path,"position_hist_relative_main_peptides.png"),
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
    png(file.path(set.path,"position_hist_relative_short.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
    hist(s.hdat$rpos[small], breaks=seq(0,1,0.1),
         main=paste0(sum(small,na.rm=TRUE),
                     " short proteins <",size.cutoff[1]," aa"),
         xlab="relative position of AAS in protein")
    par(new=TRUE)
    cdf <- ecdf(s.hdat$rpos[small])
    plot(seq(0,1,.01), cdf(seq(0,1,.01)), type="l", axes=FALSE, col=2,
         xlab=NA, ylab=NA)
    abline(a=0,b=1)
    mtext(expression(ecdf(x)), 4, 1.3)
    axis(4)
    dev.off()
    png(file.path(set.path,"position_hist_relative_long.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
    hist(s.hdat$rpos[huge], breaks=seq(0,1,0.1),
         main=paste0(sum(huge,na.rm=TRUE),
                     " long proteins >",size.cutoff[2]," aa"),
         xlab="relative position of AAS in protein")
    par(new=TRUE)
    cdf <- ecdf(s.hdat$rpos[huge])
    plot(seq(0,1,.01), cdf(seq(0,1,.01)), type="l", axes=FALSE, col=2,
         xlab=NA, ylab=NA)
    abline(a=0,b=1)
    mtext(expression(ecdf(x)), 4, 1.3)
    axis(4)
    dev.off()
    png(file.path(set.path,"position_hist_relative_mid.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
    hist(s.hdat$rpos[large], breaks=seq(0,1,0.1),
         main=paste0(sum(large,na.rm=TRUE),
                     " proteins >",size.cutoff[1]," aa, <",size.cutoff[2]," aa"),
         xlab="relative position of AAS in protein")
    par(new=TRUE)
    cdf <- ecdf(s.hdat$rpos[large])
    plot(seq(0,1,.01), cdf(seq(0,1,.01)), type="l", axes=FALSE, col=2,
         xlab=NA, ylab=NA)
    abline(a=0,b=1)
    mtext(expression(ecdf(x)), 4, 1.3)
    axis(4)
    dev.off()
    

    ## TODO: different analysis for global RAAS,
    ## e.g., position bins vs. all RAAS

    png(file.path(set.path,"position_RAAS_absolute.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(s.hdat$pos, s.hdat[,PRAAS], colf=viridis::viridis,#xlim=c(0,1e4),
            ylim=range(s.hdat[,PRAAS], na.rm=TRUE),
            xlab="absolute position of AAS in protein",
            ylab=PRAAS)
    dev.off()
    
    png(file.path(set.path,"position_RAAS_relative.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(s.hdat$rpos, s.hdat[,PRAAS], colf=viridis::viridis,
            xlim=range(s.hdat$rpos, na.rm=TRUE),
            ylim=range(s.hdat[,PRAAS], na.rm=TRUE),
            xlab="relative position of AAS in protein",
            ylab=PRAAS)
    dev.off()
    
    png(file.path(set.path,"protein_length_RAAS.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(s.hdat[,PRAAS], s.hdat$len, colf=viridis::viridis, na.rm=TRUE,
            ylim=c(0,1e4),
            xlim=range(s.hdat[,PRAAS], na.rm=TRUE),
            ylab="protein length",
            xlab=PRAAS)
    dev.off()
    
    png(file.path(set.path,"protein_length_RAAS_log.png"),
        res=300, width=4, height=3, units="in")
    par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(s.hdat[,PRAAS], log10(s.hdat$len), colf=viridis::viridis,na.rm=TRUE,
            xlim=range(s.hdat[,PRAAS], na.rm=TRUE),
            ylab=expression(log[10](protein~length)),
            xlab=PRAAS)
    dev.off()
    
}

