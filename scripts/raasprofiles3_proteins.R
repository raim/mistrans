
library(Biostrings) # for blosum62
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
source("~/work/mistrans/scripts/raasprofiles3_init.R")

pfig.path <- file.path(fig.path,"proteins")
dir.create(pfig.path, showWarnings=FALSE)
          
### START ANALYSIS


#### PROTEINS and COMPLEXES

## TODO: align use of raasProfile vs. listProfile

## MEDIAN RAAS PER UNIQUE PROTEIN SITE
source("~/work/mistrans/scripts/saap_utils.R")
site <- split(tmtf$RAAS, paste(tmtf$ensembl, tmtf$pos))
site <- listProfile(site, y=tmtf$RAAS, use.test=use.test, min=3)

if ( interactive() ) {
    ## test alternative measures of consistent RAAS
    ## CV:
    dense2d(site$cv, -log10(site$p), xlab="CV", ylab=expression(-log10(p)))
    dense2d(site$cv, site$sd, xlab="CV", ylab="SD")
    hist(site$sd)
}

## mod. column names and add protein and site info
colnames(site) <- paste0("RAAS.", colnames(site))
site$ensembl <- sub(" .*", "", rownames(site))
site$pos <- as.numeric(sub(".* ", "", rownames(site)))

## MEDIAN RAAS FOR EACH UNIQUE PROTEIN POSITION

## protein median of median RAAS
sitl <- split(site$RAAS.median, site$ensembl)
pbstat <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)


## median RAAS vs. number of sites
if ( interactive() )
    dense2d(pbstat$median, log10(pbstat$n))

## volcano
plotdev(file.path(pfig.path,paste0("proteins_volcano_sites")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
res <- volcano(pbstat, value="median",
               p.txt=-log10(0.01), v.txt=c(-2,-2), cut=5, mid=-2,
               ids=pnms, xlab=xl.site)
mtext("protein RAAS, median of site medians", 3,0)
dev.off()


## PROTEIN MEDIAN RAAS per protein w/o site-specific median first

ptl <- split(tmtf$RAAS, tmtf$ensembl)
ptstat <- listProfile(ptl, y=tmtf$RAAS, use.test=use.test, min=3)


plotdev(file.path(pfig.path,paste0("proteins_volcano_all")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
res <- volcano(ptstat, value="median",
               p.txt=12, v.txt=c(-2,-1), cut=50, mid=0,
               ids=pnms, xlab=xl.all)
mtext("protein RAAS, median of all RAAS", 3,0)
dev.off()

## compare median of all RAAS vs. median of site median RAAS
if ( interactive() ) {
    plotCor(pbstat$median, ptstat$median,
            xlab="median of site median RAAS", ylab="median of all RAAS")
}


## PROTEIN HALF-LIVES, @Mathieson2018
hlvd <- readxl::read_xlsx(math18.file)

## mean half-live over all replicates and cell types
## TODO: consider distribution
cidx <- grep("half_life", colnames(hlvd), value=TRUE)
hlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
names(hlv) <- unlist(hlvd[,1])


xl.hlf <- expression(protein~"half-life"/h)

## halflives site
plotdev(file.path(pfig.path,paste0("protein_halflives_site")),
        type="png", res=300, width=3.5,height=3.5)
idx <- match(pnms[rownames(pbstat)], names(hlv))
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
plotCor(pbstat$median, log10(hlv[idx]),
        xlab=xl.site, ylab=xl.hlf, axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
## TODO: log tick marks
##axis(2, at=rep(1:10, 4) * 10^rep(0:3, each=10), labels=FALSE)
dev.off()

## halflives all
plotdev(file.path(pfig.path,paste0("protein_halflives_all")),
        type="png", res=300, width=3.5,height=3.5)
idx <- match(pnms[rownames(ptstat)], names(hlv))
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
plotCor(ptstat$median, log10(hlv[idx]), signif = 1,
        xlab=xl.all, ylab=xl.hlf, axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
dev.off()

if ( interactive() ) {
    ## TODO: compare iupred3 and half-lives for all proteins, not just RAAS
    idx <- match(dat$name,names(hlv))
    plotCor(dat$iupred3.protein, log10(hlv[idx]))
}

for ( ctype in c("Bcells", "NK", "hepatocytes", "monocytes", "neurons") ) {
    cidx <- grep("half_life", colnames(hlvd), value=TRUE)
    cidx <- cidx[grep(ctype, cidx, ignore.case=TRUE)]
    hlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
    names(hlv) <- unlist(hlvd[,1])
    
    plotdev(file.path(pfig.path,paste0("protein_halflives_sites_",ctype)),
            type="png", res=300, width=3.5,height=3.5)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
    idx <- match(pnms[rownames(pbstat)], names(hlv))
    plotCor(pbstat$median, log10(hlv[idx]),
            xlab=xl.site, ylab=xl.hlf, axes=FALSE)
    axis(1)
    axis(2, at=1:10, labels=10^(1:10))
    figlabel(ctype, pos="bottomleft")
    dev.off()
    
    plotdev(file.path(pfig.path,paste0("protein_halflives_all_",ctype)),
            type="png", res=300, width=3.5,height=3.5)
    idx <- match(pnms[rownames(ptstat)], names(hlv))
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
    plotCor(ptstat$median, log10(hlv[idx]), signif = 1,
        xlab=xl.all, ylab=xl.hlf, axes=FALSE)
    axis(1)
    axis(2, at=1:10, labels=10^(1:10))
    figlabel(ctype, pos="bottomleft")
    dev.off()
}

## PROTEIN RAAS vs. IUPRED

## get whole protein mean iupred3 score
iu3 <- split(dat$iupred3.protein, dat$ensembl)
## QC: all protein level
table(unlist(lapply(iu3, function(x) length(unique))))
iu3 <- unlist(lapply(iu3, unique))

## NOTE: no correlation to protein-wide iupred3
plotdev(file.path(pfig.path,paste0("protein_iupred3_all")),
        type="png", res=300, width=3.5,height=3.5)
idx <- match(rownames(ptstat), names(iu3))
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
plotCor(ptstat$median, log(hlv[idx]), signif = 1,
        xlab=xl.all,
        ylab="iupred3")
dev.off()


## PROTEIN WINDOWS - 50 AA WINDOWS


## add position windows to raas table
## see cbind(1:200,round(seq(1:200)/50),
##           round((25+seq(1:200))/50))
## for window assignment of positions
tmtf$win1 <- round(tmtf$pos/50) # first is 1:25
tmtf$win2 <- round((tmtf$pos+25)/50) # first is 1:50

## search windows with high RAAS
windows <- append(split(tmtf$RAAS, paste(tmtf$ensembl, tmtf$win1)),
                  split(tmtf$RAAS, paste(tmtf$ensembl, tmtf$win2)))
hist(lengths(windows), breaks=100)

## RAAS statistics for windows
winstat <- listProfile(windows, y=tmtf$RAAS, min=3)
winstat$ensembl <- sub("\\..*","",rownames(winstat))
ids <- pnms[winstat$ensembl]
names(ids) <- rownames(winstat)

plotdev(file.path(pfig.path,paste0("windows_volcano_all")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
volcano(winstat, cut=50, v.txt=c(-2,0), p.txt=10,
        ids=ids,
        xlab=xl.all)
mtext("50 AA windows, median of all RAAS", 3,0)
dev.off()

## protein window stat based on site-specific RAAS medians

site$win1 <- round(site$pos/50) # first is 1:25
site$win2 <- round((site$pos+25)/50) # first is 1:50

## search windows with high RAAS
windows <- append(split(site$RAAS.median,
                        paste0(site$ensembl,"-w1-",site$win1)),
                  split(site$RAAS.median,
                        paste0(site$ensembl,"-w2-",site$win2)))
hist(lengths(windows), breaks=100)

winstat <- listProfile(windows, y=tmtf$RAAS, min=2)
winstat$ensembl <- sub("-w.*","",rownames(winstat))

ids <- strsplit(rownames(winstat),"-")
ids <- unlist(lapply(ids, function(x) pnms[x[1]]))
names(ids) <- rownames(winstat)

plotdev(file.path(pfig.path,paste0("windows_volcano_sites")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
volcano(winstat, value="median", cut=20, v.txt=c(-2,-2),
        p.txt=6, ids=unlist(ids),
        xlab=xl.site)
        
mtext("50 AA windows, median of site medians", 3,0)
show <- (winstat$median>0 & winstat$p<0.001) | winstat$median>1.6
shadowtext(winstat$median[show],
           -log10(winstat$p[show]),
           labels=ids[show],
           pos=4, xpd=TRUE, col=1, cex=.8)
show <- pnms[winstat$ensembl] %in% c("PGM1","PSMA1")
points(winstat$median[show],
       -log10(winstat$p[show]), col=3)
shadowtext(winstat$median[show],
           -log10(winstat$p[show]),
           labels=ids[show],
           pos=4, xpd=TRUE, col=3, font=2)
dev.off()

## PROTEIN COMPLEXES

## TODO: first boil down RAAS to median per site to avoid
## over-estimation by single sites, eg. Q->G in PSMA1.

## rename names by ensembl in position-wise RAAS medians

## distribution in protein complexes

## hu.MAP v2
if ( FALSE ) {
    humap <- read.csv(humap.file)
    hulst <- strsplit(humap[,3]," ")
    huids <- humap[,1]
}

## CORUM v4.1
humap <- read.delim(corum.file)
hulst <- strsplit(humap$subunits.UniProt.IDs.,";")
huids <- sub(" complex$","",humap[,2])

names(hulst)<- huids

## replace complex uniprot genes by ALL ensembl proteins that map to it
huens <- lapply(hulst, function(x) {
    unique(uni2ens[uni2ens[,1]%in%x,2])
})


## COMPLEXES - MEDIAN SITE RAAS
## replace complex ensembl IDs with rows in RAAS table
huraas <- lapply(huens, function(x) c(na.omit(pbstat[x,"median"])))

## RAAS statistics
hustat <- listProfile(huraas, y=tmtf$RAAS, min=3)

plotdev(file.path(pfig.path,paste0("complex_volcano_sites")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
res <- volcano(hustat, value="median",
               p.txt=3, v.txt=c(-Inf,-1.5), cut=50, mid=-1,
               xlab=xl.site)
mtext("protein complexes, median of site medians", 3,0)
dev.off()


## COMPLEXES - ALL RAAS
## replace complex ensembl IDs with all RAAS in TMT level table
huraas <- lapply(huens, function(x) tmtf$RAAS[tmtf$ensembl%in%x])

## RAAS statistics
hustat <- listProfile(huraas, y=tmtf$RAAS, min=3)

plotdev(file.path(pfig.path,paste0("complex_volcano_all")),
        type="png", res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
res <- volcano(hustat, value="median",
               p.txt=20, v.txt=c(-Inf,-1), cut=50, mid=0,
               xlab=xl.all)
mtext("protein complexes, median of all RAAS", 3,0)
dev.off()

## NOTE: EZR, MSN, 

## investigate complexes
if ( FALSE ) {
    humap[grep("Spliceosome, E", humap$ComplexName),c("ComplexName",
                                                      "Complex.comment",
                                                      "subunits.Gene.name.")]
    humap[grep("Drosha", humap$ComplexName),c("ComplexName",
                                              "Complex.comment",
                                              "subunits.Gene.name.")]
}



