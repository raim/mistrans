
library(Biostrings) # for blosum62
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

pfig.path <- file.path(fig.path,"proteins")
dir.create(pfig.path, showWarnings=FALSE)

corW <- corH <- 2.5
pmai <- c(.4,.4,.2,.2)
pmpg <- c(1,.2,0)

## overrule specific y-axis label
xl.prota <- xl.raas

### START ANALYSIS


#### PROTEINS and COMPLEXES

## TODO: align use of raasProfile vs. listProfile

## MEDIAN SITE AND PROTEIN RAAS

### median raas per unique mane protein site

sitl <- split(tmtf$RAAS, paste(tmtf$mane, tmtf$pos))
site <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## mod. column names and add protein and site info
colnames(site) <- paste0("RAAS.", colnames(site))
site$mane <- sub(" .*", "", rownames(site))
site$pos <- as.numeric(sub(".* ", "", rownames(site)))

## add gene names
site$name <- ens2nam [site$mane]

## add uniprot id
site$uniprot <- unlist(lapply(ens2u[site$mane], paste, collapse=";"))

## RAAS COLOR
site$RAAS.color <- num2col(site$RAAS.median,
                           limits=c(RAAS.MIN, RAAS.MAX), colf=arno)

## 

## protein median of site median RAAS
sitl <- split(site$RAAS.median, site$mane)
pbstat <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## protein median raas per protein w/o site-specific median first
ptl <- split(tmtf$RAAS, tmtf$mane)
ptstat <- listProfile(ptl, y=tmtf$RAAS, use.test=use.test, min=3)

if ( interactive() ) {
    ## test alternative measures of consistent RAAS
    ## CV:
    dense2d(site$RAAS.cv, -log10(site$RAAS.p.value),
            xlab="CV", ylab=expression(-log10(p)))
    dense2d(site$RAAS.cv, site$RAAS.sd, xlab="CV", ylab="SD")
    hist(site$RAAS.sd)
}

## ORDER SITES

## order proteins by RAAS
pbstat$rank <- rank(pbstat$median)

## make sure its ordered
## TODO: ORDER PROTEINS BY SIZE OR NUMBER OF AAS
site <-site[order(site$mane, site$pos),]

## order by protein RAAS rank (for hotspot plot)
site$rank <- pbstat[site$mane,"rank"]
site <-site[order(site$rank, site$pos),]

##
nsites <- as.numeric(sub(".*\\.","",tagDuplicates(site$mane)))
nsites[is.na(nsites)] <- 1
site$n <- nsites

### WRITE OUT SITE FILE
## TODO: do this upstream and make site file central protein level file
write.table(site, file=file.path(out.path, "aas_sites.tsv"), sep="\t",
            na="", row.names=FALSE, quote=FALSE)

if ( interactive() ) { # inspect some proteins
    site[grep("PSMA1", site$name),]
}


## list sites per protein

## filter here for only proteins>1, but not used!
multip <- split(site$pos, site$mane)
multip <- names(lengths(multip)[lengths(multip)> 0 ])
sites <- site[site$mane%in%multip,]
cpos <- cumsum(sites$pos)

## AAS and RAAS along all concatenated proteins

intv <- findInterval(sites$RAAS.median, vec=vbrks)
## raise to min/max
intv[intv==0] <- 1
intv[intv>length(vcols)] <- length(vcols)
raas.col <- vcols[intv]
plotdev(file.path(pfig.path,"hotspots_problem"), type=ftyp,
        width=10, height=3, res=200)
layout(t(t(1:3)), heights=c(.1,.4,.6))
par(mai=c(0.1,.5,.05,.1),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
## TODO: number of RAAS
par(mai=c(0,.5,.05,.1))
dense2d(cpos, rep(1, length(cpos)), axes=FALSE, ylab=NA, xlab=NA)
par(mai=c(0.1,.5,.01,.1))
plot(x=cpos, y=sites$RAAS.n, col=raas.col, pch=19, cex=.5,
     ylab="#RAAS per site",
     axes=FALSE, log="y", xpd=TRUE)
axis(2)
axis(3, at=cpos[sites$n==1], col=2, col.axis=2, labels=FALSE, tcl=.25)
par(mai=c(.35,.5,0.1,.1))
plot(cpos, sites$RAAS.median, type="h", col=raas.col, xpd=TRUE,
     xlab="pos. concatenated proteins with >1 AAS, sorted by protein RAAS",
     ylab=xl.site)
points(cpos, sites$RAAS.median, col=raas.col, pch=19, cex=.2, xpd=TRUE)
## TODO: protein lines
axis(3, at=cpos[sites$n==1], col=2, col.axis=2, labels=FALSE, tcl=.5)
dev.off()


### MEDIAN RAAS FOR EACH UNIQUE PROTEIN POSITION

## median RAAS vs. number of sites
if ( interactive() )
    dense2d(pbstat$median, log10(pbstat$n))

## volcano
plotdev(file.path(pfig.path,paste0("proteins_volcano_sites")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=pmpg, tcl=-.25)
res <- volcano(pbstat, value="median",
               p.txt=-log10(0.01), v.txt=c(-2,-2), cut=5, mid=-2,
               ids=pnms, xlab=xl.prots)
mtext("protein RAAS, median of site medians", 3,0)
dev.off()


## PROTEIN MEDIAN RAAS per protein w/o site-specific median first

plotdev(file.path(pfig.path,paste0("proteins_volcano_all")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=pmpg, tcl=-.25)
res <- volcano(ptstat, value="median",
               p.txt=12, v.txt=c(-2,-1), cut=50, mid=0,
               ids=pnms, xlab=xl.prota)
mtext("protein RAAS, median of all RAAS", 3,0)
dev.off()

## compare median of all RAAS vs. median of site median RAAS
if ( interactive() ) {
    plotCor(pbstat$median, ptstat$median,
            xlab="median of site median RAAS", ylab="median of all RAAS")
}

#### COMPARE PROTEIN RAAS TO VARIOUS PROTEIN LEVEL MEASURES
## TODO: collect those values for all proteins e.g. in halflives,
## and load here

### PROTEIN LENGTHS and AAS along proteins

## get protein length with saap_mapped4.tcv
plen <- split(hdat$len, hdat$ensembl)
plen <- lapply(plen, unique)
table(lengths(plen)) # check: all should be the same
plen <- unlist(lapply(plen, function(x) x[1]))

source("~/programs/segmenTools/R/plotUtils.R")
plotdev(file.path(pfig.path,paste0("protein_lengths_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(plen[rownames(ptstat)]), ptstat$median,
        ylim=range(pbstat$median), axes=FALSE,
        ylab=xl.prota, xlab="protein length", legpos="bottomright",
        legbg="#ffffff77", title=TRUE, cor.legend=FALSE)
axis(2)
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
##box()
dev.off()

## PROTEIN MELTING TEMPERATURE
therm <- as.data.frame(read_xlsx(thermo.file, sheet=2))
mlt <- therm$meltP_Jurkat
names(mlt) <- therm[,2]

## melting point vs. length
plotdev(file.path(pfig.path,paste0("protein_Tmelt_length")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(plen), mlt[match(pnms[names(plen)], names(mlt))],
        xlab="protein length",
        ylab=expression(protein~melting~point/"°C"),
        axes=FALSE)
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
axis(2)
##box()
dev.off()

## melting point all
plotdev(file.path(pfig.path,paste0("protein_Tmelt_all")),
        type=ftyp, res=300, width=corW,height=corH)
idx <- match(pnms[rownames(ptstat)], names(mlt))
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(mlt[idx], ptstat$median, ylim=range(pbstat$median),
        ylab=xl.prota, xlab=expression(protein~melting~point/"°C"),
        axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(1)
axis(2)
##box()
dev.off()

## Diff melting point with ATP
therm <- as.data.frame(read_xlsx(thatp.file, sheet=2))
mlt <- as.numeric(therm$diff_meltP_Exp1)
names(mlt) <- therm[,2]

## melting point difference with ATP
plotdev(file.path(pfig.path,paste0("protein_TmeltDelta_all")),
        type=ftyp, res=300, width=corW,height=corH)
idx <- match(pnms[rownames(ptstat)], names(mlt))
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(mlt[idx], ptstat$median, ylim=range(pbstat$median),
        ylab=xl.prota,
        xlab=expression(melting~point~difference~Delta*T[ATP-vehicle]/"°C"),
        axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(1)
axis(2)
##box()
dev.off()


## PROTEIN THERMOSTABILITY 
protstab <- read.csv(protstab.file)
## use refseq ID w/o version number as row names
rownames(protstab) <- sub("\\.[0-9]*","", protstab$id)

## get local reduced refseq mapping
## TODO: use full n<->n mapping and maximize yield
ps2ens <- refseq2ens[refseq2ens[,1] %in%rownames(pbstat),]
ps2ens <- ps2ens[ps2ens[,2] %in%rownames(protstab),]
protstab$ensembl <- ps2ens[match(rownames(protstab),ps2ens[,2]),1]
## match to RAAS table


## predicted stability 
plotdev(file.path(pfig.path,paste0("protein_protstab2_all")),
        type=ftyp, res=300, width=corW,height=corH)
idx <- match(rownames(ptstat), protstab$ensembl)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(protstab$Human_predict_Tm[idx], ptstat$median,
        ylim=range(pbstat$median),
        ylab=xl.prota, xlab="ProtStab2", axes=FALSE,
        title=TRUE, cor.legend=FALSE)
axis(1)
axis(2)
##box()
dev.off()

## TODO: thermostability and ATP/GTP from @Sridharan2019:
## Proteome-wide solubility and thermal stability profiling reveals
## distinct regulatory roles for ATP
## https://www.nature.com/articles/s41467-019-09107-y

## PROTEIN HALF-LIVES, @Mathieson2018
hlvd <- readxl::read_xlsx(math18.file)
##hlvd <- hlvd[,-grep("Mouse", colnames(hlvd))]

## mean half-live over all replicates and cell types
## TODO: consider distribution
cidx <- grep("half_life", colnames(hlvd), value=TRUE)
hlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
names(hlv) <- unlist(hlvd[,1])

plotCor(unlist(hlvd[,cidx[1]]), unlist(hlvd[,cidx[3]]))

xl.hlfm <- expression(protein~"half-life"/h)
xl.hlf <- expression(protein~"half-life"/h)

## halflives site
plotdev(file.path(pfig.path,paste0("protein_halflives_site")),
        type=ftyp, res=300, width=corW,height=corH)
idx <- match(pnms[rownames(pbstat)], names(hlv))
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(pbstat$median, log10(hlv[idx]),
        xlab=xl.prots, ylab=xl.hlfm, axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
##box()
dev.off()

## halflives all
plotdev(file.path(pfig.path,paste0("protein_halflives_all")),
        type=ftyp, res=300, width=corW,height=corH)
idx <- match(pnms[rownames(ptstat)], names(hlv))
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(hlv[idx]), ptstat$median, ylim=range(pbstat$median),
        ylab=xl.prota, xlab=xl.hlfm, axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(2)
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
##box()
dev.off()

## TODO: do this for all proteins and not just those
## in the RAAS set.
plotdev(file.path(pfig.path,paste0("protein_halflives_lengths")),
        type=ftyp, res=300, width=corW,height=corH)
idx <- match(pnms[names(plen)], names(hlv))
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(plen), log10(hlv[idx]),
        xlab="protein length", ylab=xl.hlfm, axes=FALSE)
axis(1, at=1:10, labels=10^(1:10))
axis(2, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
##box()
dev.off()

if ( interactive() ) {
    ## TODO: compare iupred3 and half-lives for all proteins, not just RAAS
    idx <- match(dat$name,names(hlv))
    plotCor(dat$iupred3.protein, log10(hlv[idx]))
}

for ( ctype in c("Bcells", "NK", "hepatocytes", "monocytes", "neurons",
                 "Mouse Neurons") ) {

    cidx <- grep("half_life", colnames(hlvd), value=TRUE)
    cidx <- cidx[grep(ctype, cidx, ignore.case=TRUE)]
    hlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
    names(hlv) <- unlist(hlvd[,1])
    
    plotdev(file.path(pfig.path,paste0("protein_halflives_sites_",
                                       sub(" ","_",ctype))),
            type=ftyp, res=300, width=corW,height=corH)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    idx <- match(pnms[rownames(pbstat)], names(hlv))
    plotCor(pbstat$median, log10(hlv[idx]),
            xlab=xl.prots, ylab=xl.hlf, axes=FALSE)
    axis(1)
    axis(2, at=1:10, labels=10^(1:10))
    figlabel(ctype, pos="bottomleft")
    dev.off()
    
    plotdev(file.path(pfig.path,paste0("protein_halflives_all_",ctype)),
            type=ftyp, res=300, width=corW,height=corH)
    idx <- match(pnms[rownames(ptstat)], names(hlv))
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plotCor(ptstat$median, log10(hlv[idx]), signif = 1,
        xlab=xl.prota, ylab=xl.hlf, axes=FALSE)
    axis(1)
    axis(2, at=1:10, labels=10^(1:10))
    figlabel(ctype, pos="bottomleft")
    dev.off()
## halflives site
    plotdev(file.path(pfig.path,paste0("protein_halflives_lengths_",ctype)),
            type=ftyp, res=300, width=corW,height=corH)
    idx <- match(pnms[names(plen)], names(hlv))
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plotCor(log10(plen), log10(hlv[idx]),
            xlab="protein length", ylab=xl.hlf, axes=FALSE)
    axis(1, at=1:10, labels=10^(1:10))
    axis(2, at=1:10, labels=10^(1:10))
    ## TODO: log tick marks
    axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)),
         tcl=-.125, labels=FALSE)
    axis(2, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)),
         tcl=-.125, labels=FALSE)
    ##box()
    figlabel(ctype, pos="bottomleft")
    dev.off()
}

### PROTEIN RAAS vs. 20S Targets - @Pepelnjak2024
## TODO: p20 data is on peptide level with many duplicated prots
    
p20 <- data.frame(readxl::read_xlsx(pepe24.file, sheet=2))

xl.20s <- expression(median~log[2]("20S"/control))
xl.20p <- expression("20S target:"~sign~"x"~log[10](p))

## MEDIAN LG2FC PER PROTEIN
p20s <- p20[!is.na(p20$Log2.FC.),]
p20l <- split(p20s[,c("Log2.FC.")],
              p20s[,"UniprotID"])
p20p <- unlist(lapply(p20l, median, na.rm=TRUE))

## median -log10(p)*sign
p20s$pscale <- sign(p20s$Log2.FC.) * -log10(p20s$p.value)
p20l <- split(p20s[,c("pscale")],
              p20s[,"UniprotID"])
p20pv <- unlist(lapply(p20l, median, na.rm=TRUE))


## get uniprot mapping
## TODO: akugb with uni2e and ens2u  in init script.
p2ens <- uni2ens
p2ens <- p2ens[p2ens[,2]%in%tmtf$mane,]
p2ens <- p2ens[p2ens[,1]%in%names(p20p),]

cat(paste(sum(duplicated(p2ens[,1])),"uniprot with multiple ensembl\n"))
cat(paste(sum(duplicated(p2ens[,2])),"uniprot with multiple ensembl\n"))
p2el <- unlist(split(p2ens[,2], p2ens[,1]))

e2u <- names(p2el)
names(e2u) <- p2el

if ( length(unique(lengths(p2el)))>1 )
    stop("non-unique uniprot 2 ensembl mapping")

##
lg2fc <- p20p
names(lg2fc) <- p2el[names(lg2fc)]
lg2fc <- lg2fc[rownames(pbstat)]
plotdev(file.path(pfig.path,paste0("protein_p20_lg2fc_site")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(pbstat$median, lg2fc, xlab=xl.prots,
        ylab=xl.20s, signif=2, round=2)
dev.off()
lg2fc <- p20p
names(lg2fc) <- p2el[names(lg2fc)]
lg2fc <- lg2fc[rownames(ptstat)]
plotdev(file.path(pfig.path,paste0("protein_p20_lg2fc_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(lg2fc, ptstat$median, ylim=range(pbstat$median),ylab=xl.prota,
        xlab=xl.20s, signif=2, round=2, legpos="bottomright",
        title=TRUE, cor.legend=FALSE)
dev.off()

lg2fc <- p20pv
names(lg2fc) <- p2el[names(lg2fc)]
lg2fc <- lg2fc[rownames(ptstat)]
plotdev(file.path(pfig.path,paste0("protein_p20_pval_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(lg2fc, ptstat$median, ylim=range(pbstat$median),ylab=xl.prota,
        xlab=xl.20p, signif=2, round=2, title=TRUE, cor.legend=FALSE)
dev.off()

## brute force
if ( FALSE ) {
    ## get uniprot mapping
    p2ens <- uni2ens
    p2ens <- p2ens[p2ens[,2]%in%tmtf$mane,]
    p2ens <- p2ens[p2ens[,1]%in%p20[,1],]
    
    cat(paste(sum(duplicated(p2ens[,1])),"uniprot with multiple ensembl\n"))
    cat(paste(sum(duplicated(p2ens[,2])),"uniprot with multiple ensembl\n"))
    p2el <- unlist(split(p2ens[,2], p2ens[,1]))
    p20$ensembl <- p2ens[match(p20[,1],p2ens[,1]),2]
    
    
    log2fc <- p20$Log2.FC.[match(rownames(pbstat), p20$ensembl)]
    logp <- -log10(p20$adj.pvalue[match(rownames(pbstat), p20$ensembl)])
    slogp <- logp*sign(log2fc)
    plotdev(file.path(pfig.path,paste0("p20_peptide_volcano")),
            type=ftyp, res=300, width=corW,height=corH)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    dense2d(log2fc, logp, xlab=xl.20s, ylab="-log10(p)")
    dev.off()
    
    plotdev(file.path(pfig.path,paste0("p20_peptide_protein_raas_pval_all")),
            type=ftyp, res=300, width=corW,height=corH)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plotCor(pbstat$median, slogp, xlab=xl.prots,
            ylab="p20/tryptic: sign(lg2fc)*-log10(p)")
    dev.off()
    
    plotdev(file.path(pfig.path,paste0("p20_peptide_protein_raas_lg2fc_all")),
            type=ftyp, res=300, width=corW,height=corH)
    par(mai=pmai, mgp=pmpg, tcl=-.25)
    plotCor(pbstat$median, log2fc, xlab=xl.prots,
            ylab=xl.20s)
    dev.off()
}



### PROTEIN RAAS vs. IUPRED

## get whole protein mean iupred3 score
iu3 <- split(hdat$iupred3.protein, hdat$ensembl)
## QC: all protein level
table(unlist(lapply(iu3, function(x) length(unique))))
iu3 <- unlist(lapply(iu3, unique))

plotdev(file.path(pfig.path,paste0("protein_iupred3_all")),
        type=ftyp, res=300, width=corW, height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(iu3[rownames(ptstat)], ptstat$median, ylim=range(pbstat$median),
        ylab=xl.prota, xlab="disordered score, IUPred3",
        axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(1)
axis(2)
dev.off()
plotdev(file.path(pfig.path,paste0("protein_iupred3_site")),
        type=ftyp, res=300, width=corW,height=corH)
plotCor(pbstat$median, iu3[rownames(pbstat)], signif = 1,
        legpos="topright", 
        xlab=xl.prota,
        ylab="disordered score, IUPred3")
dev.off()

