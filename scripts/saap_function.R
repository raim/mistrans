
## TODO:
## * enrichment analysis using ONLY avaible proteins
##   as background.
## * implement SetRank R package in segmenTools,
## https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1571-6

## * re-calculate mean RAAS from raw data where required!

library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)
library(gprofiler2)

proj.path <- "/home/raim/data/mistrans"
mam.path <- file.path("/home/raim/data/mammary")
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")
fig.path <- file.path(proj.path,"figures","saap_function")
dir.create(fig.path)

feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")
goslim.file  <- file.path(mam.path,"processedData","goslim.tsv")

## READ MTP TABLE
## output from map_peptides.R
dat <- read.delim(file.path(out.path,"saap_mapped.tsv"))

## FILTERED TABLES

dat$aacodon <- paste(dat$from,dat$codon, sep="-")

hdat <- dat[!dat$remove,]
hdat <- hdat[hdat$codon!="" & !is.na(hdat$RAAS.median),]

### REDUCE TABLE TO UNIQUE HUMAN PROTEINS

## * TODO: calculate mean or median RAAS on protein level here!
##   from all SAAP of this protein

mids <- sub("_.*", "", hdat$protein)
mids.sze <- table(mids) # sort(table(mids), decreasing=TRUE)

## NOTE: many mutations in expected, eg. ENSP00000474524
## "V region of the variable domain of immunoglobulin heavy chains
## that participates in the antigen recognition"

## analyze RAAS distribution of multiply mutated/SAAPd
raas.col <- "Mean_precursor_RAAS"
raas <- as.numeric(hdat[,raas.col])
## No Inf in current SAAP data
##raas[raas==-Inf] <- min(raas[is.finite(raas)]) -1
##raas[raas== Inf] <- max(raas[is.finite(raas)]) +1

hist(raas)

raas.lst <- split(raas, f=mids)
raas <- cbind(##n=mids.sze,
              n=unlist(lapply(raas.lst, length)),#sum(raas$n!=raas$n2)==0
              mean=unlist(lapply(raas.lst, mean)),
              sd=unlist(lapply(raas.lst, sd)))
raas <- as.data.frame(raas)

head(raas$n)

### UNIQUE PROTEIN MATRIX

## ADD GENE DESCRIPTIONS
##desc <- sub("\\[Homo sapiens\\]","",sub(".*\\|","",dat$Blast_ref_protein))
##desc.lst <- split(desc, f=mids)
##desc.lst <- lapply(desc.lst, function(x) trimws(x[x!=""]))
##desc.lst <- lapply(desc.lst, unique)

## check unique? NOTE: several with non-unique ensembl<->refseq mapping
## i.e. Blast_ref_protein vs. Top_leading_protein
##table(unlist(lapply(desc.lst, length)))
##which(unlist(lapply(desc.lst, length))==2)

##if ( any(as.numeric(names(table(unlist(lapply(desc.lst,length)))))>1) )
##    warning("non-unique mapping between ensembl and refseq")
##desc <- unlist(lapply(desc.lst, function(x) x[1]))

## TODO: collect more info here,
## * protein length

## ADD PROTEIN LENGTH
len <- hdat$len
len.lst <- split(len, f=mids)
len.lst <- lapply(len.lst, function(x) trimws(x[x!=""]))
len.lst <- lapply(len.lst, unique)

## check unique? NOTE: several with non-unique ensembl<->refseq mapping
## i.e. Blast_ref_protein vs. Top_leading_protein
table(unlist(lapply(len.lst, length)))
which(unlist(lapply(len.lst, length))==2)

if ( any(as.numeric(names(table(unlist(lapply(len.lst,length)))))>1) )
    warning("non-unique mapping between ensembl and refseq")

len <- as.numeric(unlist(len.lst))

## SAAP density
dns <- raas$n/len

## RAAS bins
brks <- c(-5,-3,-2,-1,0,1,4)
raas_bins <- cut(raas$mean, brks)

raas <- cbind.data.frame(raas, length=len, density=dns, bins=raas_bins)

## rm NA
naras <- is.na(raas$mean)
raas <- raas[!naras,]
cat(paste("removed", sum(naras), "proteins where", raas.col,"is NA, now:",nrow(raas),"\n"))

## sort by SAAP/protein
raas <- raas[order(unlist(raas$density), decreasing=TRUE),]

## plot SAAP/protein distribution
rn <- unlist(raas$n)
rn[rn>10] <- 11
## counts
rcnt <- table(rn)
## percent of total SAAPs
rprc <- rcnt/nrow(hdat)
## percent of total proteins
pprc <- rcnt/nrow(raas)

png(file.path(fig.path,"saap_per_protein.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(rcnt, xlab="unique SAAP/protein",
        ylab="# proteins")
legend("topright",c(paste0("unique SAAP: ", nrow(hdat)),
                    paste0("unique proteins: ", nrow(raas))))
dev.off()

png(file.path(fig.path,"saap_per_protein_percent_saap.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(rprc, xlab="unique SAAP/protein",
        ylab="% of total unique SAAP")
legend("topright",c(paste0("unique SAAP: ", nrow(hdat)),
                    paste0("unique proteins: ", nrow(raas))))
dev.off()
png(file.path(fig.path,"saap_per_protein_percent_protein.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(pprc, xlab="unique SAAP/protein",
        ylab="% of total unique proteins")
legend("topright",c(paste0("unique SAAP: ", nrow(hdat)),
                    paste0("unique proteins: ", nrow(raas))))
dev.off()

## plot RAAS per protein variance
## TODO: understand RAAS per protein distribution or get good value
## test whether correlating
plotdev(file.path(fig.path,"saap_raas_mean_sd"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(raas$mean, raas$sd, xlab="mean RAAS", ylab="SD RAAS", na.rm=TRUE,
        colf=viridis::viridis)
dev.off()

mx <- max(abs(c(raas$sd,raas$mean)), na.rm=TRUE)
brks <- seq(-mx-1, mx+1, 1)
hist(raas$mean, breaks=brks)
hist(raas$mean[raas$n>1], breaks=brks, col=2, add=TRUE)
hist(raas$sd[raas$n>1], xlab="standard deviation of RAAS", breaks=brks,
     border=4, add=TRUE)

hist(log10(raas$density), breaks=100,
     xlab=expression(log[10](SAAP/AA)))

## RAAS vs. SAAP DENSITY
plotdev(file.path(fig.path,"saap_raas_density"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(raas$mean, log10(raas$density), xlab="mean RAAS",
        ylab=expression(log[10](SAAP/AA)), colf=viridis::viridis)
dev.off()

plotdev(file.path(fig.path,"saap_raas_density_bins"),
        height=3.5, width=3.5, res=200)
par(mai=c(.75,.5,.5,.1), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(log10(raas$density) ~ raas$bin, ylab=expression(log[10](SAAP/AA)),
        xlab=NA, las=2, at=seq_along(levels(raas$bin)))
axis(3, at=seq_along(levels(raas$bin)), labels=table(raas$bin),las=2)
mtext("mean RAAS", 1, 2.5)
dev.off()

plotdev(file.path(fig.path,"saap_raas_number"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(raas$mean, log10(raas$n), xlab="mean RAAS", ##cor.method="spearman",
        ylab=expression(log[10](SAAP/protein)), colf=viridis::viridis)
dev.off()

plotdev(file.path(fig.path,"saap_raas_number_bins"),
        height=3.5, width=3.5, res=200)
par(mai=c(.75,.5,.5,.1), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(log10(raas$n) ~ raas$bin, ylab=expression(log[10](SAAP/protein)),
        xlab=NA, las=2, at=seq_along(levels(raas$bin)))
axis(3, at=seq_along(levels(raas$bin)), labels=table(raas$bin),las=2)
mtext("mean RAAS", 1, 2.5)
dev.off()

plotdev(file.path(fig.path,"saap_length_number_log"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(log10(raas$len), log10(raas$n),
        xlab=expression(log[10](protein~length)),
        ylab=expression(log[10](SAAP/protein)),
        colf=viridis::viridis)
dev.off()


plotdev(file.path(fig.path,"saap_raas_length"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(raas$mean, raas$len, xlab="mean RAAS",
        ylab="protein length/AA", colf=viridis::viridis)
dev.off()

plotdev(file.path(fig.path,"saap_raas_length_log"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(raas$mean, log10(raas$len), xlab="mean RAAS",
        ylab=expression(log[10](protein~length)), colf=viridis::viridis)
dev.off()

plotdev(file.path(fig.path,"saap_raas_length_bins"),
        height=3.5, width=3.5, res=200)
par(mai=c(.75,.5,.5,.1), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(log10(raas$len) ~ raas$bin, ylab=expression(log[10](protein~length)),
        xlab=NA, las=2, at=seq_along(levels(raas$bin)))
axis(3, at=seq_along(levels(raas$bin)), labels=table(raas$bin),las=2)
mtext("mean RAAS", 1, 2.5)
dev.off()


## ENRICHMENT TESTS
go.path <- file.path(fig.path,"annotation")
dir.create(go.path, showWarnings=FALSE)

## categorizations

rn <- raas$n
rn[rn>3] <- ">3" # implict cast to char
cls.mat <- cbind(number=rn)
                
rs <- raas$mean
rs[abs(rs)>10] <- 10 * sign(rs[abs(rs)>10])
brks <- c(-5,-3,-2,-1,0,1,4)
bins <- cut(raas$mean, breaks=brks)
bn <- as.character(bins)
bn[is.na(bn)] <- "na"

cls.mat <- cbind(number=rn, raas=bn)
cls.srts <- list(number=c("1","2","3",">3"), raas=levels(bins))
cls.labs <- c(number="SAAP/protein", raas="mean RAAS")
rownames(cls.mat) <- rownames(raas)

## discrete SAAP/density vs. RAAS
ovl <- clusterCluster(cls.mat[,1], cls.mat[,2],
                      cl1.srt=cls.srts[[1]], cl2.srt=cls.srts[[2]])
## calculate optimal figure height: result fields + figure margins (mai)
nh <- nrow(ovl$p.value) *.3 + 1
nw <- ncol(ovl$p.value) *.4 + 1

## RAAS vs. SAAP densitiy
plotdev(file.path(fig.path,"saap_raas_density_discrete"),
        height=nh, width=nw, res=200)
par(mai=c(.75,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovl, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=cls.labs[1])
mtext(cls.labs[2], 1, 2.5)
dev.off()

## GOslim with segmenTools

## get feature file with GO
features <- read.delim(feature.file)

## search saap proteins in features$protein and filter for MANE
gidx <- sapply(rownames(raas), function(x) grep(x, features$proteins))

## check multiples?
## table(unlist(lapply(gidx, length)))
gidx <- unlist(gidx)

## NOTE: very few (26) are non-MANE
table(features$MANE[gidx]=="")

## get ALL MANE + 26 non-MANE
midx <- which(features$MANE!="")
genes <- features[sort(union(gidx,midx)),]

gidx <- sapply(rownames(raas), function(x) grep(x, genes$proteins))
gidx <- unlist(gidx)
gen.cls <- matrix("na", ncol=ncol(cls.mat), nrow=nrow(genes))
gen.cls[gidx, ] <- cls.mat
colnames(gen.cls) <- colnames(cls.mat)


## get GOslim table for use with clusterAnnotation
got <- parseAnnotationList(genes[,c("ID","GOslim")]) 
## replace GO IDs by terms
terms <- read.delim(goslim.file)
trms <- terms[,2]
names(trms) <- terms[,1]
colnames(got) <- trms[colnames(got)]

## REDUCE AGAIN to use only available proteins with RAAS
## TODO: define proper background set!
local <- FALSE
if ( local ) {
    gen.cls <- gen.cls[gidx,]
    got <- got[gidx,]
}


## calculate enrichment!
for ( i in 1:ncol(gen.cls) ) {

    ct <- colnames(gen.cls)[i]
    ovl <- clusterAnnotation(cls=gen.cls[,i],
                             data=got, cls.srt=c(cls.srts[[i]],"na"))
    ovc <- sortOverlaps(ovl, axis=2, p.min=1e-5, cut=TRUE)
    
    ## calculate optimal figure height: result fields + figure margins (mai)
    nh <- .2; nw <- .4
    if ( nrow(ovc$p.value) < 5 ) {
        nh <- nh+.1; nw <- nw+.1
    }
    nh <- nrow(ovc$p.value) *nh + 1
    nw <- ncol(ovc$p.value) *nw + 5
    
    ## plot figure
    plotdev(file.path(go.path,paste0(ct, "_goslim")),
            height=nh, width=nw, res=200)
    par(mai=c(.75,4.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
    
    plotOverlaps(ovc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
                 xlab=NA, ylab=NA)
    mtext(cls.labs[ct], 1, 1, adj=-1.5, cex=1.2, font=2)
    dev.off()
}


### GPROFILER2

## TODO: look for manually defined functions over all annotations
my.functions <-
    list(immune=c("antigen binding","immunoglobulin","immune response"),
         extracellular=c("external", "extracellular",
                         "vesicle", "junction", "adhesion", "secretory"),
         cytoskeleton=c("actin[^g]","myosin","tubul","adhesion","fiber",
                        "polymer","filament"),
         metabolism=c("glycoly", "gluconeogene", "proteolysis", "catabolism",
                      "vacuol", "lysos", "proteaso", "translation", "ribosom"),
         organelle=c("mitochon", "organell","vacuol", "lysos"))

for ( ct in 1:ncol(cls.mat) ) {

    ## get clustering
    cid <- colnames(cls.mat)[ct]
    cls <- cls.mat[,ct]

    cl.srt <- cls.srts[[ct]]

    cat(paste("CALCULATING ANNOTATION ENRICHMENTS for", cid, "\n"))
    

    ## ENRICHMENT OVER ALL CATEGORIES in gprofiler2
    ovll <- runGost(cls, organism="hsapiens", cls.srt=cl.srt, evcodes=FALSE,
                    custom_bg=rownames(got), categories="KEGG")
###, significant=FALSE, evcodes=FALSE)
    
    ## plot enrichments
    for ( ctgy in names(ovll) ) {
    
        cat(paste("\tPLOTTING for", ctgy, "\n"))
        ovl <- ovll[[ctgy]]

        ## filter large annotations such as GO
        p.filt <- 1e-2
        cut <- TRUE
        if ( length(grep("^GO:",ctgy))>0 ) {

            ## 
            trm.srt <- ovl$num.query[,1]
            keep <- trm.srt > 10 & trm.srt <5e3
            trm.srt <- names(trm.srt)[keep]
            ovl <- sortOverlaps(ovl, srt=trm.srt)
            
            p.filt <- 1e-5
            cut <- TRUE
        } else if ( ctgy=="TF" ) {
            p.filt <- 1e-5
            cut <- TRUE
        } else if ( ctgy=="KEGG" ) {
            p.filt <- 1e-5
            cut <- TRUE
        }

        ## sort along clusters
        ovc <- sortOverlaps(ovl, p.min=p.filt, cut=cut)

        if ( nrow(ovc$p.value)==0 ) {
            cat(paste("NO SIGNIFICANT HITS REMAINING\n"))
            next
        }
        
       
        ## calculate optimal figure height: result fields + figure margins (mai)
        nh <- nrow(ovc$p.value) *.2 + 1
        nw <- ncol(ovc$p.value) *.4 + 5

        if ( nrow(ovc$p.value) < 5 ) {
            nh <- nh*2; nw <- nw*2
        }
        
        ## plot figure
        plotdev(file.path(go.path,paste0(cid,"_annot_",ctgy)),
                height=nh, width=nw, res=200)
        par(mai=c(.75,4.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        plotOverlaps(ovc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
                     xlab=NA, ylab=NA)
        mtext(cls.labs[ct], 1, 1, adj=-1.5, cex=1.2, font=2)
        axis(2, at=nrow(ovc$p.value)+1.5, label=ctgy,
             cex.axis=1.2, font=2, xpd=TRUE, las=2)
        dev.off()
    }
    
    ## ENRICHMENT OVER PRE-SELECTED TERMS in gprofiler2
    cat(paste("CALCULATING ENRICHMENTS for pre-selected functions\n"))
    ovll <- runGost(cls, organism="hsapiens", cls.srt=cl.srt,
                    terms=my.functions, evcodes=FALSE,
                    custom_bg=rownames(got)) 
    ## plot enrichments
    for ( ctgy in names(ovll) ) {
    
        cat(paste("\tPLOTTING for", ctgy, "\n"))

        ovl <- ovll[[ctgy]]

        ## filter large annotations such as GO
        p.filt <- 1e-3
        cut <- TRUE
        
        ## sort along clusters
        ovc <- sortOverlaps(ovl, p.min=p.filt, cut=cut)

        if ( nrow(ovc$p.value)==0 ) {
            cat(paste("NO SIGNIFICANT HITS REMAINING\n"))
            next
        }
        
         
        ## calculate optimal figure height: result fields + figure margins (mai)
        nh <- nrow(ovc$p.value) *.2 + 1.25
        nw <- ncol(ovc$p.value) *.4 + 5
        
        if ( nrow(ovc$p.value) < 5 ) {
            nh <- nh*2; nw <- nw*2
        }

        ## plot figure
        plotdev(file.path(go.path,paste0(cid,"_annot_",ctgy)),
                height=nh, width=nw, res=200)
        par(mai=c(.75,4.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        plotOverlaps(ovc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
                     xlab=NA, ylab=NA)
        mtext(cls.labs[ct], 1, 1, adj=-1.5, cex=1.2, font=2)
        axis(2, at=nrow(ovc$p.value)+1.5, label=ctgy,
             cex.axis=1.2, font=2, xpd=TRUE, las=2)
        dev.off()
    }
}

