
## read protein fasta, iupred3, s4pred, and AAS,
## and generate per protein plots

## TODO:
## * GET AND LOAD PFAM ANNOTATION: descriptions, clans!
## * use gene names

## TODO: why is iupred3 for ENSP00000352639 wrong?
## sequence seems to stem from a differen short protein,
## multiple annotated proteins, e.g. ENSP00000464724.1
## TODO: why is iupred3 for ENSP00000364986 wrong?
## sequence seems to stem from ENSP00000460206.1
## TODO: ENSP00000374778 et al. stop codon missing?

source("~/work/mistrans/scripts/saap_utils.R")
source("~/programs/genomeBrowser/src/genomeBrowser_utils.R")

library(viridis)
library(segmenTools)
library(seqinr)
library(Rpdb)
options(stringsAsFactors=FALSE)


#### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")
fig.path <- file.path(proj.path,"figures","proteins")
dir.create(fig.path)
dir.create(file.path(fig.path, "selected"))

mam.path <- "~/data/mammary/"

## AAS 
in.file <- file.path(out.path,"saap_mapped3.tsv")
## RAAS values
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")
## protein fasta
pfasta <- file.path(out.path,"all_proteins.fa")
## coding region fasta
tfasta <- file.path(mam.path,"processedData","coding.fa")

## protein-transcript map
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
## CDS structure
cdsmap <- file.path(mam.path,"processedData","protein_cds_structure.dat")
cdspos <- file.path(mam.path,"processedData","protein_cds_coordinates.tsv")


## structure predictions
## s4pred
s4pred <- file.path(mam.path,"processedData",
                    "Homo_sapiens.GRCh38.pep.large_s4pred.fas.gz")
## iupred3/anchor2
iupred <- file.path(mam.path,"processedData","iupred3")

## pfam hits
pfam <- file.path(mam.path,"processedData",
                  "Homo_sapiens.GRCh38.pep.large_annotations.csv")
## pfam clans
clans <- file.path(mam.path,"originalData", "pfam", "Pfam-A.clans.tsv.gz")

## genome feature file
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

## PARAMETERS

## plot
ftyp <- "png"

## only plot subset of substitutions with high RAAS
min.raas <- -1

## raas color scheme
RAAS.MIN <- -4
RAAS.MAX <-  1
colors <- "arno" # "inferno" #"rocket" # "viridis" # 

COLF <- get(colors)

## PARSE & FILTER DATA

## proteins
genes <- read.delim(feature.file)

## AAS
dat <- read.delim(in.file)

## remove SAAP/BP w/o protein
rm <- dat$ensembl==""
dat <- dat[!rm,]

## TMT Level RAAS Values
tmtf <- read.delim(tmt.file)

## exclude NA or Inf
## NOTE: since we don't to tests against the global
## distribution and only calculate protein-specific
## RAAS means, we don't need to filter RAAS data here.
rm <- is.na(tmtf$RAAS) | is.infinite(tmtf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from TMT level\n"))
tmtf <- tmtf[!rm,]

## ADD RAAS STATS

## split RAAS by BP/SAAP
tmtl <- split(tmtf$RAAS, paste(tmtf$BP, tmtf$SAAP))

## match to main data
tmt2dat <- match(paste(dat$BP, dat$SAAP), names(tmtl))
tmtl <- tmtl[tmt2dat]

## RAAS statistics
tmts <- lapply(tmtl, function(x) {
    x <- 10^x
    mn <- mean(x)
    md <- median(x)
    sd <- sd(x)
    cv <- sd/mn
    c(n=length(x),
      mean=log10(mn),
      median=log10(md),
      sd=log10(sd),
      cv=cv)
})
tmts <- as.data.frame(do.call(rbind, tmts))

## RAAS COLORS
png(file.path(fig.path,paste0("legend_raas_proteins.png")),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(tmtf$RAAS,
                          mn=RAAS.MIN, mx=RAAS.MAX,colf=COLF,
                          n=50, plot=TRUE,
                          mai=c(.5,.5,.1,.1),
                          xlab=expression(TMT~level~log[10]*RAAS))
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(colors, pos="bottomleft", cex=1)
##figlabel(LAB, pos="bottomright", cex=1)
dev.off()

## add to main data
dat <- cbind(dat, tmts)

## remove SAAP without RAAS
dat <- dat[dat$n>0,]

## add RAAS coloring
intv <- findInterval(dat$median, lraas.col$breaks)
## raise to min/max
intv[intv==0] <- 1
intv[intv>length(lraas.col$col)] <- length(lraas.col$col)
dat$color <- lraas.col$col[intv]

## remove SAAP with low RAAS
##dat <- dat[dat$median >= min.raas,]

## LIST OF AAS BY PROTEIN
aasl <- split(dat, dat$ensembl) 

## TODO: count and raas-profiles per AAS

## s4pred
s4p <- segmenTools::readFASTA(s4pred, grepID=TRUE)
## skip version number - TODO: use version numbers!
names(s4p) <- sub("\\.[0-9]+", "", names(s4p))

## iupred3
iufiles <- list.files(pattern=paste0(".*iupred3.tsv.gz"), path=iupred)
names(iufiles) <- sub("\\..*","", iufiles)

## load transcript and protein fasta
## GET ENSEMBL PROTEINS - from project mammary
pfas <- segmenTools::readFASTA(pfasta, grepID=TRUE)

## get matching transcripts
tfas <- segmenTools::readFASTA(tfasta, grepID=TRUE)

## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
## reverse map transcript-protein
pamrt <- matrix(rownames(trmap), ncol=1)
rownames(pamrt) <- trmap[,1]

## rename by protein names via trmap, for quick access
## of protein-specific transcript
names(tfas) <- pamrt[names(tfas),1]



## pfam
## fixed width format, where \s was replaced with ; by sed
## NOTE: three sets of from/to: hmm, ali, env;
## env may overlap between multiple - could be fused!
pfm <- read.csv(file=pfam, sep=";", fill=FALSE, header=FALSE,comment.char="#")
pfmh <- c(
    "target", "accession" , "tlen",            
    "query"   , "accession" , "qlen",                 
    "E-value"      , "score"     , "bias",                 
    "#"            , "of"        , "c-Evalue",             
    "i-Evalue"     , "score"     , "bias",            
    "from"         , "to"        , "from",                 
    "to"           , "FROM"      , "TO",                   
    "acc")#          , "description")
colnames(pfm) <- pfmh

## pfam clans: used to collapse domains below!
pclan <- read.delim(clans, row.names=1 , header=FALSE)
pfm$clan <- pclan[sub("\\..*","", pfm$accession),"V3"]
pfm$clan[pfm$clan==""] <- pfm$target[pfm$clan==""]

use.pclan <- TRUE# FALSE
if ( !use.pclan )
    pfm$clan <- pfm$target

pfl <- split(pfm, pfm$query)
names(pfl) <- sub("\\.[0-9]+", "", names(pfl))

## AAS per protein
app <- unlist(lapply(aasl, nrow))
appc <- app
appc   [app>30 ] <- 30

hist(appc, breaks=1:30, xlab="AAS per protein")



## proteins of interest:
## TDP-43,
## proteasome: PSM*,
## TODO: glycolysis, AA metabolism,
## glycolysis: ENO1/2/3, GAPDH,
## histones:
## globin: HBA1/2, HBB,
## Ras/Rap: RAB*, RAP1A, RAP2A,
## collagen,
## CCR4-NOT: CNOT10,
## calmodulin: CALM1, CALML3
## TODO: RAAS profiles for protein complexes, 

POI <- c(pamrt[genes[genes$name=="TARDBP","canonical"],],
         ## proteasome 20S subunit alpha 1 
         pamrt[genes[genes$name=="PSMA1","canonical"],],
         ## actins
         pamrt[genes[genes$name=="ACTA1","canonical"],],
         pamrt[genes[genes$name=="ACTB","canonical"],],
         pamrt[genes[genes$name=="ACTC1","canonical"],])

## map gene names and uniprot
gidx <-  match(trmap[names(aasl),1], genes$canonical)
pnms <- genes$name[gidx]
pnms[is.na(pnms)] <- names(aasl)[is.na(pnms)]
names(pnms)  <- names(aasl)

puni <- genes$swissprot[gidx]
puni[is.na(puni)] <- names(aasl)[is.na(puni)]
names(puni)  <- names(aasl)

## plot all proteins INCL. QC
pids <- names(aasl)#POI #
for ( pid in pids ) {

    
    ffile <- file.path(fig.path, pnms[pid])
    if ( pid%in%POI )
        sfile <- file.path(fig.path, "selected", pnms[pid])

    ##if ( file.exists(ffile) ) next

    ## collect all data for protein
    aas <- aasl[[pid]]
    pf <- pfl[[pid]]
    s4 <- gsub("H","0",gsub("E","/",gsub("C","-",s4p[[pid]]$seq)))
    iu <- read.delim(file.path(iupred,iufiles[pid]), header=FALSE)
    psq <- pfas[[pid]]$seq
    tsq <- tfas[[pid]]$seq
    ## check lengths
    tlen <- nchar(tsq)/3 -1
    plen <- nchar(psq)
    slen <- nchar(s4)
    ilen <- nrow(iu)
    if ( length(unique(c(tlen,plen,slen,ilen)))>1 )
        cat(paste("WARNING:", pid, "lengths differ",
                  paste(unique(round(c(tlen,plen,slen,ilen),1)),
                        collapse=";"), "\n"))
    ## check sequence via translate
    if ( paste0(iu[,2],collapse="")!=psq)
        cat(paste("WARNING:", pid, "wrong iupred3 seq\n"))
    if ( is.null(tsq) ) {
        cat(paste("WARNING:", pid, "no transcript\n"))
    } else {
        tpsq <- sub("\\*$","",
                    paste0(translate(strsplit(tsq,"")[[1]]),collapse=""))
        if ( tpsq!=psq )
            cat(paste("WARNING:", pid, "wrong translation\n"))
    }

    ## get domains
    pf <- pfl[[pid]]

    ## convert all to genomeBrowser style tables
    ## and use genomeBrowser plot functions.

    ## dummy plot coordinates for protein length
    coors <- c(chr=1,start=0, end=plen+1)
   
    ## fuse domains of same type
    ## using segmenTools segmentMerge via bedtools
    if ( !is.null(pf) ) {
        pf <- pf[,c("clan","FROM","TO","E-value")]
        colnames(pf) <- c("type","start","end","score")
        pf <- pf[order(pf$start),]
        pf <- cbind(pf,chr=1,strand=".")
        pf <- segmentMerge(pf, type="type", verb=0)
        pf$strand <- "."
        pf$name <- pf$type
    }

    ## PREPARE AAS PLOT
    ## order by RAAS such that higher RAAS will be on top of overlapping
    aas <-
        aas[order(aas$median),]
    aaco <- paste0(aas$from,":",aas$to)
    ## TODO: type 0,1,2,3 for closeby, not just identical
    aatp <- as.numeric(duplicated(aas$pos))
    raas <- aas$median
    aad <- data.frame(type=aaco,#aatp,
                      name=aaco,
                      start=aas$pos,
                      end=aas$pos,
                      chr=1, strand=".",
                      size=log(aas$n)+.5,
                      color=aas$color,
                      codon=aas$codon)
    aad <- aad[aas$median >= min.raas, ]

    ## skip plot if no RAAS is higher than minimum
    if ( nrow(aad)==0 & !pid%in%POI ) next

    cat(paste("PLOTTING", pnms[pid], pid, "\n"))

    plotdev(ffile, width=min(c(100,plen/30)), height=3, type=ftyp,
            res=200)
    layout(mat=rbind(1,2,3,4,5), heights=c(.5,.1,.1,.25,.1))
    par(mai=c(.05,.5,.05,.5), xaxs="i", xpd=TRUE)
    
    ## domains as arrows
    if ( !is.null(pf) ) {
        plotFeatures(pf, coors=coors, tcx=1.5, names=TRUE)
    } else plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
    text(x=plen/2, y=1, labels=paste(pnms[pid],"/",puni[pid]), cex=2)
    
    ## domains as blocks
    ##plotFeatureBlocks(pf, coors=coors)
    
    ## data as heatmap
    iud <- cbind(chr=1,coor=iu[,"V1"],vals=iu[,c("V3","V4")])
    brks <- 50:100/100
    plotHeat(iud, coors=coors, breaks=brks, colors=c(viridis(length(brks)-1)))
    
    plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
    text(1:plen, y=1, labels=strsplit(s4,"")[[1]], cex=1)
    ## indicate K|R cleavage sites - TODO: omit Keil rules
    arg <- which(strsplit(psq,"")[[1]]%in%c("R"))
    if ( length(arg) )
        arrows(x0=arg, y0=1.5, y1=1, length=.05, lwd=1)
    lys <- which(strsplit(psq,"")[[1]]%in%c("K"))
    if ( length(lys) )
        arrows(x0=lys, y0=1.5, y1=1, length=.05, lwd=1.5, col=2)

    ## indicate ALL AAS
    arrows(x0=aas$pos, y0=0, y1=1, length=.05, lwd=3)
    arrows(x0=aas$pos, y0=0, y1=1, length=.05, lwd=2.5, col=aas$color)
    ##abline(v=aas$pos, col=aas$color) ## TODO: arrows, color by RAAS etc.
    ## DETAILS for high RAAS AAS
    if (FALSE) { # custom with n~size
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        for ( tp in unique(aad$type) )
            shadowtext(x=aad$start[aad$type==tp],
                       y=1,#-rep(tp,sum(aad$type==tp))*.25,
                       labels=aaco[aad$type==tp],
                           cex=aad$size[aad$type==tp],
                       col=aad$color[aad$type==tp])
    }
    plotFeatures(aad, coors=coors, names=TRUE, tcx=2,
                 typord=TRUE, axis2=FALSE)
    axis(1, at=aad$start, labels=FALSE, col=NA)
    axis(1, tcl=-.1, mgp=c(1,.2,0))
    ## AA context
    
    ## codons
    plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
    if ( nrow(aad)>0 )
        shadowtext(x=aad$start, y=rep(1, nrow(aad)),
                   labels=aad$codon, cex=2, col=aad$color)
    dev.off()
    if ( pid%in%POI )
        file.copy(paste0(ffile,".",ftyp),
                  paste0(sfile,".",ftyp))
}
 

