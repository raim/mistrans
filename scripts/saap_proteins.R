
## read protein fasta, iupred3, s4pred, and AAS,
## and generate per protein plots

## TODO: GET AND LOAD PFAM ANNOTATION: descriptions, clans!

source("~/work/mistrans/scripts/saap_utils.R")

library(viridis)
library(segmenTools)
library(seqinr)
options(stringsAsFactors=FALSE)


#### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")

mam.path <- "~/data/mammary/"

in.file <- file.path(out.path,"saap_mapped3.tsv")
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

## genome feature file
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

## PARSE & FILTER DATA

## proteins
genes <- read.delim(feature.file)

## AAS
dat <- read.delim(in.file)
## remove SAAP/BP w/o protein
rm <- dat$ensembl==""
dat <- dat[!rm,]
aasl <- split(dat, dat$ensembl) 

## TMT Level RAAS Values
tmtf <- read.delim(tmt.file)

## s4pred
s4p <- readFASTA(s4pred, grepID=TRUE)
## skip version number - TODO: use version numbers!
names(s4p) <- sub("\\.[0-9]+", "", names(s4p))

## iupred3
iufiles <- list.files(pattern=paste0(".*iupred3.tsv.gz"), path=iupred)
names(iufiles) <- sub("\\..*","", iufiles)

## load transcript and protein fasta
## GET ENSEMBL PROTEINS - from project mammary
pfas <- readFASTA(pfasta, grepID=TRUE)

## get matching transcripts
tfas <- readFASTA(tfasta, grepID=TRUE)

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

pfl <- split(pfm, pfm$query)
names(pfl) <- sub("\\.[0-9]+", "", names(pfl))

## AAS per protein
app <- unlist(lapply(aasl, nrow))
appc <- app
appc   [app>30 ] <- 30

hist(appc, breaks=1:30, xlab"AAS per protein")

## look at TDP-43
pid <- pamrt[genes[genes$name=="TARDBP","canonical"],]

## select a nice protein
pid <- names(aasl)[which(app >= 5)[2]]


## plot all proteins INCL. QC
##for ( pid in names(aasl) ) {

    ## collect all data for protein
    aas <- aasl[[pid]]
    pf <- pfl[[pid]]
    s4 <- s4p[[pid]]$seq
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
## fuse domains of same type
## using segmenTools segmentMerge via bedtools
pf <- pf[,c("target","FROM","TO","E-value")]
colnames(pf) <- c("type","start","end","score")
pf <- pf[order(pf$start),]
pf <- cbind(pf,chr=1,strand=".")
pf <- segmentMerge(pf, type="type")
pf$strand <- "."
pf$name <- pf$type

source("~/programs/genomeBrowser/src/genomeBrowser_utils.R")

coors <- c(chr=1,start=0, end=plen+1)

layout(mat=rbind(1,2,3,4,5), heights=c(.2,.1,.1,.1,.1))
par(mai=c(.05,.5,.05,.5), xaxs="i")

## domains as arrows
plotFeatures(pf, coors=coors, names=TRUE)
axis(1)
axis(2)
## domains as blocks
##plotFeatureBlocks(pf, coors=coors)

## data as heatmap
iud <- cbind(chr=1,coor=iu[,"V1"],vals=iu[,c("V3","V4")])
plotHeat(iud, coors=coors, breaks=30:100/100)

plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE)
abline(v=aas$pos) ## TODO: arrows, color by RAAS etc.
plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE)
text(1:plen, y=1, labels=strsplit(s4,"")[[1]], cex=.5)
plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE)
text(1:plen, y=1, labels=strsplit(psq,"")[[1]], cex=.5)

##plot(1:plen, xlim=c(0,plen+1), col=NA, ylim=c(-2,2))
##rect(xleft = pf$start, xright = pf$to, ybottom = -1, ytop = 1)
##text(x=apply(cbind(pf$start,pf$to),1,mean), y=0, labels=pf$target)
}

## TODO: why is iupred3 for ENSP00000352639 wrong?
## sequence seems to stem from a differen short protein,
## multiple annotated proteins, e.g. ENSP00000464724.1
## TODO: why is iupred3 for ENSP00000364986 wrong?
## sequence seems to stem from ENSP00000460206.1
## TODO: ENSP00000374778 et al. stop codon missing?

