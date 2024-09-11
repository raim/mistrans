
### LOAD BP/SAAP MAPPING AND TMT LEVEL RAAS DATA
## for further downstream analyses, RAAS profiles by
## AA properties and codons (raasprofiles3_codons.R), and by
## protein complexes, proteins and protein windows
## (raasprofiles3_proteins.R).


## attempt to get the correct path for input data;
## these should be downloaded or generated as outlined
## in README.md
proj.path <- file.path(Sys.getenv("DECODE")) 
if ( proj.path=="" ) # author's local path
    proj.path <- "/home/raim/data/mistrans"
SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")

## CRAN packages
library(viridis) # viridis coloring scheme
library(readxl)
library(R.utils) # parsing gz files

## BioConductor packages
library(DiffLogo) # sequence diference logos

## via github, see install_packages.R
library(basicPlotteR) # for non-overlapping text
library(segmenTools) # plot utils and overlap profile sorting and plotting

## project-specific functions
source(file.path(SRC.PATH, "raas_utils.R"))

## R options
options(stringsAsFactors=FALSE) # default since R v4
options(scipen=0) # use e notation for p-values

#### PATHS AND FILES

## INPUT/OUTPUT PATHS


## input data
add.path <- file.path(proj.path, "additionalData")

## MAIN INPUT: MAPPED BP AND SAAP
## for publication with collected data
in.file <- file.path(add.path,"saap_mapped.tsv.gz")
tmt.file <- file.path(add.path,  "All_SAAP_TMTlevel_quant_df.txt.gz")

in.file <- file.path(proj.path,"processedData","saap_mapped.tsv")
##tmt.file <- file.path(proj.path, "originalData",
##                      "All_SAAP_TMTlevel_quant_df.txt")

## 20240906: positional probabilities
prob.file <- file.path(proj.path, "originalData",
                     "Supplemental_Data_3.SAAP_precursor_quant_w_fragments.xlsx")

## TODO: remove those from TMT level
## and afterwards remove all missing from  bdat
MINPP <- .9
RM.POSPROB <- FALSE # TRUE # 
MAXPP <- .5
MN.POSPROB <- TRUE # FALSE # 

if ( !file.exists(in.file) | !file.exists(tmt.file) )
    stop("INPUT FILES MISSING. MAKE SURE PATHS ARE DEFINED PROPERLY",
         "AND INPUT EXISTS, as outlined in README.md")


## output path for ALL scripts
out.path <- file.path(proj.path,"processedData")
fig.path <- file.path(proj.path,"figures")

if ( RM.POSPROB ) {
    fig.path <- file.path(fig.path, "posprob")
} else if ( MN.POSPROB ) {
    fig.path <- file.path(fig.path, "lowposprob")
}

## figures of this init script
ifig.path <- file.path(fig.path,"init")

## CREATE DIRECTORIES FOR OUTPUT
dir.create(out.path, recursive=TRUE, showWarnings=FALSE)
dir.create(ifig.path, recursive=TRUE, showWarnings=FALSE)




### SUPPLEMENTAL DATA

## @Mathieson2018: protein half-lives 
math18.file <- file.path(add.path, "41467_2018_3106_MOESM5_ESM.xlsx")

## @Wu2019: codon stability coefficient (CSC)
wu19.file <- file.path(add.path,"elife-45396-fig1-data2-v2.csv.gz")

## @Dana2014, Dana and Tuller 2014: decoding time
dana14.file <- file.path(add.path,"dana14_codons.csv")

## @McCormick2024pre: RNA pseudouridylation data by (bioRxiv) 
psi.file <- file.path(add.path, "six_cell_lines_minimal.xlsx")


### GENOME FEATURE TABLE and PROTEIN ID MAPPINGS
## via https://gitlab.com/raim/genomeBrowser/data/mammary/setup.sh

## uniprot/refseq/ensembl/name mappings
uni2ens.file <- file.path(add.path,"uniprot_ensembl.dat.gz")
uni2nam.file <- file.path(add.path,"uniprot_name.dat.gz")
refseq.file <- file.path(add.path, "ensembl_refseq_20240528.tsv.gz")

## gene name synonyms, via https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
## downloaded on 20240712
synonym.file <- file.path(add.path, "gene_synonyms.tsv.gz")

## codon frequencies, calculated in genomeBrowser
codon.file <- file.path(add.path,"coding_codons.tsv.gz")

## coding sequence fasta
tfas.file <- file.path(add.path, "coding.fa.gz")

## chromosome length index
chr.file <- file.path(add.path,"sequenceIndex.csv.gz")

## genome feature table 
feature.file <- file.path(add.path,"features_GRCh38.110.tsv.gz")
## tRNA file
trna.file <- file.path(add.path, "codons_GRCh38.tsv.gz")

## GOslim term definitions
goslim.file  <- file.path(add.path,"goslim.tsv.gz")



## DATA FILTERS

## NOTE: these labels refer to options in the original code
## leaving them here for now, because they are used in file names.
LAB <- "" # "all"
SETID <- "cancer"


### PARAMETERS

## TODO: fuse global and tight color scale?
RAAS.MIN <- -4   # broad RAAS colors vcols: codon plot
RAAS.MAX <-  1
RAAS.MINA <- -4  # tighter RAAS colors acols: all other dotplots
RAAS.MAXA <- 0

## RAAS color pallette
colors <- "arno" # "inferno" #"rocket" # "viridis" # 
COLF <- get(colors)

## p-value cutoffs for plots and filters
p.min <- 1e-10
p.txt <- 1e-5
## dotplot parameters
p.dot <- p.min # p.txt
dot.sze <- c(.3,2)

## STATISTICAL TEST TO RUN
use.test <- t.test # w.test # 

## to keep comparable, use t-test also for high pos prob set
##if ( RM.POSPROB ) ## NON-NORMAL DISTRIBUTION
##    use.test <- w.test

## SET THIS VARIABLE TO GENERATE PDFs INSTEAD OF PNGs
ftyp <- "png" # "pdf" # # 


## heatmap colors
docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))
gcols <- grey.colors(n=100, start=.9, end=0)

## RAAS axis labels
xl.raas <- expression(log[10](RAAS)) # *bar(RAAS))
xl.raaa <- expression(log[10](RAAS))
xl.raau <- expression(log[10]*bar(RAAS[unique]))
xl.site <- expression(log[10](RAAS[site]))
xl.all <- expression(log[10](RAAS[all]))
xl.prot <- expression(log[10](RAAS[protein]))
xl.prota <- expression(log[10](RAAS["protein,all"]))
xl.prots <- expression(log[10](RAAS["protein,site"]))

plab <- expression(log[10](p))



## AA PROPERTY CLASSES
## https://en.wikipedia.org/wiki/Amino_acid#/media/File:ProteinogenicAminoAcids.svg
AAPROP <- list(charged=c("R","H","K","D","E"),
               polar=c("S","T","N","Q"),
               special=c("C","U","G","P"),
               hydrophobic=c("A","V","I","L","M","F","Y","W"))
tmp <- unlist(AAPROP)
AAPROP <- sub("[0-9]+","", names(tmp))
names(AAPROP) <- tmp
aaprop <- sub("hydrophobic", "hphobic", AAPROP)

aaprop.cols <- c(charged="#00DCDC",
                 polar="#8282D2",
                 special="#33b864",
                 hphobic="#ffcc00",
                 hydrophobic="#ffcc00")
aaprop.cols["hydrophobic"] <- rgb(0.8745098,
                                  0.3254902,
                                  0.4196078)
aaprop.cols["hphobic"] <- aaprop.cols["hydrophobic"]

## copied from the Biostrings package
GENETIC_CODE <- c(
    TTT="F", TTC="F",         
    TTA="L", TTG="L",         
    TCT="S", TCC="S", TCA="S", TCG="S", AGT="S", AGC="S",         
    TAT="Y", TAC="Y",         
    TAA="*", TAG="*", TGA="*",
    TGT="C", TGC="C",    
    TGG="W",         
    CTT="L", CTC="L", CTA="L", CTG="L",         
    CCT="P", CCC="P", CCA="P", CCG="P",         
    CAT="H", CAC="H",         
    CAA="Q", CAG="Q",         
    CGT="R", CGC="R", CGA="R", CGG="R", AGA="R", AGG="R",         
    ATT="I", ATC="I", ATA="I",
    ATG="M",         
    ACT="T", ACC="T", ACA="T", ACG="T",         
    AAT="N", AAC="N",         
    AAA="K", AAG="K",         
    GTT="V", GTC="V", GTA="V", GTG="V",         
    GCT="A", GCC="A", GCA="A", GCG="A",         
    GAT="D", GAC="D",         
    GAA="E", GAG="E",         
    GGT="G", GGC="G", GGA="G", GGG="G"
)

###  CODONS
aa <- unique(GENETIC_CODE)
CODONS <- rep("", length(aa))
for ( i in seq_along(aa) )
    CODONS[i] <- paste(names(which(GENETIC_CODE==aa[i])), collapse=";")
names(CODONS) <- aa
ACODONS <- paste0(names(CODONS),": ", CODONS)
names(ACODONS) <- aa

## SORT CODONS BY AA PROPERTY
CODL <- strsplit(CODONS, ";")[names(aaprop)]
CODL <- CODL[unlist(lapply(CODL, function(x) !is.null(x)))]
CODL <- lapply(CODL, sort)
COD.SRT <- paste(Biostrings::GENETIC_CODE[unlist(CODL)], unlist(CODL), sep="-")

## CODON PWM
cod.pwm <- lapply(CODL, function(x) do.call(rbind,strsplit(x,"")))
nucs <- unique(unlist(cod.pwm))
nucm <- matrix(0,ncol=3, nrow=length(nucs))
rownames(nucm) <- nucs
CODP <- lapply(cod.pwm, function(x) {
    mat <- nucm
    for ( i in 1:3 ) {
        tb <- table(x[,i])/nrow(x)
        mat[names(tb),i] <- tb
    }
    as.data.frame(mat)
})

## 3<->1 letter code
AAC <- seqinr::aaa( seqinr::a())
names(AAC) <-  seqinr::a()
AAC <- AAC[AAC!="Stp"]

## AA colors; most divergent from
## https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
n <- length(CODL)
set.seed(1)
aa.cols <- sample(color, n)
names(aa.cols) <- names(CODL)

## clustal, via https://www.bioinformatics.nl/~berndb/aacolour.html
clustal.cols <- aa.cols
clustal.cols[] <- "gray"
clustal.cols[c("G","P","S","T")] <- "orange"
clustal.cols[c("H","K","R")] <- "red"
clustal.cols[c("F","W","Y")] <- "blue"
clustal.cols[c("I","L","M","V")] <- "green"


## lesk, via https://www.bioinformatics.nl/~berndb/aacolour.html
lesk.cols <- aa.cols
lesk.cols[c("G","A","S","T")] <- "orange" # small nonpolar
lesk.cols[c("C","V","I","L","P","F","Y","M","W")] <- "green" ## hydrophobic
lesk.cols[c("N","Q","H")] <- "magenta" # polar
lesk.cols[c("D","E")] <- "red" # negatively charged
lesk.cols[c("K","R")] <- "blue" # positively charged


## shapely/rasmol via 
shapely.cols <- aa.cols
shapely.cols[c("D","E")]   <- "#E60A0A"  # acidic
shapely.cols[c("K","R")]   <- "#145AFF"  # basic
shapely.cols[c("S","T")]   <- "#FA9600"  # S/T - P-targets
shapely.cols[c("C","M")]   <- "#E6E600"  # sulfur-containing
shapely.cols[c("P")]     <- "#DC9682"    # stiff backbone 
shapely.cols[c("F","Y")]   <- "#3232AA"  # aromatic
shapely.cols[c("W")]     <- "#B45AB4"    # large aromatic
shapely.cols[c("H")]     <- "#8282D2"    # aromatic, nitrogen
shapely.cols[c("N","Q")]   <- "#00DCDC"  # nitrogen-containing
shapely.cols[c("G")]     <- "#EBEBEB"    # just backbone
shapely.cols[c("L","V","I")] <- "#0F820F"# branched chain
shapely.cols[c("A")]     <- "#C8C8C8"    # minimal residue

## darker grays for white background
shapely.cols[c("G")]     <- "#A0A0A0"   
shapely.cols[c("A")]     <- "#909090"    # minimal residue

aa.pchs <- rep(5, length(aa.cols))
names(aa.pchs) <- names(aa.cols)
aa.pchs[c("G","P","S","T")] <- 1:4
aa.pchs[c("H","K","R")] <- 1:3
aa.pchs[c("F","W","Y")] <- 1:3
aa.pchs[c("I","L","M","V")] <- 1:4

## AA pch, based on coloring scheme
shapely.pchs <- aa.pchs
shapely.pchs[c("D","E")]   <- c(19,4)
shapely.pchs[c("C","M")]   <- c(4,19)
shapely.pchs[c("K","R")]   <- c(19,4)
shapely.pchs[c("S","T")]   <- c(19,4)
shapely.pchs[c("F","Y")]   <- c(19,4)
shapely.pchs[c("N","Q")]   <- c(19,4)
shapely.pchs[c("G")]     <- 19 #light grey
shapely.pchs[c("L","V","I")] <- c(19,2,4)
shapely.pchs[c("A")]     <- 19
shapely.pchs[c("W")]     <- 4
shapely.pchs[c("H")]     <- 19
shapely.pchs[c("P")]     <- 4


## SELECT AA COLOR SCHEME
aa.cols <- shapely.cols
aa.pchs <- shapely.pchs

## AA SORTING
aa.srt <- names(sort(aa.cols))

aa.cols <- aa.cols[aa.srt]
aa.pchs <- aa.pchs[aa.srt]

## AA sorting by property
aaprop.srt <- c("charged","polar","special","hydrophobic")
aap.srt <- names(AAPROP[order(match(AAPROP, aaprop.srt))])

## MANUAL SORTING: based on above but P as first special
aap.srt <- c("R","H","K","D","E","S","T","N","Q",
             "P","U","C","G","A","V","I","L","M","F","Y","W")

## TODO: from:to sorting by property
aap.ftsrt <- cbind(rep(aaprop.srt, length(aaprop.srt)),
                   rep(aaprop.srt, each=length(aaprop.srt)))


aap.cols <- aaprop.cols[aaprop]
names(aap.cols) <- names(aaprop)
aap.cols <- aap.cols[aap.srt]




### START

## PARSE & FILTER DATA
dat <- read.delim(in.file)
tmtf <- read.delim(tmt.file)

## convert to logical

## NOTE, 2024-08-25: for publication this script moved from using
## unfiltered BP/SAAP to the official supplement, filtered for IG
## genes; this file does not have the Keep.SAAP flag anymore, since
## since it only contains rows where Keep.SAAP was TRUE.
if ( !"Keep.SAAP"%in%colnames(tmtf) ) {
    tmtf$Keep.SAAP <- TRUE
} else {
    tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)
}

## tag albumin/hemoglobin
dat$Hemoglobin.Albumin <- dat$albumin|dat$globin 

## tag K/R substitutions
dat$KR <- dat$from%in%c("K","R") | dat$to%in%c("K","R")

dat$Keep.SAAP <- !dat$IG & !dat$KR

## 20240906: POSITIONAL PROBABILITY FILTER
pprob <- readxl::read_xlsx(prob.file)


## FILTER INDEPENDENTLY OF Keep.SAAP, but (i) filter from tmtf and
## (ii) done below(?), remove all from bdat that are not present in
## tmtf


## NOTE: with BP,SAAP,Datset,TMT/Tissue-wise pos.prob.
## we wouldn't need the split
ppl <- split(pprob$"Positional probability",
             paste(pprob$BP, pprob$SAAP,
                   pprob$Dataset, pprob$"TMT/Tissue"))
pp <- unlist(lapply(ppl, max, na.rm=TRUE))
tmtf$pp <-  pp[match(paste(tmtf$BP, tmtf$SAAP,
                           tmtf$Dataset, tmtf$TMT.Tissue),
                     names(pp))]
tmtf$pp[is.na(tmtf$pp)] <- -Inf
fal <- split(pprob$"Fragment at AAS site",
             paste(pprob$BP, pprob$SAAP,
                   pprob$Dataset, pprob$"TMT/Tissue"))
fa <- unlist(fal) #lapply(fal, any))
tmtf$fa <-  fa[match(paste(tmtf$BP, tmtf$SAAP,
                           tmtf$Dataset, tmtf$TMT.Tissue),
                     names(fa))]
tmtf$fa[is.na(tmtf$fa)] <- FALSE


plotdev(file.path(ifig.path, paste0("position_probabilities")),
        width=3, height=2, res=200, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(tmtf$pp, xlab="position probabilities", main=NA)
abline(v=MINPP, col=2)
abline(v=MAXPP, col=4)
dev.off()

## TODO: merge fragment at AAS site filter
plotdev(file.path(ifig.path, paste0("position_probabilities_bar")),
        width=3, height=3, res=200, type=ftyp)
par(mai=c(.75,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(t(table(cut(tmtf$pp, breaks=seq(0,1,.1), include.lowest = TRUE),
                tmtf$fa)), beside=TRUE, legend=TRUE, las=2,
        args.legend=list(x="top", title="AAS fragment?"))
figlabel("pos. prob.", pos="bottomleft")
dev.off()


if ( RM.POSPROB ) {
    ##dat$Keep.SAAP <- dat$Keep.SAAP & (dat$pp >= MINPP) # & dat$faas)
    tmtf <- tmtf[tmtf$pp >= MINPP & tmtf$fa,] 
} else if ( MN.POSPROB ) {
    ##dat$Keep.SAAP <- dat$Keep.SAAP & (dat$pp <= MAXPP) # & !dat$faas)
    tmtf <- tmtf[tmtf$pp <= MAXPP & !tmtf$fa,]
}

### UNIFY FILTER COLUMNS

## fuse all excluded tags for each unique SAAP
alls <- rbind(dat[,c("Keep.SAAP","SAAP")],
              tmtf[,c("Keep.SAAP","SAAP")])
alls <- split(alls$Keep.SAAP, alls$SAAP)


## remove if tagged so for any dataset
keep <- unlist(lapply(alls, all)) # all keep tags must be tree
dat$keep <- keep[dat$SAAP]
tmtf$keep <- keep[tmtf$SAAP]

## remove excluded
cat(paste("removing", sum(!dat$keep, na.rm=TRUE),
          "tagged as false positive on protein level\n"))
hdat <- dat[which(dat$keep),]

cat(paste("removing", sum(hdat$match!="good"),
          "tagged as bad blast hit\n"))
hdat <- hdat[which(hdat$match=="good"),]


## get raw RAAS data TMT level
## remove excluded
cat(paste("removing", sum(!tmtf$keep, na.rm=TRUE),
          "tagged as false positive on TMT level\n"))

tmtf <- tmtf[tmtf$keep,]
## exclude NA or Inf
rm <- is.na(tmtf$RAAS) | is.infinite(tmtf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from TMT level\n"))
tmtf <- tmtf[!rm,]

## CLEAN PROTEIN LEVEL
tmtsaap <- unique(tmtf$SAAP)
rm <- !hdat$SAAP%in%tmtsaap

cat(paste("removing", sum(rm), "SAAP without values at TMT level\n"))
hdat <- hdat[!rm,]


## first find mutated AA pairs
## (column AASin input mixes I/L)
## and split mutated AA pairs into from/to
## TODO: get this from the columns in the mapped file, once clear
## why missing
saaps <- strsplit(hdat$SAAP,"")
bases <- strsplit(hdat$BP, "")
fromtol <- lapply(1:length(saaps), function(i) {
    pos <- which(saaps[[i]]!=bases[[i]])
    c(from=bases[[i]][pos], to=saaps[[i]][pos])
})
fromto <- do.call(rbind, fromtol)

## replace from/to columns
hdat$from <- fromto[,1]
hdat$to <- fromto[,2]

## just useful for plots instead of raw codons
hdat$aacodon <- paste(hdat$from, hdat$codon, sep="-")

## fromto column
hdat$fromto <- apply(fromto, 1, paste0, collapse=":")

## generate columns by AA property
hdat$pfrom <- aaprop[fromto[,1]]
hdat$pto <- aaprop[fromto[,2]]
hdat$pfromto <- paste0(hdat$pfrom,":",hdat$pto)
hdat$frompto <- paste0(hdat$from, ":", hdat$pto)


### STRUCTURE: ADD BINS FOR NUMERIC VALUES

## secondary structure by names
ssrt <- c(E="sheet", H="helix", C="coil")
hdat$sstruc <- ssrt[hdat$s4pred]
hdat$sstruc[hdat$sstruc==""] <- "na"

## iupred3
iupred3.bins <- cut(hdat$iupred3, breaks=seq(0,1,.2), include.lowest = TRUE)
rsrt <- levels(iupred3.bins)
hdat$iupred3.bins <- as.character(iupred3.bins)
hdat$iupred3.bins[is.na(hdat$iupred3.bins)] <- "na"
## anchor2
anchor2.bins <- cut(hdat$anchor2, breaks=seq(0,1,.2), include.lowest = TRUE)
rsrt <- levels(anchor2.bins)
hdat$anchor2.bins <- as.character(anchor2.bins)
hdat$anchor2.bins[is.na(hdat$anchor2.bins)] <- "na"
## flDPnn - disordered
flDPnn.bins <- cut(hdat$flDPnn, breaks=seq(0,1,.2), include.lowest = TRUE)
rsrt <- levels(flDPnn.bins)
hdat$flDPnn.bins <- as.character(flDPnn.bins)
hdat$flDPnn.bins[is.na(hdat$flDPnn.bins)] <- "na"
## DisoRDPbind - 
DisoRDPbind.bins <- cut(hdat$DisoRDPbind, breaks=seq(0,1,.2), include.lowest = TRUE)
rsrt <- levels(DisoRDPbind.bins)
hdat$DisoRDPbind.bins <- as.character(DisoRDPbind.bins)
hdat$DisoRDPbind.bins[is.na(hdat$DisoRDPbind.bins)] <- "na"
## MMSeq2 - conservation
mmseq2 <- hdat$MMSeq2
mmseq2[mmseq2 <= 0] <- 1e-200
mmseq2[mmseq2 > 5] <- 5
MMSeq2.bins <- cut(mmseq2, breaks=seq(0,5,1), include.lowest = TRUE)
rsrt <- levels(MMSeq2.bins)
hdat$MMSeq2.bins <- as.character(MMSeq2.bins)
hdat$MMSeq2.bins[is.na(hdat$MMSeq2.bins)] <- "na"
## SCRIBER - binding
SCRIBER.bins <- cut(hdat$SCRIBER, breaks=seq(0,1,.2), include.lowest = TRUE)
rsrt <- levels(SCRIBER.bins)
hdat$SCRIBER.bins <- as.character(SCRIBER.bins)
hdat$SCRIBER.bins[is.na(hdat$SCRIBER.bins)] <- "na"
## ASAquick - surface
ASAquick.bins <- cut(hdat$ASAquick, breaks=seq(0,1,.2), include.lowest = TRUE)
rsrt <- levels(ASAquick.bins)
hdat$ASAquick.bins <- as.character(ASAquick.bins)
hdat$ASAquick.bins[is.na(hdat$ASAquick.bins)] <- "na"

## protein length bins
hdat$loglen <- log10(hdat$len)
loglen.bins <- cut(log10(hdat$len), breaks=seq(1.5,4,.5), include.lowest = TRUE)
hdat$loglen.bins <- as.character(loglen.bins)

## by ranks
rank.vals <- c("iupred3", "anchor2", "ASAquick", "SCRIBER", "MMSeq2",
               "DisoRDPbind", "flDPnn", "loglen")
val <- rank.vals[1]
rank.mat <- matrix(NA, ncol=length(rank.vals), nrow=nrow(hdat))
colnames(rank.mat) <- rank.vals
for ( val in rank.vals ) {
    ranks <- rank(hdat[,val], na.last="keep")
    rankbins <- as.character(cut(ranks/max(ranks,na.rm=TRUE), seq(0,1,.2),
                                 include.lowest = TRUE))
    rankbins[is.na(rankbins)] <- "na"
    rank.mat[,val] <- rankbins
    plotdev(file.path(ifig.path, paste0("rankbins_",val)),
            width=3, height=2, res=200, type=ftyp)
    par(mai=c(.75,.5,.1,.1), mgp=c(1.5,.3,0), tcl=-.25)
    boxplot(hdat[,val] ~ rankbins, ylab=val, las=2, xlab="")
    dev.off()
}
ranksrt <- c("na",levels(cut(1:10/10, seq(0,1,.2), include.lowest = TRUE)))
colnames(rank.mat) <- paste0(colnames(rank.mat),".rank")
hdat <- cbind(hdat, rank.mat)


### MAP PROTEIN LEVEL INFO TO TMT Data

idx.old <- match(tmtf$SAAP, hdat$SAAP)
idx <- match(paste(tmtf$BP, tmtf$SAAP), paste(hdat$BP, hdat$SAAP))


ina <- which(is.na(idx))

## REMOVE MISSING BP/SAAP
if ( length(ina)>0 ) {

    cat(paste("TODO:", length(ina), "BP/SAAP missing from BP/SAAP file.\n"))
    tmtf <- tmtf[-ina,]
    idx <- idx[-ina]    
}


## AA/codon/structure mapping
tmtf$pos <- hdat$pos[idx]
tmtf$rpos <- hdat$rpos[idx]
tmtf$from <- hdat$from[idx]
tmtf$to <- hdat$to[idx]
tmtf$fromto <- hdat$fromto[idx]
tmtf$pfrom <- hdat$pfrom[idx]
tmtf$pto <- hdat$pto[idx]
tmtf$pfromto <- hdat$pfromto[idx]
tmtf$frompto <- hdat$frompto[idx]
tmtf$codon <- hdat$codon[idx]
tmtf$aacodon <- hdat$aacodon[idx]
tmtf$iupred3 <- hdat$iupred3[idx]
tmtf$anchor2 <- hdat$anchor2[idx]
tmtf$sstruc <- hdat$sstruc[idx]

## DescribePROT
tmtf$flDPnn <- hdat$flDPnn[idx]
tmtf$DisoRDPbind <- hdat$DisoRDPbind[idx]
tmtf$MMSeq2 <- hdat$MMSeq2[idx]
tmtf$ASAquick <- hdat$ASAquick[idx]
tmtf$SCRIBER <- hdat$SCRIBER[idx]

## PFAM/CLAN
tmtf$pfam <- hdat$pfam[idx]
tmtf$clan <- hdat$clan[idx]
tmtf$pfam.ebi <- hdat$pfam.ebi[idx]
tmtf$clan.ebi <- hdat$clan.ebi[idx]

## bins of structural data
tmtf$iupred3.bins <- hdat$iupred3.bins[idx]
tmtf$anchor2.bins <- hdat$anchor2.bins[idx]
tmtf$flDPnn.bins <- hdat$flDPnn.bins[idx]
tmtf$DisoRDPbind.bins <- hdat$DisoRDPbind.bins[idx]
tmtf$MMSeq2.bins <- hdat$MMSeq2.bins[idx]
tmtf$ASAquick.bins <- hdat$ASAquick.bins[idx]
tmtf$SCRIBER.bins <- hdat$SCRIBER.bins[idx]
tmtf$loglen.bins <- hdat$loglen.bins[idx]

## ranks of structural data
tmtf$iupred3.rank <- hdat$iupred3.rank[idx]
tmtf$anchor2.rank <- hdat$anchor2.rank[idx]
tmtf$flDPnn.rank <- hdat$flDPnn.rank[idx]
tmtf$DisoRDPbind.rank <- hdat$DisoRDPbind.rank[idx]
tmtf$MMSeq2.rank <- hdat$MMSeq2.rank[idx]
tmtf$ASAquick.rank <- hdat$ASAquick.rank[idx]
tmtf$SCRIBER.rank <- hdat$SCRIBER.rank[idx]
tmtf$loglen.rank <- hdat$loglen.rank[idx]


## gene mapping
tmtf$name <- hdat$name[idx]
tmtf$gene <- hdat$gene[idx]
tmtf$transcript <- hdat$transcript[idx]
tmtf$protein <- hdat$protein[idx]
tmtf$ensembl <- hdat$ensembl[idx]
tmtf$mane <- hdat$MANE.protein[idx]
namane <- is.na(tmtf$mane)|tmtf$mane==""
tmtf$mane[namane] <- tmtf$ensembl[namane]

## tag protein sites
## NOTE: these are not unique chromosome coors, because
## for some proteins a mutated BP mapped to a different
## protein
tmtf$unique.site <- paste0(tmtf$ensembl, "_", tmtf$pos)


## gene names
tnm <- tmtf$name
tnm[tnm==""] <- tmtf$ensembl[tnm==""]
tmtf$name <- tnm

## tagging protein type
tmtf$albumin <- hdat$Hemoglobin.Albumin[idx]
tmtf$extracellular <- hdat$extracellular[idx]

## AAS position in peptide
tmtf$site <- hdat$site[idx]
tmtf$rsite <- hdat$site[idx]/nchar(hdat$BP)[idx]
tmtf$rsite.bins <- cut(tmtf$rsite, seq(0,1,.2), include.lowest = TRUE)

## coordinates in transcripts and chromosomes
tmtf$tpos <- hdat$tpos[idx]
tmtf$chr <- hdat$chr[idx]
tmtf$coor <- hdat$coor[idx]
tmtf$strand <- hdat$strand[idx]
tmtf$ssd <- hdat$ssd[idx] # distance to next splice site
tmtf$exl <- hdat$exl[idx] # exon length
tmtf$exn <- hdat$exn[idx] # exon number

## tag chromosome coordinates
tmtf$unique.coor <- paste(tmtf$chr, tmtf$coor, tmtf$strand, sep=":")


### MEAN AND MEDIAN RAAS

usaap <- paste0(tmtf$SAAP,"/",tmtf$BP,"/",tmtf$Dataset)
araas <- split(tmtf$RAAS, usaap)

raas.emedian <- unlist(lapply(araas, median))
raas.median <- unlist(lapply(araas, function(x) { log10(median(10^x)) } ))
raas.emean <- unlist(lapply(araas, mean))
raas.mean <- unlist(lapply(araas, function(x) { log10(mean(10^x)) } ))

tmtf$unique <- usaap
tmtf$RAAS.median <- raas.median[usaap]


## RAAS BINS
raas.bins <- cut(tmtf$RAAS.median, breaks=c(-6,-4,-2,-1,0,3),
                 include.lowest = TRUE)
raas.srt <- levels(raas.bins)
tmtf$raas.bins <- as.character(raas.bins)
tmtf$raas.bins[is.na(tmtf$raas.bins)] <- "na"

## FROM THIS POINT ON WE ONLY WORK WITH TMTF, WHICH
## NOW SHOULD HAVE ALL INFORMATION

## dataset sorting and labels
uds <- sort(unique(tmtf$Dataset))
uds <- c("Healthy",uds[uds!="Healthy"])
uds <- c("Cancer",uds[uds!="Cancer"])
uds <- uds[uds%in%tmtf$Dataset]

## additionally all and cancer-only
auds <- uds
auds <- c("cancer", auds)
auds <- c("all", auds)

## TODO: below should go above where amino acid sorting and colors
## is defined

## sorting of amino acid propertis
srt <- c("charged","polar","hphobic","special")
srt <- paste(rep(srt,each=4), rep(srt,4), sep=":")
## amino acid property axis labels
axex <- ftlabels(srt) # axis labels with arrows


## RAAS COLORS

## colors used for codon dotplot
plotdev(file.path(ifig.path,paste0("legend_raas_vcols")),
        res=300, type=ftyp, width=4, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(tmtf$RAAS,
                          mn=RAAS.MIN, mx=RAAS.MAX,colf=COLF,
                          n=50, plot=TRUE,
                          mai=c(.5,.5,.1,.1),
                          xlab=expression(TMT~level~log[10]*RAAS))
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(colors, pos="bottomleft", cex=1)
figlabel(LAB, pos="bottomright", cex=1)
dev.off()

## globally used RAAS colors!!
vcols <- lraas.col$col
vbrks <- lraas.col$breaks

## colors used for all other dotplots
plotdev(file.path(ifig.path,paste0("legend_raas_acols")),
        res=300, type=ftyp, width=4, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
aaprop.raas.col <- selectColors(tmtf$RAAS,
                          mn=RAAS.MINA, mx=RAAS.MAXA,colf=COLF,
                          n=50, plot=TRUE,
                          mai=c(.5,.5,.1,.1),
                          xlab=expression(TMT~level~log[10]*RAAS))
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(colors, pos="bottomleft", cex=1)
figlabel(LAB, pos="bottomright", cex=1)
dev.off()

## tight RAAS range
acols <- aaprop.raas.col$col
abrks <- aaprop.raas.col$breaks




## LEGENDS FOR DOT PLOTs

## broadest RAAS COLORS
pp <- seq(0, -log10(p.dot), length.out=3)
rs <- c(-4,-2,0,1) #seq(RAAS.MIN,RAAS.MAX, length.out=3)
pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))
colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- round(rs,1)
ovlg <- list(p.value=t(10^-pm), median=t(rm))


## slim
plotdev(file.path(ifig.path,paste0("legend_dotplot_vcols_slim")),
        height=5, width=1, res=300, type=ftyp)
layout(1:2, heights=c(.5,.4))
par(mai=c(0,0.05,.15,.6), mgp=c(1.3,.3,0), tcl=-.25)
image_matrix(y=vbrks, t(t(vbrks)), col=rev(vcols), xlab=NA, ylab=NA)
mtext(xl.raas, 4, 2, cex=1.5)
axis(4, at=-5:5, las=2, cex.axis=1.5)
ovp <- list(p.value=ovlg$p.value[,1,drop=FALSE],
            median=ovlg$median[,1,drop=FALSE])
par(mai=c(1,0.25,.15,.5))
dotprofile(ovp, value="median",
           vbrks=vbrks,
           vcols=vcols, 
           dot.sze=1.5*dot.sze, p.dot=p.dot, axis=NA,
           ylab=plab,
           xlab=NA, xpd=TRUE)
axis(4, at=nrow(ovp$p.value):1, labels=log10(ovp$p.value), las=2, col=NA,
     cex.axis=1.5)
mtext(plab, 1, 0.75, cex=1.5, adj=.4)
dev.off()


## tight RAAS range - legend for acols/abreaks
pp <- seq(0, -log10(p.dot), length.out=3)
rs <- c(-4,-2,-1,0) #seq(RAAS.MIN,RAAS.MAX, length.out=3)
pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))
colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- round(rs,1)
ovlg <- list(p.value=t(10^-pm), median=t(rm))


plotdev(file.path(ifig.path,paste0("legend_dotplot_acols_slim")),
        height=5, width=1, res=300, type=ftyp)
layout(1:2, heights=c(.5,.3))
par(mai=c(0,0.15,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
image_matrix(y=abrks, t(t(abrks)), col=rev(acols), xlab=NA, ylab=NA)
mtext(xl.raas, 4, 1.6, cex=1.2)
axis(4, at=-5:5, las=2, cex.axis=1.2)
ovp <- list(p.value=ovlg$p.value[,1,drop=FALSE],
            median=ovlg$median[,1,drop=FALSE])
par(mai=c(1,0.25,.15,.5))
dotprofile(ovp, value="median",
           vbrks=abrks,
           vcols=acols, 
           dot.sze=1.5*dot.sze, p.dot=p.dot, axis=NA,
           ylab=plab,
           xlab=NA, xpd=TRUE)
axis(4, at=nrow(ovp$p.value):1, labels=log10(ovp$p.value), las=2, col=NA,
     cex.axis=1.2)
mtext(plab, 1, 0.75, cex=1.2)
dev.off()

## legend for p.values 
plotdev(file.path(ifig.path,paste0("legend_pvals_horizontal")),
        type=ftyp, res=300, width=2, height=.6, bg=NA)
par(mai=c(.3,.05,.05,.05), mgp=c(.5,.3,0), tcl=-.25, xaxs="i")
plotOverlapsLegend(p.min=p.min, p.txt=p.txt, dir=2)
dev.off()


### GENERATE unique BP/SAAP TABLE: FILTER and ANNOTATE

bpraas <- split(tmtf$RAAS, paste(tmtf$BP, tmtf$SAAP)) #tmtf$BP)
bpraas <- listProfile(bpraas, y=tmtf$RAAS, use.test=use.test, min=3)

bdat <- hdat##[!duplicated(hdat$BP),]
rownames(bdat) <- paste(bdat$BP, bdat$SAAP)

## MISSING?
## manual inspection shows this BP had infinite RAAS
missing <- which(!rownames(bdat)%in%rownames(bpraas))
if ( length(missing) )
    cat(paste("WARNING:", length(missing), "BP missing from TMT level file\n"))

## keep only bdat for which we have RAAS values
bdat <- cbind(bdat[rownames(bpraas),], bpraas)

## use median as RAAS
bdat$RAAS <- bdat$median

## add Datasets+Tissues where each BP/SAAP appears
datasets <- split(tmtf$Dataset, paste(tmtf$BP, tmtf$SAAP))
datasets <- lapply(datasets, unique)
bdat$numDatasets <-
    lengths(datasets)[paste(bdat$BP, bdat$SAAP)]
datasets <- unlist(lapply(datasets,  paste, collapse=";"))
bdat$Datasets <- datasets[paste(bdat$BP, bdat$SAAP)]

tissues <- tmtf$TMT.Tissue
tissues[tmtf$Dataset!="Healthy"] <- "cancer"

## ADD TISSUE INFO BACK TO TMT TABLE
tmtf$tissue <- tissues

tissues <- split(tissues, paste(tmtf$BP, tmtf$SAAP))
tissues <- lapply(tissues, unique)
bdat$numTissues <-
    lengths(tissues)[paste(bdat$BP, bdat$SAAP)]
tissues <- unlist(lapply(tissues,  paste, collapse=";"))
bdat$Tissues <- tissues[paste(bdat$BP, bdat$SAAP)]




### PROTEIN PROPERTIES
## collect published protein properties, each as
## hash with ensembl protein names

### PROTEIN ID MAPPING
## generate vectors for unique mapping of ensembl protein IDs to
## gene names, uniprot and refseq IDs


## ENSEMBL <-> NAME MAPPING
## TODO: add MANE column
genes <- read.delim(feature.file)
genes <- genes[genes$proteins!="" & !is.na(genes$proteins),]
##genes <- genes[genes$name!="" & !is.na(genes$name),]
noname <- genes$name=="" | is.na(genes$name)
genes$name[noname] <- genes$ID[noname]
ptl <- strsplit(genes$proteins, ";")
ens2nam <- rep(genes$name, lengths(ptl))
names(ens2nam) <- unlist(ptl)
pnms <- ens2nam

## list of all MANE proteins
MANES <- genes$MANE.protein
MANES <- MANES[MANES!=""]

## gene name synonyms - to catch missing proteins in experimental data below
syns <- read.delim(synonym.file)
syns <- syns[syns[,2]!="",]

## UNIPROT <-> ENSEMBL MAPPING
## NOTE: ~90 duplicated ensembl IDs
uni2ens <- read.delim(uni2ens.file, header=FALSE)
uni2ens[,2] <- sub("\\..*", "", uni2ens[,2]) # remove ensembl version tag
uni2ens[,1] <- sub("-.*", "", uni2ens[,1]) # remove uniprot version tag

## filter for ensembl IDs in our data!
uni2dat <- uni2ens[uni2ens[,2]%in%c(dat$ensembl, dat$mane),]
## 1-many lists
##uni2e <- split(uni2ens[,2], uni2ens[,1])
ens2u <- split(uni2ens[,1], uni2ens[,2])

## REFSEQ <-> ENSEMBL MAPPING
refseq2ens <- read.delim(refseq.file, header=TRUE)


## PROTEIN LENGTHS via saap_mapped.tsv in hdat
plen <- split(hdat$len, hdat$ensembl)
plen <- lapply(plen, unique)
plen <- unlist(lapply(plen, function(x) x[1]))

## PROTEIN HALF-LIVES via Mathieson et al. 2018
hlvd <- readxl::read_xlsx(math18.file)
## mean half-live over all replicates and cell types
## TODO: consider halflife distributions instead of just taking mean
## NOTE: THIS includes mouse data, but they correlate.
cidx <- grep("half_life", colnames(hlvd), value=TRUE)
phlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
names(phlv) <- unlist(hlvd[,1])


## rename via synonyms
miss <- which(!names(phlv)%in%pnms)
mnms <- names(phlv)[miss]
## GET FIRST OF MULTIPLE
msyns <- sapply(mnms, function(x) syns[grep(x, syns[,2], fixed=TRUE)[1],1])
msyns <- msyns[!is.na(msyns)]
## rename
found <- miss[names(phlv)[miss]%in%names(msyns)]
names(phlv)[found] <- msyns[names(phlv)[found]]





### CALCULATE PROTEIN ABUNDANCES via Razor.protein.precursor.intensity

## remove brackets
pint <- sub("\\[","", sub("\\]","", tmtf$Razor.protein.precursor.intensity))
## split multiple values and convert to numeric
pint <- lapply(pint, function(x) as.numeric(unlist(strsplit(x, ","))))
## remove zeros and take unique value
pint <- lapply(pint, function(x) unique(x[x!=0]))
## NOTE: there are still many entries with multiple values
pint <- unlist(lapply(pint, median))

### ADD TO TMTF TABLE
tmtf$razor.intensity <- pint

## MEDIAN of PROTEINS
pint <- split(pint, tmtf$ensembl)
pint <- lapply(pint, function(x) x[!is.na(x)])

pint <- lapply(pint, as.numeric)

## median protein intensities
pint <- lapply(pint, median, na.rm=TRUE)
pint <- unlist(pint)

### ADD TO TMTF TABLE
tmtf$protein.intensity <- pint[tmtf$ensembl]


#### GENERATE COLLAPSED TABLES AT DIFFERENT LEVELS
## 1) unique chromosome coordinates,
## 2) unique protein sites,
## 3) unique protein sites X AAS types.


## add these columns to unique site tables
add.data <- c("name",
              "codon", "fromto","Dataset","tissue",
              "iupred3", "flDPnn", # disorder
              "anchor2","DisoRDPbind", # disorder-binding
              "MMSeq2", # sequence conservation
              "BP","SAAP",
              "ensembl",
              "pos","rpos", # position in protein
              "transcript",
              "tpos", "ssd", "exl", "exn", "chr", "coor", "strand") 



### GENERATE UNIQUE CHROMOSOME SITE TABLE
## median raas per unique mane protein site
## plus protein level data since this table is used
## externally for random forest et al. modeling

sitl <- split(tmtf$RAAS, paste(tmtf$unique.coor))
csite <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## mod. column names and add protein and site info
colnames(csite) <- paste0("RAAS.", colnames(csite))


## enumeration of sites within unique proteins
nsites <- as.numeric(sub(".*\\.","",tagDuplicates(csite$ensembl)))
nsites[is.na(nsites)] <- 1
csite$n <- nsites


## COLLECT DATA FOR SITES
for ( i in 1:length(add.data) ) {
    
    datl <- lapply(split(tmtf[[add.data[i]]], tmtf$unique.coor), unique)
    datl <- lapply(datl, function(x) x[x!="" & !is.na(x)])
    datl <- lapply(datl, function(x)
        ifelse(length(x)==1, x, paste(x, collapse=";")))
    datl <- unlist(datl)
    if ( any(lengths(datl)>1) )
        stop("all protein site data should be unique")
    csite[[add.data[i]]] <- datl[rownames(csite)]
}

## add gene names
csite$name <- ens2nam [csite$ensembl]

## add uniprot id
csite$uniprot <- unlist(lapply(ens2u[csite$ensembl], paste, collapse=";"))

## RAAS COLOR
csite$RAAS.color <- num2col(csite$RAAS.median,
                            limits=c(RAAS.MIN, RAAS.MAX), colf=arno)


## NOTE: where we require BP and SAAP we take the longest (but store all)
## get longest BP
bpl <- strsplit(csite$BP, ";")
bpl <- unlist(lapply(bpl, function(x) x[which.max(nchar(x))]))
csite$all.BP <- csite$BP
csite$BP <- bpl

## get longest SAAP
saapl <- strsplit(csite$SAAP, ";")
saapl <- unlist(lapply(saapl, function(x) x[which.max(nchar(x))]))
csite$all.SAAP <- csite$SAAP
csite$SAAP <- saapl

## get AA context for longest BP/SAAP pair
aal <- lapply(split(bdat$AA, paste(bdat$BP,bdat$SAAP)),unique)
if ( length(table(lengths(aal)))!=1 )
    stop("non-unique AA context per BP/SAAP pair")
aal <- lapply(split(bdat$AA, paste(bdat$BP,bdat$SAAP)),unique)
if ( any(!paste(csite$BP, csite$SAAP)%in%names(aal)) )
    warning("AA context not available for selected site BP/SAAP pair")
aal <- unlist(aal)
csite$AA <- aal[paste(csite$BP, csite$SAAP)]

## add protein data, half-live, melting point, intensity, length
csite$protein.length <- plen[csite$ensembl]
csite$protein.halflife <- phlv[pnms[csite$ensembl]]
csite$protein.intensity <- pint[csite$ensembl]


## MEDIAN RAAS as MAIN RAAS
csite$RAAS <- csite$RAAS.median


## write-out unique site-specific RAAS stats with additional data
## used for random forest modeling
csite.file <- file.path(out.path,"sites_raas_unique_genome.tsv")
##if ( !file.exists(csite.file) )
    write.table(file=csite.file, x=cbind(ID=rownames(csite),csite),
                sep="\t", row.names=FALSE, quote=FALSE, na="")



### GENERATE UNIQUE PROTEIN SITE TABLE
## median raas per unique mane protein site
## plus protein level data since this table is used
## externally for random forest et al. modeling

sitl <- split(tmtf$RAAS, paste(tmtf$unique.site))
asite <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## mod. column names and add protein and site info
colnames(asite) <- paste0("RAAS.", colnames(asite))
##asite$ensembl <- sub("_.*", "", rownames(asite))
##asite$pos <- as.numeric(sub(".*_", "", rownames(asite)))




## enumeration of sites within unique proteins
nsites <- as.numeric(sub(".*\\.","",tagDuplicates(asite$ensembl)))
nsites[is.na(nsites)] <- 1
asite$n <- nsites


## COLLECT DATA FOR SITES
for ( i in 1:length(add.data) ) {
    
    datl <- lapply(split(tmtf[[add.data[i]]], tmtf$unique.site), unique)
    datl <- lapply(datl, function(x) x[x!="" & !is.na(x)])
    datl <- lapply(datl, function(x)
        ifelse(length(x)==1, x, paste(x, collapse=";")))
    datl <- unlist(datl)
    if ( any(lengths(datl)>1) )
        stop("all protein site data should be unique")
    asite[[add.data[i]]] <- datl[rownames(asite)]
}


## add gene names
asite$name <- ens2nam [asite$ensembl]

## add uniprot id
asite$uniprot <- unlist(lapply(ens2u[asite$ensembl], paste, collapse=";"))

## RAAS COLOR
asite$RAAS.color <- num2col(asite$RAAS.median,
                           limits=c(RAAS.MIN, RAAS.MAX), colf=arno)

## NOTE: where we require BP and SAAP we take the longes (but store all)
## get longest BP
bpl <- strsplit(asite$BP, ";")
bpl <- unlist(lapply(bpl, function(x) x[which.max(nchar(x))]))
asite$all.BP <- asite$BP
asite$BP <- bpl

## get longest SAAP
saapl <- strsplit(asite$SAAP, ";")
saapl <- unlist(lapply(saapl, function(x) x[which.max(nchar(x))]))
asite$all.SAAP <- asite$SAAP
asite$SAAP <- saapl

## get AA context for longest BP/SAAP pair
aal <- lapply(split(bdat$AA, paste(bdat$BP,bdat$SAAP)),unique)
if ( length(table(lengths(aal)))!=1 )
    stop("non-unique AA context per BP/SAAP pair")
aal <- lapply(split(bdat$AA, paste(bdat$BP,bdat$SAAP)),unique)
if ( any(!paste(asite$BP, asite$SAAP)%in%names(aal)) )
    warning("AA context not available for selected site BP/SAAP pair")
aal <- unlist(aal)
asite$AA <- aal[paste(asite$BP, asite$SAAP)]

## add protein data, half-live, melting point, intensity, length
asite$protein.length <- plen[asite$ensembl]
asite$protein.halflife <- phlv[pnms[asite$ensembl]]
asite$protein.intensity <- pint[asite$ensembl]


## MEDIAN RAAS as MAIN RAAS
asite$RAAS <- asite$RAAS.median

## make sure its ordered
## TODO: ORDER PROTEINS BY SIZE OR NUMBER OF AAS
asite <- asite[order(asite$ensembl, asite$pos),]

## write-out unique site-specific RAAS stats with additional data
## used for random forest modeling
asite.file <- file.path(out.path,"sites_raas_unique.tsv")
##if ( !file.exists(asite.file) )
    write.table(file=asite.file, x=cbind(ID=rownames(asite),asite),
                sep="\t", row.names=FALSE, quote=FALSE, na="")


#### UNIQUE SITE x AAS  TABLE
## median raas per unique mane protein site, split by AAS type
## plus protein level data since this table is used
## externally for random forest et al. modeling

sitl <- split(tmtf$RAAS, paste(tmtf$unique.site, tmtf$fromto))
site <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## mod. column names and add protein and site info
colnames(site) <- paste0("RAAS.", colnames(site))


## enumeration of sites within unique proteins
nsites <- as.numeric(sub(".*\\.","",tagDuplicates(site$ensembl)))
nsites[is.na(nsites)] <- 1
site$n <- nsites


## COLLECT DATA FOR SITES
for ( i in 1:length(add.data) ) {
    
    datl <- lapply(split(tmtf[[add.data[i]]],
                         paste(tmtf$unique.site, tmtf$fromto)), unique)
    datl <- lapply(datl, function(x) x[x!="" & !is.na(x)])
    datl <- lapply(datl, function(x)
        ifelse(length(x)==1, x, paste(x, collapse=";")))
    datl <- unlist(datl)
    if ( any(lengths(datl)>1) )
        stop("all protein site data should be unique")
    site[[add.data[i]]] <- datl[rownames(site)]
}

## add gene names
site$name <- ens2nam[site$ensembl]

## add uniprot id
site$uniprot <- unlist(lapply(ens2u[site$ensembl], paste, collapse=";"))

## RAAS COLOR
site$RAAS.color <- num2col(site$RAAS.median,
                           limits=c(RAAS.MIN, RAAS.MAX), colf=arno)

## add protein data, half-live, melting point, intensity, length
site$protein.length <- plen[site$ensembl]
site$protein.halflife <- phlv[pnms[site$ensembl]]
site$protein.intensity <- pint[site$ensembl]

## make sure its ordered
## TODO: ORDER PROTEINS BY SIZE OR NUMBER OF AAS
site <- site[order(site$ensembl, site$pos),]




## write-out unique site-specific RAAS stats with additional data
## used for random forest and glm4 modeling
site.file <- file.path(out.path,"sites_raas.tsv")
##if ( !file.exists(site.file) )
    write.table(file=site.file, x=cbind(ID=rownames(site),site),
                sep="\t", row.names=FALSE, quote=FALSE, na="")

## TODO: also write out bdat or hdat, with all RAAS values added

exportCols <- c("BP","SAAP",
                "site","fromto", "codon",
                "RAAS","n", "Datasets", "Tissues",
                "name",
                "ensembl","pos",
                "transcript", "tpos",
                "gene","chr","coor","strand")

edat <- bdat[,exportCols]
rownames(edat) <- NULL
colnames(edat) <- sub("^pos$", "protein.position",
                      sub("^tpos$","transcript.position",
                          sub("^n$","RAAS.n",
                              sub("ensembl","protein",colnames(edat)))))

## NOTE: this table is used as supplemental data table 7
bpsaap.file <- file.path(out.path,"Supplemental_Data_7.SAAP_coordinates.tsv")
##if ( !file.exists(bpsaap.file) )
    write.table(file=bpsaap.file, x=edat, sep="\t",
                row.names=FALSE, quote=FALSE, na="")


### PLOT ALIGNMENT

fh <- fw <- .2 # width/height factor for invidiual fields in dotplots

CMAIL <- 1.54 ## commonly used between motif and domain figures, defined 
              ## in domain script via nchar to fit the y-axis labels
FONT <- "monospace" # font used for aligned figures in motifs and domains


##
cat(paste("KEEPING\n\t",
          nrow(bdat), "unique BP/SAAP\n\t",
          nrow(csite), "unique genomic sites\n"))
