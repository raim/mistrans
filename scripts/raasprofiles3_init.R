
### LOAD BP/SAAP MAPPING AND TMT LEVEL RAAS DATA
## for further downstream analyses, RAAS profiles by
## AA properties and codons (raasprofiles3_codons.R), and by
## protein complexes, proteins and protein windows
## (raasprofiles3_proteins.R).

library(Biostrings) # for genetic code, blosum62, etc
library(viridis)
library(segmenTools)
library(readxl)

options(stringsAsFactors=FALSE)
options(scipen=0) # use e notation for p-values

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## TODO:
## * which AA are missing from mapped file and why? likely
## because BP didnt match, generate from SAAP/BP,
## * global or local background distribution?

####!!!
## TODO 20240215:
## * rows: plot only largest contributor! Q->G, T->V,
##   TODO: find largest RAAS effect size.
## * RAAS distributions: show all and scale/color by pval
## * formal approach to missing data: map measured peptides to the
##   full peptide space (main peptides or all possible).
## * plot by protein class, e.g only Albumins.
## * CHECK HANDLING OF DUPLICATE SAAP here and in tmt retrieval


## FROM->TO BY AA PROPERTY CLASSES
## DONT FILTER FOR CODONS, use hdat

## TODO 20240222
## * plotovl: plot all function,
## * fall back on median,
## * re-blast saap, re-annotate,
## * Healthy: untangle tissues and collapse cancer.

### PARAMETERS

## TODO: use p.adjust for raasProfile and aaProfile calls
p.adjust <- "none" ## multiple hypothesis testing

## TODO: fuse global and tight color scale?
RAAS.MIN <- -4   # broad RAAS colors vcols: codon plot
RAAS.MAX <-  1
RAAS.MINA <- -4  # tighter RAAS colors acols: AA properties
RAAS.MAXA <- 0
RAAS.MINT <- -3  # tightest RAAS colors tcols: not used anymore
RAAS.MAXT <-  -1

colors <- "arno" # "inferno" #"rocket" # "viridis" # 

COLF <- get(colors)

p.min <- 1e-10
p.txt <- 1e-5

p.dot <- p.min # p.txt
dot.sze <- c(.3,2)

## STATISTICAL TEST TO RUN
use.test <- t.test # w.test # 

ftyp <- "png" # "pdf" # # 
##if ( !interactive() ) ftyp="pdf"

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


### TODO: move AA colors, codons etc. to saap_utils.R or a new
## saap_params.R

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
COD.SRT <- paste(GENETIC_CODE[unlist(CODL)], unlist(CODL), sep="-")

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

## color by codon frequency
if ( FALSE ) {
    aa.cols <-
        unlist(lapply(CODL, length))
    aa.pchs[] <- 1
}

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




#### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")

mam.path <- "~/data/mammary/"
codon.file <- file.path(mam.path,"processedData","coding_codons.tsv")
trna.file <- file.path(mam.path, "codons_GRCh38.tsv")

## external codon usage measures
## @Dana2014, Dana and Tuller 2014: decoding time
dana14.file <- file.path(dat.path,"dana14_codons.csv")
## @Wu2019: codon stability coefficient (CSC)
wu19.file <- file.path(dat.path,"elife-45396-fig1-data2-v2.csv")
gingold14.file <- file.path(dat.path, "gingold14_mmc2.xls")

in.file <- file.path(out.path,"saap_mapped4.tsv")
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")

## protein complexes
humap.file <- file.path(mam.path,"originalData","humap2_complexes_20200809.txt")
corum.file <- file.path(mam.path,"originalData","humanComplexes.txt")

## protein half-lives - @Mathieson2018
math18.file <- file.path(mam.path,"originalData",
                         "41467_2018_3106_MOESM5_ESM.xlsx")

## 20S targets - @Pepelnjak2024
pepe24.file <- file.path(mam.path,"originalData",
                         "44320_2024_15_moesm1_esm.xlsx")

##  @Watson2023 - T/osmo
## NOTE/TODO: unused since this is data from mouse cells;
## perhaps useful when mapped to human orthologs.
w23prot.file <- file.path(mam.path,"originalData",
                          "41586_2023_6626_MOESM4_ESM.xlsx")
w23pprot.file <- file.path(mam.path,"originalData",
                           "41586_2023_6626_MOESM5_ESM.xlsx")

## @Yang2022: thermal stability prediction
## https://structure-next.med.lu.se/ProTstab2/
protstab.file <- file.path(mam.path,"originalData",
                           "ProTstab2_human.csv")

## @Savitski2014: Tracking cancer drugs in living cells by thermal
## profiling of the proteome
thermo.file <- file.path(mam.path, "originalData",
                         "savitski14_tableS11.xlsx")
thatp.file <- file.path(mam.path, "originalData",
                         "savitski14_tableS3.xlsx")

## PROTEIN ID MAPPINGS

## uniprot/refseq/ensembl/name mappings
uni2ens.file <- file.path(mam.path,"originalData","uniprot_ensembl.dat")
uni2nam.file <- file.path(mam.path,"originalData","uniprot_name.dat")

refseq.file <- file.path(mam.path,"originalData",
                         "ensembl_refseq_20240528.tsv.gz")

## genome feature file
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")




## DATA FILTERS

## CANCERS OR HEALTHY TISSUES
healthy <- FALSE # TRUE #  

## TODO: extracellular is mostly album/globin - analyze
exclude.nterm <- FALSE # TRUE # 
exclude.extracellular <- FALSE # TRUE # 
exclude.albumin <- FALSE # TRUE # 
only.unique <- FALSE # TRUE # 
include.kr <- FALSE # TRUE # 

exclude.frequent <- FALSE # TRUE # 
frequent <- c("Q","W","T","S")

LAB <- "" # "all"
fig.path <- file.path(proj.path,"figures","raasprofiles3")
if (  exclude.albumin ) {
    fig.path <- paste0(fig.path,"_noalb")
    LAB <- "-Alb./Hemog."
}
if ( exclude.frequent ) {
    tmp <- paste(frequent,collapse=",")
    fig.path <- paste0(fig.path,"_", gsub(",","",tmp))
    LAB <- paste0("-", tmp)
}
if ( only.unique ) {
    fig.path <- paste0(fig.path,"_unique")
    LAB <- paste0(LAB, ", unique SAAP")
}
if ( exclude.extracellular ) {
    fig.path <- paste0(fig.path,"_no_extracellular")
    LAB <- paste0(LAB, " -extracell.")
}
if ( include.kr ) {
    fig.path <- paste0(fig.path,"_wKR")
    LAB <- paste0(LAB, "+K/R")
}
if (  exclude.nterm ) {
    fig.path <- paste0(fig.path,"_nterm")
    LAB <- paste0(LAB, "-N3")
}

SETID <- ifelse(healthy,"tissues","cancer")


## figure output paths
dir.create(fig.path, showWarnings=FALSE)
## folder for detailed distributions
dpath <- file.path(fig.path,"dists")
dir.create(dpath, showWarnings=FALSE)

### START

## PARSE & FILTER DATA
dat <- read.delim(in.file)
tmtf <- read.delim(tmt.file)

## convert to logical
tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)

if ( include.kr ) {
    tmtf$Keep.SAAP[grep("R|K", tmtf$AAS)] <- TRUE
}

## tag albumin/hemoglobin
## TODO: compare and fuse with protein level annotation
##       as.logical(dat$Hemoglobin.Albumin)
dat$Hemoglobin.Albumin <- dat$albumin|dat$globin 
dat$KR <- dat$from%in%c("K","R") | dat$to%in%c("K","R")

if ( include.kr ) {
    dat$KR <- FALSE
}

dat$Keep.SAAP <- !dat$IG & !dat$KR

#### TODO: find out which/how many TMT level SAAP are missing
## from the protein level file, and why.

### UNIFY FILTER COLUMNS


## fuse all excluded tags for each unique SAAP
alls <- rbind(dat[,c("Keep.SAAP","SAAP")],
              tmtf[,c("Keep.SAAP","SAAP")])
alls <- split(alls$Keep.SAAP, alls$SAAP)

## NOTE/TODO: better analyze conflicting tags and synchronize use with Shiri.
## tmtf has exclusion based on "Potential.uncaptured.transcript",
## Unfortunately, this excludes two nice AAS in KRAS.
allu <- lapply(alls, unique)
table(lengths(allu))


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

## STRUCTURE: add bins for numeric values
## secondary structure by names
ssrt <- c(E="sheet", H="helix", C="coil")
hdat$sstruc <- ssrt[hdat$s4pred]
hdat$sstruc[hdat$sstruc==""] <- "na"

## iupred3
iupred3.bins <- cut(hdat$iupred3, breaks=seq(0,1,.2))
rsrt <- levels(iupred3.bins)
hdat$iupred3.bins <- as.character(iupred3.bins)
hdat$iupred3.bins[is.na(hdat$iupred3.bins)] <- "na"
## anchor2
anchor2.bins <- cut(hdat$anchor2, breaks=seq(0,1,.2))
rsrt <- levels(anchor2.bins)
hdat$anchor2.bins <- as.character(anchor2.bins)
hdat$anchor2.bins[is.na(hdat$anchor2.bins)] <- "na"
## flDPnn - disordered
flDPnn.bins <- cut(hdat$flDPnn, breaks=seq(0,1,.2))
rsrt <- levels(flDPnn.bins)
hdat$flDPnn.bins <- as.character(flDPnn.bins)
hdat$flDPnn.bins[is.na(hdat$flDPnn.bins)] <- "na"
## DisoRDPbind - 
DisoRDPbind.bins <- cut(hdat$DisoRDPbind, breaks=seq(0,1,.2))
rsrt <- levels(DisoRDPbind.bins)
hdat$DisoRDPbind.bins <- as.character(DisoRDPbind.bins)
hdat$DisoRDPbind.bins[is.na(hdat$DisoRDPbind.bins)] <- "na"
## MMSeq2 - conservation
MMSeq2.bins <- cut(hdat$MMSeq2, breaks=seq(0,5,1))
rsrt <- levels(MMSeq2.bins)
hdat$MMSeq2.bins <- as.character(MMSeq2.bins)
hdat$MMSeq2.bins[is.na(hdat$MMSeq2.bins)] <- "na"
## SCRIBER - binding
SCRIBER.bins <- cut(hdat$SCRIBER, breaks=seq(0,1,.2))
rsrt <- levels(SCRIBER.bins)
hdat$SCRIBER.bins <- as.character(SCRIBER.bins)
hdat$SCRIBER.bins[is.na(hdat$SCRIBER.bins)] <- "na"
## ASAquick - surface
ASAquick.bins <- cut(hdat$ASAquick, breaks=seq(0,1,.2))
rsrt <- levels(ASAquick.bins)
hdat$ASAquick.bins <- as.character(ASAquick.bins)
hdat$ASAquick.bins[is.na(hdat$ASAquick.bins)] <- "na"

## by ranks
rank.vals <- c("iupred3", "anchor2", "ASAquick", "SCRIBER", "MMSeq2",
               "DisoRDPbind", "flDPnn")
val <- rank.vals[1]
rank.mat <- matrix(NA, ncol=length(rank.vals), nrow=nrow(hdat))
colnames(rank.mat) <- rank.vals
for ( val in rank.vals ) {
    ranks <- rank(hdat[,val])
    rankbins <- cut(ranks/max(ranks,na.rm=TRUE), seq(0,1,.2))
    rank.mat[,val] <- rankbins
    plotdev(file.path(fig.path, paste0("rankbins_",val)),
            width=3, height=2, res=200)
    par(mai=c(.75,.5,.1,.1), mgp=c(1.5,.3,0), tcl=-.25)
    boxplot(hdat[,val] ~ rankbins, ylab=val, las=2, xlab="")
    dev.off()
}
colnames(rank.mat) <- paste0(colnames(rank.mat),".rank")
hdat <- cbind(hdat, rank.mat)

### MAP PROTEIN LEVEL INFO TO TMT Data
idx.old <- match(tmtf$SAAP, hdat$SAAP)
idx <- match(paste(tmtf$BP, tmtf$SAAP), paste(hdat$BP, hdat$SAAP))

## 20240516: switched from matching only via SAAP to matching via BP/SAAP!
## compare differences between SAAP matching and BP/SAAP matching
if ( interactive() ) {
    sum(idx!=idx.old) ## 24 different matches
    df <- which(idx!=idx.old)[2]
    hdat[idx[df],]
    hdat[idx.old[df],]
}

ina <- which(is.na(idx))
if ( length(ina)>0 ) {
    cat(paste("TODO:", length(ina), "missing from unique saap file.\n"))
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

## ranks of structural data
tmtf$iupred3.rank <- hdat$iupred3.rank[idx]
tmtf$anchor2.rank <- hdat$anchor2.rank[idx]
tmtf$flDPnn.rank <- hdat$flDPnn.rank[idx]
tmtf$DisoRDPbind.rank <- hdat$DisoRDPbind.rank[idx]
tmtf$MMSeq2.rank <- hdat$MMSeq2.rank[idx]
tmtf$ASAquick.rank <- hdat$ASAquick.rank[idx]
tmtf$SCRIBER.rank <- hdat$SCRIBER.rank[idx]


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
tmtf$rsite.bins <- cut(tmtf$rsite, seq(0,1,.2))

### CHOOSE DATASETS
if ( healthy ) {
    ds <- tmtf$Dataset
    tmtf$Dataset[ds=="Healthy"] <- tmtf$TMT.Tissue[ds=="Healthy"]
    tmtf$Dataset[ds!="Healthy"] <- "Cancer"
}

### MEAN AND MEDIAN RAAS

## RAAS MEDIAN  PER DATA SET and UNIQUE BP/SAAP
##if ( only.unique ) {
    usaap <- paste0(tmtf$SAAP,"/",tmtf$BP,"/",tmtf$Dataset)
    araas <- split(tmtf$RAAS, usaap)
    
    ##araasl <- listProfile(araas, y=tmtf$RAAS, use.test=use.test, min=3)
    
    raas.emedian <- unlist(lapply(araas, median))
    raas.median <- unlist(lapply(araas, function(x) { log10(median(10^x)) } ))
    raas.emean <- unlist(lapply(araas, mean))
    raas.mean <- unlist(lapply(araas, function(x) { log10(mean(10^x)) } ))

    tmtf$unique <- usaap
    tmtf$RAAS.median <- raas.median[usaap]
    
    raas.bins <- cut(tmtf$RAAS.median, breaks=c(-6,-4,-2,-1,0,3))
    raas.srt <- levels(raas.bins)
    tmtf$raas.bins <- as.character(raas.bins)
    tmtf$raas.bins[is.na(tmtf$raas.bins)] <- "na"
##}
## TODO: albumin vs. globin
if ( FALSE ) {
    ovl <- clusterCluster(cl1=tmtf$extracellular, cl2=tmtf$albumin)
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt)
    dev.off()
}

### FILTER
if (  exclude.albumin ) {
    rmsaap <- tmtf$SAAP[tmtf$albumin]
    tmtf <- tmtf[!tmtf$SAAP%in%rmsaap,]
}
if ( exclude.frequent ) {
    rmsaap <- tmtf$SAAP[tmtf$from%in%frequent]
    tmtf <- tmtf[!tmtf$SAAP%in%rmsaap,]
}
if ( only.unique ) {
    ## just take the first of each SAAP/BP per Dataset
    tmtf <- tmtf[!duplicated(tmtf$unique),]
    tmtf$RAAS.orig <- tmtf$RAAS
    tmtf$RAAS <- tmtf$RAAS.median 
}
if (  exclude.extracellular ) {
    rmsaap <- tmtf$SAAP[tmtf$extracellular]
    tmtf <- tmtf[!tmtf$SAAP%in%rmsaap,]
}
if (  exclude.nterm ) {
    rmsaap <- tmtf$site <4
    tmtf <- tmtf[!rmsaap,]
}


## from this point on we only work with tmtf, which
## now should have all information

## sorting and labels
uds <- sort(unique(tmtf$Dataset))
uds <- c("Healthy",uds[uds!="Healthy"])
uds <- c("Cancer",uds[uds!="Cancer"])
uds <- uds[uds%in%tmtf$Dataset]

## additionally all and cancer-only
auds <- uds
if ( SETID=="cancer" )
    auds <- c("cancer", auds)
auds <- c("all", auds)

srt <- c("charged","polar","hphobic","special")
srt <- paste(rep(srt,each=4), rep(srt,4), sep=":")

axex <- ftlabels(srt) # axis labels with arrows


## RAAS COLORS
png(file.path(fig.path,paste0("legend_raas_vcols.png")),
    res=300, width=4, height=3, units="in")
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

png(file.path(fig.path,paste0("legend_raas_tcols.png")),
    res=300, width=4, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
traas.col <- selectColors(tmtf$RAAS,
                          mn=RAAS.MINT, mx=RAAS.MAXT,colf=COLF,
                          n=50, plot=TRUE,
                          mai=c(.5,.5,.1,.1),
                          xlab=expression(TMT~level~log[10]*RAAS))
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(colors, pos="bottomleft", cex=1)
figlabel(LAB, pos="bottomright", cex=1)
dev.off()


## tightest RAAS range!!
tcols <- traas.col$col
tbrks <- traas.col$breaks


## legend for dot plot
## RAAS COLORS
png(file.path(fig.path,paste0("legend_raas_acols.png")),
    res=300, width=4, height=3, units="in")
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





## legend for all two-sided statistics
png(file.path(fig.path,paste0("legend_wtests.png")),
    res=300, width=2, height=2, units="in")
par(mai=c(.6,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlapsLegend(p.min=p.min, p.txt=p.txt, type=2, col=ttcols)
dev.off()

## legend for dot plot
pp <- seq(0, -log10(p.dot), length.out=6)
rs <- seq(RAAS.MIN,RAAS.MAX, length.out=6)

pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))

colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- round(rs,1)
    
ovw <- list()
ovw$p.value <- t(10^-pm)
ovw$median <- t(rm)

plab <- expression(log[10](p))
if ( p.adjust=="q" )
    plab <- expression(log[10](q))

mai <- c(.5,.5,.1,.1)
fh <- fw <- .2
nh <- nrow(ovw$p.value) *fh + mai[1] + mai[3]
nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]

for ( colstyle in c("viridis","rocket","inferno","arno") ) {

    scols <- get(colstyle, mode="function")(length(vcols))
    
    segmenTools::plotdev(file.path(fig.path,paste0("legend_dotplot_",colstyle)),
                         height=nh, width=nw, res=300)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
    dotprofile(ovw, value="median",
               vbrks=vbrks,
               vcols=scols, 
               dot.sze=dot.sze, p.dot=p.dot, axis=1:2,
               ylab=plab,
               xlab=xl.raas)
    figlabel(colstyle, pos="bottomleft", cex=1)
    dev.off()
}
segmenTools::plotdev(file.path(fig.path,paste0("legend_dotplot_vcols")),
                     height=nh, width=nw, res=300)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovw, value="median",
           vbrks=vbrks,
           vcols=vcols, 
           dot.sze=dot.sze, p.dot=p.dot, axis=1:2,
           ylab=plab,
           xlab=xl.raas)
figlabel(colors, pos="bottomleft", cex=1)
dev.off()

## LEGENDS FOR DOT PLOTs

## broadest RAAS COLORS
pp <- seq(0, -log10(p.dot), length.out=3)
rs <- c(-4,-2,0,1) #seq(RAAS.MIN,RAAS.MAX, length.out=3)
pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))
colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- round(rs,1)
ovlg <- list(p.value=t(10^-pm), median=t(rm))

mai <- c(.4,.5,.05,.06)
fh <- fw <- .2
nh <- nrow(ovlg$p.value) *fh + mai[1] + mai[3]
nw <- ncol(ovlg$p.value) *fw + mai[2] + mai[4]

plotdev(file.path(fig.path,paste0("legend_dotplot_vcols_tight")),
                     height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovlg, value="median",
           vbrks=vbrks,
           vcols=vcols, 
           dot.sze=dot.sze, p.dot=p.dot, axis=1:2,
           ylab=plab,
           xlab=NA)
##mtext(xl.raas, 1, 1.1, adj=-.4)
text(1.5, -1, xl.raas, xpd=TRUE)
dev.off()

## slim
plotdev(file.path(fig.path,paste0("legend_dotplot_vcols_slim")),
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


## tight
mair <- mai
mair[2] <- .05
mair[3] <- .01
mair[4] <- .45
nh <- nrow(ovlg$p.value) *fh + mair[1] + mair[3]
nw <- ncol(ovlg$p.value) *fw + mair[2] + mair[4]
plotdev(file.path(fig.path,paste0("legend_dotplot_vcols_tight_right")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mair, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovlg, value="median",
           vbrks=vbrks,
           vcols=vcols, 
           dot.sze=dot.sze, p.dot=p.dot, axis=1,
           ylab=NA,
           xlab=NA)
##mtext(xl.raas, 1, 1.25, adj=2)
text(3, -1, xl.raas, xpd=TRUE)
axis(4, at=1:3, labels=c(-10,-5,0), las=2)
mtext(plab, 4, 1.5, adj=.2)
dev.off()


## tightest RAAS colors - legend for tcols/tbrks
pp <- seq(0, -log10(p.dot), length.out=3)
rs <- seq(RAAS.MINT,RAAS.MAXT, length.out=3)
pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))
colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- round(rs,1)
ovlg <- list(p.value=t(10^-pm), median=t(rm))

mai <- c(.4,.5,.05,.06)
mair <- mai
mair[2] <- .05
mair[3] <- .01
mair[4] <- .45
nh <- nrow(ovlg$p.value) *fh + mair[1] + mair[3]
nw <- ncol(ovlg$p.value) *fw + mair[2] + mair[4]
plotdev(file.path(fig.path,paste0("legend_dotplot_tcols_tight_right")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mair, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovlg, value="median",
           vbrks=tbrks,
           vcols=tcols, 
           dot.sze=dot.sze, p.dot=p.dot, axis=1,
           ylab=NA,
           xlab=NA)
##mtext(xl.raas, 1, 1.25, adj=2)
text(3, -1, xl.raas, xpd=TRUE)
axis(4, at=1:3, labels=c(-10,-5,0), las=2)
mtext(plab, 4, 1.5, adj=.2)
dev.off()

plotdev(file.path(fig.path,paste0("legend_dotplot_tcols_slim")),
        height=5, width=1, res=300, type=ftyp)
layout(1:2, heights=c(.5,.3))
par(mai=c(0,0.15,.15,.5), mgp=c(1.3,.3,0), tcl=-.25)
image_matrix(y=tbrks, t(t(abrks)), col=rev(tcols), xlab=NA, ylab=NA)
mtext(xl.raas, 4, 1.4)
axis(4, las=2)
ovp <- list(p.value=ovlg$p.value[,1,drop=FALSE],
            median=ovlg$median[,1,drop=FALSE])
par(mai=c(1,0.25,.15,.5))
dotprofile(ovp, value="median",
           vbrks=tbrks,
           vcols=tcols, 
           dot.sze=1.5*dot.sze, p.dot=p.dot, axis=NA,
           ylab=plab,
           xlab=NA, xpd=TRUE)
axis(4, at=nrow(ovp$p.value):1, labels=log10(ovp$p.value), las=2, col=NA)
mtext(plab, 1, 0.5)
dev.off()




## tigtht RAAS range - legend for acols/abreaks
pp <- seq(0, -log10(p.dot), length.out=3)
rs <- c(-4,-2,-1,0) #seq(RAAS.MIN,RAAS.MAX, length.out=3)
pm <- matrix(rep(pp, each=length(rs)), nrow=length(rs))
rm <- matrix(rep(rs, length(pp)), ncol=length(pp))
colnames(pm) <- colnames(rm) <- -pp
rownames(pm) <- rownames(rm) <- round(rs,1)
ovlg <- list(p.value=t(10^-pm), median=t(rm))

mai <- c(.4,.5,.05,.06)
fh <- fw <- .2
nh <- nrow(ovlg$p.value) *fh + mai[1] + mai[3]
nw <- ncol(ovlg$p.value) *fw + mai[2] + mai[4]

plotdev(file.path(fig.path,paste0("legend_dotplot_acols_tight")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovlg, value="median",
           vbrks=abrks,
           vcols=acols, 
           dot.sze=dot.sze, p.dot=p.dot, axis=1:2,
           ylab=plab,
           xlab=NA)
##mtext(xl.raas, 1, 1.1, adj=-.4)
text(1.5, -1, xl.raas, xpd=TRUE)
dev.off()

plotdev(file.path(fig.path,paste0("legend_dotplot_acols_slim")),
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




## UNIPROT <-> ENSEMBL MAPPING
## NOTE: ~90 duplicated ensembl IDs
uni2ens <- read.delim(uni2ens.file, header=FALSE)
uni2ens[,2] <- sub("\\..*", "", uni2ens[,2]) # remove ensembl version tag
uni2ens[,1] <- sub("-.*", "", uni2ens[,1]) # remove uniprot version tag

## filter for ensembl IDs in our data!
uni2dat <- uni2ens[uni2ens[,2]%in%c(dat$ensembl, dat$mane),]
## 1-many lists
uni2e <- split(uni2ens[,2], uni2ens[,1])
ens2u <- split(uni2ens[,1], uni2ens[,2])

refseq2ens <- read.delim(refseq.file, header=TRUE)

## remove duplicate ensembl IDs
## TODO: is the list sorted? best uniprot hit?
##uni2ens <- uni2ens[!duplicated(uni2ens[,2]),]

## UNIPROT <-> NAME MAPPING
## TODO: add MANE column
uni2nam <- read.delim(uni2nam.file, header=FALSE)
uni2nam[,1] <- sub("-.*", "", uni2nam[,1]) # remove uniprot version tag

## ENSEMBL <-> NAME MAPPING
## TODO: add MANE column
genes <- read.delim(feature.file)
genes <- genes[genes$proteins!="" & !is.na(genes$proteins),]
genes <- genes[genes$name!="" & !is.na(genes$name),]
ptl <- strsplit(genes$proteins, ";")
ens2nam <- rep(genes$name, lengths(ptl))
names(ens2nam) <- unlist(ptl)
pnms <- ens2nam

## list of all MANE proteins
MANES <- genes$MANE.protein
MANES <- MANES[MANES!=""]

## filter MANE.protein in mappings
rs2ens <- refseq2ens[refseq2ens[,1]%in%names(ens2nam),] # available
rs2mane <- refseq2ens[refseq2ens[,1]%in%MANES,] # MANE

### GLOBAL DISTRIBUTION BY CANCER TYPE
           
ylm <- range(tmtf$RAAS)
plotdev(file.path(fig.path,paste0("RAAS_distribution")),
        type=ftyp, res=300, width=3,height=2)
par(mai=c(.7,.5,.1,.1), mgp=c(1.2,.3,0), tcl=-.25, xaxs="i")
boxplot(tmtf$RAAS ~ factor(tmtf$Dataset, levels=uds), ylim=ylm, las=2,
        xlab=NA, ylab=xl.raaa, cex=.5, pch=19,
        pars=list(outcol="#00000055"), axes=FALSE)
for ( ax in c(2,4) ) axis(ax)
axis(1, at=1:length(uds), labels=uds, las=2)
figlabel(LAB, pos="bottomright", cex=.9)
dev.off()

plotdev(file.path(fig.path,paste0("RAAS_distribution_density")),
        type=ftyp, res=300, width=3,height=2)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
plot(density(tmtf$RAAS), ylim=c(0,.5), col=NA,
     xlab=xl.raaa, main=NA, xlim=ylm)
for ( i in seq_along(uds) )
    lines(density(tmtf$RAAS[tmtf$Dataset==uds[i]]), col=i, lwd=2)
legend("topright", uds, col=seq_along(uds), lty=1, seg.len=.5, lwd=2, bty="n",
       ncol=1, cex=.8, y.intersp=.9)
figlabel(LAB, pos="bottomleft", cex=.8)
dev.off()

if ( only.unique ) {

    tmtu <- tmtf[!duplicated(tmtf$unique),]
    plotdev(file.path(fig.path,paste0("RAAS_distribution_unique")),
            type=ftyp, res=300, width=3,height=2)
    par(mai=c(.7,.5,.1,.1), mgp=c(1.2,.3,0), tcl=-.25, xaxs="i")
    boxplot(tmtu$RAAS.median ~ factor(tmtu$Dataset, levels=uds),
            ylim=ylm, las=2, xlab=NA, ylab=xl.raau,
            cex=.5, pch=19, pars=list(outcol="#00000055"), axes=FALSE)
    for ( ax in c(2,4) ) axis(ax)
    axis(1, at=1:length(uds), labels=uds, las=2)
    figlabel(LAB, pos="bottomright", cex=.9)
    dev.off()
    
    plotdev(file.path(fig.path,paste0("RAAS_distribution_unique_density")),
            type=ftyp, res=300, width=3,height=2)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
    plot(density(tmtu$RAAS.median), ylim=c(0,.5), col=NA,
         xlab=NA, main=NA, xlim=ylm)
    mtext(xl.raau, 1, 1.5)
    for ( i in seq_along(uds) )
        lines(density(tmtu$RAAS.median[tmtu$Dataset==uds[i]]), col=i, lwd=2)
    legend("topright", uds, col=seq_along(uds),
           lty=1, seg.len=.5, lwd=2, bty="n",
           ncol=1, cex=.8, y.intersp=.9)
    figlabel(LAB, pos="bottomleft", cex=.8)
    dev.off()
}

## TODO: pairs of same SAAP
##tmtl <- split(tmtf$Dataset, tmtf$SAAP)


### MOTIFS & KRAQ
## TODO: make sure this doesnt overwrite any of the above
## it should mainly providebdat

## additional input
pep.file <- file.path(dat.path, "All_main_tryptic_peptide_list.txt")
non.file <- file.path(dat.path, "All_main_nontryptic_peptide_list.txt")
nkr.file <- file.path(dat.path, "All_main_nontryptic_noKR_peptide_list.txt")
akr.file <- file.path(dat.path, "All_main_ArgC_LysC_peptide_list.txt")


## FILTER UNIQUE BP - since those have the same AA CONTEXT
## NOTE: that this looses different SAAPs for the same BP

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

### NOTE: check use of bdat vs. hdat, where bdat has one less
## row, due to missing RAAS values.


### PLOT ALIGNMENT

CMAIL <- 1.54 ## commonly used between motif and domain figures, defined 
              ## in domain script via nchar to fit the y-axis labels
FONT <- "monospace" # font used for aligned figures in motifs and domains

