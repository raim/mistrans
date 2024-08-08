
## comparing AAS with mRNA modifications

library(segmenTools)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

rfig.path <- file.path(out.path,"rna")
dir.create(rfig.path, showWarnings=FALSE)

## additional data

## chromosome length index
chr.file <- file.path(mam.path,"chromosomes","sequenceIndex.csv")

chrMap <- read.delim(chr.file)
chrS <- c(0,cumsum(as.numeric(chrMap[,3])))
chrIdx <- chrMap[,1]
names(chrIdx) <- chrMap[,2]

## RNA pseudouridylation data sent by
## Oleksandra Fanari <fanari.o@northeastern.edu>
psi.file <- file.path(dat.path, "six_cell_lines_minimal(A549).csv")

psi <- read.csv(psi.file)

colnames(psi) <- sub("position","coor", colnames(psi))

## convert chromosomes ID to index
psi$chr <-  chrIdx[sub("^M$","MT",sub("chr","",psi$chr))]
bdat$chr <- chrIdx[bdat$chr]

psidx <- coor2index(psi[,c("chr","coor")], chrS=chrS)
bdidx <- coor2index(bdat[,c("chr","coor")], chrS=chrS)

bdat[which(bdidx[,2]%in%intersect(psidx[,2], bdidx[,2])),"name"]

##dstm <- dist(psidx[,2], bdidx[,2])
bid <- which(bdidx[,2]%in%intersect(psidx[,2], bdidx[,2]))
pid <- which(psidx[,2]%in%intersect(psidx[,2], bdidx[,2]))

unique(bdat[bid,"name"])
unique(psi[pid,"Annotation"])
table(bdat[bid,"fromto"])

bdat[bid,c("chr","coor")]
psi[pid,c("chr","coor")]
