
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

hist(psi$mm.Direct)

## convert chromosomes ID to index via segmenTools'
## chromosome coordinate indexing system.
psi$chr <-  chrIdx[sub("^M$","MT",sub("chr","",psi$chr))]
bdat$chr <- chrIdx[bdat$chr]

psidx <- coor2index(psi[,c("chr","coor")], chrS=chrS)
bd2dx <- coor2index(bdat[,c("chr","coor")], chrS=chrS)
bd1 <- cbind(chr=bdat[,c("chr")],
             coor=bdat[,c("coor")] + ifelse(bdat[,"strand"]=="+",-1,1))
bd1dx <- coor2index(bd1, chrS=chrS)
bd3 <- cbind(chr=bdat[,c("chr")],
             coor=bdat[,c("coor")] - ifelse(bdat[,"strand"]=="+",-1,1))
bd3dx <- coor2index(bd3, chrS=chrS)

##bdat[which(bd2dx[,2]%in%intersect(psidx[,2], bd2dx[,2])),"name"]

##dstm <- dist(psidx[,2], bd2dx[,2])
bid <- which(bd3dx[,2]%in%intersect(psidx[,2], bd3dx[,2]))
pid <- which(psidx[,2]%in%intersect(psidx[,2], bd3dx[,2]))

unique(bdat[bid,"name"])
bdat[bid,"RAAS"]

unique(psi[pid,"Annotation"])
psi[pid,"mm.Direct"]

table(bdat[bid,"fromto"])

##bdat[bid,c("chr","coor","strand")]
##psi[pid,c("chr","coor")]

length(bid)
length(pid)
