
## comparing AAS with mRNA modifications

library(segmenTools)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

rfig.path <- file.path(fig.path,"rna")
dir.create(rfig.path, showWarnings=FALSE)

## additional data
## protein lengths
plen.file <- file.path(mam.path, "processedData", "protein_length.tsv")

## chromosome length index
chr.file <- file.path(mam.path,"chromosomes","sequenceIndex.csv")

## get all protein lengths
plengths <- read.delim(plen.file, header=FALSE, row.names=1)
plengths <- setNames(plengths[,1], rownames(plengths))
names(plengths) <- sub("\\.[0-9]*","", names(plengths))

## get chromosome index
chrMap <- read.delim(chr.file)
chrS <- c(0,cumsum(as.numeric(chrMap[,3])))
chrIdx <- chrMap[,1]
names(chrIdx) <- chrMap[,2]

## RNA pseudouridylation data sent by
## Oleksandra Fanari <fanari.o@northeastern.edu>
psi.file <- file.path(dat.path, "six_cell_lines_minimal(A549).csv")

psi <- read.csv(psi.file)
colnames(psi) <- sub("position","coor", colnames(psi))

## match gene name
## TODO: find  missing 440
psi$proteinID <- names(ens2nam[match(psi$Annotation,  ens2nam)]) 
psi$len <- plengths[psi$proteinID]

## distribution of psi percentages
hist(psi$mm.Direct)

## convert chromosomes ID to index via segmenTools'
## chromosome coordinate indexing system.
psi$chr <-  chrIdx[sub("^M$","MT",sub("chr","",psi$chr))]
site$chr <- chrIdx[site$chr]

psidx <- coor2index(psi[,c("chr","coor")], chrS=chrS)

## AAS 2nd codon position
bd2dx <- coor2index(site[,c("chr","coor")], chrS=chrS)
## AAS 1st codon position
bd1 <- cbind(chr=site[,c("chr")],
             coor=site[,c("coor")] + ifelse(site[,"strand"]=="+",-1,1))
bd1dx <- coor2index(bd1, chrS=chrS)
## AAS 3rd codon position
bd3 <- cbind(chr=site[,c("chr")],
             coor=site[,c("coor")] - ifelse(site[,"strand"]=="+",-1,1))
bd3dx <- coor2index(bd3, chrS=chrS)

psite <- site
for ( i in 1:3 ) {

    bdidx <- get(paste0("bd",i,"dx"))
                 
    ## find same positions
    bid <- which(bdidx[,2]%in%intersect(psidx[,2], bdidx[,2]))
    pid <- which(psidx[,2]%in%intersect(psidx[,2], bdidx[,2]))
    
    psi2site <- match(bdidx[,2], psidx[,2])
    
    sitepsi <- psi[psi2site, c("mm.Direct", "p.value.Direct")]
    colnames(sitepsi) <- paste0("codon",i, "_", c("psi","psi.p"))
    psite <- cbind(psite, sitepsi)

    unique(site[bid,"name"])
    site[bid,"RAAS"]
    
    unique(psi[pid,"Annotation"])
    psi[pid,"mm.Direct"]
    
    table(site[bid,"fromto"])
    
    ##site[bid,c("chr","coor","strand")]
    ##psi[pid,c("chr","coor")]
    
    length(bid)
    length(pid)
}
## collapse all psi
## NOTE: using median as a helper
psite$psi <- apply(psite[,paste0("codon",1:3,"_psi")], 1, median, na.rm=TRUE)

## correlation of psi % to RAAS?

for ( i in 1:3 ) {
    ccol <- paste0("codon",i,"_psi")
    plotdev(file.path(rfig.path,paste0("raas_psi_codon",i)), type=ftyp,
            height=2.5, width=2.5, res=200)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(psite$RAAS.median, psite[,ccol], density=FALSE,
        xlab=xl.raas, ylab=paste0("codon position ",i,", psi/%"), col=i)
    dev.off()
}



## correlation of RAAS to psi %
plotdev(file.path(rfig.path,"raas_psi"), type=ftyp,
        height=2.5, width=2.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(psite$RAAS.median, psite$psi, density=FALSE,
        xlab=xl.raas, ylab="psi/%", pch=NA)
points(psite$RAAS.median, psite$codon1_psi, col=1, pch=19, cex=.5)
points(psite$RAAS.median, psite$codon2_psi, col=2, pch=19, cex=.5)
points(psite$RAAS.median, psite$codon3_psi, col=3, pch=19, cex=.5)
legend("topright",paste(1:3), title="codon position", col=1:3, cex=.5,
       pch=19, pt.cex=.75, bty="n")
dev.off()

## HYPERGEO, USING NUMBER OF NUCLEOTIDES AS BACKGROUND

## TODO: remove those psi w/o protein name match?
totaa <- sum(plengths[unique(c(psite$ensembl, psi$proteinID))],na.rm=TRUE)

## total balls
tot <- totaa*3 # all nucleotides coding for the psi+AAS subset of proteins

## white balls: psi
m <- length(unique(paste(psi$chr, psi$coor)))
n <- tot-m # black balls

## balls drawn: AAS * 3 for codons
k <- length(unique(paste(psite$chr, psite$coor, psite$strand)))*3

## q: white balls drawn number of AAS with psi
q <- sum(!is.na(psite$psi[!duplicated(paste(psite$chr,
                                            psite$coor, psite$strand))]))

##m: the number of white balls in the urn.
##n: the number of black balls in the urn.
##k: the number of balls drawn from the urn

phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)
