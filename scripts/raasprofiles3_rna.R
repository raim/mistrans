
## comparing AAS with mRNA modifications

library(segmenTools)
library(readxl)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

rfig.path <- file.path(fig.path,"rna")
dir.create(rfig.path, showWarnings=FALSE)

## use only AAS where we have a codon
usite <- site[!is.na(site$codon),]


### ADDITIONAL DATA

## coding sequence fasta
tfas.file <- file.path(mam.path, "processedData","coding.fa")
## protein-transcript map
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
## protein lengths
plen.file <- file.path(mam.path, "processedData", "protein_length.tsv")
## chromosome length index
chr.file <- file.path(mam.path,"chromosomes","sequenceIndex.csv")

## RNA pseudouridylation data sent by
## Oleksandra Fanari <fanari.o@northeastern.edu>
psi.file <- file.path(dat.path, "six_cell_lines_minimal(A549).csv")

## @Zhang2023 - PRAISE to find psi sites
zh23.file <- file.path(dat.path, "zhang23_stables_1-6.xlsx")
## @Dai2023 - BID-seq
dai.file <- file.path(dat.path, "dai23_stables_1-23.xlsx")

## prepare AAS data and coordinates

## get all protein lengths
plengths <- read.delim(plen.file, header=FALSE, row.names=1)
plengths <- setNames(plengths[,1], rownames(plengths))
names(plengths) <- sub("\\.[0-9]*","", names(plengths))

## get chromosome index
chrMap <- read.delim(chr.file)
chrS <- c(0,cumsum(as.numeric(chrMap[,3])))
chrIdx <- chrMap[,1]
names(chrIdx) <- chrMap[,2]

## convert AAS chromosome to index
usite$chr <- chrIdx[usite$chr]

## AAS 2nd codon position
bd2dx <- coor2index(usite[,c("chr","coor")], chrS=chrS)[,2]
## AAS 1st codon position
bd1 <- cbind(chr=usite[,c("chr")],
             coor=usite[,c("coor")] + ifelse(usite[,"strand"]=="+",-1,1))
bd1dx <- coor2index(bd1, chrS=chrS)[,2]
## AAS 3rd codon position
bd3 <- cbind(chr=usite[,c("chr")],
             coor=usite[,c("coor")] - ifelse(usite[,"strand"]=="+",-1,1))
bd3dx <- coor2index(bd3, chrS=chrS)[,2]

### COLLECT VARIOUS PSI DATA
PSI <- NULL

### RNA pseudouridylation data by Zhang et al.
## NOTE: sheet 2 has only 1 overlap (codons 2 and 3 in one psi site),
##       and sheet 3 has none.

## @Zhang2023 Supp. Dataset 2
## sheet 2: Ψ identified by PRAISE in HEK293T mRNAs and ncRNAs
psiz <- as.data.frame(read_xlsx(zh23.file, sheet=2, skip=2))

## calculate mean deletion ratio of two replications
psiz$deletion_ratio <- apply(psiz[,c("rep1_deletion_ratio",
                                     "rep2_deletion_ratio")],1,mean)


## parse chromosoe coordiantes - "PRAISE-tools firstly got all possible
## mapped positions to reference transcriptome (GRCh38)"

## first, through away all w/o chromosome mapping
## TODO: alternatively use transcript mapping
psiz <- psiz[psiz$chr_site!="NONE",]

coors <- strsplit(psiz$chr_site, "_")
coors <- lapply(lapply(coors, strsplit, "-"), unlist)
## add coordinate for single coor positions
coors <- lapply(coors, function(x) {if(length(x)==2) x<-c(x,x[2]); x})
## handling location tags of the form "chr17_KV575245v1_fix" 
coors <- lapply(coors, function(x) {if(length(x)==4) x<-c(x,x[4]); x})
fixtag <- unlist(lapply(coors, length))==5
coors[fixtag] <- lapply(coors[fixtag], function(x) x[c(1,4,5)])
coors <- as.data.frame(cbind(do.call(rbind, coors), fix=fixtag))
coors[,1] <- sub("^chr","", coors[,1])
coors[,2:3] <- apply(coors[,2:3], 2, as.numeric)

colnames(coors) <- c("chr","start","end","fix")

## add strand info: testing some coors, it seems that end<start
## implies negative strand
coors$strand <- ifelse(coors[,3]>=coors[,2], "+","-")

## collect relevant columns
psi <- cbind(coors[,c("chr","start","end","strand")],
             info="",
             gene=psiz$gene_name,
             psi=psiz$deletion_ratio, source="zhang23_stable2")
PSI <- rbind(PSI, psi)
      

## @Zhang2023 Supp. Dataset 3
## sheet 3: PUS-dependent Ψ list in HEK293T mRNAs and ncRNAs
psiz <- as.data.frame(read_xlsx(zh23.file, sheet=3, skip=2))

## calculate mean deletion ratio of two replications
psiz$deletion_ratio <- apply(psiz[,c("ko_rep1_deletion_ratio",
                                     "ko_rep2_deletion_ratio")],1,mean)


## parse chromosoe coordiantes - "PRAISE-tools firstly got all possible
## mapped positions to reference transcriptome (GRCh38)"

## first, through away all w/o chromosome mapping
## TODO: alternatively use transcript mapping
psiz <- psiz[psiz$chr_site!="NONE",]

coors <- strsplit(psiz$chr_site, "_")
coors <- lapply(lapply(coors, strsplit, "-"), unlist)
## add coordinate for single coor positions
coors <- lapply(coors, function(x) {if(length(x)==2) x<-c(x,x[2]); x})
## handling location tags of the form "chr17_KV575245v1_fix" 
coors <- lapply(coors, function(x) {if(length(x)==4) x<-c(x,x[4]); x})
fixtag <- unlist(lapply(coors, length))==5
coors[fixtag] <- lapply(coors[fixtag], function(x) x[c(1,4,5)])
coors <- as.data.frame(cbind(do.call(rbind, coors), fix=fixtag))
coors[,1] <- sub("^chr","", coors[,1])
coors[,2:3] <- apply(coors[,2:3], 2, as.numeric)

colnames(coors) <- c("chr","start","end","fix")

## add strand info: testing some coors, it seems that end<start
## implies negative strand
coors$strand <- ifelse(coors[,3]>=coors[,2], "+","-")


## collect relevant columns
psi <- cbind(coors[,c("chr","start","end","strand")],
             info=psiz$enzyme_dependency,
             gene=psiz$gene_name,
             psi=psiz$deletion_ratio, source="zhang23_stable3")
PSI <- rbind(PSI, psi)
      

### @Dai2023: BID-seq
## sheet/STable 5: Ψ on wild-type HeLa mRNA (>10% fraction)
## sheet/Stable 6: Ψ on wild-type HEK293T mRNA (>10% fraction)
## sheet/STable 7: Ψ on wild-type A549 mRNA (>10% fraction)

ids <- c("HeLa", "HEK293T", "A549")
tcol <- "Frac_Ave %" # TODO: use deletion average or fraction?
for ( i in 1:3 ) {
    psiz <- as.data.frame(read_xlsx(dai.file, sheet=i+4, skip=3))

    psiz$start <- psiz$end <- psiz$pos
    psiz$chr <- sub("^chr","", psiz$chr)

    ## only take CDS
    psiz <- psiz[psiz$seg=="CDS",]
    
    psi <- cbind(psiz[,c("chr","start","end","strand")],
                 info=ids[i],
                 gene=psiz$name,
                 psi=psiz[,tcol]/100,
                 source=paste0("dai23_stable",i+4))
    PSI <- rbind(PSI, psi)
}

### RNA pseudouridylation data sent by
### Oleksandra Fanari <fanari.o@northeastern.edu>

psi <- read.csv(psi.file)
colnames(psi) <- sub("position","coor", colnames(psi))

## match gene name
## TODO: find  missing 440
psi$proteinID <- names(ens2nam[match(psi$Annotation,  ens2nam)]) 
psi$len <- plengths[psi$proteinID]


## convert chromosomes ID to index via segmenTools'
## chromosome coordinate indexing system.
psi$chr <-  sub("^M$","MT",sub("chr","",psi$chr))
psi$start <- psi$end <- psi$coor

psiz <- cbind(psi[,c("chr","start","end")],
              strand=".",
              info="",
              gene=psi$Annotation,
              psi=psi[,"mm.Direct"]/100, # % -> fraction
              source="fanari24pre")

PSI <- rbind(PSI, psiz)

## get ensembl protein ID
PSI$ensembl <- names(ens2nam[match(PSI$gene,  ens2nam)]) 


## take mean position - only affects @Zhang which provides ranges,
## and we found only one overlap (2nd and 3rd position of one AAS)
PSI$coor <- round(apply(PSI[,c("start","end")],1,mean))

## chromosome index
PSI$chr <-  chrIdx[PSI$chr]

psidx <- coor2index(PSI[,c("chr","coor")], chrS=chrS)[,2]

## collect overlapping sites
psite <- usite
for ( i in 1:3 ) {

    bdidx <- get(paste0("bd",i,"dx"))
                 
    ## find same positions
    bid <- which(bdidx%in%intersect(psidx, bdidx))
    ##pid <- which(psidx%in%intersect(psidx, bdidx))
    cat(paste("found", length(bid), "overlapping sites", i, "\n"))
    
    psi2site <- match(bdidx, psidx)
    
    sitepsi <- PSI[psi2site, ]
    colnames(sitepsi) <- paste0("codon",i, "_", colnames(sitepsi))
    psite <- cbind(psite, sitepsi)
}
## collapse all psi
## NOTE: using median as a helper
psite$psi <- apply(psite[,paste0("codon",1:3,"_psi")], 1, median, na.rm=TRUE)
psite$psi.source <- apply(psite[,paste0("codon",1:3,"_source")], 1,
                          function(x) paste0(x[!is.na(x)],collapse=";"))

## CORRELATION OF PSI % TO RAAS

sources <- unique(psite$psi.source)
sources <- sources[sources!=""]
pch.source <- 1:length(sources) 
names(pch.source) <- sources

for ( i in 1:3 ) {
    ccol <- paste0("codon",i,"_psi")
    scol <- paste0("codon",i,"_source")
    plotdev(file.path(rfig.path,paste0("raas_psi_codon",i)), type=ftyp,
            height=2.5, width=2.5, res=200)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(psite$RAAS.median, psite[,ccol], density=FALSE,
#            pch=pch.source[psite[,scol]],
            xlab=xl.raas, ylab=paste0("codon position ",i,", psi/%"), col=NA)
    points(psite$RAAS.median, psite[,ccol], pch=pch.source[psite[,scol]],
           cex=.75, col=i)
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

## HYPERGEO, USING NUMBER OF Us AS BACKGROUND


## COUNT ONLY Us in CODING SEQUENCE
tfas <- readFASTA(tfas.file, grepID=TRUE)
## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
## reverse map transcript-protein
pamrt <- matrix(rownames(trmap), ncol=1)
rownames(pamrt) <- trmap[,1]
## rename by protein names via trmap, for quick access
## of protein-specific transcript
names(tfas) <- pamrt[names(tfas),1]

pids <- unique(c(psite$ensembl, PSI$ensembl))
tfas <- tfas[pids[pids%in%names(tfas)]]

## total balls: all Us in coding sequences of the total gene set
totaa <- sum(unlist(lapply(tfas,
                           function(x) sum(unlist(strsplit(x$seq,""))=="T"))))

## white balls: all detected psi
m <- length(unique(paste(psi$chr, psi$coor)))
n <- tot-m # black balls

## balls drawn: all U in AAS codons
## take only unique site here!
filt <- !duplicated(paste(psite$chr,
                          psite$coor, psite$strand))]
k <- sum(unlist(strsplit(psite$codon[filt],""))=="T")

## q: white balls drawn number of AAS with psi
q <- sum(!is.na(psite$psi[!duplicated(paste(psite$chr,
                                            psite$coor, psite$strand))]))

##m: the number of white balls in the urn.
##n: the number of black balls in the urn.
##k: the number of balls drawn from the urn

phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)
