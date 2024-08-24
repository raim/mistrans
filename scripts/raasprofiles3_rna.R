
## comparing AAS with mRNA modifications

library(segmenTools)
library(readxl)
library(basicPlotteR) # for non-overlapping text

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

rfig.path <- file.path(fig.path,"rna")
dir.create(rfig.path, showWarnings=FALSE)

## use only unique AAS sites where we have a codon
usite <- asite[!is.na(asite$codon),]

## use only AAS mapped to a MANE transcripts
usite <- usite[usite$transcript%in%genes$MANE,]

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
psi.file <- file.path(dat.path, "six_cell_lines_minimal.xlsx")#(A549).csv")


## RNAMD - mA6 modifications
rmd.file <- file.path(mam.path, "originalData", "human_all_logFC.txt.gz")
rma.file <- file.path(mam.path, "originalData", "hg38_CL_Tech.txt.gz")

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


### RNA pseudouridylation data sent by
### Oleksandra Fanari <fanari.o@northeastern.edu>
psi <- NULL
for ( i in excel_sheets(psi.file) ) {
    
    tmp <-as.data.frame(read_xlsx(psi.file, sheet=i))

    ## collapse duplicate sites
    ## check whether
    upos <- paste(tmp$chr, tmp$position)
    dpos <-upos[which(duplicated(upos))]

    ## all annotations as list
    danl <- sapply(dpos, function(x) tmp$Annotation[which(upos==x)])
    ## filter those that are present in protein name vector
    danl <- lapply(danl, function(x) x[x%in%pnms])

    ## TODO: further reduce to a single annotation
    ## e.g. by length

    dann <- sapply(dpos, function(x)
        paste(tmp$Annotation[which(upos==x)], collapse=";"))
    fpos <- sapply(dpos, function(x)
        which(upos==x)[1])

    ## add annotation lists to be filtered later
    tmp$Annotation.all <- tmp$Annotation.best <- tmp$Annotation.first <-
        tmp$Annotation
    tmp$Annotation.all[fpos] <- dann
    tmp$Annotation.best[fpos] <- unlist(lapply(danl, paste, collapse=";"))
    tmp$Annotation.first[fpos] <- unlist(lapply(danl, function(x) x[1]))

    ## remove all duplicated
    tmp <- tmp[!duplicated(upos),]
    
    psi <- rbind(psi, cbind(tmp, cell=i))
}
colnames(psi) <- sub("position","coor", colnames(psi))



## match gene name
## TODO: find  missing 440
##psi$proteinID <- names(ens2nam[match(psi$Annotation,  ens2nam)]) 
##psi$len <- plengths[psi$proteinID]


## convert chromosomes ID to index via segmenTools'
## chromosome coordinate indexing system.
psi$chr <-  sub("^M$","MT",sub("chr","",psi$chr))
psi$start <- psi$end <- psi$coor

psiz <- cbind(psi[,c("chr","start","end")],
              strand=".", info="",
              gene=psi$Annotation.first,
              psi=psi[,"mm.DirectMINUSmm.IVT"]/100, # % -> fraction
              source=paste0("fanari24pre","_",psi$cell))

## TODO: track different annotations!

##  e.g. psia["NA.46990",] annotated as RP1-120G22.12
## but this has coors chr1 6264900-6265840 while psi site
## is annotated as chr1 6197653

## line 26694: RP1-120G22.12	chr1	6197653	TTTTG	4013	13	2249	6	4062	2282	0.266075388026608	0.322901142573274	0.266075388026608	0	1	0.056825754546666	4026	2255	1	RP1-120G22.12chr16197653	1	pIVT	chr16197653

## line 27997: RPL22	chr1	6197653	TTTTG	4013	13	2249	6	4062	2282	0.266075388026608	0.322901142573274	0.266075388026608	0	1	0.056825754546666	4026	2255	1	RPL22chr16197653	1	pIVT	chr16197653

## FILTER
##psiz <- psiz[psiz$psi>.1,]

PSI <- rbind(PSI, psiz)


### PROCESS PSI DATA

## get ensembl protein ID
PSI$ensembl <- names(ens2nam[match(PSI$gene,  ens2nam)]) 

## TODO: find in synonyms, trace where NA come from (find synonyms)
## already while parsing and filtering multiple hits.
PSI$gene[is.na(PSI$ensembl)]

## take mean position - only affects @Zhang which provides ranges,
## and we found only one overlap (2nd and 3rd position of one AAS)
PSI$coor <- round(apply(PSI[,c("start","end")],1,mean))

## chromosome index
PSI$chr <-  chrIdx[PSI$chr]


### PSI DISTRIBUTIONS

plotdev(file.path(rfig.path, "psi_boxplots"),
        width=2.5, height=2.5, res=200, type=ftyp)
par(mai=c(.8,.5,.5,.25), mgp=c(1.4,.3,0), tcl=-.25)
bp <- boxplot(PSI$psi ~ sub(".*_","",PSI$source),
              las=2, ylim=c(0,1), xlab=NA, ylab=NA, cex=.5)
mtext(expression(psi~fraction), 2, 1.4)
axis(3, at=1:length(bp$n), labels=bp$n, las=2, cex.axis=.8)
abline(h=.3, col=2, lty=1)
axis(4, at=.3, label="suggested cutoff", col.axis=2, col=2, mgp=c(1.3,.2,0))
##abline(h=1, col=2, lty=3)
dev.off()


### MAP PSI and AAS SITES
psidx <- coor2index(PSI[,c("chr","coor")], chrS=chrS)[,2]

## REANNOTATE PSI, EACH TO ENSEMBL TRANSCRIPTS
cds.file <- "/home/raim/data/mammary/originalData/transcript_coordinates.tsv"
cds <- read.delim(cds.file)
cds <- cds[cds$transcript%in%genes$MANE,]
cdsi <- cds[,c("transcript","chr","start","end")] # SKIP STRAND!
cdsi$chr <- chrIdx[cdsi$chr]
cdsi <- coor2index(cdsi, chrS=chrS)

## find PSI sites in ensembl MANE transcript coordinates
## (ignoring strand)
## TODO: more efficient via cut or findInterval,see
## cut(x, breaks = interval.vector, include.lowest = TRUE)
psiu <- unique(psidx) # TODO: search only unique sites and map back
psi2ens <- sapply(psiu, function(x) which(cdsi$start <= x &  cdsi$end >= x))

psi2trans <- lapply(psi2ens,
                    function(x) {
                        if ( length(x)>0 ) {
                            x <- cdsi[x,"transcript"]
                        } else { x <- NA }
                        x
                    })

## investigate dual hits
if ( FALSE ) {
    table(lengths(psi2ens))
    unlist(psi2trans[which(lengths(psi2ens)==2)])%in%genes$MANE
}

## expand back to site
PSI$transcripts <- unlist(lapply(psi2trans[match(psidx,psiu)],
                                function(x) paste(x, collapse=";")))

## take only first
PSI$transcript <- unlist(lapply(psi2trans[match(psidx,psiu)],
                                function(x) x[1]))

PSI$transcript[PSI$transcript=="NA"] <- NA

## TODO: find transcripts for already mapped proteins
## and find out why we don't find them directly

## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
## reverse map transcript-protein
pamrt <- matrix(rownames(trmap), ncol=1)
rownames(pamrt) <- trmap[,1]


## add missing transcripts via proteins
## TODO: this assumes transcripts that may not cover the psi site
PSI$transcript.missing <- is.na(PSI$transcript)
PSI$transcript.annotated <- trmap[PSI$ensembl,1]

### ONLY CONSIDER PSI SITES THAT WE COULD MAP TO
### MANE ENSEMBL TRANSCRIPTS
cat(paste("NOTE: REMOVING", sum(is.na(PSI$transcript)), "of", nrow(PSI),
          "PSI SITES w/o MANE TRANSCRIPT MATCH\n"))

PSI <- PSI[!is.na(PSI$transcript),]
### RE-MAP PSI SITES COORDINATES
psidx <- coor2index(PSI[,c("chr","coor")], chrS=chrS)[,2]

## collect overlapping sites
psite <- usite
psia <- PSI
for ( i in 1:3 ) {

    bdidx <- get(paste0("bd",i,"dx"))
                 
    ## find same positions
    bid <- which(bdidx%in%intersect(psidx, bdidx))
    ##pid <- which(psidx%in%intersect(psidx, bdidx))
    cat(paste("found", length(bid), "overlapping sites", i, "\n"))
    
    ## annotate AAS
    psi2site <- match(bdidx, psidx)
    sitepsi <- PSI[psi2site, ]
    sitepsi$chr <-
        gsub("chrNA","",paste0("chr", chrMap[sitepsi$chr,2]))
    sitepsi$chr <- gsub(";","",sitepsi$chr)
    sitepsi$chr[sitepsi$chr==""] <- NA
    
    colnames(sitepsi) <- paste0("codon",i, "_", colnames(sitepsi))
    psite <- cbind(psite, sitepsi)

    ## annotate psi sites
    aas2psi <- match(psidx, bdidx)
    aasite <- psite[aas2psi,c("transcript","fromto","codon","RAAS.median")]
    colnames(aasite) <- paste0("codon",i, "_", colnames(aasite))
    psia <- cbind(psia, aasite)
}

## collapse all psi
## NOTE: using median as a helper

pcodons <- function(x) {
    y <- which(!is.na(x))
    if (length(y)==0) y<-NA;
    paste(y[!is.na(y)],collapse=";")
}
pvals <- function(x) paste0(x[!is.na(x)],collapse=";")
pmax <- function(x) {
    y <- x[!is.na(x)]
    ifelse(length(y)==0, NA, max(y))
}

## add summary to AAS table
psite$psi.n <- apply(psite[,paste0("codon",1:3,"_psi")], 1,
                     function(x) sum(!is.na(x)))
psite$psi.max <- apply(psite[,paste0("codon",1:3,"_psi")], 1, pmax)
psite$psi <- apply(psite[,paste0("codon",1:3,"_psi")], 1, pvals)
psite$psi.codonpos <- apply(psite[,paste0("codon",1:3,"_psi")], 1,  pcodons)
psite$psi.source <- apply(psite[,paste0("codon",1:3,"_source")], 1, pvals)
psite$psi.chr <- apply(psite[,paste0("codon",1:3,"_chr")], 1, pvals)
psite$psi.coor <- apply(psite[,paste0("codon",1:3,"_coor")], 1, pvals)

## add summary to PSI table
psia$aas.n <- apply(psia[,paste0("codon",1:3,"_RAAS.median")], 1,
                    function(x) sum(!is.na(x)))
psia$aas.transcript <- apply(psia[,paste0("codon",1:3,"_transcript")], 1, pvals)

### NOTE: mitochondrial transcripts missing


## export site for Sara and Sasha
## TODO: rm separate codon columns, and instead add affected
## codon position:
ecols <- c("chr", "coor", "strand", "RAAS.median", "name",
           "ensembl", "transcript",
           "pos", "codon",
           "fromto", "psi",
           "psi.codonpos", "psi.chr", "psi.coor", "psi.source")

expsite <- psite[grep("fanari24pre",psite$psi.source),ecols]
expsite$chr <- paste0("chr", chrMap[expsite$chr,2])
#expsite$psi.chr <- paste0("chr", chrMap[expsite$psi.chr,2])

expsite$psi.source <- gsub("fanari24pre_","",expsite$psi.source)

write.table(file=file.path(rfig.path,"aas_vs_psi.tsv"),
            x=expsite,
            sep="\t", quote=FALSE, row.names=FALSE)


## CORRELATION OF PSI % TO RAAS

sources <- unique(na.omit(unlist(psite[,paste0("codon",1:3,"_source")])))
pch.source <- 1:length(sources) 
names(pch.source) <- sources

## for each source/cell line separately
cors <- list()
for ( src in sources ) {
    ##ssite <- psite[grep(src,psite$psi.source),]
    src <- sub(".*_","",src)
    cors[[src]] <- list()
    for ( i in 1:3 ) {

        
        ccol <- paste0("codon",i,"_psi")
        scol <- paste0("codon",i,"_source")

        ssite <- psite[grep(src,psite[,scol]),]
        nvals <- sum(!is.na(ssite[,ccol]))

        pc <- NA

        plotdev(file.path(rfig.path,paste0("psi_raas_codon",i,"_",src)),
                type=ftyp, height=2.5, width=2.5, res=200)
        par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
        if (nvals==0) {
            plot(1, col=NA, xlab=NA, ylab=NA, axes=FALSE)
            box()
        } else if ( nvals>2) {
            pc <- plotCor(ssite$RAAS.median, ssite[,ccol], density=FALSE,
                          xlab=NA, ylab=NA, signif=2, col=NA,
                          legcex=.8)$p
            
        } else {
            plot(ssite$RAAS.median, ssite[,ccol], col=NA,
                 xlab=NA, ylab=NA)
        }
        if ( nvals>0 ) {
            points(ssite$RAAS.median, ssite[,ccol],
                   pch=1,#pch.source[ssite[,scol]],
                   cex=.75, col=i)
            basicPlotteR::addTextLabels(ssite$RAAS.median, ssite[,ccol],
                                        labels=ssite$fromto, #xpd=TRUE,
                                        col.label=1, col.line=i,
                                        cex.label=.5)#, pos=4)
        }
        mtext(paste0(src, ", codon pos. ", i) , 3, 0)
        mtext(xl.raas, 1, 1.3)
        mtext(bquote(psi~fraction), 2, 1.3)
        dev.off()
        cors[[src]][[i]] <- pc
        
    }
}

for ( i in 1:3 ) {
    ccol <- paste0("codon",i,"_psi")
    scol <- paste0("codon",i,"_source")
    plotdev(file.path(rfig.path,paste0("psi_raas_codon",i)), type=ftyp,
            height=2.5, width=2.5, res=200)
    par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(psite$RAAS.median, psite[,ccol], density=FALSE,
#            pch=pch.source[psite[,scol]],
            xlab=xl.raas, ylab=bquote("codon pos"~.(i)*","~ psi~fraction),
            col=NA)
    points(psite$RAAS.median, psite[,ccol], pch=pch.source[psite[,scol]],
           cex=.75, col=i)
    text(psite$RAAS.median, psite[,ccol],labels=psite$fromto, xpd=TRUE,
         col=1, cex=.5, pos=4)
    dev.off()
}


## split multiple psi per AAS
psil <- lapply(strsplit(psite$psi, ";"), as.numeric)
xy <- cbind(rep(psite$RAAS.median, lengths(psil)), unlist(psil))

plotdev(file.path(rfig.path,"psi_raas_codons"), type=ftyp,
        height=2.5, width=2.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(xy[,1], xy[,2], density=FALSE,
        xlab=xl.raas, ylab=expression(psi~fraction), pch=NA)
points(psite$RAAS.median, psite$codon1_psi, col=1, pch=19, cex=.5)
points(psite$RAAS.median, psite$codon2_psi, col=2, pch=19, cex=.5)
points(psite$RAAS.median, psite$codon3_psi, col=3, pch=19, cex=.5)
legend("right",paste(1:3), title="codon position", col=1:3, cex=.5,
       pch=19, pt.cex=.75, bty="n")
dev.off()

## HYPERGEO, USING NUMBER OF Us AS BACKGROUND


## COUNT ONLY Us in CODING SEQUENCE
if ( !exists("tfas")  ) {
    tfas <- readFASTA(tfas.file, grepID=TRUE)
    #### protein-transcript map
    ##trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
    #### reverse map transcript-protein
    ##pamrt <- matrix(rownames(trmap), ncol=1)
    ##rownames(pamrt) <- trmap[,1]
    #### rename by protein names via trmap, for quick access
    #### of protein-specific transcript
    ##names(tfas) <- pamrt[names(tfas),1]
}
pids <- unique(c(psite$transcript, PSI$transcript))
sfas <- tfas[pids[pids%in%names(tfas)]]

## total balls: all Us in coding sequences of the total gene set
totaa <- sum(unlist(lapply(sfas,
                           function(x) sum(unlist(strsplit(x$seq,""))=="T"))))

## white balls: all UNIQUE detected psi
m <- length(unique(paste(PSI$chr, PSI$coor)))
n <- totaa-m # black balls

## balls drawn: all U in AAS codons
## take only unique site here!
FILT <- !duplicated(paste(psite$chr,
                          psite$coor, psite$strand))
k <- sum(unlist(strsplit(psite$codon[FILT],""))=="T")

## q: white balls drawn number of AAS with psi
## TODO: account for AAS with multiple hits!!
q <- sum(psite$psi.n[FILT])

##m: the number of white balls in the urn.
##n: the number of black balls in the urn.
##k: the number of balls drawn from the urn

p <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)

## pdist

plotdev(file.path(rfig.path,"psi_hypergeotest"), res=200,
        width=2.5, height=2.5, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
qs <- 1:(2*q)
plot(qs, phyper(q=qs-1, m=m, n=n, k=k, lower.tail=FALSE),
     xlab="overlap count", ylab="probability", type="l", lwd=2)
legend("topright", c(paste("transcripts=", round(length(pids)/1000,1),"k"),
                     paste("total U=", round(totaa/1000000,1),"M"),
                     paste("psi sites=", round(m/1000,1),"k"),
                     paste("U in AAS=", round(k/1000,1),"k"),
                     paste("overlap=", q)), cex=.8,
       box.col=NA, bg=NA)
points(q, p, pch=4, cex=1.2, col=2)
shadowtext(q, p,label=paste0("p=",signif(p,1)), pos=3, col=2, font=2)
##figlabel("union",pos="bottomleft", font=2)
dev.off()
plotdev(file.path(rfig.path,"psi_hypergeotest_log"), res=200,
        width=2.5, height=2.5, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
qs <- 1:700
plot(qs, log10(phyper(q=qs-1, m=m, n=n, k=k, lower.tail=FALSE)),
     xlab="overlap count", ylab=expression(log[10](p)), type="l", lwd=2)
points(q, log10(p), pch=4, cex=1.2, col=2)
shadowtext(q, log10(p), label=paste0("p=",signif(p,1)), pos=4, col=2, font=2)
dev.off()

### order of magnitude test
## 10k proteins of 300 AA length
tot <- 300*10e3*3/4
black <- 80e3
draw <- 7e3*3/4
dblack <- 150

### USING ONLY INTERSECT

## TODO: find out which psi are NOT the intersect version and why
## find those via synonym list.

pids <- unique(psite$transcript[psite$psi.n>0]) 
sfas <- tfas[pids[pids%in%names(tfas)]]

## total balls: all Us in coding sequences of the total gene set
totaa <- sum(unlist(lapply(sfas,
                           function(x) sum(unlist(strsplit(x$seq,""))=="T"))))

## white balls: all detected psi in these transcripts
## TODO: unclear whether we loose sites mapped to a different but
## overlapping transcript
m <- length(unique(paste(PSI$chr, PSI$coor)[PSI$transcript%in%pids]))
n <- totaa-m # black balls

## balls drawn: all U in AAS codons
## take only unique site here!
FILT <- !duplicated(paste(psite$chr,
                          psite$coor, psite$strand)) & psite$transcript%in%pids
k <- sum(unlist(strsplit(psite$codon[FILT],""))=="T")

## q: white balls drawn number of AAS with psi
## TODO: account for AAS with multiple hits!!
q <- sum(psite$psi.n[FILT])

##m: the number of white balls in the urn.
##n: the number of black balls in the urn.
##k: the number of balls drawn from the urn

p <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)

## pdist

plotdev(file.path(rfig.path,"psi_hypergeotest_intersect"), res=200,
        width=2.5, height=2.5, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
qs <- 1:(2*q)
plot(qs, phyper(q=qs-1, m=m, n=n, k=k, lower.tail=FALSE),
     xlab="overlap count", ylab="probability", type="l", lwd=2)
legend("topright", c(paste("transcripts=", length(pids)),
                     paste("total U=", round(totaa/1000,1),"k"),
                     paste("psi sites=", round(m),""),
                     paste("U in AAS=", round(k),""),
                     paste("overlap=", q)), cex=.8,
       box.col=NA, bg=NA)
points(q, p, pch=4, cex=1.2, col=2)
shadowtext(q, p,label=paste0("p=",signif(p,1)), pos=3, col=2, font=2)
##figlabel("intersect",pos="bottomleft", font=2)
dev.off()
plotdev(file.path(rfig.path,"psi_hypergeotest_intersect_log"), res=200,
        width=2.5, height=2.5, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
qs <- 1:700
plot(qs, log10(phyper(q=qs-1, m=m, n=n, k=k, lower.tail=FALSE)),
     xlab="overlap count", ylab=expression(log[10](p)), type="l", lwd=2)
points(q, log10(p), pch=4, cex=1.2, col=2)
shadowtext(q, log10(p), label=paste0("p=",signif(p,1)), pos=4, col=2, font=2)
dev.off()



### m6A MODIFICATIONS

M6A <- NULL

## m6A quantification from rnamd database
rmd <- read.delim(rmd.file)

## m6A site annotation
rma <- read.delim(rma.file)

## TODO: understand data structure,
## log2FC in different cell lines?
##apply(rmd,2, function(x) sum(!is.na(x)))

rmd$chr <- sub("^chr", "", rmd$seqnames)

## take median log2FC over all experiments
## TODO: instead select a single column with most sites?
expCols <- c(grep("^GSE", colnames(rmd)),
             grep("PRJ",  colnames(rmd)))

##apply(rmd[,expCols], 2, function(x) sum(!is.na(x)))

rmd$lg2fc <- apply(rmd[,expCols], 1, median, na.rm=TRUE)
rmd$lg2fc.n <- apply(rmd[,expCols], 1, function(x) sum(!is.na(x)))

## NOTE: distribution of values -> filter?
plotdev(file.path(rfig.path, "m6a_values_per_site"),
        width=2.5, height=2.5, res=200, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
hist(rmd$lg2fc.n, xlab="values per site", main=NA)
dev.off()

## FILTER by number of observations
if ( FALSE ) 
    rmd <- rmd[rmd$lg2fc.n > median(rmd$lg2fc.n),]

## add gene name info
d2a <- match(rmd$ID, rma$ID)
rmd$gene <- rma[d2a,"Gene_Name"]
rmd$annotation <- rma[d2a, "annotation"]

##  REDUCE TO EXONS
rmd <- rmd[rmd$annotation=="Exon",]


## manually use this to test code below
rmdz <- cbind(rmd[,c("chr","start","end","strand")],
              info="m6A",
              gene=rmd$gene, # TODO: load via site info file!
              lg2fc=rmd$lg2fc, 
              source="RNAMD")
M6A <- rbind(M6A, rmdz)


plotdev(file.path(rfig.path, "m6a_boxplots"),
        width=2.5, height=2.5, res=200, type=ftyp)
par(mai=c(.1,.5,.5,.25), mgp=c(1.4,.3,0), tcl=-.25)
bp <- boxplot(M6A$lg2fc ~ sub(".*_","",M6A$source),
              las=2,  xlab=NA, ylab=NA, cex=.5)
mtext(expression(log[2](fold~change)), 2, 1.4)
axis(3, at=1:length(bp$n), labels=bp$n, las=2, cex.axis=.8)
##abline(h=1, col=2, lty=3)
dev.off()



## get ensembl protein ID
M6A$ensembl <- names(ens2nam[match(M6A$gene,  ens2nam)]) 


## take mean position - only affects @Zhang which provides ranges,
## and we found only one overlap (2nd and 3rd position of one AAS)
M6A$coor <- round(apply(M6A[,c("start","end")],1,mean))

## chromosome index
M6A$chr <-  chrIdx[M6A$chr]


### MAP M6A and AAS SITES
m6adx <- coor2index(M6A[,c("chr","coor")], chrS=chrS)[,2]

## collect overlapping sites
psite <- usite
for ( i in 1:3 ) {

    bdidx <- get(paste0("bd",i,"dx"))
                 
    ## find same positions
    bid <- which(bdidx%in%intersect(m6adx, bdidx))
    ##pid <- which(m6adx%in%intersect(m6adx, bdidx))
    cat(paste("found", length(bid), "overlapping sites", i, "\n"))
    
    m6a2site <- match(bdidx, m6adx)
    
    sitem6a <- M6A[m6a2site, ]
    colnames(sitem6a) <- paste0("codon",i, "_", colnames(sitem6a))
    psite <- cbind(psite, sitem6a)
}
## collapse all m6a
## NOTE: using median as a helper

psite$lg2fc.n <- apply(psite[,paste0("codon",1:3,"_lg2fc")], 1, function(x) sum(!is.na(x)))
psite$lg2fc.max <- apply(psite[,paste0("codon",1:3,"_lg2fc")], 1, pmax)
psite$lg2fc <- apply(psite[,paste0("codon",1:3,"_lg2fc")], 1, pvals)
psite$lg2fc.codonpos <- apply(psite[,paste0("codon",1:3,"_lg2fc")], 1,  pcodons)
psite$lg2fc.source <- apply(psite[,paste0("codon",1:3,"_source")], 1, pvals)
psite$lg2fc.chr <- apply(psite[,paste0("codon",1:3,"_chr")], 1, pvals)
psite$lg2fc.coor <- apply(psite[,paste0("codon",1:3,"_coor")], 1, pvals)

## CORRELATION OF M6A % TO RAAS

sources <- unique(psite$lg2fc.source)
sources <- sources[sources!=""]
pch.source <- 1:length(sources) 
names(pch.source) <- sources

for ( i in 1:3 ) {
    ccol <- paste0("codon",i,"_lg2fc")
    scol <- paste0("codon",i,"_source")
    plotdev(file.path(rfig.path,paste0("m6a_raas_codon",i)), type=ftyp,
            height=2.5, width=2.5, res=200)
    par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(psite$RAAS.median, psite[,ccol], density=FALSE,
#            pch=pch.source[psite[,scol]],
            xlab=xl.raas, ylab=bquote("codon pos"~.(i)*","~ log[2](FC)),
            col=NA)
    points(psite$RAAS.median, psite[,ccol], pch=pch.source[psite[,scol]],
           cex=.75, col=i)
    text(psite$RAAS.median, psite[,ccol],labels=psite$fromto, xpd=TRUE,
         col=1, cex=.5, pos=4)
    dev.off()
}


## correlation of RAAS to m6a %
## split multiple psi per AAS
lg2fcl <- lapply(strsplit(psite$lg2fc, ";"), as.numeric)
xy <- cbind(rep(psite$RAAS.median, lengths(lg2fcl)), unlist(lg2fcl))

plotdev(file.path(rfig.path,"m6a_raas_codons"), type=ftyp,
        height=2.5, width=2.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(xy[,1], xy[,2], density=FALSE,
        xlab=xl.raas, ylab=expression(log[2](FC)), pch=NA)
points(psite$RAAS.median, psite$codon1_lg2fc, col=1, pch=19, cex=.5)
points(psite$RAAS.median, psite$codon2_lg2fc, col=2, pch=19, cex=.5)
points(psite$RAAS.median, psite$codon3_lg2fc, col=3, pch=19, cex=.5)
legend("topright",paste(1:3), title="codon position", col=1:3, cex=.5,
       pch=19, pt.cex=.75, bty="n")
dev.off()

## HYPERGEO, USING NUMBER OF Us AS BACKGROUND


## COUNT ONLY Us in CODING SEQUENCE
if ( !exists("tfas")  ) {

    tfas <- readFASTA(tfas.file, grepID=TRUE)
    ## protein-transcript map
    trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
    ## reverse map transcript-protein
    pamrt <- matrix(rownames(trmap), ncol=1)
    rownames(pamrt) <- trmap[,1]
    ## rename by protein names via trmap, for quick access
    ## of protein-specific transcript
    names(tfas) <- pamrt[names(tfas),1]

}

pids <- unique(c(psite$ensembl, M6A$ensembl))
sfas <- tfas[pids[pids%in%names(tfas)]]

## total balls: all Us in coding sequences of the total gene set
totaa <- sum(unlist(lapply(sfas,
                           function(x) sum(unlist(strsplit(x$seq,""))=="A"))))

## white balls: all detected m6a
m <- length(unique(paste(M6A$chr, M6A$coor, M6A$strand)))
n <- totaa-m # black balls

## balls drawn: all U in AAS codons
## take only unique site here!
FILT <- !duplicated(paste(psite$chr,
                          psite$coor, psite$strand))

## investigate duplicated, shouldnt be there?
## TODO: why are their equal sites matched to different proteins??
## ENSP00000365998_74 and ENSP00000366005_74
if ( FALSE ) {
    tmp <- paste(psite$chr,
                 psite$coor, psite$strand)
    idx <- which(tmp==tmp[!FILT][2])
    psite[idx,]
}

k <- sum(unlist(strsplit(psite$codon[FILT],""))=="A")

## q: white balls drawn number of AAS with m6a
q <- sum(psite$lg2fc.n[FILT])

##m: the number of white balls in the urn.
##n: the number of black balls in the urn.
##k: the number of balls drawn from the urn

p <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)

## pdist

plotdev(file.path(rfig.path,"m6a_hypergeotest"), res=200,
        width=2.5, height=2.5, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
qs <- 1:(2*q)
plot(qs, phyper(q=qs-1, m=m, n=n, k=k, lower.tail=FALSE),
     xlab="overlap count", ylab="probability", type="l", lwd=2)
points(q, p, pch=4, cex=1.2, col=2)
shadowtext(q, p, label=paste0("p=",signif(p,1)), pos=3, col=2, font=2)
legend("topright", c(paste("total A=", round(totaa/1000000,1),"M"),
                     paste("psi sites=", round(m/1000,1),"k"),
                     paste("A in AAS=", round(k/1000,1),"k"),
                     paste("overlap=", q)), cex=.8,
       box.col=NA, bg=NA)
dev.off()
plotdev(file.path(rfig.path,"m6a_hypergeotest_log"), res=200,
        width=2.5, height=2.5, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
qs <- 1:900
plot(qs, log10(phyper(q=qs-1, m=m, n=n, k=k, lower.tail=FALSE)),
     xlab="overlap count", ylab=expression(log[10](p)), type="l", lwd=2)
points(q, log10(p), pch=4, cex=1.2, col=2)
shadowtext(q, log10(p), label=paste0("p=",signif(p,1)), pos=4, col=2, font=2)
dev.off()

## export site list
ecols <- c("chr", "coor", "strand", "RAAS.median", "name","ensembl",
           "pos", "codon",
           "fromto", "lg2fc",
           "lg2fc.codonpos", "lg2fc.chr", "lg2fc.coor", "lg2fc.source")
m6asite <- psite[psite$lg2fc.n>0,ecols]
