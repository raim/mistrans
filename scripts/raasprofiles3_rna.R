
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
## TODO: use unique sites?
usite <- asite[!is.na(asite$codon),]


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

## @Zhang2023 - PRAISE to find psi sites
zh23.file <- file.path(mam.path, "originalData", "zhang23_stables_1-6.xlsx")
## @Dai2023 - BID-seq
dai.file <- file.path(mam.path, "originalData", "dai23_stables_1-23.xlsx")

## @Song2020: piano database of psi sites
pia.file <- file.path(mam.path, "originalData", "piano_psi_human.txt")

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
##PSI <- rbind(PSI, psi)
      

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
##PSI <- rbind(PSI, psi)
      

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
    ##PSI <- rbind(PSI, psi)
}


### RNA pseudouridylation data sent by
### Oleksandra Fanari <fanari.o@northeastern.edu>
psi <- NULL
for ( i in excel_sheets(psi.file) ) {
    tmp <-as.data.frame(read_xlsx(psi.file, sheet=i))
    psi <- rbind(psi, cbind(tmp, cell=i))
}
##psi <- read.csv(psi.file)
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
              strand=".",
              info="",
              gene=psi$Annotation,
              psi=psi[,"mm.DirectMINUSmm.IVT"]/100, # % -> fraction
              source=paste0("fanari24pre","_",psi$cell))

## TODO: track different annotations!

## QC: duplicated sites annotated for different genes!
## also see below
tmp <- paste(psiz$chr,psiz$start,psiz$source)
psiz[which(tmp==tmp[which(duplicated(tmp))[2]]),]

##  e.g. psia["NA.46990",] annotated as RP1-120G22.12
## but this has coors chr1 6264900-6265840 while psi site
## is annotated as chr1 6197653

## line 26694: RP1-120G22.12	chr1	6197653	TTTTG	4013	13	2249	6	4062	2282	0.266075388026608	0.322901142573274	0.266075388026608	0	1	0.056825754546666	4026	2255	1	RP1-120G22.12chr16197653	1	pIVT	chr16197653

## line 27997: RPL22	chr1	6197653	TTTTG	4013	13	2249	6	4062	2282	0.266075388026608	0.322901142573274	0.266075388026608	0	1	0.056825754546666	4026	2255	1	RPL22chr16197653	1	pIVT	chr16197653

## FILTER
##psiz <- psiz[psiz$psi>.1,]

PSI <- rbind(PSI, psiz)

## @Song2020 - PIANO
## NOTE: unclear to which genome release the current version refers to,
## and all loglikelihood ratios=100
pia <- read.delim(pia.file)

pia$chr <- sub("^chr","", pia$Chromosome)
pia$start <- pia$end <- pia[,"Ψ_Site_POS"]

## only take CDS
pia <- pia[pia$Gene_Region=="CDS",]

psiz <- cbind(pia[,c("chr","start","end","strand")],
              info=pia[,"Ψ_Source"],
              gene=pia$Gene, ## TODO: use ensembl gene ID?
              psi=pia[,"likelihood_ratio"]/100, # % -> fraction
              source="piano")
##PSI <- rbind(PSI, psiz)



### PROCESS PSI DATA

## get ensembl protein ID
PSI$ensembl <- names(ens2nam[match(PSI$gene,  ens2nam)]) 


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
    aasite <- psite[aas2psi,c("ensembl","fromto","codon","RAAS.median")]
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

psite$psi.n <- apply(psite[,paste0("codon",1:3,"_psi")], 1,
                     function(x) sum(!is.na(x)))
psite$psi.max <- apply(psite[,paste0("codon",1:3,"_psi")], 1, pmax)
psite$psi <- apply(psite[,paste0("codon",1:3,"_psi")], 1, pvals)
psite$psi.codonpos <- apply(psite[,paste0("codon",1:3,"_psi")], 1,  pcodons)
psite$psi.source <- apply(psite[,paste0("codon",1:3,"_source")], 1, pvals)
psite$psi.chr <- apply(psite[,paste0("codon",1:3,"_chr")], 1, pvals)
psite$psi.coor <- apply(psite[,paste0("codon",1:3,"_coor")], 1, pvals)

## TODO: use psia to track redundant annotations, see above

psia$aas.n <- apply(psia[,paste0("codon",1:3,"_RAAS.median")], 1,
                    function(x) sum(!is.na(x)))
psia$aas.gene <- apply(psia[,paste0("codon",1:3,"_ensembl")], 1, pvals)


## export site for Sara and Sasha
## TODO: rm separate codon columns, and instead add affected
## codon position:
ecols <- c("chr", "coor", "strand", "RAAS.median", "name","ensembl",
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

sources <-
    unique(na.omit(unlist(psite[,paste0("codon",1:3,"_source")])))
pch.source <- 1:length(sources) 
names(pch.source) <- sources

## for each source/cell line separately
for ( src in sources ) {
    ssite <- psite[grep(src,psite$psi.source),]
    src <- sub(".*_","",src)
    for ( i in 1:3 ) {
        ccol <- paste0("codon",i,"_psi")
        scol <- paste0("codon",i,"_source")
        nvals <- sum(!is.na(ssite[,ccol]))
        plotdev(file.path(rfig.path,paste0("psi_raas_codon",i,"_",src)),
                type=ftyp, height=2.5, width=2.5, res=200)
        par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
        if (nvals==0) {
            plot(1, col=NA, xlab=NA, ylab=NA, axes=FALSE)
            box()
        } else if ( nvals>2) {
            plotCor(ssite$RAAS.median, ssite[,ccol], density=FALSE,
                    xlab=NA, ylab=NA, signif=2, col=NA)
        } else {
            plot(ssite$RAAS.median, ssite[,ccol], col=NA,
                 xlab=NA, ylab=NA)
        }
        if ( nvals>0 ) {
            points(ssite$RAAS.median, ssite[,ccol], pch=pch.source[ssite[,scol]],
                   cex=.75, col=i)
            text(ssite$RAAS.median, ssite[,ccol],labels=ssite$fromto, xpd=TRUE,
                 col=1, cex=.5, pos=4)
        }
        mtext(src, 3, 0)
        mtext(xl.raas, 1, 1.3)
        mtext(bquote("codon pos"~.(i)*","~ psi~fraction), 2, 1.3)
        dev.off()
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


## correlation of RAAS to psi %



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
    ## protein-transcript map
    trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
    ## reverse map transcript-protein
    pamrt <- matrix(rownames(trmap), ncol=1)
    rownames(pamrt) <- trmap[,1]
    ## rename by protein names via trmap, for quick access
    ## of protein-specific transcript
    names(tfas) <- pamrt[names(tfas),1]
}
pids <- unique(c(psite$ensembl, PSI$ensembl))
sfas <- tfas[pids[pids%in%names(tfas)]]

## total balls: all Us in coding sequences of the total gene set
totaa <- sum(unlist(lapply(sfas,
                           function(x) sum(unlist(strsplit(x$seq,""))=="T"))))

## white balls: all detected psi
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
     xlab="overlap count", ylab=expression(p), type="l", lwd=2)
legend("topright", c(paste("total U=", round(totaa/1000000,1),"M"),
                     paste("psi sites=", round(m/1000,1),"k"),
                     paste("U in AAS=", round(k/1000,1),"k"),
                     paste("overlap=", q)), cex=.8,
       box.col=NA, bg=NA)
points(q, p, pch=4, cex=1.2, col=2)
shadowtext(q, p,label=paste0("p=",signif(p,1)), pos=3, col=2, font=2)
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
pids <- unique(intersect(psite$ensembl, PSI$ensembl))
sfas <- tfas[pids[pids%in%names(tfas)]]

## total balls: all Us in coding sequences of the total gene set
totaa <- sum(unlist(lapply(sfas,
                           function(x) sum(unlist(strsplit(x$seq,""))=="T"))))

## white balls: all detected psi
m <- length(unique(paste(PSI$chr, PSI$coor)[PSI$ensembl%in%pids]))
n <- totaa-m # black balls

## balls drawn: all U in AAS codons
## take only unique site here!
FILT <- !duplicated(paste(psite$chr,
                          psite$coor, psite$strand)) & psite$ensembl%in%pids
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
     xlab="overlap count", ylab=expression(p), type="l", lwd=2)
legend("topright", c(paste("total U=", round(totaa/1000000,1),"M"),
                     paste("psi sites=", round(m/1000,1),"k"),
                     paste("U in AAS=", round(k/1000,1),"k"),
                     paste("overlap=", q)), cex=.8,
       box.col=NA, bg=NA)
points(q, p, pch=4, cex=1.2, col=2)
shadowtext(q, p,label=paste0("p=",signif(p,1)), pos=3, col=2, font=2)
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
     xlab="overlap count", ylab=expression(p), type="l", lwd=2)
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
