
## load AAS mapping, 
## load protein fasta,
## load transcript fasta,
## get sequences surrounding the AAS and analyze here,
## and export for sequence motif analysis.

## TODO:
## * filter unique proteins?
## * align at internal K,
## * grep

## TODO:
## * still contains some genes that may be filtered, eg. IGLL5,
## *

## TODO:

## * PCA of AA FREQUENCIES in -7:7

source("~/work/mistrans/scripts/saap_utils.R")

histn <- function(...) hist(main=NA, ...)

library(segmenTools)
require(ggplot2)
library(readxl)
require(ggseqlogo)
library(Peptides)
options(stringsAsFactors=FALSE)

library(Biostrings)
AAS <- sort(unique(GENETIC_CODE))
AAT <- AAS[AAS!="*"]
diAAT <- sort(paste(AAT, rep(AAT,each=length(AAT)),sep=""))

##AAT <- c(AAT,"-")

## shapely
aa.cols <- character()
aa.cols[c("D","E")]   <- "#E60A0A"
aa.cols[c("C","M")]   <- "#E6E600"
aa.cols[c("K","R")]   <- "#145AFF"
aa.cols[c("S","T")]   <- "#FA9600"
aa.cols[c("F","Y")]   <- "#3232AA"
aa.cols[c("N","Q")]   <- "#00DCDC"
aa.cols[c("G")]     <- "#EBEBEB" #light grey
aa.cols[c("L","V","I")] <- "#0F820F"
aa.cols[c("A")]     <- "#C8C8C8"
aa.cols[c("W")]     <- "#B45AB4"
aa.cols[c("H")]     <- "#8282D2"
aa.cols[c("P")]     <- "#DC9682"
aa.pchs <- numeric()
aa.pchs[c("D","E")]   <- c(19,4)
aa.pchs[c("C","M")]   <- c(4,19)
aa.pchs[c("K","R")]   <- c(19,4)
aa.pchs[c("S","T")]   <- c(19,4)
aa.pchs[c("F","Y")]   <- c(19,4)
aa.pchs[c("N","Q")]   <- c(19,4)
aa.pchs[c("G")]     <- 19 #light grey
aa.pchs[c("L","V","I")] <- c(19,2,4)
aa.pchs[c("A")]     <- 19
aa.pchs[c("W")]     <- 4
aa.pchs[c("H")]     <- 19
aa.pchs[c("P")]     <- 4

## randomize protein names in FASTA to
## test randomly picked sequences
RANDOMIZE <- FALSE # TRUE #  

## DISTANCE AROUND AAS TO ANALYZE
dst <- 25  # left/right distance for frequency plots
sdst <- 5 # left/right distance for fasta file
mdst <- 3# left/right distance for PTM search with momo

## do diAA enrichment for all filters
do.dimer <- TRUE #FALSE
## PCA analysis of AA frequencies
do.PCA <- FALSE
## peptide property analysis
do.props <- FALSE

## threshold p-values
p.min <- 1e-10
p.txt <- 1e-5

## filtering with bonferoni correction
p.sig <- .05
l.sig <- 3
p.flt <- p.sig/((l.sig*2+1)*20) # bonferoni-adjusted p.value filter

## p-value colors
docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))

### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"processedData")
out.path <- file.path(proj.path,"processedData","motifs")
fig.path <- file.path(proj.path,"figures", "saap_context")

if ( RANDOMIZE ) {
    fig.path <- file.path(fig.path, "random")
    out.path <- file.path(out.path, "random")
}
    
dir.create(out.path, showWarnings=FALSE)
dir.create(fig.path, showWarnings=FALSE)

## INPUT FILES
## list of SAAPs with coordinates
in.file <- file.path(dat.path,"saap_mapped3.tsv")
## protein fasta
fasta <- file.path(dat.path,"all_proteins.fa")
## coding region fasta
##transcr <- file.path(mam.path,"processedData","coding.fa")

degron.file <- file.path(proj.path,"originalData",
                         "DEGRONOPEDIA_degron_dataset.xlsx")

## analyze degron patterns
degrons <- as.data.frame(read_xlsx(degron.file, sheet=2))
grep("W",degrons$Degron_regex, value=TRUE)
degl <- strsplit(degrons$Degron,"")
dlen <- unlist(lapply(degl, function(x) sum(x%in%AAT)))
hist(dlen)

## filter nonAA chars
degl <- lapply(degl, function(x) x[x%in%AAT])

## remove short and (partially regex
degl <- degl[dlen>5]

## GET SAAP/BP DATA
dat <- read.delim(in.file)
cat(paste("loaded", nrow(dat), "SAAP/BP\n"))

## mahlon's degron candidates
bpids <- c("LLIYTNNQRPSGVPDR", "EIASTLMESEMMEILSVLAK", "MDQPAGLQVDYVFR")
dat[which(dat$BP%in%bpids),"name"]

## remove by exclude tag - IG/Albumin/globin
cat(paste("removing", sum(dat$IG,na.rm=TRUE), "immunoglobulins\n"))
dat <- dat[!dat$IG ,]

## remove those w/o protein info
cat(paste("removing", sum(is.na(dat$pos)), "w/o protein info\n"))
dat <- dat[!is.na(dat$pos), ]

## remove protein with substitutions at the same site
## TODO: more complex handling of this
dupls <- duplicated(paste(dat$protein, dat$pos))
cat(paste("removing", sum(dupls),
          "AAS with equal position in proteins\n"))
dat <- dat[!dupls,]

## rm AAS w/o from info
## TODO: where do these come from?
cat(paste("removing", sum(dat$from==""), "w/o from AA info\n"))
dat <- dat[dat$from!="",]

## ONLY USE GOOD BLAST MATCHES
cat(paste("removing", sum(dat$match!="good"), "with bad blast match\n"))
dat <-dat[dat$match=="good",]

## add new tags
## tag trypsin 
dat$KR <- dat$from%in%c("K","R")

## tag miscleavage
ncut <- 2
ccut <- 2
bps <- substr(dat$BP,1+ncut,nchar(dat$BP)-ccut)
dat$miscleavage <- rep(FALSE, nrow(dat))
dat$miscleavage[grep("[KR]", bps)] <- TRUE

## FILTER ONLY UNIQUE BP
dupls <- duplicated(dat$BP)
cat(paste("removing", sum(dupls), "duplicated BP\n"))
bpd <- dat[!dupls,]


### ANALYZE MUTATION WITHIN PEPTIDE

## AA vs. position bins
site.bins <- cut(bpd$site/nchar(bpd$SAAP), seq(0,1,.2))
cls <- clusterCluster(bpd$from, site.bins, cl2.srt=levels(site.bins),
                      alternative="two.sided")
##cls <- sortOverlaps(cls, p.min=.1)
plotdev(file.path(fig.path,paste0("peptide_AA_overlap")),
        height=5, width=3, res=300)
par(mai=c(1,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(cls, p.min=1e-10, p.txt=1e-5, ylab="encoded AA at AAS", xlab=NA,
             show.total=TRUE, show.sig=FALSE)#, col=ttcols)
mtext("relative position of AAS in peptide", 1, 3.5)
dev.off()


## Categorize by N+i and C-i
mx <- 9
nt <- bpd$site
ct <- nchar(bpd$BP)-bpd$site +1
## fuse central AAS
nt[nt>mx] <- ct[ct>mx] <- mx+1
## take smaller
at <- paste0("N",nt)
at[ct<nt] <- paste0("C",ct[ct<nt])
## re-name central AAS
at[at==paste0("N",mx+1)] <- paste0(">",mx)
at.srt <- c(
    paste0("N",1:mx),
    paste0(">",mx),
    paste0("C",mx:1))

asite.bins <- at
asite.srt <- at.srt
cls <- clusterCluster(bpd$from, asite.bins, cl2.srt=asite.srt,
                      alternative="two.sided")
##cls <- sortOverlaps(cls, p.min=.1)
plotdev(file.path(fig.path,paste0("peptide_AA_overlap_absolute")),
        height=5, width=mx/2+1, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=1e-10, p.txt=1e-5, ylab="encoded AA at AAS", xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()


## TODO: use p.value for color and ratio for size
plotdev(file.path(fig.path,
                  paste0("peptide_AA_overlap_absolute_enrichment")),
        height=5, width=mx/2+1, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
plotOverlaps2(cls, max.val=5, xlab=NA, ylab="encoded AA at AAS",
              p.min=1e-10, dot.sze=c(.2,2), bg=TRUE)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

source("~/work/mistrans/scripts/saap_utils.R")
plotdev(file.path(fig.path,paste0("peptide_AA_overlap_absolute_dotplot")),
        height=5, width=mx/2+1, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
dotprofile(cls, value="ratio", vcols=ttcols, xlab=NA,
           p.dot=1e-10, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2),test=FALSE,
           ylab="encoded AA at AAS", show.total=TRUE)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()


## STORE to compare with miscleavage pattern below
OVL.AA <- cls

## TODO: plot frequency ratios at site==1 vs. miscleavage K|R +1
## negative correlation!


##bpd <- bpd[bpd$from!="Q",] # no difference
plotdev(file.path(fig.path,paste0("peptide_AAS")),
        height=3, width=3, res=300)
par(mfcol=c(1,1), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
hist(bpd$site/nchar(bpd$SAAP), breaks=seq(0,1,.1),
     xlab="relative position of AAS in peptide", main=NA)
legend("topright", paste(nrow(bpd), "unique BP"))
dev.off()

plotdev(file.path(fig.path,paste0("peptide_AAS_absolute")),
        height=3, width=mx/2+1, res=300)
par(mfcol=c(1,1), mai=c(.5,.5,.1,.1), mgp=c(1.5,.05,0), tcl=.25, yaxs="i")
barplot(table(asite.bins)[asite.srt], xlab="position of AAS in peptide",las=2)
legend("topright", paste(nrow(bpd), "unique BP"))
dev.off()

##bpd <- bpd[bpd$from!="Q",] # no difference
plotdev(file.path(fig.path,paste0("peptide_AAS_n2")),
        height=3, width=6, res=300)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
barplot(table(bpd$site[bpd$site>2]), xlab="position of AAS in peptide",las=2)
hist(bpd$site[bpd$site>2]/nchar(bpd$SAAP[bpd$site>2]), breaks=seq(0,1,.1),
     xlab="relative position of AAS in peptide", main=NA)
legend("topright", paste(sum(bpd$site>2,na.rm=TRUE), "unique BP"))
dev.off()

hist(nchar(bpd$SAAP), xlab="peptide length", main=NA, breaks=seq(1,60,1))

## only for unique BP - w/o from Q
bpnq <- bpd[bpd$from!="Q",] # no difference
plotdev(file.path(fig.path,paste0("peptide_AAS_noQ")),
        height=3, width=6, res=300)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
barplot(table(bpnq$site), xlab="position of AAS in peptide",las=2)
hist(bpnq$site/nchar(bpnq$SAAP), breaks=seq(0,1,.1),
     xlab="relative position of AAS in peptide", main=NA)
legend("topright", paste(nrow(bpnq), "w/o 'from:Q'"))
dev.off()


## POSITION OF K/R IN PEPTIDE - MIS-CLEAVAGE
ncut <- 0
ccut <- 0
bps <- substr(bpd$BP,1+ncut,nchar(bpd$BP)-ccut)
des <- gregexpr("[KR]", bps) #bpd$BP)
des <- lapply(des, function(x) as.numeric(x[x>0]+ncut))
der <- des
for ( i in 1:length(bps) ) {
    der[[i]] <- (der[[i]]+ncut)/nchar(bps[i])
}
des <- unlist(des)
der <- unlist(der)
plotdev(file.path(fig.path,paste0("peptide_KR")),
        height=3, width=6, res=300)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")

tb1 <- table(bpd$qlen)
tb2 <- table(des[des>0])
tbn <- as.character(sort(as.numeric(unique(c(names(tb1), names(tb2))))))
tb <- rbind(tb2[tbn], tb1[tbn])
barplot(tb, beside=FALSE,
        xlab="position of K|R in peptide",las=2,
        col=c("#ff000099","#77777777"))
legend("topright", paste(length(des), "K|R"))
hist(der[der>0], #breaks=seq(0,1,.1),
     xlab="relative position of K|R in peptide", main=NA)
dev.off()

## POSITION OF K/R IN PEPTIDE - CUT ENDS
ncut <- 2
ccut <- 2
bps <- substr(bpd$BP,1+ncut,nchar(bpd$BP)-ccut)
des <- gregexpr("[KR]", bps) #bpd$BP)
des <- lapply(des, function(x) as.numeric(x[x>0]+ncut))
der <- des
for ( i in 1:length(bps) ) {
    der[[i]] <- (der[[i]]+ncut)/nchar(bps[i])
}
des <- unlist(des)
der <- unlist(der)
plotdev(file.path(fig.path,paste0("peptide_KR_cut")),
        height=3, width=6, res=300)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")

tb1 <- table(bpd$qlen)
tb2 <- table(des[des>0])
tbn <- as.character(sort(as.numeric(unique(c(names(tb1), names(tb2))))))
tb <- rbind(tb2[tbn], tb1[tbn])
barplot(tb, beside=FALSE,
        xlab="position of K|R in peptide",las=2,
        col=c("#ff000099","#77777777"))
legend("topright", paste(length(des), "K|R"))
hist(der[der>0], #breaks=seq(0,1,.1),
     xlab="relative position of K|R in peptide", main=NA)
dev.off()

## POSITION OF D/E IN PEPTIDE
des <- gregexpr("[DE]", bpd$BP)
der <- des
for ( i in 1:nrow(bpd) ) {
    der[[i]] <- der[[i]]/nchar(bpd$BP[i])
}
des <- unlist(des)
der <- unlist(der)
plotdev(file.path(fig.path,paste0("peptide_DE")),
        height=3, width=6, res=300)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
tb1 <- table(bpd$qlen)
tb2 <- table(des[des>0])
tbn <- as.character(sort(as.numeric(unique(c(names(tb1), names(tb2))))))
tb <- rbind(tb2[tbn], tb1[tbn])
barplot(tb, beside=FALSE,
        xlab="position of D|E in peptide",las=2,
        col=c("#ff000099","#77777777"))
##barplot(table(des[des>0]), xlab="position of D|E in peptide",las=2)
legend("topright", paste(length(des), "D|E"))
hist(der[der>0], breaks=seq(0,1,.1),
     xlab="relative position of D|E in peptide", main=NA)
dev.off()


### MISCLEAVAGE
## remove all where K|R are the mutation
bpk <- bpd[!bpd$KR,]
## cut FIRST and last AA since this is mostly K|R cleavage site
## NOTE: why is above diagonal missing by cutting first?
ncut <- 2
ccut <- 2
bps <- substr(bpk$BP, 1+ncut, nchar(bpk$BP)-ccut)
krs <- gregexpr("[KR]", bps)
krs <- lapply(krs, function(x) as.numeric(x[x>0]+ncut))
## get mean (first) KR site
krf <- unlist(lapply(krs, function(x) mean(x)))
plotdev(file.path(fig.path,paste0("miscleavage_KR_vs_AAS")),
        height=2, width=6, res=300)
par(mfcol=c(1,3), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(bpk$site[krf>0], krf[krf>0],  ##colf=arno,
        ylab="mean position of K|R in peptide",
        xlab="position of AAS in peptide", xlim=c(0,40), ylim=c(0,40), cex=.5)
abline(a=0, b=1)
hist((krf[krf>0]- bpk$site[krf>0]), main=NA,
     xlab="abs. dist. (K|R - AAS)", xlim=c(-10,10),
     breaks=-100:100-.5)
legend("topleft", paste(sum(krf>1,na.rm=TRUE),"w K|R"), bty="n")
hist((krf[krf>0]- bpk$site[krf>0])/nchar(bpk$BP)[krf>0], main=NA,
     xlab="(K|R - AAS)/peptide length")
dev.off()

## distance to D|E
## all where K|R are the mutation
bpk <- bpd[!bpd$KR,]
## don't cut last two AA 
bps <- bpk$BP #substr(bpk$BP,1,nchar(bpk$BP)-2)
head(grep("[DE]", bps, value=TRUE))
krs <- gregexpr("[DE]", bps)
krs <- lapply(krs, function(x) as.numeric(x))
## get mean (first)
krf <- unlist(lapply(krs, function(x) mean(x)))
plotdev(file.path(fig.path,paste0("miscleavage_DE_vs_AAS")),
        height=2, width=6, res=300)
par(mfcol=c(1,3), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(bpk$site[krf>0], krf[krf>0], 
        ylab="mean position of D|E in peptide",
        xlab="position of AAS in peptide", xlim=c(0,40), ylim=c(0,40), cex=.5)
abline(a=0, b=1)
hist((krf[krf>0]- bpk$site[krf>0]), main=NA,
     xlab="abs. dist. (D|E - AAS)", xlim=c(-30,30),
     breaks=-100:100)
legend("topright", paste(sum(krf>1),"w D|E"), bty="n")
hist((krf[krf>0]- bpk$site[krf>0])/nchar(bpk$BP)[krf>0], main=NA,
     xlab="(D|E - AAS)/peptide length")
dev.off()

### ANALYZE SEQUENCE CONTEXT OF MUTATION
## load fasta files, analyze frequencies, and write out
## fasta for sequence analysis.
pfas <- readFASTA(fasta, grepID=TRUE)

## inspect some known degrons, eg.
## @Poirson2024: PRR20A, ENSP00000367164, residues 168–221
pid <- "ENSP00000367164"
unlist(strsplit(pfas[[pid]]$seq,""))[168:221]
## @Zhang2015: EID1, ENSP00000431162, residues 160–172
pid <- "ENSP00000431162"
unlist(strsplit(pfas[[pid]]$seq,""))[160:172]
## "A" "F" "I" "E" "E" "L" "F" "S" "L" "M" "V" "V" "N"


## get only the required protein sequences
## in the correct order!

pidx <- match(bpd$protein, names(pfas))
### RANDOMIZE PROTEIN NAMES as control
if ( RANDOMIZE ) {
    pidx <- sample(1:length(pfas))[1:nrow(bpd)]
    ## random positions WITHIN the randomly selected protein
    bpd$pos <- unlist(lapply(pfas[pidx], function(x)
        sample(1:nchar(x$seq))[1]))
    
}

pfas <- pfas[pidx]


pctx <- sapply(1:nrow(bpd), function(i) {
    ##for ( i in 1:nrow(dat) ){
    fsq <- pfas[[i]]$seq # protein sequence
    nsq <- nchar(fsq)
    
    rrng <- seq(-dst,dst,1)
    sq <- rep("-",length(rrng)) ## GAP
    names(sq) <- rrng

    ## GET RANGE AROUND AAS
    rng <- (bpd$pos[i]-dst):(bpd$pos[i]+dst)
    ## cut range to available
    rrng <- rrng[rng>0 & rng<=nsq]
    rng <- rng[rng>0 & rng<=nsq]
    sq[as.character(rrng)] <- unlist(strsplit(fsq,""))[rng]
    list(sq)
    ##}
})
## as matrix
pctx <- do.call(rbind, pctx)

## gene names
rownames(pctx) <-
    paste0(bpd$name,"_",bpd$pos-dst,"_",bpd$pos+dst)

sum(pctx=="-")

## convert AA matrix to diAA matrix
## TODO: grep specifically the found diAA patterns
## M[S|T], W[ILV], AQ, KA, RA
dctx <- character()
for ( i in 1:(ncol(pctx)-1) )
    dctx <- cbind(dctx, apply(pctx[,i:(i+1)],1,paste,collapse=""))
colnames(dctx) <- head(colnames(pctx), ncol(pctx)-1)
rownames(dctx) <- rownames(pctx)

if ( interactive() ) {
    x <-c(table(c(pctx[,20:25]))[AAT])
    y <-c(table(c(pctx[,26:31]))[AAT])
    plot(x, y, col=NA)
    abline(a=0, b=1)
    text(x, y, names(x),xpd=TRUE)
}

### PCA OF FREQUENCIES
## TODO: do by substitution type
if ( do.PCA ) {

pcalab <- expression(bold(AAS%+-%10))

pcax <- pctx[,as.character(-10:10)]
pcax <- apply(pcax,1,function(x) table(x)[AAT])
rownames(pcax) <- AAT

## frequencies
pcax <- apply(pcax,2, function(x) x/sum(x,na.rm=TRUE))
pcax[is.na(pcax)] <- 0

## cluster?
if ( FALSE ) {
    kcl <- kmeans(t(pcax), 5)
    boxplot(pcax["A",] ~ kcl$cluster)
    boxplot(pcax["G",] ~ kcl$cluster)
    boxplot(pcax["P",] ~ kcl$cluster)
    plot(pcax["G",], pcax["S",])
}

## enrichments instead of frequencies?
if ( FALSE ) {
    
    ## get frequencies in background
    pcab <- pctx[,as.character(c(-25:-11),c(11:25))]
    pcab <- apply(pcab,1,function(x) table(x)[AAT])
    rownames(pcab) <- AAT
    
    ##pcab <- apply(pcab,2, function(x) x/sum(x,na.rm=TRUE))
    ##pcab[is.na(pcab)] <- 0
    
    
    ## make sure columns and rows are the same
    ##pcab <- pcab[rownames(pcax),]
    ##pcab <- pcab[,colnames(pcax)]
    
    ## background over all
    bg <- apply(pcab,1, sum, na.rm=TRUE)/sum(pcab, na.rm=TRUE)

    ## ENRICHMENT OVER LOCAL BACKGROUND
    pcax <- pcax/bg
}

pcax <- t(pcax)

## TODO: why positive correlations for randomized and
## negative for AAS and degrons??
tmp <- cor(pcax, use="pairwise.complete")
tmp[upper.tri(tmp)] <- NA
plotdev(file.path(fig.path,paste0("PCA_AA_correlation")),
        type="png", res=300, width=4,height=4)
par(mai=c(.25,.25,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
image_matrix(tmp, col=ttcols, breaks=seq(-.5,.5,length=length(ttcols)+1),
             axis=1:2, xlab=NA, ylab=NA)
dev.off()

pcax <- pcax - apply(pcax,1,mean) # center rows (peptides)
pcx <- prcomp(pcax, scale=TRUE)

values <- pcx$sdev^2
vectors <- pcx$rotation
var <- values/sum(values) # % of variance explained by PCn

## add % variance explained to colnames, for axis labels
colnames(pcx$x) <-
    paste0(colnames(pcx$x), " (",round(100*var), "%)")


plotdev(file.path(fig.path,paste0("PCA_AA")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcx$rotation[,1],pcx$rotation[,2],
     col=NA,#aa.cols[rownames(pcx$rotation)],
     pch=aa.pchs[rownames(pcx$rotation)],
     xlab=colnames(pcx$x)[1],
     ylab=colnames(pcx$x)[2])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcx$rotation[,1],pcx$rotation[,2],labels=rownames(pcx$rotation),
     col=aa.cols[rownames(pcx$rotation)], font=2, cex=1.2)
figlabel(pcalab, pos="bottomright", font=2, cex=1.2)
dev.off()

plotdev(file.path(fig.path,paste0("PCA_AA_13")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcx$rotation[,1],pcx$rotation[,3],
     col=NA,#aa.cols[rownames(pcx$rotation)],
     pch=aa.pchs[rownames(pcx$rotation)],
     xlab=colnames(pcx$x)[1],
     ylab=colnames(pcx$x)[3])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcx$rotation[,1],pcx$rotation[,3],labels=rownames(pcx$rotation),
     col=aa.cols[rownames(pcx$rotation)], font=2, cex=1.2)
figlabel(pcalab, pos="bottomright", font=2, cex=1.2)
dev.off()
plotdev(file.path(fig.path,paste0("PCA_AA_24")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcx$rotation[,2],pcx$rotation[,4],
     col=NA,#aa.cols[rownames(pcx$rotation)],
     pch=aa.pchs[rownames(pcx$rotation)],
     xlab=colnames(pcx$x)[2],
     ylab=colnames(pcx$x)[4])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcx$rotation[,2],pcx$rotation[,4],labels=rownames(pcx$rotation),
     col=aa.cols[rownames(pcx$rotation)], font=2, cex=1.2)
figlabel(pcalab, pos="bottomright", font=2, cex=1.2)
dev.off()
plotdev(file.path(fig.path,paste0("PCA_AA_23")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcx$rotation[,2],pcx$rotation[,3],
     col=NA,#aa.cols[rownames(pcx$rotation)],
     pch=aa.pchs[rownames(pcx$rotation)],
     xlab=colnames(pcx$x)[2],
     ylab=colnames(pcx$x)[3])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcx$rotation[,2],pcx$rotation[,3],labels=rownames(pcx$rotation),
     col=aa.cols[rownames(pcx$rotation)], font=2, cex=1.2)
figlabel(pcalab, pos="bottomright", font=2, cex=1.2)
dev.off()

## CONTROL

pcclab <- expression(bold(AAS+15))

pcac <- pctx[,as.character(15:25)]
pcac <- apply(pcac,1,function(x) table(x)[AAT])
rownames(pcac) <- AAT

## frequencies
pcac <- apply(pcac,2, function(x) x/sum(x,na.rm=TRUE))
pcac[is.na(pcac)] <- 0


pcac <- t(pcac)

tmp <- cor(pcac, use="pairwise.complete")
tmp[upper.tri(tmp)] <- NA
plotdev(file.path(fig.path,paste0("PCA_AA_correlation_control")),
        type="png", res=300, width=4,height=4)
par(mai=c(.25,.25,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
image_matrix(tmp, col=ttcols, breaks=seq(-.5,.5,length=length(ttcols)+1),
             axis=1:2, xlab=NA, ylab=NA)
dev.off()

pcac <- pcac - apply(pcac,1,mean) # center rows (peptides)
pcc <- prcomp(pcac, scale=TRUE)

values <- pcc$sdev^2
vectors <- pcc$rotation
var <- values/sum(values) # % of variance explained by PCn

## add % variance explained to colnames, for axis labels
colnames(pcc$x) <-
    paste0(colnames(pcc$x), " (",round(100*var), "%)")


plotdev(file.path(fig.path,paste0("PCA_AA_control")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcc$rotation[,1],pcc$rotation[,2],
     col=NA,#aa.cols[rownames(pcc$rotation)],
     pch=aa.pchs[rownames(pcc$rotation)],
     xlab=colnames(pcc$x)[1],
     ylab=colnames(pcc$x)[2])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcc$rotation[,1],pcc$rotation[,2],labels=rownames(pcc$rotation),
     col=aa.cols[rownames(pcc$rotation)], font=2, cex=1.2)
figlabel(pcclab, pos="bottomright", font=2, cex=1.2)
dev.off()


## PCA of diAA content

pcdx <- dctx[,as.character(-7:7)]
pcdx <- apply(pcdx,1,function(x) table(x)[diAAT])
rownames(pcdx) <- diAAT

pcdx <- apply(pcdx,2, function(x) x/sum(x,na.rm=TRUE))
pcdx[is.na(pcdx)] <- 0

pcdx <- t(pcdx)

image_matrix(cor(pcdx), col=ttcols, breaks=seq(-1,1,length=length(ttcols)+1),
             axis=1:2)

pcdx <- pcdx - apply(pcdx,1,mean) # center genes
pcdx <- prcomp(pcdx, scale=TRUE)

values <- pcdx$sdev^2
vectors <- pcdx$rotation
var <- values/sum(values) # % of variance explained by PCn

## add % variance explained to colnames, for axis labels
colnames(pcdx$x) <-
    paste0(colnames(pcdx$x), " (",round(100*var), "%)")

dense2d(pcdx$rotation[,1],pcdx$rotation[,2],
             xlab=colnames(pcdx$x)[1],
             ylab=colnames(pcdx$x)[2])
plotdev(file.path(fig.path,paste0("PCA_diAA")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcdx$rotation[,1],pcdx$rotation[,2],
     col=NA,#aa.cols[rownames(pcdx$rotation)],
     pch=aa.pchs[rownames(pcdx$rotation)],
     xlab=colnames(pcdx$x)[1],
     ylab=colnames(pcdx$x)[2])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
cols <- unlist(lapply(strsplit(rownames(pcdx$rotation),""),function(x) x[1]))
shadowtext(pcdx$rotation[,1],pcdx$rotation[,2],labels=rownames(pcdx$rotation),
     col=aa.cols[cols],cex=.7)
figlabel(pcalab, pos="bottomright", font=2, cex=1.2)
dev.off()

## PCA OF DEGRONS
pcad <- lapply(degl,function(x) table(x)[AAT])
pcad <- do.call(cbind, pcad)
rownames(pcad) <- AAT

pcad <- apply(pcad,2, function(x) x/sum(x,na.rm=TRUE))
pcad[is.na(pcad)] <- 0


pcad <- t(pcad)

tmp <- cor(pcad, use="pairwise.complete")
tmp[upper.tri(tmp)] <- NA
plotdev(file.path(fig.path,paste0("PCA_AA_correlation_degrons")),
        type="png", res=300, width=4,height=4)
par(mai=c(.25,.25,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
image_matrix(tmp, col=ttcols, breaks=seq(-.5,.5,length=length(ttcols)+1),
             axis=1:2, xlab=NA, ylab=NA)
dev.off()


## TODO: mostly negative correlation of raw frequencies
## -> instead correlate enrichments?
image_matrix(cor(pcad), col=ttcols, breaks=seq(-1,1,length=length(ttcols)+1),
             axis=1:2)

pcad <- pcad - apply(pcad,1,mean) # center genes
pcd <- prcomp(pcad, scale=TRUE)

values <- pcd$sdev^2
vectors <- pcd$rotation
var <- values/sum(values) # % of variance explained by PCn

## add % variance explained to colnames, for axis labels
colnames(pcd$x) <-
    paste0(colnames(pcd$x), " (",round(100*var), "%)")

dense2d(pcd$rotation[,1],pcd$rotation[,2],
             xlab=colnames(pcd$x)[1],
             ylab=colnames(pcd$x)[2])

plotdev(file.path(fig.path,paste0("PCA_AA_degrons")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcd$rotation[,1],pcd$rotation[,2],
     col=NA,#aa.cols[rownames(pcd$rotation)],
     pch=aa.pchs[rownames(pcd$rotation)],
     xlab=colnames(pcd$x)[1],
     ylab=colnames(pcd$x)[2])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcd$rotation[,1],pcd$rotation[,2],labels=rownames(pcd$rotation),
     col=aa.cols[rownames(pcd$rotation)], font=2, cex=1.2)
figlabel("degrons", pos="bottomright", font=2, cex=1.2)
dev.off()
plotdev(file.path(fig.path,paste0("PCA_AA_degrons_13")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcd$rotation[,1],pcd$rotation[,3],
     col=NA,#aa.cols[rownames(pcd$rotation)],
     pch=aa.pchs[rownames(pcd$rotation)],
     xlab=colnames(pcd$x)[1],
     ylab=colnames(pcd$x)[3])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcd$rotation[,1],pcd$rotation[,3],labels=rownames(pcd$rotation),
     col=aa.cols[rownames(pcd$rotation)], font=2, cex=1.2)
figlabel("degrons", pos="bottomright", font=2, cex=1.2)
dev.off()
plotdev(file.path(fig.path,paste0("PCA_AA_degrons_24")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcd$rotation[,2],pcd$rotation[,4],
     col=NA,#aa.cols[rownames(pcd$rotation)],
     pch=aa.pchs[rownames(pcd$rotation)],
     xlab=colnames(pcd$x)[2],
     ylab=colnames(pcd$x)[4])
abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcd$rotation[,2],pcd$rotation[,4],labels=rownames(pcd$rotation),
           col=aa.cols[rownames(pcd$rotation)], font=2, cex=1.2)
figlabel("degrons", pos="bottomright", font=2, cex=1.2)
dev.off()

## correlation of pc components
mpx <- pcx$rotation
colnames(mpx)<- sub("PC","", colnames(mpx))
mpd <- pcd$rotation
colnames(mpd)<- sub("PC","", colnames(mpd))
mpc <- cor(mpd, mpx)
mpt <- mpc
mpt[] <- sub("0\\.",".",round(c(mpc),1))
mpt[mpt=="0"] <- ""
plotdev(file.path(fig.path,paste0("PCA_AA_crosscorrelation")),
        type="png", res=300, width=4,height=4)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
image_matrix(mpc,
             col=ttcols, breaks=seq(-1,1,length=length(ttcols)+1),
             axis=1:2, ylab="degron PC", xlab="AAS PC",
             text=mpt, text.cex=.7)
dev.off()

## shouldn't happen
if ( any(rownames(pcx$rotation)!=rownames(pcd$rotation)) )
    pcd$rotation <- pcd$rotation[rownames(pcx$rotation),]

for ( i in 1:3 ) {
    for ( j in 1:3 ) {
        plotdev(file.path(fig.path,paste0("PCA_AA_crosscorrelation_",i,"_",j)),
                type="png", res=300, width=3.5, height=3.5)
        par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
        plotCor(pcx$rotation[,i], pcd$rotation[,j], density=FALSE, col=NA,
                xlab=paste("AAS -",colnames(pcx$x)[i]),
                ylab=paste("degrons -",colnames(pcd$x)[j]))
        abline(v=0,lwd=.5);abline(h=0,lwd=.5)
        shadowtext(pcx$rotation[,i],pcd$rotation[,j],
                   labels=rownames(pcd$rotation),
                   col=aa.cols[rownames(pcd$rotation)], xpd=TRUE,
                   font=2, cex=1.2)
        dev.off()
    }
}

plotdev(file.path(fig.path,paste0("PCA_AA_crosscorrelation_4_7")),
        type="png", res=300, width=3.5, height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(pcx$rotation[,4], pcd$rotation[,7], density=FALSE, col=NA,
     xlab=paste("AAS -",colnames(pcx$x)[4]),
     ylab=paste("degrons -",colnames(pcd$x)[7]))
    abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcx$rotation[,4],pcd$rotation[,7],labels=rownames(pcd$rotation),
           col=aa.cols[rownames(pcd$rotation)], xpd=TRUE, font=2, cex=1.2)
dev.off()

plotdev(file.path(fig.path,paste0("PCA_AA_crosscorrelation_17_17")),
        type="png", res=300, width=3.5, height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(pcx$rotation[,17], pcd$rotation[,17], density=FALSE, col=NA,
     xlab=paste("AAS -",colnames(pcx$x)[17]),
     ylab=paste("degrons -",colnames(pcd$x)[17]))
    abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcx$rotation[,17],pcd$rotation[,17],labels=rownames(pcd$rotation),
           col=aa.cols[rownames(pcd$rotation)], xpd=TRUE, font=2, cex=1.2)
dev.off()

## cluster PCs ?

xcl <- kmeans(pcx$rotation, 4)
dcl <- kmeans(pcd$rotation, 4)
table(dcl$cluster, xcl$cluster)

## coPCA of degrons and AAS

pcaa <- rbind(pcax,pcad)

    
    types <- c(rep(1,nrow(pcax)),
               rep(2, nrow(pcad)))
    cpca <- multigroup::FCPCA(pcaa, types, Scale = TRUE, graph = TRUE)
    
    plotdev(file.path(fig.path,paste0("PCA_AA_both_FCPCA")),
            type="png", res=300, width=3.5,height=3.5)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plot(cpca$loadings.common[,1:2], col=NA)
    abline(v=0,lwd=.5);abline(h=0,lwd=.5)
    shadowtext(cpca$loadings.common[,1], cpca$loadings.common[,2],
               labels=rownames(cpca$loadings.common),
               col=aa.cols[rownames(cpca$loadings.common)], font=2, cex=1.2)
    dev.off()

pcb <- prcomp(pcaa, scale=TRUE)

values <- pcb$sdev^2
vectors <- pcb$rotation
var <- values/sum(values) # % of variance explained by PCn

## add % variance explained to colnames, for axis labels
colnames(pcb$x) <-
    paste0(colnames(pcb$x), " (",round(100*var), "%)")

dense2d(pcb$rotation[,1],pcb$rotation[,2],
             xlab=colnames(pcb$x)[1],
             ylab=colnames(pcb$x)[2])

plotdev(file.path(fig.path,paste0("PCA_AA_both")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcb$rotation[,1],pcb$rotation[,2],
     col=NA,#aa.cols[rownames(pcb$rotation)],
     pch=aa.pchs[rownames(pcb$rotation)],
     xlab=colnames(pcb$x)[1],
     ylab=colnames(pcb$x)[2])
    abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcb$rotation[,1],pcb$rotation[,2],labels=rownames(pcb$rotation),
     col=aa.cols[rownames(pcb$rotation)])
figlabel("degrons+AAS", pos="bottomright", font=2, cex=1.2)
dev.off()
plotdev(file.path(fig.path,paste0("PCA_AA_both_13")),
        type="png", res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(pcb$rotation[,1],pcb$rotation[,3],
     col=NA,#aa.cols[rownames(pcb$rotation)],
     pch=aa.pchs[rownames(pcb$rotation)],
     xlab=colnames(pcb$x)[1],
     ylab=colnames(pcb$x)[3])
    abline(v=0,lwd=.5);abline(h=0,lwd=.5)
shadowtext(pcb$rotation[,1],pcb$rotation[,3],labels=rownames(pcb$rotation),
     col=aa.cols[rownames(pcb$rotation)], font=2, cex=1.2)
figlabel("degrons+AAS", pos="bottomright", font=2, cex=1.2)
dev.off()

barplot(var*100)


## TODO: how to use coPCA?
if ( FALSE ) {
    ump <- umap::umap(pcaa)
    kcl <- kmeans(ump$layout, 5)
    dense2d(ump$layout[,1], ump$layout[,2])
    plot(ump$layout[,1], ump$layout[,2], col=kcl$cluster)
    ##pcl <- cluster::pam(ump$layout, 5)
    ##plot(ump$layout[,1], ump$layout[,2], col=pcl$clustering)
}
}
### PEPTIDE PROPERTIES
## also do for Main peptides
if ( do.props ) {    
seqs <- apply(pctx[,as.character(-10:10)],1,paste, collapse="") # AAS
deqs <- unlist(lapply(degl, paste, collapse="")) # DEGRONS
aeqs <- apply(pctx[,as.character(-25:-15)],1,paste, collapse="") # AAS left
beqs <- apply(pctx[,as.character(15:25)],1,paste, collapse="") # AAS right

scols <- c("#0000ff","#ff0000", gray.colors(2))
names(scols) <- c("AAS","degrons", "AAS-15", "AAS+15")

plotdev(file.path(fig.path,paste0("peptide_boman")),
        type="png", res=300, width=3.5,height=3.5)
par(mfrow=c(2,1), mai=c(.15,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
histn(boman(seqs), col=paste0(scols[1],77),
      freq=FALSE, breaks=-100:100, xlab="Boman index (protein interaction)",
      xlim=c(-6,12))
hist(boman(deqs), col=paste0(scols[2],77),
     freq=FALSE, breaks=-100:100, add=TRUE)
legend("right", names(scols[1:2]), col=scols, pch=19)
par(mai=c(.1,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
yrev <- par("usr")[4:3]
histn(boman(aeqs), col=paste0(scols[3],77),
      freq=FALSE, breaks=-100:100, add=FALSE, ylim=yrev, axes=FALSE,
      xlim=c(-6,12))
hist(boman(beqs), col=paste0(scols[4],77),
     freq=FALSE, breaks=-100:100, add=TRUE)
legend("right", names(scols[3:4]), col=scols[3:4], pch=19)
axis(3, labels=FALSE)
axis(2)
figlabel("Boman index (protein interaction)", pos="bottomleft")
dev.off()

msh <- do.call(rbind, mswhimScores(seqs))

## NOTE: protein structure class, simply taking the most frequent per peptide
memp <- membpos(seqs)
memd <- membpos(deqs)
mema <- membpos(aeqs)
memb <- membpos(beqs)
cols <- c("Globular","Surface","Transmembrane")
memt <- 
    rbind(AAS=table(unlist(lapply(memp, function(x)
        tail(names(sort(table(x[,"MembPos"]))),1))))[cols],
        degrons=table(unlist(lapply(memd, function(x)
            tail(names(sort(table(x[,"MembPos"]))),1))))[cols],
        "AAS-15"=table(unlist(lapply(mema, function(x)
            tail(names(sort(table(x[,"MembPos"]))),1))))[cols],
        "AAS+15"=table(unlist(lapply(memb, function(x)
            tail(names(sort(table(x[,"MembPos"]))),1))))[cols])
memt <- t(apply(memt,1,function(x) x/sum(x)))

plotdev(file.path(fig.path,paste0("peptide_membpos")),
        type="png", res=300, width=4,height=4)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(memt, beside=TRUE, legend=TRUE, col=scols)
mtext("mempbos, Eisenberg (1984)",1,1.5)
dev.off()

plotdev(file.path(fig.path,paste0("peptide_pI")),
        type="png", res=300, width=3.5,height=3.5)
par(mfrow=c(2,1), mai=c(.15,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
shst <- histn(pI(seqs), col=paste0(scols[1],77),
             breaks=seq(0,14,1),
             freq=FALSE, xlab="peptide isoelectric point")
dhst <- hist(pI(deqs), col=paste0(scols[2],77),
             add=TRUE, breaks=seq(0,14,1), freq=FALSE)
legend("topright", names(scols[1:2]), col=scols, pch=19)
par(mai=c(.1,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
yrev <- par("usr")[4:3]
shst <- histn(pI(aeqs), col=paste0(scols[3],77),
             breaks=seq(0,14,1),axes=FALSE,
             freq=FALSE, ylim=yrev)
dhst <- hist(pI(beqs), col=paste0(scols[4],77),
             add=TRUE, breaks=seq(0,14,1), freq=FALSE)
legend("bottomright", names(scols[3:4]), col=scols[3:4], pch=19)
axis(3, labels=FALSE)
axis(2)
figlabel("peptide isoelectric point", pos="bottom")
dev.off()

plotdev(file.path(fig.path,paste0("peptide_charge")),
        type="png", res=300, width=3.5,height=3.5)
par(mfrow=c(2,1), mai=c(.15,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
shst <- histn(charge(seqs), col=paste0(scols[1],77),
              breaks=-100:100, freq=FALSE, xlab="peptide net charge",
              ylim=c(0,.35), xlim=c(-16,12))
dhst <- hist(charge(deqs), col=paste0(scols[2],77),
             add=TRUE, breaks=-100:100, freq=FALSE)
legend("topright", names(scols), col=scols, pch=19)
par(mai=c(.1,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
yrev <- par("usr")[4:3]
shst <- histn(charge(aeqs), col=paste0(scols[3],77),
             breaks=-100:100,axes=FALSE,
             freq=FALSE, ylim=yrev, xlim=c(-16,12))
dhst <- hist(charge(beqs), col=paste0(scols[4],77),
             add=TRUE, breaks=-100:100, freq=FALSE)
##legend("bottomright", names(scols[3:4]), col=scols[3:4], pch=19)
axis(3, labels=FALSE)
axis(2)
figlabel("peptide net charge", pos="bottom")
dev.off()

plotdev(file.path(fig.path,paste0("peptide_hydrophobicity")),
        type="png", res=300, width=3.5,height=3.5)
par(mfrow=c(2,1), mai=c(.15,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
shst <- histn(hydrophobicity(seqs), col=paste0(scols[1],77),
             breaks=seq(-6,6,.5), freq=FALSE, xlab="peptide hydrophobicity")
dhst <- hist(hydrophobicity(deqs), col=paste0(scols[2],77),
             add=TRUE, breaks=seq(-6,6,.5), freq=FALSE)
legend("topright", names(scols), col=scols, pch=19)
par(mai=c(.1,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
yrev <- par("usr")[4:3]
shst <- histn(hydrophobicity(aeqs), col=paste0(scols[3],77),
             breaks=seq(-6,6,.5),axes=FALSE,
             freq=FALSE, ylim=yrev)
dhst <- hist(hydrophobicity(beqs), col=paste0(scols[4],77),
             add=TRUE, breaks=seq(-6,6,.5), freq=FALSE)
##legend("bottomright", names(scols[3:4]), col=scols[3:4], pch=19)
axis(3, labels=FALSE)
axis(2)
figlabel("peptide hydrophobicity", pos="bottom")
dev.off()

scrp <- do.call(rbind, crucianiProperties(seqs))
dcrp <- do.call(rbind, crucianiProperties(deqs))
acrp <- do.call(rbind, crucianiProperties(aeqs))
bcrp <- do.call(rbind, crucianiProperties(beqs))

plotdev(file.path(fig.path,paste0("peptide_cruciani_pp1")),
        type="png", res=300, width=3.5,height=3.5)
par(mfrow=c(2,1), mai=c(.15,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
shst <- histn(scrp[,1], col=paste0(scols[1],77),
              breaks=seq(-1,1,.1), freq=FALSE, xlab="Cruciano, polarity")
dhst <- hist(dcrp[,1], col=paste0(scols[2],77),
             add=TRUE, breaks=seq(-1,1,.1), freq=FALSE)
legend("topright", names(scols), col=scols, pch=19)
par(mai=c(.1,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
yrev <- par("usr")[4:3]
shst <- histn(acrp[,1], col=paste0(scols[3],77),
             breaks=seq(-1,1,.1),axes=FALSE,
             freq=FALSE, ylim=yrev)
dhst <- hist(bcrp[,1], col=paste0(scols[4],77),
             add=TRUE, breaks=seq(-1,1,.1), freq=FALSE)
##legend("bottomright", names(scols[3:4]), col=scols[3:4], pch=19)
axis(3, labels=FALSE)
axis(2)
figlabel("Cruciano, polarity", pos="bottom")
dev.off()

plotdev(file.path(fig.path,paste0("peptide_cruciani_pp2")),
        type="png", res=300, width=3.5,height=3.5)
par(mfrow=c(2,1), mai=c(.15,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
shst <- histn(scrp[,2], col=paste0(scols[1],77),
              breaks=seq(-1,1,.1), freq=FALSE, xlab="Cruciano, hydrophobicity")
dhst <- hist(dcrp[,2], col=paste0(scols[2],77),
             add=TRUE, breaks=seq(-1,1,.1), freq=FALSE)
legend("topright", names(scols), col=scols, pch=19)
par(mai=c(.1,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
yrev <- par("usr")[4:3]
shst <- histn(acrp[,2], col=paste0(scols[3],77),
             breaks=seq(-1,1,.1),axes=FALSE,
             freq=FALSE, ylim=yrev)
dhst <- hist(bcrp[,2], col=paste0(scols[4],77),
             add=TRUE, breaks=seq(-1,1,.1), freq=FALSE)
##legend("bottomright", names(scols[3:4]), col=scols[3:4], pch=19)
axis(3, labels=FALSE)
axis(2)
figlabel("Cruciano, hydrophobicity", pos="bottom")
dev.off()

plotdev(file.path(fig.path,paste0("peptide_cruciani_pp3")),
        type="png", res=300, width=3.5,height=3.5)
par(mfrow=c(2,1), mai=c(.15,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
shst <- histn(scrp[,3], col=paste0(scols[1],77),
              breaks=seq(-1,1,.1), freq=FALSE, xlab="Cruciano, H-bonding")
dhst <- hist(dcrp[,3], col=paste0(scols[2],77),
             add=TRUE, breaks=seq(-1,1,.1), freq=FALSE)
legend("topright", names(scols), col=scols, pch=19)
par(mai=c(.1,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
yrev <- par("usr")[4:3]
shst <- histn(acrp[,3], col=paste0(scols[3],77),
             breaks=seq(-1,1,.1),axes=FALSE,
             freq=FALSE, ylim=yrev)
dhst <- hist(bcrp[,3], col=paste0(scols[4],77),
             add=TRUE, breaks=seq(-1,1,.1), freq=FALSE)
##legend("bottomright", names(scols[3:4]), col=scols[3:4], pch=19)
axis(3, labels=FALSE)
axis(2)
figlabel("Cruciano, H-bonding", pos="bottom")
dev.off()

}

### INTRODUCE REQUIRED TAGS

## after some WQRK we found a KRAQ in our data, and
## many other patterns, WIVL that soMe of those are
## rather i'NST'ing

## TODO: make useful previous/next filters
## for M and W enrichments
## from I|V|L : W context        - WIVL that soMe are
## from S|T|P|N:M|E| : M context - i'NST'ing
bpd$fromto <- paste0(bpd$from, ":", bpd$to)
bpd$prev <- pctx[,"-1"]
bpd$nxt <- pctx[,"1"]

## tag KRAQ motif specifically

## TODO: grep only in correct position to avoid
## spurious with different AAS?

tctx <- pctx[,as.character(-5:5)]
txsq <- apply(tctx,1, paste, collapse="")
## grep K|RA.Q motif
#kraqidx <- grep("[KR]A[A-Z]{0,1}Q",txsq)
## grep K|RAQ motif
kraqidx <- grep("[KR]AQ",txsq)
## filter only those where Q is the AAS
kraqidx <- kraqidx[bpd$from[kraqidx]=="Q"]
bpd$KRAQ <- rep(FALSE, nrow(bpd))
bpd$KRAQ[kraqidx] <- TRUE

## tag M?TM? motif specifically
tctx <- pctx[,as.character(-5:5)]
txsq <- apply(tctx,1, paste, collapse="")
## grep MTM motif
kraqidx <- unique(c(grep("MT",txsq),
                    grep("TM",txsq)))
## filter only those where T is the AAS
kraqidx <- kraqidx[bpd$from[kraqidx]=="T"]
bpd$TM <- rep(FALSE, nrow(bpd))
bpd$TM[kraqidx] <- TRUE

## filter methionine left and right
bpd$methionine <- rep(FALSE, nrow(bpd))
midx <- grep("M", apply(tctx[,c("-2","-1","1","2")],1, paste, collapse=""))
bpd$methionine[midx] <- TRUE
## filter tryptophane left and right
bpd$tryptophane <- rep(FALSE, nrow(bpd))
midx <- grep("W", apply(tctx[,c("-2","-1","1","2")],1, paste, collapse=""))
bpd$tryptophane[midx] <- TRUE

## branched chain AA encoded
bpd$ftype <- bpd$from
bpd$ftype[bpd$from%in%c("I","L","V")] <- "branched"

## custom filters
filters <- rbind(
    c(column="site", pattern="1"),
    c(column="site", pattern="2"),
    c(column="methionine", pattern="TRUE"),
    c(column="tryptophane", pattern="TRUE"),
    c(column="ftype", pattern="branched"),
    c(column="KRAQ",pattern="TRUE"),
    c(column="miscleavage",pattern="TRUE"),
    c(column="all",pattern="")
    )
apats <- list(c(),
              c())

### all fromto, from and to
if ( !interactive() ) {
    
    ft <- unique(bpd$fromto)
    ft <- cbind(column=rep("fromto", length(ft)),
                pattern=ft)
    filters <- rbind(filters,ft)
    ft <- unique(bpd$from)
    ft <- cbind(column=rep("from", length(ft)),
            pattern=ft)
    filters <- rbind(filters,ft)
    ft <- unique(bpd$to)
    ft <- cbind(column=rep("to", length(ft)),
                pattern=ft)
    filters <- rbind(filters,ft)
}

##TODO:
## 1. align filtered sequences with clustalw via msa(seq)
## 2. split into well aligned subsets
## 3. plot alignments of subsets using AA colors after
##https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3526-6

## analyze AA frequencies
for ( i in 1:nrow(filters) ) {

    column <- filters[i,"column"]
    pattern <- filters[i,"pattern"]

    cat(paste("calculating", pattern, "in", column, "\n"))

    if ( column=="all" ) {
        filt <- rep(TRUE, nrow(bpd))
    } else {
        filt <- as.character(bpd[,column])==pattern
    }
    
    ctx <- pctx[filt,,drop=FALSE]
    rownames(ctx) <-
        paste0(bpd$name[filt],"_",bpd$pos[filt]-dst,"_",bpd$pos[filt]+dst)


    aa1 <- pattern
    if ( length(grep(":",aa1)) ) 
        aa1 <- strsplit(aa1,":")[[1]][1]

    ## generate frequencies at each position
    ctl <- apply(ctx, 2, table)
    ## aaids
    aaids <- unique(unlist(lapply(ctl, names)))
    ctl <- do.call(cbind, lapply(ctl, function(x) x[aaids]))
    rownames(ctl) <- aaids
    ctl[is.na(ctl)] <- 0

    ## remove gaps
    if ( "-" %in% rownames(ctl) ) {
        empty <- ctl["-",]
        ctl <- ctl[rownames(ctl)!="-",]
    }
    ## get frequency per position
    ##ctl <- t(t(ctl)/apply(ctl,2,sum))

    ## log2 mean ratio over all positions
    ctm <- log2(ctl/apply(ctl,1,mean))

    ## filter rare unusual U and X
    ctm <- ctm[rownames(ctm)%in%AAS,]


    ## HYPERGEO TESTS - tigther context
    ovl <- aaProfile(ctx, abc=AAT)
    ## set 0 count to p=1
    ##ovl$p.value[ovl$count==0] <- 1

    ## ADJUST P-VALUES
    if ( FALSE )  # removes a lot of informative structure from the plot
        ovl$p.value[] <- p.adjust(ovl$p.value, method="hochberg")

    ## sort along N- to C-terminus
    ## TODO: implement two-sided sort
    ##ovl <- sortOverlaps(ovl, axis=2, p.min=p.min)

    
    plotdev(file.path(fig.path,paste0("seqcontext_",
                                      column,"_",pattern, "_overlap")),
            type="png", res=300, width=15,height=5)
    par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovl, p.min=1e-10, p.txt=1e-5, xlab=NA,
                 ylab="amino acid", col=ttcols, show.total=TRUE, text.cex=.8)
    mtext("position relative to AAS", 1, 1.3)
    figlabel(paste0(column,"==",pattern), pos="bottomleft", font=2, cex=1.2)
    figlabel(bquote(n[seq]==.(nrow(ctx))), pos="bottomright", font=2, cex=1.2)
    dev.off()
    plotdev(file.path(fig.path,paste0("seqcontext_",
                                      column,"_",pattern, "_dotplot")),
            type="png", res=300, width=15,height=5)
    par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
    dotprofile(ovl, value="ratio", vcols=ttcols, xlab=NA,
               p.dot=1e-10, lg2=TRUE, mxr=2,
               dot.sze=c(.3,2), ylab="amino acid", axis=1:2, show.total=TRUE)
    mtext("position relative to AAS", 1, 1.3)
    figlabel(paste0(column,"==",pattern), pos="bottomleft", font=2, cex=1.2)
    figlabel(bquote(n[seq]==.(nrow(ctx))), pos="bottomright", font=2, cex=1.2)
    dev.off()

    

    ## FILTER BY P-VALUE
    
    ## sort w/o AAS
    ovs <- sortOverlaps(ovl, axis=1, srt=as.character(-5:5))

    ## only use enriched at AAS 
    ## set p.values at AAS site to NA
    ## TODO: base on used filter!
    aasp <- ovs$sign #[,"0"]
    ovs$p.value[aasp==-1] <- 1

    ## FILTER
    ovs <- sortOverlaps(ovs, axis=2, p.min=p.txt, cut=TRUE)

    ## re sort with AAS
    ovc <- sortOverlaps(ovl, axis=1, srt=as.character(-7:7))
    ovc <- sortOverlaps(ovc, axis=2, srt=rownames(ovs$p.value))

    #ovs <- sortOverlaps(ovl, axis=1, srt=as.character(-7:7))
   

    if ( nrow(ovc$p.value)>0 ) {
       
        plotdev(file.path(fig.path,paste0("seqcontext_",
                                          column,"_",pattern,"_overlap_tight")),
                type="png", res=300,
                width=.3*ncol(ovc$p.value)+1,,height=.25*nrow(ovc$p.value)+1)
        par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        plotOverlaps(ovc, p.min=1e-10, p.txt=1e-5, xlab=NA,
                     ylab="amino acid", col=ttcols, show.total=TRUE,text.cex=.8)
        mtext("position relative to AAS", 1, 1.3)
        figlabel(paste0(column,"==",pattern), pos="bottomleft", font=2, cex=1.2)
        figlabel(bquote(n[seq]==.(nrow(ctx))),pos="bottomright",font=2, cex=1.2)
        dev.off()

        plotdev(file.path(fig.path,paste0("seqcontext_",
                                          column,"_",pattern,"_dotplot_tight")),
                type="png", res=300, 
                width=.3*ncol(ovc$p.value)+1,,height=.25*nrow(ovc$p.value)+1)
        par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        dotprofile(ovc, value="ratio", vcols=ttcols, xlab=NA,
                   p.dot=1e-10, lg2=TRUE, mxr=2,
                   dot.sze=c(.3,2),
                   ylab="amino acid", axis=1:2, show.total=TRUE)
        mtext("position relative to AAS", 1, 1.3)
        figlabel(paste0(column,"==",pattern), pos="bottomleft", font=2, cex=1.2)
        figlabel(bquote(n[seq]==.(nrow(ctx))),pos="bottomright",font=2, cex=1.2)
        dev.off()
    }

    ## analyze dimer frequencies - TAKES LONG
    ## TODO: select diAA suggested by monoAA analysis?

    if ( !do.dimer ) next
    
    ctx <- dctx[filt,,drop=FALSE]
    rownames(ctx) <-
        paste0(bpd$name[filt],"_",bpd$pos[filt]-dst,"_",bpd$pos[filt]+dst)
    
    
    dovl <- aaProfile(ctx, verb=0, p.min=.005, abc=diAAT)
    dovs <- sortOverlaps(dovl, axis=1, srt=as.character(-7:7))
    ## TODO: fix sortOverlap2 such that pos and neg are sorted
    ## separately
    dovs <- sortOverlaps2(dovs, axis=2, p.min=1e-10, cut=TRUE)
    
    if ( nrow(dovs$p.value)>0 ) {
        plotdev(file.path(fig.path,paste0("seqcontext_",
                                          column,"_",pattern,
                                          "_diAA")),
                type="png", res=300, 
                width=.3*ncol(dovs$p.value)+1,,height=.25*nrow(dovs$p.value)+1)
        par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        plotOverlaps(dovs, p.min=1e-50, p.txt=1e-25, xlab=NA,
                     ylab="di-amino acid", col=ttcols, show.total=TRUE,
                     text.cex=.8)
        figlabel(paste0(column,"==",pattern), pos="bottomleft", font=2, cex=1.2)
        figlabel(bquote(n[seq]==.(nrow(ctx))),pos="bottomright",font=2, cex=1.2)
        dev.off()
        plotdev(file.path(fig.path,paste0("seqcontext_",
                                          column,"_",pattern,
                                          "_diAA_dotplot")),
                type="png", res=300, 
                width=.3*ncol(dovs$p.value)+1,,height=.25*nrow(dovs$p.value)+1)
        par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        dotprofile(dovs, value="ratio", vcols=ttcols, xlab=NA,
           p.dot=1e-50, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), axis=1:2,
           ylab="di-amino acid", show.total=TRUE)
        figlabel(paste0(column,"==",pattern), pos="bottomleft", font=2, cex=1.2)
        figlabel(bquote(n[seq]==.(nrow(ctx))),pos="bottomright",font=2, cex=1.2)
        dev.off()
    }
}

### WRITE OUT FASTA FILES

fidx <- which(filters[,"column"] %in% c("methionine","tryptophane"))
for ( i in 1:which(filters[,"column"]=="all") ) {
     
    column <- filters[i,"column"]
    pattern <- filters[i,"pattern"]

    if ( column=="all" ) {
        filt <- rep(TRUE, nrow(bpd))
    } else {
        filt <- as.character(bpd[,column])==pattern
    }
    
    ctx <- pctx[filt,]


    sqs <- apply(ctx,1, paste, collapse="")
    fname <- file.path(out.path,paste0("seqcontext_",
                                       column,"_",pattern,".fa"))
    out.rng <- as.character(-sdst:sdst)

    
    sqnms <- tagDuplicates(rownames(ctx))
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) {
        osq <- gsub("-","",paste0(ctx[j,out.rng], collapse=""))
        cat(paste0(">", sqnms[j],"\n", osq, "\n"),
            file=fname, append=TRUE)
    }


    ## required for meme
    fname <- file.path(out.path,paste0("seqcontext_",
                                       column,"_",pattern,"_long.fa"))
    out.rng <- as.character(-sdst:sdst)
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) 
        if ( sum(ctx[j,out.rng]!="-")>7 )
            cat(paste0(">",sqnms[j],"\n",
                       gsub("-","",paste0(ctx[j,out.rng],collapse="")),"\n"),
                file=fname, append=TRUE)

    ## required for MoMo
    fname <- file.path(out.path,paste0("seqcontext_",
                                       column,"_",pattern,"_motif.fa"))
    out.rng <- as.character(-mdst:mdst)
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) 
        if ( sum(ctx[j,out.rng]!="-")==7 )
            cat(paste0(">",sqnms[j],"\n",
                       gsub("-","",paste0(ctx[j,out.rng],collapse="")),"\n"),
                file=fname, append=TRUE)
}
sum(duplicated(apply(ctx,1, paste, collapse="")))

## WRITE OUT FASTA FILES FOR IDENTIFIED PATTERNS ONLY

## cut to tighter context
tctx <- pctx[,as.character(-10:10)]
txsq <- apply(tctx,1, paste, collapse="")
## grep K|RA.Q motif
#kraqidx <- grep("[KR]A[A-Z]{0,1}Q",txsq)
## grep K|RAQ motif
kraqidx <- grep("[KR]AQ",txsq)
## filter only those where Q is the AAS
kraqidx <- kraqidx[bpd$from[kraqidx]=="Q"]
kraqsq <- txsq[kraqidx]
fname <- file.path(out.path,paste0("seqcontext_KRAQ.fa"))
if ( file.exists(fname) ) unlink(fname)
for ( j in 1:length(kraqsq) ) 
    cat(paste0(">",names(kraqsq)[j],"\n", 
               kraqsq[j],"\n"), file=fname, append=TRUE)

## TODO: tag SAAP WITH the identified motif and split Q:G by this,
## GO analysis of motif-containing genes.


if (FALSE) { ## TODO: filter useful
    ## analyze dimer frequencies
    if ( pattern%in%c("Q:G", "T:V") ) { # takes long, only do for selected
        dictx <- character()
        for ( i in 1:(ncol(ctx)-1) )
            dictx <- cbind(dictx, apply(ctx[,i:(i+1)],1,paste,collapse=""))
        colnames(dictx) <- head(colnames(ctx), ncol(ctx)-1)
        rownames(dictx) <- rownames(ctx)

        dovl <- aaProfile(dictx, verb=0, p.min=.005, abc=diAAT)
        dovs <- sortOverlaps(dovl, axis=1, srt=as.character(-7:7))
        dovs <- sortOverlaps(dovs, p.min=1e-5, cut=TRUE)
        
        plotdev(file.path(fig.path,paste0("seqcontext_",
                                          column,"_",pattern,
                                          "_overlap_dimer")),
                type="png", res=300, width=5,
                height=.25*nrow(dovs$p.value)+1)
        par(mai=c(.5,.5,.5,.5), mgp=c(1.5,.3,0), tcl=-.25)
        plotOverlaps(dovs, p.min=1e-10, p.txt=1e-5, xlab=NA,
                     ylab="di-amino acid", col=ttcols, show.total=TRUE,
                     text.cex=.8)
        mtext("position relative to AAS", 1, 1.3)
        figlabel(paste0(column,"==",pattern), pos="bottomleft", font=2, cex=1.2)
        figlabel(bquote(n[seq]==.(nrow(ctx))), pos="bottomright",
                 font=2, cex=1.2)
        dev.off()
    }
    
    plotdev(file.path(fig.path,paste0("seqcontext_",
                                      column,"_",pattern, "_all")),
            height=3.5, width=5, res=300)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    par(xpd=TRUE)
    matplot(colnames(ctm),t(ctm), type="p", lty=1, pch=rownames(ctm),
            xlab="distance from AAS",
            ylab=expression(log[2](AA[i]/bar(AA[51]))))
    figlabel(paste0(column,"==",pattern), pos="topleft",
             region="plot", font=2, cex=1.2)
    dev.off()

    ## CUSTOM
    ## also plot most frequent at -1 and +1
    aa3 <- c("K","R")#names(which.max(table(ctx[,"1"])))
    aa2 <- c("E","D")#names(which.max(table(ctx[,"-1"])))
    aa4 <- c("M")#names(which.max(table(ctx[,"1"])))
    aa5 <- c("A")#names(which.max(table(ctx[,"1"])))
   
    plotdev(file.path(fig.path,paste0("seqcontext_",
                                      column,"_",pattern)),
            height=3.5, width=5, res=300)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
    plot(-dst:dst,apply(ctx,2, function(x) mean(x%in%aa1)), ylim=c(0,1),
         type="h", lwd=2, xlim=c(-15,15),
         xlab="distance from AAS", ylab="AA frequency")
    axis(1, at=-1000:1000, tcl=par("tcl")/2, labels=FALSE)
    abline(h=mean(pctx%in%aa1))
    lines(-dst:dst-.2,apply(ctx,2, function(x) mean(x%in%aa2)),
          col=2, lwd=2, type="h")
    abline(h=mean(pctx%in%aa2),col=2)
    lines(-dst:dst-.1,apply(ctx,2, function(x) mean(x%in%aa3)),
          col=3, lwd=2, type="h")
    abline(h=mean(pctx%in%aa3),col=3)
    lines(-dst:dst+.1,apply(ctx,2, function(x) mean(x%in%aa4)),
          col=4, lwd=2, type="h")
    abline(h=mean(pctx%in%aa4),col=4)
    lines(-dst:dst+.2,apply(ctx,2, function(x) mean(x%in%aa5)),
          col=5, lwd=2, type="h")
    abline(h=mean(pctx%in%aa5),col=5)
    legend("topright", c(paste(aa1,collapse="|"),
                         paste(aa2,collapse="|"),
                         paste(aa3,collapse="|"),
                         paste(aa4,collapse="|"),
                         paste(aa5,collapse="|")), col=c(1,2,3,4,5), lty=1,
           title=paste0("n=",sum(filt)))
    figlabel(paste0(column,"==",pattern), pos="topleft",
             region="plot", font=2, cex=1.2)
    ##    figlabel(paste0(column,"==",pattern), pos="bottomleft")
    dev.off()
    
    
    ## SEQUENCE LOGO
    ## generate position weight matrix
    
    
    aaa <- unique(c(ctx))
    aaa <- aaa[aaa!=""]
    pwm <- matrix(0, ncol=ncol(ctx), nrow=length(aaa))
    colnames(pwm) <- colnames(ctx)
    rownames(pwm) <- aaa
    
    for ( j in 1:ncol(ctx) ) {
        tb <- table(ctx[,j])
        tb <- tb[names(tb)!=""]
        pwm[names(tb),j] <- tb/sum(tb)
    }
    
    plotdev(file.path(fig.path,paste0("seqcontext_",
                                      column,"_",pattern,"_logo_bits")),
            height=3.5, width=5, res=300)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
    lg <- ggseqlogo(pwm[,as.character(-10:10)], method="bits")
    print(lg)
    dev.off()
    
    plotdev(file.path(fig.path,paste0("seqcontext_",
                                      column,"_",pattern,"_logo_prob")),
            height=3.5, width=5, res=300)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
    lg <- ggseqlogo(pwm[,as.character(-10:10)], method="probability")
    print(lg)
    dev.off()
}


## CLEAVAGE BIAS?
## remove all where K|R are the mutation
bpk <- bpd[!bpd$KR,]
## align at K|R
ncut <- 2
ccut <- 2
bps <- substr(bpk$BP,1+ncut,nchar(bpk$BP)-ccut)

## NOTE: 1 [KR]AQ vs. 0 [KR]AG
grep("[KR]A[QG]",dat$BP, value=TRUE)
grep("[KR]A[QG]",bps, value=TRUE)

krs <- gregexpr("[KR]", bps) #bpd$BP)
krs <- lapply(krs, function(x) as.numeric(x[x>0]+ncut))

## only data where K|R appears
haskr <- unlist(lapply(krs, length))>0
krd <- bpk[haskr,]
krs <- krs[haskr]
## take only first!
krs <- unlist(lapply(krs, function(x) x[1]))
kfas <- pfas[krd$protein]

##krs <- krs+krd$pos
dense2d(krd$site, krs)
abline(a=0, b=1)

## TODO: plot to file - AAS accumulate left of miscleavage K
## but this is likely an effect of short K-K peptides not being
## mapped or lost already in MS.
dff <- krd$site -krs
par(yaxs="i")
hist(dff, breaks=-100:100-.5, xlim=range(dff),
     xlab="relative position AAS - K", col="#00000055",
     main="distance of non-cleaved K from AAS")
abline(v=0, col=2, lwd=3)
dev.off()

kctx <- sapply(1:nrow(krd), function(i) {
    ##for ( i in 1:nrow(dat) ){
    fsq <- kfas[[i]]$seq # protein sequence
    nsq <- nchar(fsq)
    
    rrng <- seq(-dst,dst,1)
    sq <- rep("-",length(rrng))
    names(sq) <- rrng

    ## K|R site
    site <- krd$pos[i]-krd$site[i] + krs[i] 
    ##unlist(strsplit(fsq,""))[site]
    
    ## cut range to available
    rng <- (site-dst):(site+dst) 
    rrng <- rrng[rng>0 & rng<=nsq]
    rng <- rng[rng>0 & rng<=nsq]
    asq <-  unlist(strsplit(fsq,""))
    sq[as.character(rrng)] <- asq[rng]
    list(sq)
    ##}
})
##
kctx <- do.call(rbind, kctx)

## CALCULATE POSITION-WISE ENRICHMENTS
ovk <- aaProfile(kctx, abc=AAT)

## set 0 count to p=1
##ovk$p.value[ovk$count==0] <- 1

ovs <- sortOverlaps(ovk, axis=1, srt=as.character(-7:7))

plotdev(file.path(fig.path,paste0("miscleavage_AA_profile")),
        type="png", res=300, width=5,height=5)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovs, p.min=1e-10, p.txt=1e-5, xlab=NA,
             ylab="amino acid", col=ttcols, show.total=TRUE, text.cex=.8)
mtext("position relative to internal K|R", 1, 1.3)
figlabel(paste0("miscleavage"), pos="bottomleft", font=2, cex=1.2)
figlabel(bquote(n[seq]==.(nrow(kctx))), pos="bottomright", font=2, cex=1.2)
dev.off()

source("~/work/mistrans/scripts/saap_utils.R")
plotdev(file.path(fig.path,paste0("miscleavage_AA_profile_dotplot")),
        type="png", res=300, width=5,height=5)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(ovs, value="ratio", vcols=ttcols, xlab=NA,
           p.dot=1e-10, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2),test=TRUE,
           ylab="amino acid", axis=1:2, show.total=TRUE)
mtext("position relative to internal K|R", 1, 1.3)
figlabel(paste0("miscleavage"), pos="bottomleft", font=2, cex=1.2)
figlabel(bquote(n[seq]==.(nrow(kctx))), pos="bottomright", font=2, cex=1.2)
dev.off()

plotdev(file.path(fig.path,paste0("miscleavage_AA_profile_volcano")),
        type="png", res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
volcano(ovs, value="ratio", cut=50, lg2=TRUE, mxr=3, xlim=c(-4,4),
        xlab="log2 frequency ratio")#, adjust="bonferroni")
dev.off()

OVL.MIS <- ovk

aa <- rownames(OVL.MIS$ratio)
aa <- aa[!aa%in%c("K","R")]
mis <- OVL.MIS$ratio[aa,"1"]
aas <- OVL.AA$ratio[aa,"N1"]
##mis <- -log10(mis)* OVL.MIS$sign[aa,"1"]
##aas <- -log10(aas)* OVL.AA$sign[aa,"N1"]
lg2 <- TRUE
if ( lg2 ) {
    mis <- log2(mis)
    aas <- log2(aas)
}
plotdev(file.path(fig.path,paste0("miscleavage_AA_site1")),
        type="png", res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(mis, aas, density=FALSE, col=NA,
        xlab="lg2 ratio, miscleavage site +1",
        ylab="lg2 ratio, AAS at N1 of peptide")
abline(v=ifelse(lg2,0,1), lwd=.5)
abline(h=ifelse(lg2,0,1), lwd=.5)
points(mis, aas, col=aa.cols[aa], pch=aa.pchs[aa], lwd=2)
dev.off()
plotdev(file.path(fig.path,paste0("miscleavage_AA_site1_labels")),
        type="png", res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(mis, aas, density=FALSE, col=NA,
        xlab="lg2 ratio, miscleavage site +1",
        ylab="lg2 ratio, AAS at N1 of peptide")
abline(v=ifelse(lg2,0,1), lwd=.5)
abline(h=ifelse(lg2,0,1), lwd=.5)
points(mis, aas, col=aa.cols[aa], pch=aa.pchs[aa], lwd=2)
shadowtext(mis, aas, labels=aa, pos=4, xpd=TRUE, col="black")
dev.off()



### dimer ANALYSIS

## convert AA matrix to diAA matrix
dctx <- character()
for ( i in 1:(ncol(kctx)-1) )
    dctx <- cbind(dctx, apply(kctx[,i:(i+1)],1,paste,collapse=""))
colnames(dctx) <- head(colnames(kctx), ncol(kctx)-1)
rownames(dctx) <- rownames(kctx)

## analyze dimer frequencies
dovl <- aaProfile(dctx, verb=0, p.min=.005, abc=diAAT)
dovs <- sortOverlaps(dovl, axis=1, srt=as.character(-7:7))

plotdev(file.path(fig.path,paste0("miscleavage_diAA_profile")),
        type="png", res=300, width=5,height=45)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(dovs, p.min=1e-10, p.txt=1e-5, xlab=NA,
             ylab="di-amino acid", col=ttcols, show.total=TRUE, text.cex=.8)
dev.off()
plotdev(file.path(fig.path,paste0("miscleavage_diAA_dotplot")),
        type="png", res=300, width=5,height=45)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(dovs, value="ratio", vcols=ttcols, xlab=NA,
           p.dot=1e-10, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), axis=1:2,
           ylab="di-amino acid", show.total=TRUE)
dev.off()

dovc <- sortOverlaps2(dovs, axis=2, p.min=p.min, cut=TRUE)
mai <- c(.5,.5,.5,.5)
fh <- fw <- .2
nh <- nrow(dovc$p.value) *fh + mai[1] + mai[3]
nw <- ncol(dovc$p.value) *fw + mai[2] + mai[4]

plotdev(file.path(fig.path,paste0("miscleavage_diAA_profile_cut")),
        type="png", res=300, width=nw,height=nh)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(dovc, p.min=1e-10, p.txt=1e-5, xlab=NA,
             ylab="di-amino acid", col=ttcols, show.total=TRUE, text.cex=.8)
dev.off()
plotdev(file.path(fig.path,paste0("miscleavage_diAA_dotplot_cut")),
        type="png", res=300, width=nw,height=nh)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
dotprofile(dovc, value="ratio", vcols=ttcols, xlab=NA,
           p.dot=1e-10, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), axis=1:2,
           ylab="di-amino acid", show.total=TRUE)
dev.off()

## CLEAVAGE BIASES?
aaf <- matrix(0, nrow=length(AAS), ncol=30)
rownames(aaf) <- AAS
for ( i in 1:30 ) {
    tb <- table(unlist(lapply(strsplit(bpd$BP[bpd$from!="Q"],""),
                              function(x) ifelse(i<=length(x),x[i],NA))))
    aaf[names(tb),i] <- tb
}
aaf <- t(t(aaf)/apply(aaf,2,sum,na.rm=TRUE))
matplot(t(aaf), type="l")
points(aaf["Q",])
points(aaf["A",],pch=4)
points(aaf["K",],pch=3, col=2)

### KMER ANALYSIS
library(kmer)

kbins <- as.AAbin(kctx)
kmr.bg <- kcount(kbins, k=2)
kmr <- table(apply(kctx[,as.character(0:1)],1,paste,collapse=""))
sum(kmr.bg[,"KR"])

## MANUAL: clustalx and phylo tree

library(phylogram)
library(dendextend)
fname <- file.path(out.path, paste0("seqcontext_","fromto","_","Q:G",".phb"))
x <- read.dendrogram(file=fname)
labels(x) <- sub("_.*", "", labels(x))
par(cex=.7)
plot(x, yaxt = "n", horiz=TRUE) #, type="triangle")


##library(seqinr)
##virusaln  <- read.alignment(file = fname, format = "phylip")
##library(phangorn)
##plotBS(x)

library(ape)
x <- ape::read.tree(file=fname)
labels(x) <- sub("_.*", "", labels(x))
par(cex=.8)
plot(as.phylo(x), type="fan")
