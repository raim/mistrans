
### ANALYZE SEQUENCE CONTEXT around AAS

## 1: all -> just do hypergeo enrichment tests,
##    as in previous get_aas_context.R
## 2: analyze position of AAS in peptide.

## Motivated by above analysis, select subsets of BP/sequences,
## and do a DiffLogo and a RAAS profile.

## 1: M or W in -2/-1/+1/+2
## 2: higher than threshold E/D
## 3: AAS in peptide position 1-3.

## ... and finally loop over all of certain type:
## 1: from AA,
## 2: to AA,
## 3: from/to AA,
## 4: from/to by AA prop.


library(stringr)
library(DiffLogo)
library(Biostrings)
AAS <- sort(unique(GENETIC_CODE))
AAT <- AAS[AAS!="*"]
diAAT <- sort(paste(AAT, rep(AAT,each=length(AAT)),sep=""))

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")


mfig.path <- file.path(fig.path,"kraq")
dir.create(mfig.path, showWarnings=FALSE)

### POSITION OF AAS IN PEPTIDES

## Categorize by N+i and C-i
mx <- 10

nt <- bdat$site
ct <- nchar(bdat$BP)-bdat$site +1

## fuse central AAS
nt[nt>mx] <- ct[ct>mx] <- mx+1

## take smaller
## TODO: does this introduce a bias?
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

## RAAS vs. POSITION
## TODO: dotprofile per cancer type and site bins
plotdev(file.path(kfig.path,paste0("peptides_sites_RAAS_cor")),
        height=3.5, width=3.5, res=300)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05)
plotCor(bdat$site, bdat$median, xlab="AAS position in peptide",
        ylab=xl.raas, legpos="topright")
dev.off()

plotdev(file.path(kfig.path,paste0("peptides_sites_RAAS")),
        height=3.5, width=7, res=300)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.05, xaxs="i")
boxplot(bdat$median ~ factor(asite.bins, levels=at.srt),
        ylab=expression(log[10](RAAS)), xlab="Position in peptide", las=2)
dev.off()

## AAS vs. POSITION

cls <- clusterCluster(bdat$from, asite.bins, cl2.srt=asite.srt,
                      alternative="two.sided")
cls <- sortOverlaps2(cls, p.min=p.txt)
plotdev(file.path(kfig.path,paste0("peptides_AAS_encoded")),
        height=5, width=mx/2+1, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=p.min, p.txt=p.txt, ylab="Encoded AA at AAS", xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS, encoded", pos="bottomleft",cex=1.2, font=2)
dev.off()

cls <- clusterCluster(bdat$to, asite.bins, cl2.srt=asite.srt,
                      alternative="two.sided")
cls <- sortOverlaps2(cls, p.min=p.txt)
plotdev(file.path(kfig.path,paste0("peptides_AAS_incorporated")),
        height=5, width=mx/2+1, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=p.min, p.txt=p.txt, ylab="Incorporated AA at AAS",
             xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS, incorporated", pos="bottomleft",cex=1.2, font=2)
dev.off()

cls <- clusterCluster(bdat$fromto, asite.bins, cl2.srt=asite.srt,
                      alternative="two.sided")
cls <- sortOverlaps(cls, p.min=1e-3)
plotdev(file.path(kfig.path,paste0("peptides_AAS_AAStype")),
        height=.2*nrow(cls$p.value)+1, width=mx/2+1, res=300)
par(mai=c(.5,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=p.min, p.txt=p.txt, ylab="AAS type", xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()
## sort by only enriched: TODO: implement in sortOverlaps
tmp <- cls
tmp$p.value[tmp$sign<0] <- 1
tmp <- sortOverlaps(tmp, p.min=1e-3, cut=TRUE)
cls <- sortOverlaps(cls, srt=rownames(tmp$p.value))

plotdev(file.path(kfig.path,paste0("peptides_AAS_AAStype_cut")),
        height=.2*nrow(cls$p.value)+1, width=mx/2+1, res=300)
par(mai=c(.5,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=p.min, p.txt=p.txt, ylab="AAS type", xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

## * GENERAL PEPTIDE ANALYSIS

## TODO
## * all main peptides,
## * kmer analysis, e.g. in all peptides vs. proteome,
##   show only for K|R|A|Q involving kmers,
## * proteome scan for K|RAQ with permutations.

pp.min <- 1e-20
pp.txt <- 1e-10

## from N-term
## TODO: implement min dist from N and C
paa <- strsplit(bdat$BP, "")
paa <- do.call(rbind,lapply(paa, function(x) x[1:7]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps2(ovl, p.min=p.txt)
plotdev(file.path(kfig.path,paste0("peptides_AA_Nterm_BP")),
        height=5, width=5, res=300, type=ftyp)
par(mai=rep(.5,4), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=pp.min, p.txt=pp.txt, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA",
             show.sig=FALSE)
figlabel("BP  w/o Keil", pos="bottomleft", font=2, cex=1.2)
dev.off()

## from N-term - WITH keil rule AA
paa <- strsplit(bdat$BP, "")
paa <- do.call(rbind,lapply(paa, function(x) x[1:7]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT)) #[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps2(ovl, p.min=p.txt)
plotdev(file.path(kfig.path,paste0("peptides_AA_Nterm_BP_withKeil")),
        height=5, width=5, res=300, type=ftyp)
par(mai=rep(.5,4), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=pp.min, p.txt=pp.txt, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA",
             show.sig=FALSE)
figlabel("BP", pos="bottomleft", font=2, cex=1.2)
dev.off()

## from N-term w/o Q->G
paa <- strsplit(bdat$BP[bdat$fromto!="Q:G"], "")
paa <- do.call(rbind,lapply(paa, function(x) x[1:7]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT))#[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps2(ovl, p.min=p.txt)
plotdev(file.path(kfig.path,paste0("peptides_AA_Nterm_BP_noQG")),
        height=5, width=5, res=300, type=ftyp)
par(mai=rep(.5,4), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=pp.min, p.txt=pp.txt, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA",
             show.sig=FALSE)
figlabel("BP w/o Q:G", pos="bottomleft", font=2, cex=1.2)
dev.off()

## from C-term
paa <- strsplit(bdat$BP, "")
paa <- do.call(rbind,lapply(paa, function(x) rev(x)[1:7]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps(ovl, p.min=p.txt)
plotdev(file.path(kfig.path,paste0("peptides_AA_Cterm")),
        height=5, width=5, res=300, type=ftyp)
par(mai=rep(.5,4), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=pp.min, p.txt=pp.txt, show.total=TRUE,
             xlab="Distance from peptide C-terminus", ylab="AA",
             show.sig=FALSE)
dev.off()

## ALL MAIN PEPTIDES

tryptic <- readLines(pep.file)
hist(nchar(tryptic))
## remove all that are in BP
tryptic <- tryptic[!tryptic%in%bdat$BP]
paa <- strsplit(tryptic, "")
paa <- do.call(rbind,lapply(paa, function(x) x[1:7]))
paa[is.na(paa)] <- "-"

ovl <- aaProfile(paa, abc=c(AAT))##[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps(ovl, p.min=1e-50)
plotdev(file.path(kfig.path,paste0("peptides_AA_Nterm_tryptic_withKeil")),
        height=5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.6,.6), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=1e-100, p.txt=1e-50, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA",
             show.sig=FALSE)
figlabel("tryptic w Keil", pos="bottomleft", font=2, cex=1.2)
dev.off()
ovl <- aaProfile(paa, abc=c(AAT[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps(ovl, p.min=1e-50)
plotdev(file.path(kfig.path,paste0("peptides_AA_Nterm_tryptic")),
        height=5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.6,.6), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=1e-100, p.txt=1e-50, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA",
             show.sig=FALSE)
figlabel("tryptic", pos="bottomleft", font=2, cex=1.2)
dev.off()

##
nontry <- readLines(non.file)
hist(nchar(nontry))
## remove all that are in BP
nontry <- nontry[!nontry%in%bdat$BP]

paa <- strsplit(nontry, "")
paa <- do.call(rbind,lapply(paa, function(x) x[1:7]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps(ovl, p.min=1e-50)
plotdev(file.path(kfig.path,paste0("peptides_AA_Nterm_nontry")),
        height=5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.6,.6), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=1e-100, p.txt=1e-50, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA",
             show.sig=FALSE)
figlabel("no trypsin", pos="bottomleft", font=2, cex=1.2)
dev.off()


##
nontry <- readLines(nkr.file)
hist(nchar(nontry))
## remove all that are in BP
nontry <- nontry[!nontry%in%bdat$BP]

paa <- strsplit(nontry, "")
paa <- do.call(rbind,lapply(paa, function(x) x[1:7]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps(ovl, p.min=1e-50)
plotdev(file.path(kfig.path,paste0("peptides_AA_Nterm_noKR")),
        height=5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.6,.6), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=1e-100, p.txt=1e-50, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA",
             show.sig=FALSE)
figlabel("no K/R cleavage", pos="bottomleft", font=2, cex=1.2)
dev.off()

##
nontry <- readLines(akr.file)
hist(nchar(nontry))
## remove all that are in BP
nontry <- nontry[!nontry%in%bdat$BP]

paa <- strsplit(nontry, "")
paa <- do.call(rbind,lapply(paa, function(x) x[1:7]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT[!AAT%in%c("K","R","P","D","E")]))
ovl <- sortOverlaps(ovl, p.min=1e-50)
plotdev(file.path(kfig.path,paste0("peptides_AA_Nterm_ArgC_LysC")),
        height=5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.6,.6), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=1e-100, p.txt=1e-50, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA",
             show.sig=FALSE)
figlabel("ArgC/LysC", pos="bottomleft", font=2, cex=1.2)
dev.off()


### WHOLE PROTEOME KMER ANALYSIS

## ? K|RA, AQ, K|RAQ, compared to other K|R 

library(ape)
library(kmer)

## proteome kmer count
pfasta <- file.path(out.path,"all_proteins.fa")

fas <- readFASTA(pfasta, grepID=TRUE)
fas <- fas[names(fas)%in%MANES]
seql <- lapply(fas, function(x) x$seq)
conc <- paste0(unlist(seql), collapse="")
kproteome <- kcount(as.AAbin(conc), k=2)
kproteome <- apply(kproteome,2,sum)

## TODO: K|R.x

## tryptic peptide kmer count
tryptic <- readLines(pep.file)
## remove all that are in BP
tryptic <- tryptic[!tryptic%in%bdat$BP]
conc <- paste0(tryptic, collapse="")
ktryptic <-kcount(as.AAbin(conc), k=2)

## base peptide kmer count
conc <- paste0(bdat$BP, collapse="")
kbase <-kcount(as.AAbin(conc), k=2)

## cleavage site kmer count
paa <- strsplit(bdat$AA, "")
paa <- do.call(rbind, paa)
DST <- (ncol(paa)-1)/2
colnames(paa) <- -DST:DST

## get cleavage site from AA,
## analyze upstream site, it is NOT always K|R

## diAA dat cleavge site
cuts <- sapply(bdat$site, function(x)
    list(paste0("-",as.character(c(x,x-1)))))
pal <- rep("", nrow(paa))
for ( i in seq_along(cuts) )
    if ( all(cuts[[i]]%in%colnames(paa)) ) ## NOTE: loosing sites >DST!
        pal[i] <-paste0(paa[i,cuts[[i]]], collapse="")
pal[pal==""] <- ">dst" # TODO: get longer DST in map_peptides3
kcuts <- table(pal)

## NON-TRYPTIC CLEAVAGE BIAS?
## * by RAAS 
## * TODO: overlaps with AAS type?
bdat$cleavage2 <- pal
bdat$cleavage <- unlist(lapply(strsplit(pal,""), function(x) x[1]))

plotdev(file.path(kfig.path,paste0("cleavage_diAA_RAAS")),
        height=3.5, width=10, res=300, type=ftyp)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
bp <- boxplot(bdat$median ~ bdat$cleavage2, las=2,
              ylab=expression(log[10](bar(RAAS))),
              xlab="non-tryptic cleavage sites")
axis(3, at=1:length(bp$names), labels=table(bdat$cleavage2)[bp$names], las=2)
axis(4)
dev.off()

plotdev(file.path(kfig.path,paste0("cleavage_monoAA_RAAS")),
        height=3.5, width=5, res=300, type=ftyp)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
bp <- boxplot(bdat$median ~ bdat$cleavage, las=2,
              ylab=expression(log[10](bar(RAAS))),
              xlab="non-tryptic cleavage sites")
axis(3, at=1:length(bp$names), labels=table(bdat$cleavage)[bp$names], las=2)
stripchart(median ~ cleavage, vertical = TRUE, data = bdat, 
           method = "jitter", add = TRUE, pch = 20,
           col="#00000077", cex=1,
           axes=FALSE)
dev.off()

## CLEAVAG SITE FREQUENCIES
## total frequencies at cleavage site vs. proteomic
## NOTE: compare only proper K|R sites!
knms <- grep("^[KR]", diAAT, value=TRUE)
fcuts <- rbind(BP=kcuts[knms],
               proteome=kproteome[knms])
fcuts <- fcuts[,order(fcuts[2,], decreasing=TRUE)]
fcuts <- fcuts/apply(fcuts, 1, sum)

plotdev(file.path(kfig.path,paste0("cleavage_frequencies_cor")),
        height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, family="monospace")
plotCor(fcuts[2,], fcuts[1,], ##outliers=fcuts[1,]<.001,
        density=FALSE, xlim=c(0,.08), ylim=c(0,.08),
        xlab="proteomic frequency", ylab="BP frequency")
text(fcuts[2,], fcuts[1,], labels=colnames(fcuts), pos=4)
abline(a=0, b=1)
dev.off()

plotdev(file.path(kfig.path,paste0("cleavage_frequencies_bar")),
        height=3.5, width=10, res=300, type=ftyp)
par(mai=c(.5,1,.1,.1), mgp=c(2.75,.3,0), tcl=-.25, family="monospace", xaxs="i")
barplot(fcuts, las=2, beside=TRUE, legend=TRUE, ylab="rel. frequency")
dev.off()

## overall frequencies

kmers <- rbind(proteome=kproteome,
               tryptic=ktryptic[1,names(kproteome)],
               base=kbase[1,names(kproteome)])

## normalize each per total kmer count
## TODO: normalize below, for R|K only.
fmers <- kmers/apply(kmers,1,sum,na.rm=TRUE)

plotdev(file.path(kfig.path,paste0("cleavage_KR_kmers")),
        height=7, width=10, res=300, type=ftyp)
par(mai=c(.5,1,.1,.1), mgp=c(2.75,.3,0), tcl=-.25, family="monospace", xaxs="i")
par(mfcol=c(2,1))
barplot(fmers[,grep("^R",colnames(kmers))], ylab="rel. frequency",
        beside=TRUE, col=1:nrow(rmers), legend=TRUE)
barplot(fmers[,grep("^K",colnames(kmers))], ylab="rel. frequency",
        beside=TRUE, col=1:nrow(fmers), legend=TRUE)
dev.off()

## add K and R frequencies
kkmer <- fmers[,grep("^K",colnames(kmers))]
rkmer <- fmers[,grep("^R",colnames(kmers))]
colnames(kkmer) <- sub("^K", "", colnames(kkmer))
colnames(rkmer) <- sub("^R", "", colnames(rkmer))

rkmer <- rkmer + kkmer[,colnames(rkmer)]

plotdev(file.path(kfig.path,paste0("cleavage_KR_fused_kmers")),
        height=3.5, width=10, res=300, type=ftyp)
par(mai=c(.5,1,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, family="monospace", xaxs="i")
barplot(rkmer, beside=TRUE, col=1:nrow(rkmer), legend=TRUE,
        xlab="K|R. Dimers")
mtext("rel. frequency", 2, 2.75)
dev.off()

## COMPARE CLEAVAGE SITES DIRECTLY
## add first nucleotide of tryptic digest: just also have R|K

## first AA positions in BP, SAAP and Main 
bp1 <- table(unlist(lapply(strsplit(bdat$BP,""), function(x) x[1])))
saap1 <- table(unlist(lapply(strsplit(bdat$SAAP,""), function(x) x[1])))
trp1 <- table(unlist(lapply(strsplit(tryptic,""), function(x) x[1])))
bp2 <- table(unlist(lapply(strsplit(bdat$BP,""), function(x) x[2])))
saap2 <- table(unlist(lapply(strsplit(bdat$SAAP,""), function(x) x[2])))
trp2 <- table(unlist(lapply(strsplit(tryptic,""), function(x) x[2])))
bp3 <- table(unlist(lapply(strsplit(bdat$BP,""), function(x) x[3])))
saap3 <- table(unlist(lapply(strsplit(bdat$SAAP,""), function(x) x[3])))
trp3 <- table(unlist(lapply(strsplit(tryptic,""), function(x) x[3])))

## proteome-wide frequency of AA next to K|R
k.cnt <- kproteome[grep("^[K]", names(kproteome))]
r.cnt <- kproteome[grep("^[R]", names(kproteome))]
names(k.cnt) <-  substr(names(k.cnt), 2, 2)
names(r.cnt) <-  substr(names(r.cnt), 2, 2)
kr.cnt <- k.cnt + r.cnt[names(k.cnt)]
#k.frq <- kr.cnt/sum(kr.cnt)

### TODO
## * add pattern match K|R.x counts via regex.
pat <- "[KR]"
krxq <- gregexpr(pattern=pat, seql)
kcnt1 <- table(unlist(sapply(seq_along(seql), function(x) {
    unlist(strsplit(seql[[x]],""))[krxq[[x]]+1]
})))
kcnt2 <- table(unlist(sapply(seq_along(seql), function(x) {
    unlist(strsplit(seql[[x]],""))[krxq[[x]]+2]
})))
kcnt3 <- table(unlist(sapply(seq_along(seql), function(x) {
    unlist(strsplit(seql[[x]],""))[krxq[[x]]+3]
})))

cmer3 <- rbind(uniprot=kcnt3[AAT],
               tryptic=trp3[AAT],
               BP=bp3[AAT],
               SAAP=saap3[AAT])
cmer3 <- cmer3/apply(cmer3,1,sum,na.rm=TRUE)

cmer2 <- rbind(uniprot=kcnt2[AAT],
               tryptic=trp2[AAT],
               BP=bp2[AAT],
               SAAP=saap2[AAT])
cmer2 <- cmer2/apply(cmer2,1,sum,na.rm=TRUE)

cmer1 <- rbind(uniprot=kcnt1[AAT],
               tryptic=trp1[AAT],
               BP=bp1[AAT],
               SAAP=saap1[AAT])
cmer1 <- cmer1/apply(cmer1,1,sum,na.rm=TRUE)


## sort by frequency in BP
cmer3 <- cmer3[,order(cmer1["uniprot",], decreasing=TRUE)]
cmer2 <- cmer2[,order(cmer1["uniprot",], decreasing=TRUE)]
cmer1 <- cmer1[,order(cmer1["uniprot",], decreasing=TRUE)]

## R4 colors for gas exchange rates
col4 = c("#000000", "#DF536B", "#61D04F", "#2297E6", "#28E2E5")
cols <- col4[-3]

          

plotdev(file.path(kfig.path,paste0("cleavage_KR_positions")),
        height=5, width=5, res=300, type=ftyp)
par(mai=c(.3,.3,.05,.1),
    mgp=c(1.3,.3,0), tcl=-.25, family="monospace", xaxs="i")
par(mfcol=c(3,1))
barplot(cmer1, beside=TRUE, col=cols, legend=TRUE,
        xlab=expression("proteomic K|R vs. 1st position in peptides"))
barplot(cmer2, beside=TRUE, col=cols, legend=FALSE,
        xlab=expression("proteomic K|R vs. 2nd position in peptides"),
        ylab="rel. frequencies")
barplot(cmer3, beside=TRUE, col=cols, legend=FALSE,
        xlab=expression("proteomic K|R vs. 3rd position in peptides"))
dev.off()
## scatter plot first position

plotdev(file.path(kfig.path,paste0("cleavage_KR_frequencies_cor")),
        height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, family="monospace")
plotCor(cmer1["uniprot",], cmer1["BP",], ##outliers=fcuts[1,]<.001,
        density=FALSE, ##xlim=c(0,.08), ylim=c(0,.08),
        xlab="proteomic frequency", ylab="BP frequency")
text(cmer1["uniprot",], cmer1["BP",], labels=colnames(cmer1), pos=4, xpd=TRUE)
abline(a=0, b=1)
dev.off()


## PLOT AGAIN w/o KRP Frequencies


AATkeil <- AAT[!AAT%in%c("K","R","P")]

cmer3 <- rbind(uniprot=kcnt3[AATkeil],
               tryptic=trp3[AATkeil],
               BP=bp3[AATkeil],
               SAAP=saap3[AATkeil])
cmer3 <- cmer3/apply(cmer3,1,sum,na.rm=TRUE)

cmer2 <- rbind(uniprot=kcnt2[AATkeil],
               tryptic=trp2[AATkeil],
               BP=bp2[AATkeil],
               SAAP=saap2[AATkeil])
cmer2 <- cmer2/apply(cmer2,1,sum,na.rm=TRUE)

cmer1 <- rbind(uniprot=kcnt1[AATkeil],
               tryptic=trp1[AATkeil],
               BP=bp1[AATkeil],
               SAAP=saap1[AATkeil])
cmer1 <- cmer1/apply(cmer1,1,sum,na.rm=TRUE)


## sort by frequency in BP
cmer3 <- cmer3[,order(cmer1["uniprot",], decreasing=TRUE)]
cmer2 <- cmer2[,order(cmer1["uniprot",], decreasing=TRUE)]
cmer1 <- cmer1[,order(cmer1["uniprot",], decreasing=TRUE)]

## R4 colors for gas exchange rates
col4 = c("#000000", "#DF536B", "#61D04F", "#2297E6", "#28E2E5")
cols <- col4[-3]

          

plotdev(file.path(kfig.path,paste0("cleavage_KR_positions_noKeil")),
        height=5, width=5, res=300, type=ftyp)
par(mai=c(.3,.3,.05,.1),
    mgp=c(1.3,.3,0), tcl=-.25, family="monospace", xaxs="i")
par(mfcol=c(3,1))
barplot(cmer1, beside=TRUE, col=cols, legend=TRUE,
        xlab=expression("proteomic K|R vs. 1st position in peptides"))
barplot(cmer2, beside=TRUE, col=cols, legend=FALSE,
        xlab=expression("proteomic K|R vs. 2nd position in peptides"),
        ylab="rel. frequencies")
barplot(cmer3, beside=TRUE, col=cols, legend=FALSE,
        xlab=expression("proteomic K|R vs. 3rd position in peptides"))
dev.off()

## scatter plot first position

plotdev(file.path(kfig.path,paste0("cleavage_KR_frequencies_cor_noKeil")),
        height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, family="monospace")
plotCor(cmer1["uniprot",], cmer1["BP",], ##outliers=fcuts[1,]<.001,
        density=FALSE, ##xlim=c(0,.08), ylim=c(0,.08),
        xlab="proteomic frequency", ylab="BP frequency")
text(cmer1["uniprot",], cmer1["BP",], labels=colnames(cmer1), pos=4, xpd=TRUE)
abline(a=0, b=1)
dev.off()

## K|RAQ proteome-wide STATISTICS: more then expected?

## TODO:
## * account for protein ends in total number of possible motifs (L-m+1)?
## * just shuffle sequences with 2mer preservation, and count occurence
##   of K|RAQ, using universalmotif::shuffle_sequences,
## * get proteome-wide combinations of K|R .. Q
##   and count the AA between.
## * ASAquick: is the KRAQ buried?
##   count all KRAQ in proteome,
##   and compare to ASAquick prediction.
boxplot(bdat$ASAquick ~ bdat$fromto)

## TODO:
## * read the Introduction to sequence motifs by Benjamin Jean-Marie Tremblay 
## * use universalmotifs to shuffle with 2mer preservation.

## DEFINE MOTIF
pat <- "[KR]AQ"
m <- 3

## SIMPLE MOTIF PROB.

## count in base peptides
bcnt <- stringr::str_count(bdat$AA, pat)
bpcnt <- table(unlist(strsplit(bdat$AA,"")))
bpfrq <- bpcnt/sum(bpcnt)
Lbp <- nchar(bdat$AA)
p <- (bpfrq["K"]+bpfrq["R"]) * bpfrq["A"] * bpfrq["Q"]
expected <- sum(Lbp-m+1)*p

cat(paste("found", sum(bcnt), "K|RAQ in the base peptides, expected:",
          round(expected), "\n"))

## COUNT IN PROTEOME
conc <- paste0(unlist(seql), collapse="")
aacnt <- table(unlist(strsplit(conc,"")))
aafrq <- aacnt/sum(aacnt)
L <- nchar(unlist(seql))

p <- (aafrq["K"]+aafrq["R"]) * aafrq["A"] * aafrq["Q"]
expected <- sum(L-m+1)*p

## count pattern in all MANE proteins
gcnt <- sum(stringr::str_count(seql, pat))

cat(paste("found", gcnt, "K|RAQ in the proteome, expected:",
          round(expected), "\n"))


## TODO: cluster sequences kmer content, as eg.
## shown in univseralmotifs package.
library(universalmotif)
library(Biostrings)

numperm <- 1000
nums <- rep(0, length(numperm))
pcnt <- 0
if ( !interactive() )
    for ( i in 1:numperm ) {
        cat(paste(i,","))
        aaset <- AAStringSet(x=unlist(seql))
        seqr <- shuffle_sequences(aaset, k=2)
        ##vcnt <- vcountPattern(pat, seqr)
        rcnt <- sum(stringr::str_count(as.character(seqr), pat))
        nums[i] <- rcnt
        if ( rcnt >= gcnt )
            pcnt <- pcnt+1
    }

## NOTE: KRAQ as expected also when 2mer shuffling, p~0.4
## e.g.:
presult <-  467
if ( !interactive() ) presult <- pcnt

cat(paste("found", sum(gcnt), "K|RAQ in the proteome, expected:",
          round(expected), "\n",
          "p-value with 2-mer conserving shuffled sequences:",
          presult/numperm,"\n"))
