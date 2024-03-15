
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

source("~/work/mistrans/scripts/saap_utils.R")

library(segmenTools)
require(ggplot2)
require(ggseqlogo)
options(stringsAsFactors=FALSE)

library(Biostrings)
AAS <- sort(unique(GENETIC_CODE))
AAT <- AAS[AAS!="*"]
diAAT <- sort(paste(AAT, rep(AAT,each=length(AAT)),sep=""))

## threshold p-values
p.min <- 1e-10
p.txt <- 1e-5

## p-value colors
docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))

### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"processedData")
out.path <- file.path(proj.path,"processedData","motifs")
fig.path <- file.path(proj.path,"figures", "saap_context")
dir.create(out.path, showWarnings=FALSE)
dir.create(fig.path, showWarnings=FALSE)

## INPUT FILES
## list of SAAPs with coordinates
in.file <- file.path(dat.path,"saap_mapped3.tsv")
## protein fasta
fasta <- file.path(dat.path,"all_proteins.fa")
## coding region fasta
##transcr <- file.path(mam.path,"processedData","coding.fa")


## GET SAAP/BP DATA
dat <- read.delim(in.file)

## remove by exclude tag - IG/Albumin/globin
dat <- dat[!dat$IG ,]
## remove those w/o protein info
dat <- dat[!is.na(dat$pos), ]

## remove protein with substitutions at the same site
## TODO: more complex handling of this
dupls <- duplicated(paste(dat$protein, dat$pos))
cat(paste("removing", sum(dupls),
          "BP with equal position in proteins\n"))
dat <- dat[!dupls,]

## rm AAS w/o from info
## TODO: where do these come from?
dat <- dat[dat$from!="",]

## ONLY USE GOOD BLAST MATCHES
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

##bpd <- bpd[bpd$from!="Q",] # no difference
plotdev(file.path(fig.path,paste0("peptide_AAS")),
        height=3, width=6, res=300)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
barplot(table(bpd$site), xlab="position of AAS in peptide",las=2)
hist(bpd$site/nchar(bpd$SAAP), breaks=seq(0,1,.1),
     xlab="relative position of AAS in peptide", main=NA)
legend("topright", paste(nrow(bpd), "unique BP"))
dev.off()

##bpd <- bpd[bpd$from!="Q",] # no difference
plotdev(file.path(fig.path,paste0("peptide_AAS_n2")),
        height=3, width=6, res=300)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
barplot(table(bpd$site[bpd$site>2]), xlab="position of AAS in peptide",las=2)
hist(bpd$site[bpd$site>2]/nchar(bpd$SAAP[bpd$site>2]), breaks=seq(0,1,.1),
     xlab="relative position of AAS in peptide", main=NA)
legend("topright", paste(nrow(bpd), "unique BP"))
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

asite.bins <- as.character(bpd$site)
asite.bins[bpd$site>10] <- ">10"
##asite.bins[bpd$site>20] <- ">20"
asite.srt <- c(1:10,">10") #,">20")
cls <- clusterCluster(bpd$from, asite.bins, cl2.srt=asite.srt,
                      alternative="two.sided")
##cls <- sortOverlaps(cls, p.min=.1)
plotdev(file.path(fig.path,paste0("peptide_AA_overlap_absolute")),
        height=5, width=5, res=300)
par(mai=c(1,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(cls, p.min=1e-10, p.txt=1e-5, ylab="encoded AA at AAS", xlab=NA,
             show.total=TRUE, show.sig=FALSE)
mtext("absolute position of AAS in peptide", 1, 3.5)
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

## get only the required
pfas <- pfas[bpd$protein]

dst <- 25  # left/right distance for frequency plots
sdst <- 5 # left/right distance for fasta file
mdst <- 3# left/right distance for PTM search with momo

pctx <- sapply(1:nrow(bpd), function(i) {
    ##for ( i in 1:nrow(dat) ){
    fsq <- pfas[[i]]$seq # protein sequence
    nsq <- nchar(fsq)
    
    rrng <- seq(-dst,dst,1)
    sq <- rep("-",length(rrng)) ## GAP
    names(sq) <- rrng

    ## cut range to available
    rng <- (bpd$pos[i]-dst):(bpd$pos[i]+dst)
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

if ( interactive() ) {
    x <-c(table(c(pctx[,20:25]))[AAT])
    y <-c(table(c(pctx[,26:31]))[AAT])
    plot(x, y, col=NA)
    abline(a=0, b=1)
    text(x, y, names(x),xpd=TRUE)
}

## introduce required tags

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

## branched chain AA encoded
bpd$ftype <- bpd$from
bpd$ftype[bpd$from%in%c("I","L","V")] <- "branched"

filters <- rbind(
    c(column="site", pattern="1"),
    c(column="site", pattern="2"),
    c(column="ftype", pattern="branched"),
    c(column="from", pattern="S"),
    c(column="from", pattern="T"),
    c(column="to", pattern="N"),
    c(column="TM",pattern="TRUE"),
    c(column="nxt", pattern="M"),  
    c(column="nxt", pattern="W"),
    c(column="prev", pattern="M"),
    c(column="prev", pattern="W"),
    c(column="from", pattern="H"),
    c(column="from", pattern="W"),
    c(column="from", pattern="Q"),
    c(column="fromto", pattern="S:G"),
    c(column="fromto", pattern="G:S"),
    c(column="fromto", pattern="T:V"),
    c(column="fromto", pattern="Q:G"),
    c(column="KRAQ",pattern="TRUE"),
    c(column="miscleavage",pattern="TRUE"),
    c(column="all",pattern="")
    )
apats <- list(c(),
              c())

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

    ## FILTER BY P-VALUE
    
    ## sort w/o AAS
    ovs <- sortOverlaps(ovl, axis=1, srt=as.character(-7:7))

    ## only use enriched at AAS 
    ## set p.values at AAS site to NA
    ## TODO: base on used filter!
    aasp <- ovs$sign #[,"0"]
    ovs$p.value[aasp==-1] <- 1
    
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
    }
}
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

### WRITE OUT FASTA FILES

for ( i in 1:nrow(filters) ) {
     
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
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) {
        osq <- gsub("-","",paste0(ctx[j,out.rng], collapse=""))
        cat(paste0(">", rownames(ctx)[j],"\n", osq, "\n"),
            file=fname, append=TRUE)
    }


    ## required for meme
    fname <- file.path(out.path,paste0("seqcontext_",
                                       column,"_",pattern,"_long.fa"))
    out.rng <- as.character(-sdst:sdst)
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) 
        if ( sum(ctx[j,out.rng]!="-")>7 )
            cat(paste0(">",rownames(ctx)[j],"\n",
                       gsub("-","",paste0(ctx[j,out.rng],collapse="")),"\n"),
                file=fname, append=TRUE)

    ## required for MoMo
    fname <- file.path(out.path,paste0("seqcontext_",
                                       column,"_",pattern,"_motif.fa"))
    out.rng <- as.character(-mdst:mdst)
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) 
        if ( sum(ctx[j,out.rng]!="-")==7 )
            cat(paste0(">",rownames(ctx)[j],"\n",
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
figlabel(paste0("Keil rules"), pos="bottomleft", font=2, cex=1.2)
figlabel(bquote(n[seq]==.(nrow(kctx))), pos="bottomright", font=2, cex=1.2)
dev.off()

### KMER ANALYSIS
library(kmer)

kbins <- as.AAbin(kctx)
kmr.bg <- kcount(kbins, k=2)
kmr <- table(apply(kctx[,as.character(0:1)],1,paste,collapse=""))
sum(kmr.bg[,"KR"])

## dimer ANALYSIS

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
