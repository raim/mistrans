
## load AAS mapping, 
## load protein fasta,
## load transcript fasta,
## get sequences surrounding the AAS and analyze here,
## and export for sequence motif analysis.

## TODO:
## * filter unique proteins?

library(segmenTools)
require(ggplot2)
require(ggseqlogo)


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

dat <- read.delim(in.file)

## remove by exclude tag
dat <- dat[!dat$exclude ,]
## remove those w/o protein info
dat <- dat[!is.na(dat$pos), ]

## tag trypsin 
dat$KR <- dat$from%in%c("K","R")

## remove protein with substitutions at the same site
## TODO: more complex handling of this
dat <- dat[!duplicated(paste(dat$protein,dat$pos)),]
           
## load fasta files, analyze frequencies, and write out
## fasta for sequence analysis.
pfas <- readFASTA(fasta, grepID=TRUE)

## get only the required
pfas <- pfas[dat$protein]

dst <- 25  # left/right distance for frequency plots
sdst <- 5 # left/right distance for fasta file
mdst <- 3# left/right distance for PTM search with momo

pctx <- sapply(1:nrow(dat), function(i) {
    ##for ( i in 1:nrow(dat) ){
    fsq <- pfas[[i]]$seq # protein sequence
    nsq <- nchar(fsq)
    
    rrng <- seq(-dst,dst,1)
    sq <- rep("_",length(rrng))
    names(sq) <- rrng

    ## cut range to available
    rng <- (dat$pos[i]-dst):(dat$pos[i]+dst)
    rrng <- rrng[rng>0 & rng<=nsq]
    rng <- rng[rng>0 & rng<=nsq]
    sq[as.character(rrng)] <- unlist(strsplit(fsq,""))[rng]
    list(sq)
    ##}
})
##
pctx <- do.call(rbind, pctx)

dat$fromto <- paste0(dat$from, ":", dat$to)

filters <- rbind(c(column="from", pattern="Q"),
                 c(column="fromto", pattern="S:G"),
                 c(column="fromto", pattern="T:V"),
                 c(column="fromto", pattern="Q:G"))
apats <- list(c(),
              c())
##TODO:
## 1. align with clustalw via msa(seq)
## 2. split into well aligned subsets
## 3. plot alignments of subsets using AA colors after
## https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3526-6
for ( i in 1:nrow(filters) ) {

    column <- filters[i,"column"]
    pattern <- filters[i,"pattern"]
    
    filt <- dat[,column]==pattern
    ctx <- pctx[filt,]
    rownames(ctx) <-
        paste0(dat$ensembl[filt],"_",dat$pos[filt]-dst,"_",dat$pos[filt]+dst)
    saabp <- paste0(dat$SAAP[filt], "/", dat$BP[filt])

    aa1 <- pattern
    if ( length(grep(":",aa1)) ) 
        aa1 <- strsplit(aa1,":")[[1]][1]

    ## also plot most frequent at -1 and +1
    aa3 <- c("K","R")#names(which.max(table(ctx[,"1"])))
    aa2 <- c("E","D")#names(which.max(table(ctx[,"-1"])))
    aa4 <- c("M")#names(which.max(table(ctx[,"1"])))
    aa5 <- c("A")#names(which.max(table(ctx[,"1"])))

    
    segmenTools::plotdev(file.path(fig.path,paste0("seqcontext_",
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
    figlabel(paste0(column,"==",pattern), pos="bottomleft")
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
    ggseqlogo(pwm[,as.character(-10:10)], method="bits")
    dev.off()
    
    plotdev(file.path(fig.path,paste0("seqcontext_",
                                      column,"_",pattern,"_logo_prob")),
            height=3.5, width=5, res=300)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
    lg <- ggseqlogo(pwm[,as.character(-10:10)], method="probability")
    print(lg)
    dev.off()

    sqs <- apply(ctx,1, paste, collapse="")
    fname <- file.path(out.path,paste0("seqcontext_",
                                       column,"_",pattern,".fa"))
    out.rng <- as.character(-sdst:sdst)
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) {
        osq <- gsub("_","",paste0(ctx[j,out.rng], collapse=""))
        cat(paste0(">", osq,"_", rownames(ctx)[j],"\n", osq, "\n"),
            file=fname, append=TRUE)
    }

    ## special subsets, post-discover

    ## required for meme
    fname <- file.path(out.path,paste0("seqcontext_",
                                       column,"_",pattern,"_long.fa"))
    out.rng <- as.character(-sdst:sdst)
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) 
        if ( sum(ctx[j,out.rng]!="_")>7 )
            cat(paste0(">",rownames(ctx)[j],"\n",
                       gsub("_","",paste0(ctx[j,out.rng],collapse="")),"\n"),
                file=fname, append=TRUE)

    ## required for MoMo
    fname <- file.path(out.path,paste0("seqcontext_",
                                       column,"_",pattern,"_motif.fa"))
    out.rng <- as.character(-mdst:mdst)
    if ( file.exists(fname) ) unlink(fname)
    for ( j in 1:nrow(ctx) ) 
        if ( sum(ctx[j,out.rng]!="_")==7 )
            cat(paste0(">",rownames(ctx)[j],"\n",
                       gsub("_","",paste0(ctx[j,out.rng],collapse="")),"\n"),
                file=fname, append=TRUE)
}
sum(duplicated(apply(ctx,1, paste, collapse="")))

## TODO: tag SAAP WITH the identified motif and split Q:G by this,
## GO analysis of motif-containing genes.

## CLEAVAGE BIASES?
library(Biostrings)
aas <- sort(unique(GENETIC_CODE))
aaf <- matrix(0, nrow=length(aas), ncol=30)
rownames(aaf) <- aas
for ( i in 1:30 ) {
    tb <- table(unlist(lapply(strsplit(dat$BP[dat$from!="Q"],""),
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

plot(as.hclust(x))

##library(seqinr)
##virusaln  <- read.alignment(file = fname, format = "phylip")
##library(phangorn)
##plotBS(x)

library(ape)
x <- ape::read.tree(file=fname)
labels(x) <- sub("_.*", "", labels(x))
par(cex=.8)
plot(as.phylo(x), type="fan")
