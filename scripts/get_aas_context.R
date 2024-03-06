
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
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")
fig.path <- file.path(proj.path,"figures", "saap_context")
dir.create(fig.path, showWarnings=FALSE)

## INPUT FILES
## list of SAAPs with coordinates
in.file <- file.path(out.path,"saap_mapped3.tsv")
## protein fasta
fasta <- file.path(out.path,"all_proteins.fa")
## coding region fasta
##transcr <- file.path(mam.path,"processedData","coding.fa")

dat <- read.delim(in.file)

## remove by exclude tag
dat <- dat[!dat$exclude,]
## remove those w/o protein info
dat <- dat[!is.na(dat$pos), ]

## remove protein with substitutions at the same site
## TODO: more complex handling of this
dat <- dat[!duplicated(paste(dat$protein,dat$pos)),]
           
## load fasta files, analyze frequencies, and write out
## fasta for sequence analysis.
pfas <- readFASTA(fasta, grepID=TRUE)

## get only the required
pfas <- pfas[dat$protein]

dst <- 25
pctx <- sapply(1:nrow(dat), function(i) {
    ##for ( i in 1:nrow(dat) ){
    fsq <- pfas[[i]]$seq # protein sequence
    nsq <- nchar(fsq)
    
    rrng <- seq(-dst,dst,1)
    sq <- rep("",length(rrng))
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

filters <- rbind(c(column="fromto", pattern="S:G"),
                 c(column="fromto", pattern="T:V"),
                 c(column="fromto", pattern="Q:G"))
                       
for ( i in 1:nrow(filters) ) {

    column <- filters[i,"column"]
    pattern <- filters[i,"pattern"]
    
    filt <- dat[,column]==pattern
    ctx <- pctx[filt,]

    aa1 <- pattern
    if ( length(grep(":",aa1)) ) 
        aa1 <- strsplit(aa1,":")[[1]][1]

    ## plot most frequent at -1 and +1
    aa3 <-
        names(which.max(table(ctx[,"1"])))
    aa2 <-
        names(which.max(table(ctx[,"-1"])))
    
    
    segmenTools::plotdev(file.path(fig.path,paste0("seqcontext_",
                                                   column,"_",pattern)),
                         height=3.5, width=5, res=300)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
    plot(-dst:dst,apply(ctx,2, function(x) mean(x==aa1)), ylim=c(0,1), type="l",
         xlab="distance from AAS", ylab="AA frequency")
    axis(1, at=-1000:1000, tcl=par("tcl")/2, labels=FALSE)
    abline(h=mean(pctx==aa1))
    lines(-dst:dst,apply(ctx,2, function(x) mean(x==aa2)), col=2)
    abline(h=mean(pctx==aa2),col=2)
    lines(-dst:dst,apply(ctx,2, function(x) mean(x==aa3)), col=3)
    abline(h=mean(pctx==aa3),col=3)
    legend("topright", c(aa1,aa2, aa3), col=c(1,2,3), lty=1,
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
    ggseqlogo(pwm[,as.character(-10:10)], method="probability")
    dev.off()
    
}

apply(ctx[,23:29],1, paste, collapse="")
