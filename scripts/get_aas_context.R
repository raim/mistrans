
## load AAS mapping, 
## load protein fasta,
## load transcript fasta,
## get sequences surrounding the AAS and analyze here,
## and export for sequence motif analysis.

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
transcr <- file.path(mam.path,"processedData","coding.fa")

dat <- read.delim(in.file)

## remove by exclude tag
dat <- dat[!dat$exclude,]
## remove those w/o protein info
dat <- dat[!is.na(dat$pos), ]

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

filt <- dat$fromto=="Q:G"
segmenTools::plotdev(file.path(fig.path,paste0("seqcontext_","Q:G")),
                     height=3.5, width=3.5, res=300)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(-dst:dst,apply(pctx[filt,],2,
                    function(x) mean(x=="Q")), ylim=c(0,.5), type="l",
     xlab="distance from AAS", ylab="AA frequency")
abline(h=mean(pctx=="Q"))
lines(-dst:dst,apply(pctx[filt,],2,
                     function(x) mean(x=="N")), col=2)
abline(h=mean(pctx=="N"),col=2)
lines(-dst:dst,apply(pctx[filt,],2,
                     function(x) mean(x=="A")), col=3)
abline(h=mean(pctx=="A"),col=3)
legend("topright", c("Q","N", "A"), col=c(1,2,3), lty=1)
dev.off()

table(pctx[filt,"-1"])
