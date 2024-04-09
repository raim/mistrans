
## load AAS mapping, 
## load protein fasta,
## load transcript fasta,
## get sequences surrounding the AAS and analyze here,
## and export for sequence motif analysis.


source("~/work/mistrans/scripts/saap_utils.R")
options(stringsAsFactors=FALSE)
library(segmenTools)
library(seqinr)

## PATHS AND FILES

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"processedData")

## DISTANCE AROUND AAS TO ANALYZE
dst <- 25  # left/right distance for frequency plots

## INPUT FILES
## list of SAAPs with coordinates
in.file <- file.path(dat.path,"saap_mapped3.tsv")
## RAAS values, TMT level
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")
## protein fasta
pfasta <- file.path(dat.path,"all_proteins.fa")
## coding region fasta and protein-transcript mapping
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
tfasta <- file.path(mam.path,"processedData","coding.fa")

## GET SAAP/BP DATA
dat <- read.delim(in.file)
cat(paste("loaded", nrow(dat), "SAAP/BP\n"))

## get RAAS values, TMT level
tmtf <- read.delim(tmt.file)

## load protein fasta file
pfas <- readFASTA(pfasta, grepID=TRUE)
pidx <- match(dat$protein, names(pfas))

## load transcript fasta file
## NOTE: using ensembl match instead of original proteins with mutations!
tfas <- readFASTA(tfasta, grepID=TRUE)
tidx <- match(dat$transcript, names(tfas))

## SKIP ALL W/O protein or transcript match
if ( interactive() ) 
    table(is.na(tidx), is.na(pidx))

seqna <- is.na(tidx) | is.na(pidx)

cat(paste("removing", sum(seqna),
          "sequences w/o protein or transcript match\n"))
pidx <- pidx[!seqna]
tidx <- tidx[!seqna]
dat <- dat[!seqna,]

pfas <- pfas[pidx]
tfas <- tfas[tidx]

## QC: seqinr::translate transcripts and compare to proteins

## translate to proteins
t2p <- unlist(lapply(tfas, function(x) {
    sub("\\*$","",paste0(translate(unlist(strsplit(x$seq,""))),collapse=""))
}))
## proteins as vectors
p2p <- unlist(lapply(pfas, function(x) x$seq))

wrongp <- which(p2p!=t2p)
cat(paste(length(wrongp),
          "translated sequences do not match protein sequence\n"))
mutp <- grep("_", names(wrongp), invert=TRUE)

cat(paste(length(mutp),
          "non-matching sequences DO NOT contain a _ mutation tag!\n"))

## inspect some
if ( interactive() ) {
    head(wrongp[mutp])
    p2p[wrongp[mutp][1]]
    t2p[wrongp[mutp][1]]
}

## non IG genes
isig <- dat$IG[wrongp[mutp]]

cat(paste(sum(!isig),
          "non-matching sequences ARE NOT immunoglobulins!\n"))
isu <- grep("U",p2p[wrongp[mutp[!isig]]]) 
isu <- wrongp[mutp[!isig]][isu]
cat(paste(length(isu),
          "non-matching sequences contain a selenocysteine\n"))
isx <- grep("X",p2p[wrongp[mutp[!isig]]])
isx <- wrongp[mutp[!isig]][isx]
cat(paste(length(isx),
          "non-matching sequences contain an X\n"))

## unexplained differences - ONLY THREE
unex <- wrongp[mutp[!isig]]
unex <- unex[!unex%in% c(isu,isx)]
## all 3 unexplained have a non AUG start codon,
## 2 CTG and 1 GTG
if ( interactive() ) {
     tfas[[unex[1]]]$seq
     tfas[[unex[2]]]$seq
     tfas[[unex[3]]]$seq
}
cat(paste(length(unex),
          "non-matching sequences have a non-standard START codon (CTG,GTG)",
          "as indicated by manual inspecting\n"))
if ( length(unex)!=3 )
    stop("manual inspection of non-matching transcript and protein sequences",
         "is not valid anymore - RECHECK")

### GET SEQUENCE CONTEXT
pctx <- sapply(1:nrow(dat), function(i) {
    ##for ( i in 1:nrow(dat) ){
    fsq <- pfas[[i]]$seq # protein sequence
    nsq <- nchar(fsq)
    
    rrng <- seq(-dst,dst,1)
    sq <- rep("-",length(rrng)) ## GAP
    names(sq) <- rrng

    ## GET RANGE AROUND AAS
    rng <- (dat$pos[i]-dst):(dat$pos[i]+dst)
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
    paste0(dat$name,"_",dat$pos-dst,"_",dat$pos+dst)

sum(pctx=="-")


### WRITE OUT FASTA FILES

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

