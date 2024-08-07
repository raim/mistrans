library(segmenTools)

## TEST GENOMIC AND TRANSCRIPT COORDINATES OF AAS from map_peptides.R

dat.path <- "/home/raim/data/"
mam.path <- file.path(dat.path, "mammary")
mis.path <- file.path(dat.path, "mistrans")


## BP/SAAP - AAS mapping by map_peptides.R
in.file <- file.path(mis.path, "processedData", "saap_mapped.tsv")

## protein fasta
pfas <- file.path(mis.path, "processedData", "all_proteins.fa")
## transcript fasta
tfas <- file.path(mam.path, "originalData", "Homo_sapiens.GRCh38.cdna.all.fa.gz")
## chromosome fasta
gfas <- file.path(mam.path, "originalData", "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")

## cds coordinates
cds.file <- file.path(mam.path, "processedData", "protein_cds_coordinates.tsv")
## utr lengths
utr.file <- file.path(mam.path, "processedData", "transcript_utr_lengths.tsv")

## BP/SAAP
dat <- read.delim(in.file)

## only consider BP/SAAP with an identified protein
## with a blast hit tagged as "good"
dat <- dat[dat$match=="good",]

## proteins
prt <- readFASTA(pfas, grepID=TRUE)
prt <- lapply(prt, function(x) unlist(strsplit(x$seq,"")))

## compare to "from" AA in BP
ptest <- unlist(sapply(1:nrow(dat), function(i) 
    prt[[dat$protein[i]]][dat$pos[i]] == dat$from[i]))
paa <- unlist(sapply(1:nrow(dat), function(i) 
    prt[[dat$protein[i]]][dat$pos[i]]))

cat(paste(sum(ptest), "correct AAs;", sum(!ptest), "wrong AAs\n"))

### NOTE: all non-matching BP/SAAP are tagged as bad/wrong blast match!

## TEST TRANSCRIPT COORS
trn <- readFASTA(tfas, grepID=TRUE)
trn <- lapply(trn, function(x) unlist(strsplit(x$seq,"")))

## remove version tag
names(trn) <- sub("\\.[0-9]*$", "", names(trn))

## test only were codon is available!
cdat <- dat[dat$codon!="",]

ttest <- unlist(sapply(1:nrow(cdat), function(i) {
    res <- NA
    if ( !is.na(cdat$tpos[i]) ) {
        tps <- cdat$tpos[i]
        tps <- (tps-1):(tps+1) ## NOTE: position is 2nd codon pos
        res <- paste(trn[[cdat$transcript[i]]][tps],collapse="")==cdat$codon[i]
    }
    res
}))
cat(paste(sum( ttest), "correct codons;",
          sum(!ttest), "wrong codons\n"))

### NOTE: all non-matching BP/SAAP are tagged as bad/wrong blast match!

## TEST GENOME COORS
gen <- readFASTA(gfas, grepID=TRUE)
## reduce to sequences we have
gen <- gen[names(gen)%in%dat$chr]
gen <- lapply(gen, function(x) unlist(strsplit(x$seq,"")))

## reduce to BP/SAAP where we have both genome coors and a codon
cdat <- dat[dat$chr!="" & dat$codon!="",]

## CDS structure
cds <- read.delim(cds.file, header=FALSE, sep=" ")
cdl <- split(cds[,2:ncol(cds)], cds[,1])
##cdl <- lapply(strsplit(cds[,1], ";"), as.numeric)
##names(cdl) <- rownames(cds)

utr <- read.delim(utr.file, row.names=1)


gtest <- rep(FALSE, nrow(cdat))

for ( i in 1:nrow(cdat) ) {

    ## get recorded genome coordinates
    chr <- cdat$chr[i]
    cps <- cdat$coor[i] 
    str <- cdat$strand[i]

    ## account for codons that span splice sites 
    ## expand CDS to vector of coding sequences
    cd <- cdl[[cdat$ensembl[i]]]
    coors <- unlist(apply(cd,1,function(x) x[2]:x[3]))
    ## find coordinate in CDS coordinates
    dst <- abs(coors-cps)

    ## recorded position must be a CDS coordinate
    if ( min(dst)!=0 ) 
        stop("not in CDS", min(dst),"\t",
             i,cdat$name[i],":\t",str,sq,cdat$codon[i])

    idx <- which(dst==0)
    cpss <- coors[(idx-1):(idx+1)] # get coordinates of codon

    ## get codon and compare to recorded codon (via transcript)
    sq <- paste(gen[[chr]][cpss],collapse="")
    if ( str=="-" ) 
        sq <- revcomp(sq)
    if ( sq!=cdat$codon[i] ) 
        cat(paste(i,cdat$name[i],":\t",str,sq,cdat$codon[i],"\n"))
    gtest[i] <- sq==cdat$codon[i]
}

cat(paste(sum( gtest), "correct codons;",
          sum(!gtest), "wrong codons\n"))



