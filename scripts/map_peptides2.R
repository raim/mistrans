
## TODO 20240216
## * unify mapped proteins, eg. AMGIMNSFVNDIFER maps to
##   five different histone H2B genes (different for each dataset).

## Q for shiri:
## * handle sublists with ; in Leading razor,
## * ENSP00000496574_I90L: no I at position 90,
## * ENSP00000426880_R28 : missing target AA
## * 

### TODO:
## * mean of mean of means: better use median?
## * take mean/median RAAS of duplicate SAAP before removing them.

## BACKGROUND: main peptides
## * background frequencies: map all "main peptides",
## * background frequencies: map all genomic mutations,
## STATISTICS: of missing hits (peptides not in protein),
## BACKGROUND: codons
## * get codons from all BP, all "main peptides",
##   all mapped proteins or ALL proteins.

library(readxl)
library(segmenTools)
library(Biostrings)
options(stringsAsFactors=FALSE)

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"

gen.path <- file.path(mam.path, "originalData")
dat.path <- file.path(proj.path,"originalData")
fig.path <- file.path(proj.path,"figures","saap_mapping")
out.path <- file.path(proj.path,"processedData")

## global human feature file
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")
## protein fasta
fasta <- file.path(gen.path,"Homo_sapiens.GRCh38.pep.all.fa.gz")
## coding region fasta
transcr <- file.path(mam.path,"processedData","coding.fa")#"Homo_sapiens.GRCh38.pep.all_cds.fa.gz")
## protein length file
ptlen <- file.path(gen.path,"Homo_sapiens.GRCh38.pep.all_lengths.tsv")
## protein-transcript map
tpmap <- file.path(gen.path,"protein_transcript_map.tsv")
## transcript structure
exmap <- file.path(gen.path,"transcript_coordinates.tsv")

## BP re-blast results
bpmap <- file.path(out.path,"bp_mapped.tsv")

## s4pred
s4pred <- file.path(mam.path,"processedData",
                    "Homo_sapiens.GRCh38.pep.large_s4pred.fas.gz")
iupred <- file.path(mam.path,"processedData","iupred3")

## figures
dir.create(fig.path, showWarnings=FALSE)

##dat <- read_xlsx("All_MTP_BP_sharedPep_quant.xlsx")
##colnames(dat) <- gsub(" ","_", colnames(dat))

## NOTE: input has changed significantly, proceed in a map_peptides2.R

##saap.file <- "All_MTP_BP_sharedPep_quant_03Oct23.xlsx"))
saap.file <- file.path(dat.path, "All_SAAP_protein_filter_df.txt")
tmt.file <- file.path(dat.path, "All_SAAP_TMTlevel_quant_df.txt")
pat.file <- file.path(dat.path, "All_SAAP_patient_level_quant_df.txt")

## COLUMN NAMES
PMATCH <- "Leading.razor.proteins..all.TMT." # column with protein matches
PRAAS <- "Mean.precursor.RAAS" # precursor ion/experiment level
RRAAS <- "Mean.reporter.RAAS"  # reporter ion/patient level

dat <- read.delim(saap.file) #data.frame(read_xlsx(saap.file))
##colnames(dat) <- gsub("\\.","_", colnames(dat))

### FILTER
filtercols <- c("Potential.contaminant","Immunoglobulin",
                "K.R.AAS","Potential.uncaptured.transcript","Trypsin")
##filtercols <- c(filtercols, "Hemoglobin.Albumin")
logic.cols <- c("Hemoglobin.Albumin","Keep.SAAP")
## convert to logic
for ( f in c(filtercols,logic.cols) ) {
    dat[,f] <- as.logical(dat[,f])
    dat[is.na(dat[,f]),f] <- FALSE
}

## TODO: carry over to proteins for function analysis
rm <- apply(dat[,filtercols], 1, any) 
dat$remove <- rm

cat(paste("tagged", sum(rm),
          "potentially false positives; now:",nrow(dat),"\n"))


## get list of mean RAAS values for unique SAAP
## NOTE: keeping this for now, but in fact, we
## are calculating a mean/median RAAS at each
## level of analysis.
slst <- split(10^unlist(dat[,PRAAS]), f=dat$SAAP)
slst <- slst[dat$SAAP] # map to main data
slen <- unlist(lapply(slst, length))
smds <- unlist(lapply(slst, median, na.rm=TRUE))
smns <- unlist(lapply(slst, mean, na.rm=TRUE))
smns.nrm <- unlist(lapply(slst, function(x) mean(log10(x),na.rm=TRUE)))
sdev <- unlist(lapply(slst, sd, na.rm=TRUE))
scvs <- sdev/smns

## log of means vs. mean of logs
png(file.path(fig.path,"raas_means_duplicate_saap_delogged.png"),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(smns.nrm[slen>1], log10(smns[slen>1]),
        ylab=expression(log[10](mean(RAAS))),
        xlab=expression(mean(log[10]*RAAS)))
abline(a=0, b=1)
dev.off()

## duplicate SAAP/datatype
if ( FALSE ) {
    idx <- which(duplicated(paste(dat$SAAP, dat$Dataset, dat$BP)))
    head(dat[dat$SAAP==dat$SAAP[idx[1]],1:10])
}

## tag duplicates
dups <- slen > 1
cat(paste("tagged", sum(dups), "duplicated SAAP\n"))

## attach mean and CV RAAS for each unique SAAP
dat <- cbind(dat, RAAS.median=log10(smds),
             RAAS.mean=log10(smns), RAAS.cv=scvs, RAAS.n=slen)
## TODO: calculate RAAS means for unique SAAP for raw data!

### RESOLVE MAPPED PROTEINS

## get as much as possible longest ensembl proteins

if ( FALSE ) {

### THIS APPROACH TO GET THE BEST ENSEMBL HIT
### LOOSES MORE THAN IT GAINS, since we loose some
### base peptides that only match when accounting for
### the mutation suffixes in ensembl protein matches

### TODO: carry over resolved columns of all razor proteins
### to data file for later consideration.

    
    ## first split by ;
    plst <- strsplit(dat[,PMATCH],  ";")
    ## clean up strings
    plst <- lapply(plst, function(x) gsub("\\['","",gsub("\\']","",x)))
    ## split by ', '
    plst <- lapply(plst, function(x) unlist(strsplit(x, "', '")))
    
    ## SPLIT BY MAPPING PIPELINE/PROTEIN ANNOTATION
    unip.lst <- lapply(plst, function(x)
        unique(sub("-1","",sub("^CON__","",x[grep("^CON",x)]))))
    strg.lst <- lapply(plst, function(x) x[grep("^STRG",x)])
    ensp.lst <- lapply(plst, function(x) x[grep("^ENSP",x)])
    
    ## get unique ensembl proteins w/o mutation tags
    ensu.lst <- lapply(ensp.lst, function(x) unique(sub("_.*","",x)))
    
    
    ## find uniprot proteins in features and replace by ensembl proteins
    unips <- unique(unlist(unip.lst[unlist(lapply(ensu.lst, length))==0]))
    features <- read.delim(feature.file)


    ## SELECT LARGEST OF ENSEBML
    plen <- read.delim(ptlen, header=FALSE, row.names=1)
    rownames(plen) <- sub("\\.[0-9]*", "", rownames(plen))
    
    ensl.lst <- lapply(ensu.lst, function(x) {
        y <- x[which.max(plen[x,1])]
        if (length(y)==0) y <- x
        y
    })
    ensl.lst[unlist(lapply(ensl.lst,length))==0] <- list("")
    bpid <- unlist(ensl.lst) ## LONGEST ENSEMBL PROTEIN
    bpid[bpid==""] <- NA
}


## old simple approach of just taking the first ensembl
pids <- dat[,PMATCH]
pids <- gsub("\\[","", gsub("\\]","", pids))
pids <- lapply(strsplit(pids, "', '"), trimws)
## remove leading and trailing '
pids <- lapply(pids, function(x) gsub("'","", x))
## split other type of list present in some cases
pids <- lapply(pids, function(x) unlist(strsplit(x,";")))

## remove all non ENSEMBL proteins
pids <- lapply(pids, function(x) x[grep("ENSP",x)])

## take FIRST of multiple
pids <- lapply(pids, function(x) x[1])

## fill empty
pids[unlist(lapply(pids, length))==0] <- list(NA)

pids <- unlist(pids)

## remove mutation tag
eids <- sub("_.*", "", pids)

dat$protein <- pids
dat$ensembl <- eids

## LOAD RE-BLAST
bpm <- read.delim(bpmap)

bp2dat <- match(dat$BP, bpm$BP)
bpm <- bpm[bp2dat, c("ensembl", "sstart","send",
                     "gene","MANE","match",
                     "IG","albumin","globin","extracellular","exclude")]
colnames(bpm) <- paste0("reblast.", colnames(bpm))
dat <- cbind(dat, bpm)


## TODO: UNIFY ENSEMBL FOR IDENTICAL BP! 
## eg. collect all identical SAAP (before taking the first),
## look up if one is MANE and take this. Inherit the currently
## inconsistent tagging in TMT and Protein level files.

### TODO: mix approaches
##dat$protein <- pids  # first ensembl -mutation tags
##dat$ensembl <- bpid  # longest ensembl


#### FIND POSITIONS of SAAP in protein and transcript

## GET ENSEMBL PROTEINS - from project mammary
fas <- readFASTA(fasta, grepID=TRUE)
names(fas) <- sub("\\.[0-9]+", "", names(fas))



## only take proteins where Top_leading_protein is present in fasta
keep <- sub("_.*","",dat$protein)%in%names(fas)

dat$no_protein <- !keep

cat(paste("tagged", sum(!keep),
          "proteins where we don't have a sequence; now:", nrow(dat), "\n"))


## get mutations
notmutated <- grep("_", dat$protein, invert=TRUE)
muts <- sub(".*_","", dat$protein)

## nucleotide level mutations
nmuts <- grep(">", muts)
indels <- rep("",nrow(dat))
indels[nmuts] <- sub(".*_","", dat$protein[nmuts])

## protein level mutations
muts[notmutated] <- ""
muts[nmuts] <- ""
muts <- strsplit(muts, ",")


## get matching transcripts
trfas <- readFASTA(transcr)
desc <- lapply(trfas, function(x) {unlist(strsplit(x$desc, " "))})
trids <- unlist(lapply(desc, function(x) sub("\\..*","",x[1])))
trids <- sub(",$","", trids)
names(trfas) <- trids

## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)

## reverse map transcript-protein
pamrt <- matrix(rownames(trmap), ncol=1)
rownames(pamrt) <- trmap[,1]

## re-order transcripts for available proteins
missing <- which(!names(fas)%in%rownames(trmap))
## TODO: required or do we valuable loose protein info
cat(paste("removing", length(missing), "proteins w/o matching transcript\n"))
fas <- fas[-missing]

## rename by protein names via trmap, for quick access
## of protein-specific transcript
trfas <- trfas[trmap[names(fas),1]]
names(trfas) <- names(fas)
## remember transcript ID
trids <- trmap[names(fas),1]
names(trids) <-
    rownames(trmap[names(fas),,drop=FALSE])

## transcript exons genomic coordinates
trexo <- read.delim(exmap)
## only keep those where we have a protein
trexo <- trexo[trexo$transcript%in%trmap[names(fas),1],]
trexo$protein <- pamrt[trexo$transcript,1]

## s4pred
s4p <- readFASTA(s4pred, grepID=TRUE)
## skip version number - TODO: use version numbers!
names(s4p) <- sub("\\.[0-9]+", "", names(s4p))
## reorder
cat(paste("s4pred:", sum(!names(fas)%in%names(s4p)), "proteins not found\n"))
s4p <- s4p[names(fas)]

## map iupred files here! 
iufiles <- list.files(pattern=paste0(".*iupred3.tsv.gz"), path=iupred)
names(iufiles) <- sub("\\..*","", iufiles)
iufiles <- iufiles[names(fas)]

cat(paste("iupred3", sum(is.na(iufiles)), "proteins not found\n"))

## grep each base peptide (dat$BP) in the corresponding fasta
## brute force, each peptide
pos <- len <- cdn <- aaf <- aat <-  aas <-
    sss <- anc <- iup <- iubg <- anbg <- rep(NA, nrow(dat))
sssbg <- matrix(NA, nrow=nrow(dat), ncol=3)
colnames(sssbg) <- c("C","E","H")
mpos <- list() ## main peptide positions
testit <- TRUE # TEST whether the mutation position is correct

## count errors
mform <- merr <- nomatch <- noens <- noprt <-
    noseq <- wcodon <- wsss <- wiup <- character()
errors <- matrix(0, nrow=nrow(dat), ncol=8)
colnames(errors) <- c("noens", "noprt",
                      "mform", "merr", "mreblast",
                      "nomatch",
                      "noseq", "wcodon")
for ( i in 1:nrow(dat) ) {

    ## skip duplicated that have been calculated before
    if ( dat$RAAS.n[i]>1 ) {
        sidx <- which(dat$SAAP==dat$SAAP[i])
        if ( any(sidx<i) ) {
            cat(paste0("skipping duplicated SAAP: ", i , 
                       ", previous:", sidx[1], "\n"))
            next
        }
    }
    
    oid <- dat$protein[i] # ID with mutation index
    gid <- dat$ensembl[i] # original gene ID

    if ( is.na(oid) ) { # no mapped ensembl protein, check re-blast
        oid <- gid <- dat$reblast.ensembl[i]
        cat(paste("WARNING:",i, "using re-blast ensembl hit\n"))
        errors[i,"mreblast"] <- 1
    }
    if ( is.na(oid) ) { # no mapped ensembl protein
        cat(paste("WARNING:",i, "no protein identified\n"))
        noens <- c(noens, paste(i, gid, muts[[i]]))
        errors[i,"noens"] <- 1
        next
    }
    
    j <- grep(gid, names(fas),value=FALSE, fixed=TRUE)
    if ( length(j)==0 ) { # atm not appearing
        cat(paste("WARNING:",i, oid, "no protein sequence found\n"))
        noprt <- c(noprt, paste(i, gid, muts[[i]]))
        errors[i,"noprt"] <- 1
        next
    }
    if ( length(j)>1 ) { # atm not appearing
        stop(paste("WARNING:",i, oid, "more than one hit\n"))
        next
    }
    query <- unlist(dat[i,"BP"])
    target <- fas[[j]]$seq

    ## introduce genomic mutations
    if ( length(muts[[i]])>0 ) {
        if ( !interactive() )
            cat(paste(i,"replacing",
                      length(muts[[i]]), "mutations in",
                      j, names(fas)[j], "\n"))
        ## sanitize mutation
        ft <- do.call(cbind,strsplit(muts[[i]], "[0-9]+"))
        if ( nrow(ft)==1 ) {
            cat(paste("PROBLEM:", i ,"wrong mutation format",muts[[i]],
                      "using unmutated\n"))
            mform <- c(mform, paste(i, gid, muts[[i]]))
            errors[i,"mform"] <- 1
        } else {
            ntar <- try(mutatePositions(target, muts[[i]]))
            if ( class(ntar)=="try-error" ) {
                cat(paste("PROBLEM:", i ,"introducing mutations failed;",
                          "using unmutated\n"))
                merr <- c(merr, paste(i, gid, muts[[i]]))
                errors[i,"merr"] <- 1
            } else target <- ntar
        } 
    }

    ## GET POSITION OF MUTATION
    ## 1: MAP PETPTIDE 
    res <- gregexpr(query, target)
    if ( res[[1]][1]==-1 ) {
        cat(paste("WARNING:",i, oid, "no match in", j, names(fas)[j],"\n"))
        nomatch <- c(nomatch, paste(i, gid, query))
        errors[i,"nomatch"] <- 1
        next
    }
    if ( length(res)>1 ) {
        stop(paste("PROBLEM:", i, oid, 
                  "more than two hits in",j, names(fas)[j],
                  "; taking first\n"))
    }

    ## 2: AAS within peptide
    saap <- unlist(dat[i,"SAAP"])
    mut <- which(strsplit(saap,"")[[1]]!=strsplit(query,"")[[1]])
    pos[i] <- res[[1]][1] + mut -1
    len[i] <- length(unlist(strsplit(target,"")))

    ## test location of mutation
    if ( testit ) { # atm not happening!
        ntarg <- sub(query,saap,target)
        if ( strsplit(target,"")[[1]][pos[i]] ==
             strsplit(ntarg,"")[[1]][pos[i]] ) {
            stop(i, gid, j, names(fas)[j], "error: no mutation detected\n")
        }
    }

    ## GET s4pred
    s4s <- s4p[[gid]]
    if ( is.null(s4s) ) {
        cat(paste("WARNING:", i, "no s4pred found for", oid, "\n"))
        wsss <- c(wsss, paste(i, gid, "not found"))
    } else {
        ## test length
        if ( nchar(s4s$seq)!=len[i] ) {
            cat(paste("WARNING:", i, "s4 prediction has wrong length for",
                      oid, "\n"))
            wsss <- c(wsss, paste(i, gid, "wrong length"))
        } else {
            sss[i] <- substr(s4s$seq, pos[i], pos[i])
            ## background: complete protein sequence
            tb <- table(unlist(strsplit(s4s$seq,"")))
            sssbg[i,names(tb)] <- tb
        }
    }

    ## GET IUPRED3
    ##TODO: remove version tag from iupred files to avoid searching!
    ## or search once above
    iufile <- iufiles[gid]
    if ( is.na(iufile) )  {
        cat(paste("WARNING:", i, "no iupred3 found for", oid, "\n"))
        wiup <- c(wiup, paste(i, gid, "not found"))
    } else {
        iud <- read.delim(file.path(iupred,iufile), header=FALSE)
        if ( nrow(iud)!=len[i] ) {
            cat(paste("WARNING:", i, "iupred3 prediction has wrong length for",
                      oid, "\n"))
            wiup <- c(wiup, paste(i, gid, "wrong length"))
        } else {
            anc[i] <- iud[pos[i], 4]
            iup[i] <- iud[pos[i], 3]
            ## whole protein mean
            anbg[i] <- mean(iud[, 4])
            iubg[i] <- mean(iud[, 3])
        }
    }
    
    ## GET CODON
    ## protein:transcript mapping
    ## load transcript, load UTR file, get codon
    nt <- trfas[[gid]]
    if ( is.null(nt) ) {
        cat(paste("WARNING:", i, "no transcript found for", oid, "\n"))
        noseq <- c(noseq, paste(i, gid))
        errors[i,"noseq"] <- 1
        next
    }

    ## position of codon
    npos <- (pos[i]-1)*3+1

    codon <- substr(nt$seq, npos, npos+2)

    ## store codon and from/to AA
    ##codon <- codons(DNAString(nt$seq))[pos[i]]

    ## don't store wrong codon!
    aas[i] <- strsplit(target,"")[[1]][pos[i]]
    aaf[i] <- strsplit(query,"")[[1]][mut]
    aat[i] <- strsplit(saap,"")[[1]][mut]
    if ( GENETIC_CODE[codon] != aaf[i] ) {
        cat(paste("WARNING:",i,"wrong codon", aaf[i],
                  "vs", codon, GENETIC_CODE[codon],"in", gid, "\n"))
        wcodon <- c(wcodon, paste(i, gid, aaf[i],
                                  "vs", codon, GENETIC_CODE[codon]))
        errors[i,"wcodon"] <- 1
        next
    }
    
    cdn[i] <- codon

    ## get genome coor-based data

    gmap <- trexo[trexo$protein==gid,]

    ## TODO:
    ## load chromosome index chrS
    ## load phastcons genomic data, convert with coor2index
    ## do this externally, write fasta like file, but with numbers
    ## for riboseq, phastcons, etc.
 
     
    ## report
    if ( !interactive() ) {
        tmp <- ifelse(length(muts[[i]])==0, "", paste(muts[[i]], collapse=";"))
        tcodons <- paste(names(which(GENETIC_CODE%in%aat[i])),collapse=";")
        cat(paste("DONE:", i, codon, aaf[i], GENETIC_CODE[codon], "->", aat[i],
                  tcodons, tmp, "\n"))
    }
    
    if ( aas[i] != aaf[i] ) # doesnt happen since search is exact
        stop(paste("WARNING: BP pos", aaf[i], "vs",
                  aas[i], GENETIC_CODE[codon], "\n"))
}


if ( !interactive() ) {
    cat(paste("PEPTIDE NOT MATCHED,",length(nomatch),":",
              paste(nomatch, collapse=" ; "),"\n"))
    cat(paste("WRONG MUTATION FORMAT,",length(mform),":",
              paste(mform, collapse=" ; "),"\n"))
    cat(paste("MUTATION NOT FOUND,",length(merr),":",
              paste(merr, collapse=" ; "),"\n"))
    cat(paste("NO TRANSCRIPT SEQUENCE,",length(noseq),":",
              paste(noseq, collapse=" ; "),"\n"))
    cat(paste("WRONG CODON,",length(wcodon),":",
              paste(wcodon, collapse=" ; "),"\n"))
}



## relative position
rpos <- pos/len
## bind to data frame
colnames(sssbg) <- paste0(colnames(sssbg), ".protein")
dat <- cbind(dat, pos=pos, len=len, rpos=pos/len, from=aaf, to=aat, codon=cdn,
             s4pred=sss, sssbg,
             iupred3=iup, iupred3.protein=iubg,
             anchor2=anc, anchor2.protein=anbg)

## PLOT ERRORS
## TODO: legend!
png(file.path(fig.path,"mapping_loss.png"),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(.75,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(apply(errors, 2, sum, na.rm=TRUE), las=2,
        main=paste(nrow(errors),"total SAAP"))
dev.off()

png(file.path(fig.path,"mapping_loss_filtered.png"),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(.75,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(apply(errors[!dat$remove,], 2, sum, na.rm=TRUE), las=2,
        main=paste(sum(!dat$remove),"included SAAP"))
dev.off()

## check: reporter vs. precursor RAAS

png(file.path(fig.path,paste0("raas_precursor_reporter.png")),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(dat$Mean.reporter.RAAS, dat$Mean.precursor.RAAS, na.rm=TRUE,
        xlab="mean reporter RAAS", ylab="mean precursor RAAS",
        colf=viridis::viridis)
abline(a=0, b=1, col=1, lty=2, lwd=2)
dev.off()

## FILTER DATA

## missing position - no match found from peptide in protein!!
## TODO: perhaps test LONGEST of the leading razor proteins.
rm <- which(is.na(dat$pos)) #dat$pos<0)
if ( length(rm)>0 ) {
    ##dat <- dat[-rm,]
    ##errors <- errors[-rm,]
    cat(paste("WARNING:", length(rm), "SAAP without protein match\n"))
}

## remove where length or position are NA - same as "no match"!
rm <- which(is.na(dat$len))
if ( length(rm)>0 ) {
    ##dat <- dat[-rm,]
    ##errors <- errors[-rm,]
    cat(paste("WARNING:", length(rm), "SAAP no protein length\n"))
}

## remove where no codon was assigned
## no transcript found or wrong codon
rm <- which(is.na(dat$codon))
if ( length(rm)>0 ) {
    ##dat <- dat[-rm,]
    cat(paste("WARNING:", length(rm), "SAAP w/o codon\n"))
}

dat$transcript <- trids[dat$ensembl]


### WRITE OUT TABLE with positions for downstream analysis
sdat <- dat[,c(colnames(dat)[1:5],grep("RAAS",colnames(dat),value=TRUE),
               "ensembl","protein","transcript",
               "len","pos","rpos","from","to",
               "codon","s4pred","C.protein","E.protein","H.protein",
               "iupred3","iupred3.protein",
               "anchor2","anchor2.protein",
               "remove","Hemoglobin.Albumin","Keep.SAAP")]
write.table(sdat, file=file.path(out.path,"saap_mapped.tsv"),
            sep="\t", quote=FALSE, na="", row.names=FALSE)



### FIND MAIN PEPTIDE POSITIONS
## TODO: move this to a separate file, takes long
## TODO: add codon etc. frequencies as a background set.

## FIND main peptide positions
mpos <- rep(list(NA), nrow(dat)) ## main peptide positions
miss <- mpos
for ( i in 1:nrow(dat) ) {

    oid <- dat$protein[i]
    gid <- dat$ensembl[i]
    j <- grep(gid, names(fas),value=FALSE, fixed=TRUE)
    
    if ( length(j)==0 ) {
        cat(paste(i, oid, "not found\n"))
        next
    }
    if ( length(j)>1 ) {
        cat(paste(i, oid, "more than one hit\n"))
        next
    }
    target <- fas[[j]]$seq

    ## introduce genomic mutations
    if ( length(muts[[i]])>0 ) {
        cat(paste(i,"replacing",
                  length(muts[[i]]), "mutations in", j,
                  names(fas)[j], oid, "\n"))
        ntar <- try(mutatePositions(target, muts[[i]]))
        if ( class(ntar)=="try-error" )
            cat(paste("error introducing mutations; using unmutated\n"))
        else target <- ntar
    }

    ## MAP ALL MAIN PETPTIDES
    queries <- strsplit(dat$Peptides.associated.with.leading.razor.protein[i],
                        ",")
    queries <- unlist(lapply(queries, function(x) gsub("[^[:alnum:]]","",x)))
    apos <- rep(NA,length(queries))
    mss <- c()
    for ( k in 1:length(queries) ) {
        res <- gregexpr(queries[k], target)
        if ( length(res)>1 ) {
            warning(paste(i, oid, k,  "more than two hits in",j, names(fas)[j],
                          "; taking first\n"))
        }
        if ( res[[1]][1]==-1 ) {
            ##warning(paste(i, oid, "no match in", j, names(fas)[j]))
            cat(paste("MAIN PEPTIDE WARNING:",
                      i, oid, "no match for", queries[k],
                      "in", j, names(fas)[j], "\n"))
            mss <- c(mss, k)
        } else apos[k] <- res[[1]][1]        
    }
    miss[[i]] <- mss
    mpos[[i]] <- apos/len[i]
}

## also save background distribution of main peptide positions
save(mpos, file=file.path(out.path, "mapped_peptides.rda"))
