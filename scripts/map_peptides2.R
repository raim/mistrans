
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
fig.path <- file.path(proj.path,"figures","saap_analysis")
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

## s4pred
s4pred <- file.path(mam.path,"processedData","Homo_sapiens.GRCh38.pep.large_s4pred.fas.gz")
iupred <- file.path(mam.path,"processedData","iupred3")

## figures
dir.create(fig.path, showWarnings=FALSE)

##dat <- read_xlsx("All_MTP_BP_sharedPep_quant.xlsx")
##colnames(dat) <- gsub(" ","_", colnames(dat))

## NOTE: input has changed significantly, proceed in a map_peptides2.R

##saap.file <- "All_MTP_BP_sharedPep_quant_03Oct23.xlsx"))
saap.file <- "All_SAAP_protein_filter_df_20240111.xlsx"
PMATCH <- "Leading_razor_proteins__all_TMT_" # column with protein matches
PRAAS <- "Mean_precursor_RAAS" # precursor ion/experiment level
RRAAS <- "Mean_reporter_RAAS"  # reporter ion/patient level

dat <- data.frame(read_xlsx(file.path(dat.path, saap.file)))
colnames(dat) <- gsub("\\.","_", colnames(dat))

### FILTER
## TODO: carry over to proteins for function analysis
rm <- (dat$Potential_contaminant | dat$Immunoglobulin) |
    (dat$"K_R_AAS" | dat$Potential_uncaptured_transcript) |
    dat$"Hemoglobin_Albumin"
dat$remove <- rm

cat(paste("tagged", sum(rm),
          "potentially false positives; now:",nrow(dat),"\n"))

## FILTER unique mistranslated or base peptides
## and attach mean and CV RAAS,
## SAAP reflect individual mutations at potentially different locations

## first get list of RAAS values for duplicates
slst <- split(10^unlist(dat[,PRAAS]), f=dat$SAAP)
slst <- slst[dat$SAAP] # map to main data

slen <- unlist(lapply(slst, length))
smds <- unlist(lapply(slst, median, na.rm=TRUE))
smns <- unlist(lapply(slst, mean, na.rm=TRUE))
sdev <- unlist(lapply(slst, sd, na.rm=TRUE))
scvs <- sdev/smns

## tag duplicates
dups <- slen > 1
cat(paste("tagged", sum(dups), "duplicated SAAP\n"))

## attach mean and CV RAAS to duplicate-filtered list
dat <- cbind(dat, RAAS.median=log10(smds),
             RAAS.mean=log10(smns), RAAS.cv=scvs, RAAS.n=slen)

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
pids <- gsub("'","",
             sub("\\[","",
                 sub("\\]","", dat[,PMATCH])))
pids <- lapply(strsplit(pids, ","), trimws)
## remove all non ENSEMBL proteins
pids <- lapply(pids, function(x) x[grep("ENSP",x)])

## take FIRST of multiple
pids <- unlist(lapply(pids, function(x) x[1]))

## sublists: split and take first again
## TODO: handle sublists betters
pids <- strsplit(pids, ";")
pids <- unlist(lapply(pids, function(x) x[1]))

## remove mutation tag
eids <- sub("_.*", "", pids)

dat$protein <- pids
dat$ensembl <- eids

### TODO: mix approaches
##dat$protein <- pids  # first ensembl -mutation tags
##dat$ensembl <- bpid  # longest ensembl


#### FIND POSITIONS of SAAP in protein and transcript

## GET ENSEMBL PROTEINS - from project mammary
fas <- readFASTA(fasta, grepID=TRUE)
names(fas) <- sub("\\.[0-9]+", "", names(fas))



## only take proteins where Top_leading_protein is present in fasta
mids <- sub("_.*","",dat$protein)
keep <- mids%in%names(fas)

dat$no_protein <- !keep

cat(paste("tagged", sum(!keep),
          "proteins where we don't have a sequence; now:", nrow(dat), "\n"))

## refresh mids: remove mutation annotation
mids <- sub("_.*","", dat$protein)

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

## re-order transcripts for available proteins
missing <- which(!names(fas)%in%rownames(trmap))
## TODO: required or do we valuable loose protein info
cat(paste("removing", length(missing), "proteins w/o matching transcript\n"))
fas <- fas[-missing]
## order by protein names via trmap
trfas <- trfas[trmap[names(fas),1]]
names(trfas) <- names(fas)


## s4pred
s4p <- readFASTA(s4pred, grepID=TRUE)
## skip version number - TODO: use version numbers!
names(s4p) <- sub("\\.[0-9]+", "", names(s4p))
## reorder
cat(paste("s4pred:", sum(!names(fas)%in%names(s4p)), "proteins not found\n"))
s4p <- s4p[names(fas)]

## map iupred files here! 
iufiles <- list.files(pattern=paste0(".*iupred3.tsv"), path=iupred)
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
errors <- matrix(0, nrow=nrow(dat), ncol=7)
colnames(errors) <- c("noens", "noprt",
                      "mform", "merr", 
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
    gid <- mids[i] # original gene ID

    if ( is.na(oid) ) { # no mapped ensembl protein
        cat(paste("WARNING:",i, "no protein identified\n"))
        noens <- c(noens, paste(i, gid, muts[[i]]))
        errors[i,"noens"] <- 1
        next
    }
    
    j <- grep(mids[i], names(fas),value=FALSE, fixed=TRUE)
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
            stop(i, mids[i], j, names(fas)[j], "error: no mutation detected\n")
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
    ## GENETIC_CODE[codon]
    tcodons <- paste(names(which(GENETIC_CODE%in%aat[i])),collapse=";")

    ## TODO: get secondary structure, disordered probability
    
    ## report
    tmp <- ifelse(length(muts[[i]])==0, "", paste(muts[[i]], collapse=";"))
    if ( !interactive() )
        cat(paste("DONE:", i, codon, aaf[i], GENETIC_CODE[codon], "->", aat[i],
                  tcodons, tmp, "\n"))

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

## FILTER DATA
## TODO: find sources!
## TODO: don't remove but carry over as background

## TODO: move the following interactive exploration to saap_analysis.R
if ( interactive() ) {
    ## explore s4pred results
    ss.tot <- apply(sssbg,2, sum,na.rm=TRUE)
    ss.tot <- ss.tot/sum(ss.tot)
    ss.aas <- table(sss)
    ss.aas <- ss.aas/sum(ss.aas)
    ss.tab <- rbind(AAS=ss.aas, total=ss.tot)
    ## slight enrichment of beta/alpha
    barplot(ss.tab, beside=TRUE, legend=TRUE)

    ## TODO: why negative iup?



    ## slight positive trend of unstructured/anchor vs RAAS
    plotCor(iup, iubg) # positive correlation to background

    hist(iup, border=NA, col=NA)
    hist(iubg, col=1, border=1, add=TRUE)
    hist(iup, border=2, col=paste0(rgb(t(col2rgb(2))/255),77), add=TRUE)
    legend("topright", paste0("p=",signif(wilcox.test(iup, iubg)$p.value,2)))
    
    layout(t(1:2), widths=c(1,.25))
    par(mai=c(.5,.5,.1,.1),yaxs="i")
    plotCor(dat$Mean_precursor_RAAS, iup, na.rm=TRUE)
    iubgh <- hist(iubg, breaks=0:10/10, plot=FALSE)
    par(mai=c(.5,0,.1,.2),yaxs="i")
    barplot(iubgh$counts,names.arg=iubgh$mids, horiz=TRUE, las=2, space=0)

    plotCor(anc, anbg) # positive correlation to background

    hist(anc, border=NA, col=NA)
    hist(anbg, col=1, border=1, add=TRUE)
    hist(anc, border=2, col=paste0(rgb(t(col2rgb(2))/255),77), add=TRUE)
    legend("topright", paste0("p=",signif(wilcox.test(anc, anbg)$p.value,2)))

    layout(t(1:2), widths=c(1,.25))
    par(mai=c(.5,.5,.1,.1))
    plotCor(dat$Mean_precursor_RAAS, anc, na.rm=TRUE)
    anbgh <- hist(anbg, breaks=0:10/10, plot=FALSE)
    par(mai=c(.5,0,.1,.2),yaxs="i")
    barplot(anbgh$counts,names.arg=anbgh$mids, horiz=TRUE, las=2, space=0)
    axis(4)

    plotCor(dat$Mean_precursor_RAAS, anc/anbg, na.rm=TRUE)
    plotCor(dat$Mean_precursor_RAAS, iup/iubg, na.rm=TRUE)

 }

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


### WRITE OUT TABLE with positions for downstream analysis
sdat <- dat[,c(colnames(dat)[1:5],grep("RAAS",colnames(dat),value=TRUE),
               "ensembl","protein","len","pos","rpos","from","to",
               "codon","s4pred",
               "iupred3","iupred3.protein",
               "anchor2","anchor2.protein",
               "remove")]
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
    j <- grep(mids[i], names(fas),value=FALSE, fixed=TRUE)
    
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
    queries <- strsplit(dat$Peptides_associated_with_leading_razor_protein[i],
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
