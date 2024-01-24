
## Q for shiri:
## * handle sublists with ; in Leading razor,
## * ENSP00000496574_I90L: no I at position 90,
## * ENSP00000426880_R28 : missing target AA
## * 

### TODO:
## BACKGROUND: main peptides
## * background frequencies: map all "main peptides",
## * background frequencies: map all genomic mutations,
## STATISTICS: of missing hits (peptides not in protein),

library(readxl)
library(segmenTools)
options(stringsAsFactors=FALSE)

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"

gen.path <- file.path(mam.path, "originalData")
dat.path <- file.path(proj.path,"originalData")
fig.path <- file.path(proj.path,"figures","saap_analysis")
out.path <- file.path(proj.path,"processedData")

## protein fasta
fasta <- file.path(gen.path,"Homo_sapiens.GRCh38.pep.all.fa.gz")
## coding region fasta
transcr <- file.path(mam.path,"processedData","coding.fa")
## protein-transcript map
tpmap <- file.path(gen.path,"protein_transcript_map.tsv")

## 
dir.create(fig.path, showWarnings=FALSE)

##dat <- read_xlsx("All_MTP_BP_sharedPep_quant.xlsx")
##colnames(dat) <- gsub(" ","_", colnames(dat))

## NOTE: input has changed significantly, proceed in a map_peptides2.R

##saap.file <- "All_MTP_BP_sharedPep_quant_03Oct23.xlsx"))
saap.file <- "All_SAAP_protein_filter_df_20240111.xlsx"

dat <- read_xlsx(file.path(dat.path, saap.file))
colnames(dat) <- gsub(" ","_", colnames(dat))

### FILTER
rm <-
    (dat$Potential_contaminant | dat$Immunoglobulin) |
    (dat$"Hemoglobin/Albumin" | dat$Potential_uncaptured_transcript)
dat <- dat[!rm,]

cat(paste("removed", sum(rm),
          "potentially false positives; now:",nrow(dat),"\n"))

## FILTER unique mistranslated or base peptides
## SAAP reflect individual mutations at potentially different locations
dups <- duplicated(dat$SAAP)
cat(paste("removing", sum(dups), "duplicated SAAP\n"))
dat <- dat[!dups,]

cat(paste("removed", sum(dups), "duplicated SAAP; now:",nrow(dat),"\n"))

## get mapped proteins
pids <- gsub("'","",
             sub("\\[","",
                 sub("\\]","", dat$"Leading_razor_proteins_(all_TMT)")))
pids <- lapply(strsplit(pids, ","), trimws)
## remove all non ENSEMBL proteins
pids <- lapply(pids, function(x) x[grep("ENSP",x)])

## take FIRST of multiple
pids <- unlist(lapply(pids, function(x) x[1]))

## sublists: split and take first again
## TODO: handle sublists betters
pids <- strsplit(pids, ";")
pids <- unlist(lapply(pids, function(x) x[1]))

dat$protein <- pids




#### FIND POSITIONS of SAAP in protein and transcript


## GET ENSEMBL PROTEINS - from project mammary
fas <- readFASTA(fasta)

desc <- lapply(fas, function(x) {unlist(strsplit(x$desc, " "))})
ids <- unlist(lapply(desc, function(x) sub("\\..*","",x[1])))
names(fas) <- ids

## only take proteins where Top_leading_protein is present in fasta
mids <- sub("_.*","",dat$protein)
keep <- mids%in%ids

dat <- dat[keep,]

cat(paste("removed", sum(!keep),
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
cat(paste("removing", length(missing), "proteins w/o matching transcript\n"))
fas <- fas[-missing]

trfas <- trfas[trmap[names(fas),1]]
names(trfas) <- names(fas)



## grep each base peptide (dat$BP) in the corresponding fasta
## brute force, each peptide
pos <- len <- cdn <- aaf <- aat <-  aas <- rep(NA, nrow(dat))
mpos <- list() ## main peptide positions
testit <- TRUE # TEST whether the mutation position is correct
## count errors
mform <- merr <- nomatch <- noseq <- wcodon <- character()
errors <- matrix(NA, nrow=nrow(dat), ncol=5)
colnames(errors) <- c("mform", "merr", "nomatch", "noseq", "wcodon")
for ( i in 1:nrow(dat) ) {

    oid <- dat$protein[i] # ID with mutation index
    gid <- mids[i] # original gene ID
    
    j <- grep(mids[i], names(fas),value=FALSE, fixed=TRUE)
    if ( length(j)==0 ) { # atm not appearing
        stop(paste("WARNING:",i, oid, "not found\n"))
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

    ## MAP PETPTIDE 
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
    ## get position of mutation
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
    
    ## GET CODON
    ## protein:transcript mapping
    ## load transcript, load UTR file, get codon
    nt <- trfas[[gid]]
    if ( is.null(nt) ) {
        cat(paste("WARNING:", i, "no sequence found for", oid, "\n"))
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

    
    ## report
    tmp <- ifelse(length(muts[[i]])==0, "", paste(muts[[i]], collapse=";"))
    if ( !interactive() )
        cat(paste(i, codon, aaf[i], GENETIC_CODE[codon], "->", aat[i],
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
            warning(paste(i, oid, "no match in", j, names(fas)[j]))
            cat(paste0(oid, " no match for ", queries[k], "\n"))
            mss <- c(mss, k)
        } else apos[k] <- res[[1]][1]        
    }
    miss[[i]] <- mss
    mpos[[i]] <- apos/len[i]
}


## relative position
rpos <- pos/len
## bind to data frame
dat <- cbind(dat, pos=pos, len=len, rpos=pos/len, from=aaf, to=aat, codon=cdn)

apply(errors, 2, sum, na.rm=TRUE)

## FILTER DATA
## TODO: find sources!

## missing position - no match found from peptide in protein!!
## TODO: perhaps test LONGEST of the leading razor proteins.
rm <- which(is.na(dat$pos)) #dat$pos<0)
if ( length(rm)>0 ) {
    dat <- dat[-rm,]
    errors <- errors[-rm,]
    cat(paste("removed", length(rm), "SAAP without protein match, now:",
              nrow(dat), "\n"))
}

## remove where length or position are NA - same as "no match"!
rm <- which(is.na(dat$len))
if ( length(rm)>0 ) {
    dat <- dat[-rm,]
    errors <- errors[-rm,]
    cat(paste("removed", length(rm), "SAAP no protein length, now:",
              nrow(dat), "\n"))
}

## remove where no codon was assigned
## no transcript found or wrong codon
rm <- which(is.na(dat$codon))
if ( length(rm)>0 ) {
    dat <- dat[-rm,]
    cat(paste("removed", length(rm), "SAAP w/o codon, now:",
              nrow(dat), "\n"))
}


### WRITE OUT TABLE with positions for downstream analysis
sdat <- dat[,c(colnames(dat)[1:5],grep("RAAS",colnames(dat),value=TRUE),
               "protein","len","pos","rpos","from","to","codon")]
write.table(sdat, file=file.path(out.path,"saap_mapped.tsv"),
            sep="\t", quote=FALSE, na="", row.names=FALSE)

## also save background distribution of main peptide positions
save(mpos, file=file.path(out.path, "mapped_peptides.rda"))
