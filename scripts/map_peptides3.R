
## MAP AAS TO ENSEMBL PROTEINS & TRANSCRIPTS
## and collect sequence & structure data at AAS.

## TODO: also load full list of unique SAAP/BP and do
## statistics on blast, AA and codon matching.

library(segmenTools)
library(Biostrings)
options(stringsAsFactors=FALSE)

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"

dat.path <- file.path(proj.path,"originalData")
fig.path <- file.path(proj.path,"figures","saap_mapping3")
out.path <- file.path(proj.path,"processedData")

## SAAP/BP pairs
saapf <- file.path(out.path,"unique_saap.tsv")

## BP:protein blast results
bpmap <- file.path(out.path,"bp_mapped.tsv")


## protein fasta
fasta <- file.path(out.path,"all_proteins.fa")

## coding region fasta
transcr <- file.path(mam.path,"processedData","coding.fa")

## protein-transcript map
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
## transcript structure
exmap <- file.path(mam.path,"originalData","transcript_coordinates.tsv")

## structure predictions
## s4pred
s4pred <- file.path(mam.path,"processedData",
                    "Homo_sapiens.GRCh38.pep.large_s4pred.fas.gz")
## iupred3/anchor2
iupred <- file.path(mam.path,"processedData","iupred3")

## figures
dir.create(fig.path, showWarnings=FALSE)



#### FIND POSITIONS of SAAP in protein and transcript

## GET ENSEMBL PROTEINS - from project mammary
fas <- readFASTA(fasta, grepID=TRUE)

## get matching transcripts
trfas <- readFASTA(transcr, grepID=TRUE)

## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
## reverse map transcript-protein
pamrt <- matrix(rownames(trmap), ncol=1)
rownames(pamrt) <- trmap[,1]

## rename by protein names via trmap, for quick access
## of protein-specific transcript
names(trfas) <- pamrt[names(trfas),1]

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


### READ IN SAAP/BP and blast results
dat <- read.delim(saapf, header=FALSE)
colnames(dat) <- c("SAAP","BP")
bmap <- read.delim(bpmap)
dat <- merge(dat, bmap, by="BP", all=TRUE)


## grep each base peptide (dat$BP) in the corresponding fasta
## brute force, each peptide
pos <- len <- cdn <- aaf <- aat <-  aas <-
    sss <- anc <- iup <- iubg <- anbg <- rep(NA, nrow(dat))

## secondary structure frequencies in whole protein
sssbg <- matrix(NA, nrow=nrow(dat), ncol=3)
colnames(sssbg) <- c("C","E","H")

## genomic coordinate of AAS
gcoor <- matrix(NA, nrow=nrow(dat), ncol=3)
colnames(gcoor) <- c("chr","coor","strand")

testit <- TRUE # TEST whether the mutation position is correct
use.regex <- FALSE # use TRUE to test blast results vs. direct regex

## count errors
nomatch <- noprt <-
    noseq <- wcodon <- wsss <- wiup <- character()
errors <- matrix(0, nrow=nrow(dat), ncol=6)
colnames(errors) <- c("noprt",
                      "nomatch","nomut",
                      "noseq", "wcodon",
                      "wcoor")
for ( i in 1:nrow(dat) ) {
   
    oid <- dat$protein[i] # ID with mutation index
    gid <- dat$ensembl[i] # original gene ID

    j <- which(names(fas)==oid)
    if ( length(j)==0 ) { # 
        cat(paste("WARNING:",i, oid, "no protein sequence found\n"))
        noprt <- c(noprt, paste(i, gid))
        errors[i,"noprt"] <- 1
        next
    }
    if ( length(j)>1 ) { # atm not appearing
        stop(paste("WARNING:",i, oid, "more than one hit\n"))
        next
    }
    query <- unlist(dat[i,"BP"])
    target <- fas[[j]]$seq

    ## GET POSITION OF MUTATION
    ## 1: MAP PETPTIDE - TODO: cross-check or replace by blast result
    if ( use.regex ) {
        res <- gregexpr(query, target)
        if ( res[[1]][1]==-1 ) {
            cat(paste("WARNING:",i, oid, "no match in", j, names(fas)[j],"\n"))
            nomatch <- c(nomatch, paste(i, gid, query))
            errors[i,"nomatch"] <- 1
            next
        }
        if ( length(res)>1 ) 
            stop(paste("PROBLEM:", i, oid, 
                       "more than two hits in",j, names(fas)[j],
                       "; taking first\n"))
        if ( res[[1]][1] != dat[i,"sstart"] ) 
            ## if this never occurs we can use blast
            stop(paste("PROBLEM:", i, oid, 
                       "regex doesn't match blast\n"))
        AAS <- res[1][1]
    } else AAS <- dat[i,"sstart"]

    ## 2: AAS within protein
    saap <- unlist(dat[i,"SAAP"])
    mut <- which(strsplit(saap,"")[[1]]!=strsplit(query,"")[[1]])
    pos[i] <- AAS + mut -1
    len[i] <- length(unlist(strsplit(target,"")))

    ## test mutation
    if ( testit ) { # 
        ntarg <- sub(query,saap,target)
        if ( strsplit(target,"")[[1]][pos[i]] ==
             strsplit(ntarg,"")[[1]][pos[i]] ) {
            cat(paste("WARNING:", i, gid, j, names(fas)[j],
                       "no mutation detected\n"))
            nomatch <- c(nomatch, paste(i, gid, query))
            errors[i,"nomut"] <- 1
            ##next            
        }
    }

    ## test consistency with blast

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
    } else {    
        cdn[i] <- codon
    }
    
    ## get genome coor-based data and add genomic position of AAS

    gmap <- trexo[trexo$protein==gid,]
    chr <- unique(gmap[,"chr"])
    strand <- unique(gmap[,"strand"])

    ## order CDS/cons
    gmap <- gmap[order(gmap[,"start"], decreasing=strand!="-"),]

    ## CDS length
    rlen <- gmap[,"end"]-gmap[,"start"]+1
    
    cds <- c(1,cumsum(rlen))
    
    if ( max(cds)/3-1 != len[i] ) {
        cat(paste("WARNING:",i,"wrong CDS sum length", max(cds)/3-1,
                  "vs", len[i],"in", gid, "on strand",strand,"\n"))
        errors[i,"wcoor"] <- 1
        ##stop(i, gid, ": wrong length of CDS map")
    } else {
        cds <- tail(which(cds<=npos), 1) # which exon in gmap?
        ## genomic coordinate of AAS; 2nd codon position!
        dst <- npos+1
        coor <- gmap[cds, ifelse(strand=="-", "end","start")] +
            ifelse(strand=="-", -1, 1)*dst
        gcoor[i,] <- c(chr, coor, strand)
    }
    
    ## TODO:
    ## load chromosome index chrS
    ## load phastcons genomic data, convert with coor2index
    ## do this externally, write fasta like file, but with numbers
    ## for riboseq, phastcons, etc.
 
     
    ## report
    if ( !interactive() ) {        
        tcodons <- paste(names(which(GENETIC_CODE%in%aat[i])),collapse=";")
        cat(paste("DONE:", i, codon, aaf[i], GENETIC_CODE[codon], "->", aat[i],
                  tcodons,  "\n"))
    }
    
    if ( aas[i] != aaf[i] ) 
        warning(paste("WARNING: BP pos from AA:", aaf[i], "vs",
                  aas[i], GENETIC_CODE[codon], gid, "\n"))
}


if ( !interactive() ) {
    cat(paste("PEPTIDE NOT MATCHED,",length(nomatch),":",
              paste(nomatch, collapse=" ; "),"\n"))
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
             gcoor,
             s4pred=sss, sssbg, iupred3=iup, iupred3.protein=iubg,
             anchor2=anc, anchor2.protein=anbg)

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

## TODO: bar plot with errors


### WRITE OUT TABLE with positions for downstream analysis
write.table(dat, file=file.path(out.path,"saap_mapped3.tsv"),
            sep="\t", quote=FALSE, na="", row.names=FALSE)


