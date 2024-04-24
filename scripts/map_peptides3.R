
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
bpmap <- file.path(out.path,"bp_mapped2.tsv")


## protein fasta
pfasta <- file.path(out.path,"all_proteins.fa")

## coding region fasta
tfasta <- file.path(mam.path,"processedData","coding.fa")

## protein-transcript map
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
## CDS structure
cdsmap <- file.path(mam.path,"processedData","protein_cds_structure.dat")
cdspos <- file.path(mam.path,"processedData","protein_cds_coordinates.tsv")


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
fas <- readFASTA(pfasta, grepID=TRUE)

## get matching transcripts
trfas <- readFASTA(tfasta, grepID=TRUE)

## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
## reverse map transcript-protein
pamrt <- matrix(rownames(trmap), ncol=1)
rownames(pamrt) <- trmap[,1]

## rename by protein names via trmap, for quick access
## of protein-specific transcript
names(trfas) <- pamrt[names(trfas),1]


## list of splice sites for each CDS
cds <- read.delim(cdsmap,header=FALSE, row.names=1)
cdl <- strsplit(cds[,1],";")
names(cdl) <- rownames(cds)
cdl <- lapply(cdl, as.numeric)
## genomic position of CDS
cpos <- read.delim(cdspos, header=FALSE, sep=" ")
colnames(cpos) <- c("ID","chr","start","end","strand")
posl <- split(cpos[,2:5], cpos[,1])
## NOTE: reverse order for minus strand
posl <- lapply(posl, function(x) x[order(x$start,
                                         decreasing=x$strand=="-"),])

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
mut <- pos <- len <- cdn <- aaf <- aat <-  aas <-
    sss <- anc <- iup <- iubg <- anbg <- rep(NA, nrow(dat))

## secondary structure frequencies in whole protein
sssbg <- matrix(NA, nrow=nrow(dat), ncol=3)
colnames(sssbg) <- c("C","E","H")

## genomic coordinate of AAS
gcoor <- matrix(NA, nrow=nrow(dat), ncol=3)
colnames(gcoor) <- c("chr","coor","strand")

testit <- TRUE # TEST whether the mutation position is correct
use.regex <- TRUE #FALSE # use TRUE to test blast results vs. direct regex

## count errors
wiup <- character()
errors <- matrix(0, nrow=nrow(dat), ncol=8)
colnames(errors) <- c("no protein",
                      "using blast","no BP",
                      "wrong s4pred len","wrong iupred3 len",
                      "no transcript", "wrong codon",
                      "AAS > CDS len")
if ( !use.regex ) errors[,"using blast"] <- 1


## trouble shooting 20240424
pid="ENSP00000354876" # no transcript even though its in protein_transcript_map
## it's not in coding.fa, why??


for ( i in 1:nrow(dat) ) {
   
    oid <- dat$protein[i] # ID with mutation index
    gid <- dat$ensembl[i] # original gene ID

    j <- which(names(fas)==oid)
    if ( length(j)==0 ) { # 
        cat(paste("WARNING:",i, oid, "no protein sequence found\n"))
        errors[i,"no protein"] <- 1
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
        if ( length(res)>1 ) {
            stop(paste("PROBLEM:", i, oid, 
                       "more than two hits in",j, names(fas)[j], "\n"))
        } else if ( res[[1]][1]==-1 ) {
            cat(paste("WARNING:",i, oid, "no match in", j, names(fas)[j],
                      ", using blast result.\n"))
            errors[i,"using blast"] <- 1
            AAS <- dat[i,"sstart"]
        } else if ( res[[1]][1] != dat[i,"sstart"] ) {
            ## if this never occurs we can use blast results!
            stop(paste("PROBLEM:", i, oid, 
                       "regex doesn't match blast\n"))
        } else  AAS <- res[[1]][1]
    } else AAS <- dat[i,"sstart"]

    ## 2: AAS within protein
    saap <- unlist(dat[i,"SAAP"])
    mut[i] <- which(strsplit(saap,"")[[1]]!=strsplit(query,"")[[1]])
    pos[i] <- AAS + mut[i] -1
    len[i] <- length(unlist(strsplit(target,"")))

    ## test mutation
    if ( testit ) { #
        if ( length(grep(query,target))==0 ) {
            cat(paste("WARNING:", i, gid, j, names(fas)[j],
                      "BP is not in target protein\n"))
            errors[i,"no BP"] <- 1
        } else {
            ntarg <- sub(query,saap,target)
            if ( strsplit(target,"")[[1]][pos[i]] ==
                 strsplit(ntarg,"")[[1]][pos[i]] ) {
                stop(paste("WARNING:", i, gid, j, names(fas)[j],
                           "no mutation detected\n"))
                
            }
        }
    }

    ## test consistency with blast

    ## GET s4pred
    s4s <- s4p[[gid]]
    if ( is.null(s4s) ) {
        cat(paste("WARNING:", i, "no s4pred found for", oid, "\n"))
        errors[,"wrong s4pred len"] <- "no s4pred result"
    } else {
        ## test length
        if ( nchar(s4s$seq)!=len[i] ) {
            cat(paste("WARNING:", i, "s4 prediction has wrong length for",
                      oid, "\n"))
            errors[,"wrong s4pred len"] <- 1
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
        errors[i,"wrong iupred3 len"] <- "no iupred3 result"
    } else {
        iud <- read.delim(file.path(iupred,iufile), header=FALSE)
        if ( nrow(iud)!=len[i] ) {
            cat(paste("WARNING:", i, "iupred3 prediction has wrong length for",
                      oid, "iupred:",nrow(iud), "vs. protein:", len[i],"\n"))
            errors[i,"wrong iupred3 len"] <- 1
        } else {
            anc[i] <- iud[pos[i], 4] #anchor2
            iup[i] <- iud[pos[i], 3] #iupred3
            ## whole protein mean
            anbg[i] <- mean(iud[, 4])
            iubg[i] <- mean(iud[, 3])
        }

        ## NOTE: COPY TO SELECTED FOLDER FOR TRANSFER TO LAPTOP
        tof <- file.path(paste0(iupred,"_selected"),iufile)
        if ( Sys.info() ["nodename"]=="exon" ) 
            file.copy(from=file.path(iupred,iufile), to=tof, overwrite = FALSE)
    }
    
    ## GET CODON
    ## protein:transcript mapping
    ## load transcript, load UTR file, get codon
    nt <- trfas[[gid]]
    if ( is.null(nt) ) {
        cat(paste("WARNING:", i, "no transcript found for", oid, "\n"))
        errors[i,"no transcript"] <- 1
        next
    }

    ## position of codon
    npos <- (pos[i]-1)*3+1

    codon <- substr(nt$seq, npos, npos+2)

    ## store codon and from/to AA
    ##codon <- codons(DNAString(nt$seq))[pos[i]]

    ## don't store wrong codon!
    aas[i] <- strsplit(target,"")[[1]][pos[i]]
    aaf[i] <- strsplit(query,"")[[1]][mut[i]]
    aat[i] <- strsplit(saap,"")[[1]][mut[i]]
    if ( GENETIC_CODE[codon] != aaf[i] ) {
        cat(paste("WARNING:",i,"wrong codon", aaf[i],
                  "vs", codon, GENETIC_CODE[codon],"in", gid, "\n"))
        errors[i,"wrong codon"] <- 1
    } else {    
        cdn[i] <- codon
    }
    
    ## get genomic position via CDS list and posl
    cds <- cdl[[gid]]
    coors <- posl[[gid]]

    exon <- which(cds>=npos)[1]
    if ( is.na(exon) ) {
        cat(paste("WARNING:",i,"wrong genomic position", npos, "vs", max(cds),
                  "in", gid, "\n"))
         errors[i,"AAS > CDS len"] <- 1
    } else {
    
        chr <- unique(coors$chr)
        strand <- unique(coors$strand)
        ## get genomic position of AAS 2nd codon
        rpos <- npos
        if ( exon>1 ) # subtract prior exons/CDS
            rpos <- npos-cds[exon-1]
        if ( strand=="-") {
            coor <- coors[exon, "end"] - rpos 
        } else {
            coor <- coors[exon, "start"] + rpos  
        }
        gcoor[i,] <-
            c(chr, coor, strand)
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





## relative position
rpos <- pos/len
## bind to data frame
bgsss <- sssbg
colnames(bgsss) <- paste0(colnames(bgsss), ".protein")
dat <- cbind(dat,
             site=mut,  # mutated site within peptide
             pos=pos,   # position of AAS within protein
             len=len,
             rpos=pos/len,
             from=aaf, to=aat,
             codon=cdn, gcoor, 
             s4pred=sss, bgsss, iupred3=iup, iupred3.protein=iubg,
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

## remove where no codon was assigned
## no transcript found or wrong codon
rm <- which(is.na(dat$chr))
if ( length(rm)>0 ) {
    ##dat <- dat[-rm,]
    cat(paste("WARNING:", length(rm), "SAAP w/o genome coordinates\n"))
}

### ERROR AND RECOVERY STATS

errs <- errors#[,apply(errors,2,sum)>0]
er <- c(SAAP=nrow(dat),
        table(dat$match)[c("good","bad","wrong")],
        "with mutation"=sum(dat[,"protein"]!=dat[,"ensembl"],na.rm=TRUE),
        apply(errs,2,sum))
png(file.path(fig.path,"mapping_errors.png"),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(1.5,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25)
bp <- barplot(er,las=2)
text(bp, er, er, pos=3, xpd=TRUE, cex=.8)
dev.off()

## TODO: better characterize matches to mutated proteins
recovery <- c("SAAP"=nrow(dat),
              "protein"=sum(!is.na(dat[,"protein"])),
              ##"length"=sum(!is.na(dat[,"len"])),
              "s4pred"=sum(!is.na(dat[,"s4pred"])),
              "iupred3"=sum(!is.na(dat[,"iupred3"])),
              "genome"=sum(!is.na(dat[,"chr"])),
              "codon"=sum(!is.na(dat[,"codon"])))
png(file.path(fig.path,"mapping_recovery.png"),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(.75,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25)
bp <- barplot(recovery,las=2,ylab="")
text(bp, recovery, recovery, pos=3, xpd=TRUE, cex=.8)
dev.off()

types <- c("SAAP"=nrow(dat),
           "MANE"=sum(!is.na(dat[,"MANE.protein"])),
           "extracellular"=sum(dat[,"extracellular"], na.rm=TRUE),
           "IG"=sum(dat[,"IG"], na.rm=TRUE),
           "albumin"=sum(dat[,"albumin"], na.rm=TRUE),
           "globin"=sum(dat[,"globin"], na.rm=TRUE))

types <- sort(types, decreasing=TRUE)
png(file.path(fig.path,"mapping_types.png"),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(.75,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25, cex=.8)
bp <- barplot(types,las=2,ylab="")
text(bp, types, types, pos=3, xpd=TRUE)
dev.off()

## inspect some overlaps
table(dat[,"extracellular"], dat[,"globin"])
table(dat[,"extracellular"], dat[,"albumin"])
table(dat[,"extracellular"], dat[,"IG"])


### WRITE OUT TABLE with positions for downstream analysis
if ( !interactive() ) {
    write.table(dat, file=file.path(out.path,"saap_mapped3.tsv"),
                sep="\t", quote=FALSE, na="", row.names=FALSE)
}

