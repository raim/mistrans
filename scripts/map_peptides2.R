
library(readxl)
library(segmenTools)
options(stringsAsFactors=FALSE)

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
fig.path <- file.path(proj.path,"figures")
out.path <- file.path(proj.path,"processedData")

##dat <- read_xlsx("All_MTP_BP_sharedPep_quant.xlsx")
##colnames(dat) <- gsub(" ","_", colnames(dat))

## NOTE: input has changed significantly, proceed in a map_peptides2.R

saap.file <- "All_MTP_BP_sharedPep_quant_03Oct23.xlsx"))
##saap.file <- "All_SAAP_protein_filter_df.xlsx"

dat <- read_xlsx(file.path(dat.path, saap.file))
colnames(dat) <- gsub(" ","_", colnames(dat))



## get ALL proteins - from project mammary
genes <- read.delim(file.path(mam.path,"features_GRCh38.110.tsv"))
##genes <- genes[genes$type=="protein_coding",]

## GET ENSEMBL PROTEINS - from project mammary
fasta.url <- "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
fasta <- file.path(mam.path,"originalData","Homo_sapiens.GRCh38.pep.all.fa.gz")
if ( !file.exists(fasta) )
        download.file(inurl, infile)
fas <- readFASTA(fasta)

desc <- lapply(fas, function(x) {unlist(strsplit(x$desc, " "))})
ids <- unlist(lapply(desc, function(x) sub("\\..*","",x[1])))
names(fas) <- ids

## only take proteins where Top_leading_protein is present in fasta
mids <- sub("_.*","",dat$Top_leading_protein)
dat <- dat[mids%in%ids,]

## FILTER unique mistranslated or base peptides
## SAAP reflect individual mutations at potentially different locations
dat <- dat[!duplicated(dat$SAAP_sequence),]

## refresh mids: remove mutation annotation
mids <- sub("_.*","", dat$Top_leading_protein)

## get mutations
notmutated <- grep("_", dat$Top_leading_protein, invert=TRUE)
muts <- sub(".*_","", dat$Top_leading_protein)
muts[notmutated] <- ""
muts <- strsplit(muts, ",")

## grep each base peptide (dat$BP) in the corresponding fasta
## brute force, each peptide
pos <- len <- rep(NA, nrow(dat))
mpos <- list() ## main peptide positions
testit <- TRUE # TEST whether the mutation position is correct
for ( i in 1:nrow(dat) ) {
    oid <- dat$Top_leading_protein[i]
    j <- grep(mids[i], names(fas),value=FALSE, fixed=TRUE)
    if ( length(j)==0 ) {
        cat(paste(i, oid, "not found\n"))
        next
    }
    if ( length(j)>1 ) {
        cat(paste(i, oid, "more than one hit\n"))
        next
    }
    query <- unlist(dat[i,3])
    target <- fas[[j]]$seq

    ## introduce genomic mutations
    if ( length(muts[[i]])>0 ) {
        cat(paste(i,"replacing",
                  length(muts[[i]]), "mutations in", j, names(fas)[j], "\n"))
        ntar <- try(mutatePositions(target, muts[[i]]))
        if ( class(ntar)=="try-error" )
            cat(paste("error introducing mutations; using unmutated\n"))
        else target <- ntar
    }

    ## MAP PETPTIDE 
    res <- gregexpr(query, target)
    if ( res[[1]][1]==-1 ) {
        cat(paste(i, oid, 
                  "no match in", j, names(fas)[j],"\n"))
        next
    }
    if ( length(res)>1 ) {
        warning(paste(i, oid, 
                      "more than two hits in",j, names(fas)[j],
                      "; taking first\n"))
    }
    ## get position of mutation
    saap <- unlist(dat[i,2])
    mut <- which(strsplit(saap,"")[[1]]!=strsplit(query,"")[[1]])
    pos[i] <- res[[1]][1] + mut -1
    len[i] <- length(unlist(strsplit(target,"")))

    ## test location of mutation
    if ( testit ) {
        ntarg <- sub(query,saap,target)
        if ( strsplit(target,"")[[1]][pos[i]] ==
             strsplit(ntarg,"")[[1]][pos[i]] ) {
            stop(i, mids[i], j, names(fas)[j], "error: no mutation detected\n")
        }
    }   
}

mpos <- rep(list(NA), nrow(dat)) ## main peptide positions
miss <- mpos
for ( i in 1:nrow(dat) ) {

    oid <- dat$Top_leading_protein[i]
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
                  length(muts[[i]]), "mutations in", j, names(fas)[j], "\n"))
        ntar <- try(mutatePositions(target, muts[[i]]))
        if ( class(ntar)=="try-error" )
            cat(paste("error introducing mutations; using unmutated\n"))
        else target <- ntar
    }

    ## MAP ALL MAIN PETPTIDES
    queries <- strsplit(dat$Main_peptides[i],",")
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
            warning(i, oid, " no match in ", j, names(fas)[j])
            cat(paste0(oid, " no match for ", queries[k], "\n"))
            mss <- c(mss, k)
        } else apos[k] <- res[[1]][1]        
    }
    miss[[k]] <- mss
    mpos[[i]] <- apos/len[i]
}

## TODO
## 2012: Jonathan Weissman, Ignolia, Cell paper: embryonic stem cells,
## ribosome density at mutated position.
##  https://doi.org/10.1016/j.cell.2011.10.002
## Bacterial mistranslation v ribosome density.

## ribosome density: ribosome density per codon divided
## by the mean ribosome density of transcript.

## relative position
rpos <- pos/len


## bind to data frame
dat <- cbind(dat, proteinID=mids, pos=pos, len=len, rpos=pos/len)

## filter
rm <- which(pos<0)
if ( length(rm)>0 ) {

    cat(paste("removing", length(rm), "SAAP without exact match\n"))

    pos <- pos[-rm]
    len <- len[-rm]
    dat <- dat[-rm,]
}

### TODO:
## BACKGROUND: main peptides
## * background frequencies: map all "main peptides",
## * background frequencies: map all genomic mutations,
## STATISTICS: of missing hits (peptides not in protein),

### WRITE OUT TABLE with positions for downstream analysis
write.table(dat, file=file.path(out.path,"saap_mapped.tsv"),
            sep="\t", quote=FALSE, na="")


### HOTSPOTS
## look for proteins with more than one mutation
## in close vicinity
hotspots <- split(dat$pos, f=dat$proteinID)
hotspots <- lapply(hotspots, sort)

mutn <- unlist(lapply(hotspots, length))

idx <-  which.max(mutn)
#idx <- which(names(hotspots)=="ENSP00000506126")
id <- names(hotspots)[idx]
name <- genes$name[grep(id, genes$protein)]
if ( length(name)==0 ) name <- id
pos <- table(hotspots[[idx]])

png(file.path(fig.path,"hotspots_example.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(x=as.numeric(names(pos)), y=pos, type="h",
     main=name, xlab="position in protein",
     ylab="# SAAP")
axis(2)
dev.off()

#### PLOTS

### POSITION v LENGTH

size.cutoff <- 3e2
small <- len<size.cutoff
large <- len>size.cutoff




hist(pos, breaks=seq(0,4e4,100), xlim=c(0,3000), xlab="absolute position",
     main=NA)
hist(len, breaks=seq(0,4e4,100), border=2, add=TRUE)
legend("topright", c("BP position","protein length"),
       col=c(1,2), lty=1)


png(file.path(fig.path,"absolute_position_length.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i",
    xpd=TRUE)
dense2d(len, pos, xlim=c(0,1e3), ylim=c(0,1e3),
        ylab=expression(SAAP~position),
        xlab=expression(protein~length), nbin=512)
par(xpd=FALSE)
abline(a=0, b=1)
abline(a=0, b=.4)
abline(v=size.cutoff)
dev.off()

dense2d(log10(len), log10(pos), ylim=c(0,5), xlim=c(0,5),
        ylab=expression(SAAP~position/log[10]),
        xlab=expression(protein~length/log[10]))
abline(a=0,b=1)

png(file.path(fig.path,"relative_position_hist.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(rpos, breaks=seq(0,1,.1),
     main=paste(sum(!is.na(rpos)), "unique SAAP"),
     xlab="relative position of peptide in protein")
loc.sze <- table(rpos>.5)
dev.off()

## BACKGROUND: main peptides
png(file.path(fig.path,"relative_position_hist_main_peptides.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(unlist(mpos), breaks=seq(0,1,.1),
     main=paste(sum(!is.na(unlist(mpos))), "main peptides"),
     xlab="relative position of peptide in protein")
dev.off()

## SMALL v LARGE PROTEINS
png(file.path(fig.path,"relative_position_hist_small.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(rpos[small], breaks=seq(0,1,.1),
     main=paste0(sum(small), " small proteins <",size.cutoff," aa"),
     xlab="relative position of peptide in protein")
dev.off()
png(file.path(fig.path,"relative_position_hist_large.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(rpos[large], breaks=seq(0,1,.1),
     main=paste0(sum(large), " large proteins >",size.cutoff," aa"),
     xlab="relative position of peptide in protein")
dev.off()

png(file.path(fig.path,"protein_length.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(log10(len),
     xlab=expression(protein~length/log[10]), breaks=50)
abline(v=log10(size.cutoff))
dev.off()


png(file.path(fig.path,"relative_position_length.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(rpos, log10(len),
        ##len, log="y",
        xlab="relative position of peptide in protein",
        ylab=expression(protein~length/log[10]))
abline(h=log10(size.cutoff))
dev.off()

png(file.path(fig.path,"relative_position_intensity.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(dat$rpos, dat[,"SAAP_intensity"]/dat[,"BP_intensity"], log="y",
        xlab="relative position of peptide in protein",
        ylab="relative intensity SAAP/BP")
dev.off()

png(file.path(fig.path,"absolute_position_RAAS.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(dat$pos, dat[,"RAAS"], xlim=c(0,1e3),
        xlab="absolute position of peptide in protein/log10",
        ylab="RAAS")
dev.off()

png(file.path(fig.path,"relative_position_RAAS.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(dat$rpos, dat[,"RAAS"], 
        xlab="relative position of peptide in protein",
        ylab="RAAS")
dev.off()



