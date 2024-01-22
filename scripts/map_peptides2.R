
library(readxl)
library(segmenTools)
library(Biostrings) # for blosum62
data(BLOSUM62)
data(PAM250)
library(vioplot)
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

cat(paste("removing", sum(rm),
          "proteins that are potentially false positives\n"))
dat <- dat[!rm,]

## get mapped proteins
pids <- gsub("'","",
             sub("\\[","",
                 sub("\\]","", dat$"Leading_razor_proteins_(all_TMT)")))
pids <- lapply(strsplit(pids, ","), trimws)
## remove all non ENSEMBL proteins
pids <- lapply(pids, function(x) x[grep("ENSP",x)])
##hist(unlist(lapply(pids, length)))

## take FIRST of multiple
pids <- unlist(lapply(pids, function(x) x[1]))
dat$protein <- pids

## FILTER unique mistranslated or base peptides
## SAAP reflect individual mutations at potentially different locations
##cat(paste("removing", sum(duplicated(dat$SAAP)), "duplicated SAAP\n"))
##dat <- dat[!duplicated(dat$SAAP),]


#### REPLACED AA SIMILARITIES

## first find mutated AA pairs
## (column AASin input mixes I/L)
## and split mutated AA pairs into from/to

saaps <- strsplit(dat$SAAP,"")
bases <- strsplit(dat$BP, "")
fromto <- lapply(1:length(saaps), function(i) {
    pos <- which(saaps[[i]]!=bases[[i]])
    c(from=bases[[i]][pos], to=saaps[[i]][pos])
})


## analyze AA similarity by different matrices
## TODO: remove extra columns
simmats <- c("BLOSUM62", "PAM250")

## reporter vs. precursor RAAS

png(file.path(fig.path,paste0("raas_precursor_reporter.png")),
    res=300, width=3.5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(dat$Mean_reporter_RAAS, dat$Mean_precursor_RAAS, na.rm=TRUE,
        xlab="mean reporter RAAS", ylab="mean precursor RAAS")
abline(a=0, b=1, col=1, lty=2, lwd=2)
dev.off()

for ( i in seq_along(simmats) ) {

    mid <- simmats[i]
    MAT <- get(mid)[AA_STANDARD,AA_STANDARD]

    ylab <- "similarity between replaced AA" #"BLOSUM62 similarity"

    ## similarity of replaced AA
    sim <- unlist(lapply(fromto, function(x) MAT[x[1],x[2]]))

    raasr <- dat$Mean_reporter_RAAS # 9965 NA
    raasp <- dat$Mean_precursor_RAAS # 582 NA
    
    rid <- "Mean_reporter_RAAS"
    raas <- unlist(dat[,rid])
    
    raas_bins <- cut(raas, c(seq(-6,4,1)))
    raas_dual <- rep("RAAS>0", length(raas))
    raas_dual[raas < 0] <- "RAAS<0"
    
    df <- data.frame(raas=raas,
                     bins=raas_bins,
                     dual=raas_dual,
                     sim=sim)
    df <- df[!is.na(raas),]# & !is.infinite(raas),]
    dff <- df[!is.infinite(df$raas),]
    
    
    
    ## correlation
    ##ct <- cor.test(dff$raas, dff$sim)
    ct <- cor.test(dff$raas, dff$sim, method="spearman")
    fit <- lm(dff$sim ~ dff$raas)
    ##fit2 <- tls::tls(sim ~ raas, data=df)
    
    png(file.path(fig.path,paste0(mid,"_raas_dense.png")),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(dff$raas, dff$sim,
            ylab=ylab, xlab=rid)
    ##abline(fit)
    ##legend("topright", paste0("p=", signif(ct$p.value,1)),
    ##       seg.len=0, y.intersp=0, bty="n")
    mtext(paste0("p=", signif(ct$p.value,0)), 3,0, adj=1)
    dev.off()
    png(file.path(fig.path,paste0(mid,"_raas_violin.png")),
        res=300, width=5, height=3.5, units="in")
    par(mai=c(.75,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
    vioplot(sim ~ bins, data=df, ylab=ylab,
            xlab=NA, las=2)
    mtext(rid, 1, 2.5)
    dev.off()
    png(file.path(fig.path,paste0(mid,"_raas_boxplot.png")),
        res=300, width=5, height=3.5, units="in")
    par(mai=c(.75,.5,.5,.1), mgp=c(1.3,.3,0), tcl=-.25)
    boxplot(sim ~ bins, data=df, ylab=ylab,
            xlab=NA, las=2, at=seq_along(levels(df$bins)))
    axis(3, at=seq_along(levels(df$bins)), labels=table(df$bins),las=2)
    mtext(rid, 1, 2.5)
    dev.off()
    

    png(file.path(fig.path,paste0(mid,"_raas_boxplot_dual.png")),
        res=300, width=3.5, height=3.5, units="in")
    b62 <- cbind.data.frame(MAT[lower.tri(MAT)],
                            mid)
    ddff <- data.frame(sim=c(dff[,"sim"], b62[,1]),
                       dual=c(dff[,"dual"],b62[,2]))
    ddff[,2] <- as.factor(ddff[,2])
    ylim <- range(ddff[,1])
    ylim[2]  <- ylim[2]*1.1
    par(mai=c(.5,.5,.25,0), mgp=c(1.3,.3,0), tcl=-.25)
    boxplot(sim ~ dual, data=ddff, xlab=NA,
            ylab=ylab, axes=FALSE, at=1:3, box=FALSE, ylim=ylim)
    axis(2)
    axis(1, at=1:3, labels=paste0(levels(ddff[,2]), "\n", table(ddff[,2])),
         mgp=c(0,1.5,0))
    ## TODO: significance bars
    tt <- wilcox.test(ddff$sim[ddff$dual=="RAAS<0"],
                      ddff$sim[ddff$dual=="RAAS>0"])
    text(2.5, par("usr")[4]*.9, paste("p =",signif(tt$p.value,1)),
         xpd=TRUE, pos=3)
    arrows(x0=2, y0=par("usr")[4]*.9, x1=3, angle=90,
           code=3, length=.05, xpd=TRUE)
    tt <- wilcox.test(ddff$sim[ddff$dual==mid],
                      ddff$sim[ddff$dual=="RAAS<0"])
    text(1.5, par("usr")[4], paste("p =",signif(tt$p.value,1)), xpd=TRUE, pos=3)
    arrows(x0=1, y0=par("usr")[4], x1=2, angle=90, code=3, length=.05, xpd=TRUE)
    dev.off()
    
    png(file.path(fig.path,paste0(mid,"_similarities_hist.png")),
        res=300, width=5, height=3.5, units="in")
    hist(MAT[lower.tri(MAT)])
    dev.off()
}


#### FIND POSITIONS of SAAP in protein and transcript


## GET ENSEMBL PROTEINS - from project mammary
fas <- readFASTA(fasta)

desc <- lapply(fas, function(x) {unlist(strsplit(x$desc, " "))})
ids <- unlist(lapply(desc, function(x) sub("\\..*","",x[1])))
names(fas) <- ids

## only take proteins where Top_leading_protein is present in fasta
mids <- sub("_.*","",dat$protein)
keep <- mids%in%ids
cat(paste("removing", sum(!keep), "proteins where we don't have a sequence\n"))

dat <- dat[keep,]


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
pos <- len <- rep(NA, nrow(dat))
mpos <- list() ## main peptide positions
testit <- TRUE # TEST whether the mutation position is correct
for ( i in 1:nrow(dat) ) {
    oid <- dat$protein[i]
    gid <- mids[i] # original gene ID
    
    j <- grep(mids[i], names(fas),value=FALSE, fixed=TRUE)
    if ( length(j)==0 ) {
        cat(paste(i, oid, "not found\n"))
        next
    }
    if ( length(j)>1 ) {
        cat(paste(i, oid, "more than one hit\n"))
        next
    }
    query <- unlist(dat[i,"BP"])
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
    saap <- unlist(dat[i,"SAAP"])
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
    
    ## GET CODON
    ## protein:transcript mapping
    ## load transcript, load UTR file, get codon
    nt <- trfas[[gid]]
    if ( is.null(nt) ) {
        cat(paste("no sequence found for", oid, "\n"))
        next
    }

    ## position of codon
    npos <- (pos[i]-1)*3+1

    codon <- substr(nt$seq, npos, npos+2)
    ##codon <- codons(DNAString(nt$seq))[pos[i]]
    aaf <- strsplit(target,"")[[1]][pos[i]]
    aaf <- strsplit(query,"")[[1]][mut]
    aat <- strsplit(saap,"")[[1]][mut]
    ## GENETIC_CODE[codon]
    tcodons <- paste(names(which(GENETIC_CODE==aat)),collapse=";")

    ## report
    
    tmp <- ifelse(length(muts[[i]])==0, "", paste(muts[[i]], collapse=";"))
    cat(paste(i, codon, aaf, GENETIC_CODE[codon], "->", aat,
              tcodons, tmp, "\n"))
}

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
                  length(muts[[i]]), "mutations in", j, names(fas)[j], "\n"))
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
dat <- cbind(dat, pos=pos, len=len, rpos=pos/len)

## filter
rm <- which(pos<0)
if ( length(rm)>0 ) {

    cat(paste("removing", length(rm), "SAAP without exact match\n"))

    dat <- dat[-rm,]
}

## remove where NA
## TODO: find source
nas <- which(is.na(len) | is.na(rpos))

if ( length(nas)>0 ) {

    cat(paste("removing", length(nas), "SAAP with NA in length or position\n"))

    dat <- dat[-nas,]
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
## get ALL proteins - from project mammary
genes <- read.delim(file.path(mam.path,"features_GRCh38.110.tsv"))
##genes <- genes[genes$type=="protein_coding",]

## look for proteins with more than one mutation
## in close vicinity
hotspots <- split(dat$pos, f=dat$protein)
hotspots <- lapply(hotspots, sort)

mutn <- unlist(lapply(hotspots, length))

idx <-  which.max(mutn)
#idx <- which(names(hotspots)=="ENSP00000506126")
id <- names(hotspots)[idx]
name <- genes$name[grep(id, genes$protein)]
if ( length(name)==0 ) name <- id
hpos <- table(hotspots[[idx]])

png(file.path(fig.path,"hotspots_example.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(x=as.numeric(names(hpos)), y=hpos, type="h",
     main=name, xlab="position in protein",
     ylab="# SAAP")
axis(2)
dev.off()

#### PLOTS

### POSITION v LENGTH



## size cut off
size.cutoff <- 3e2
small <- dat$len<size.cutoff
large <- dat$len>size.cutoff



hist(dat$pos, breaks=seq(0,4e4,100), xlim=c(0,3000), xlab="absolute position",
     main=NA)
hist(dat$len, breaks=seq(0,4e4,100), border=2, add=TRUE)
legend("topright", c("BP position","protein length"), col=c(1,2), lty=1)


png(file.path(fig.path,"absolute_position_length.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i", yaxs="i",
    xpd=TRUE)
dense2d(dat$len, dat$pos, xlim=c(0,1e3), ylim=c(0,1e3),
        ylab=expression(SAAP~position),
        xlab=expression(protein~length), nbin=512)
par(xpd=FALSE)
abline(a=0, b=1)
abline(a=0, b=.5)
abline(v=size.cutoff)
dev.off()

dense2d(log10(dat$len), log10(dat$pos), ylim=c(0,5), xlim=c(0,5),
        ylab=expression(SAAP~position/log[10]),
        xlab=expression(protein~length/log[10]))
abline(a=0,b=1)

png(file.path(fig.path,"relative_position_hist.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(dat$rpos, breaks=seq(0,1,.1),
     main=paste(sum(!is.na(dat$rpos)), "unique SAAP"),
     xlab="relative position of peptide in protein")
loc.sze <- table(dat$rpos>.5)
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
hist(dat$rpos[small], breaks=seq(0,1,.1),
     main=paste0(sum(small,na.rm=TRUE), " small proteins <",size.cutoff," aa"),
     xlab="relative position of peptide in protein")
dev.off()
png(file.path(fig.path,"relative_position_hist_large.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(dat$rpos[large], breaks=seq(0,1,.1),
     main=paste0(sum(large,na.rm=TRUE), " large proteins >",size.cutoff," aa"),
     xlab="relative position of peptide in protein")
dev.off()

png(file.path(fig.path,"protein_length.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(log10(dat$len),
     xlab=expression(protein~length/log[10]), breaks=50)
abline(v=log10(size.cutoff))
dev.off()


png(file.path(fig.path,"relative_position_length.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(dat$rpos, log10(dat$len),
        ##len, log="y",
        xlab="relative position of peptide in protein",
        ylab=expression(protein~length/log[10]))
abline(h=log10(size.cutoff))
dev.off()

raas.col <- "Mean_precursor_RAAS"
png(file.path(fig.path,"absolute_position_RAAS.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(dat$pos, dat[,raas.col], xlim=c(0,1e3),
        ylim=range(dat[,raas.col], na.rm=TRUE),
        xlab="absolute position of peptide in protein/log10",
        ylab=raas.col)
dev.off()

png(file.path(fig.path,"relative_position_RAAS.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(dat$rpos, dat[,raas.col],
        xlim=range(dat$rpos, na.rm=TRUE),
        ylim=range(dat[,raas.col], na.rm=TRUE),
        xlab="relative position of peptide in protein",
        ylab=raas.col)
dev.off()



