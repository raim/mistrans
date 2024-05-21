
### ANALYZE SEQUENCE CONTEXT around AAS

## 1: all -> just do hypergeo enrichment tests,
##    as in previous get_aas_context.R
## 2: analyze position of AAS in peptide.

## Motivated by above analysis, select subsets of BP/sequences,
## and do a DiffLogo and a RAAS profile.

## 1: M or W in -2/-1/+1/+2
## 2: higher than threshold E/D
## 3: AAS in peptide position 1-3.

## ... and finally loop over all of certain type:
## 1: from AA,
## 2: to AA,
## 3: from/to AA,
## 4: from/to by AA prop.

library(DiffLogo)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
source("~/work/mistrans/scripts/raasprofiles3_init.R")

mfig.path <- file.path(fig.path,"motifs")
dir.create(mfig.path, showWarnings=FALSE)

## FILTER UNIQUE BP - since those have the same AA CONTEXT
## NOTE: that this looses different SAAPs for the same BP

bpraas <- split(tmtf$RAAS, tmtf$BP)
bpraas <- listProfile(bpraas, y=tmtf$RAAS, use.test=use.test, min=3)

bdat <- hdat[!duplicated(hdat$BP),]
rownames(bdat) <- bdat$BP

## match rows

## MISSING?
## manual inspection shows this BP had infinite RAAS
missing <- which(!rownames(bdat)%in%rownames(bpraas))
if ( length(missing) )
    cat(paste("WARNING:", length(missing), "BP missing from TMT level file\n"))

## keep only bdat for which we have RAAS values
bdat <- bdat[rownames(bpraas),]
bdat <- cbind(bdat, bpraas)

### POSITION OF AAS IN PEPTIDES

## Categorize by N+i and C-i
mx <- 20
nt <- bdat$site
ct <- nchar(bdat$BP)-bdat$site +1
## fuse central AAS
nt[nt>mx] <- ct[ct>mx] <- mx+1
## take smaller
at <- paste0("N",nt)
at[ct<nt] <- paste0("C",ct[ct<nt])
## re-name central AAS
at[at==paste0("N",mx+1)] <- paste0(">",mx)
at.srt <- c(
    paste0("N",1:mx),
    paste0(">",mx),
    paste0("C",mx:1))

asite.bins <- at
asite.srt <- at.srt

cls <- clusterCluster(bdat$from, asite.bins, cl2.srt=asite.srt,
                      alternative="two.sided")
##cls <- sortOverlaps(cls, p.min=.1)
plotdev(file.path(mfig.path,paste0("peptide_position_encoded")),
        height=5, width=mx/2+1, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=p.min, p.txt=p.txt, ylab="Encoded AA at AAS", xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

cls <- clusterCluster(bdat$to, asite.bins, cl2.srt=asite.srt,
                      alternative="two.sided")
plotdev(file.path(mfig.path,paste0("peptide_position_incorporated")),
        height=5, width=mx/2+1, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=p.min, p.txt=p.txt, ylab="Incorporated AA at AAS",
             xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()

## get
?createDiffLogoObject

### ALL SEQUENCES HYPERGEO
aam <- do.call(rbind, strsplit(bdat$AA,""))
nc <- (ncol(aam)-1)/2
colnames(aam) <- -nc:nc

## HYPERGEO TESTS - tigther context
ovl <- aaProfile(aam, abc=AAT)


plotOverlaps(ovl, p.min=p.min, p.txt=p.txt)

dotprofile(ovl, value="ratio", vcols=ttcols, xlab=NA,
           p.dot=p.min, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), ylab="amino acid", axis=1:2, show.total=TRUE)

## DEFINE CLASSES

kr <- apply(aam[,as.character(-3:-1)], 1,
            function(x) any(x%in%c("K","R")))
met <- apply(aam[,as.character(c(-2,-1,1,2))], 1,
             function(x) sum(x%in%c("M")))
barplot(table(met))
met <- met>0
trp <- apply(aam[,as.character(c(-2,-1,1,2))], 1,
             function(x) any(x%in%c("W")))
acd <-  apply(aam[,as.character(c(-6:-1,1:6))], 1,
              function(x) sum(x%in%c("E","D")))
barplot(table(acd))
acd <- acd>3

boxplot(bdat$median ~ met)

## CREATE POSITION WEIGHT MATRICES
## for certain sequence classes

filt <- met # trp # acd #bdat$site==2 #kr

## generate frequencies at each position
getPFM <- function(aa, alphabet=ASN$chars) {
    ## column-wise table of AA
    ctl <- apply(aa, 2, table)
    ## aaids
    aaids <- unique(unlist(lapply(ctl, names)))
    ctl <- do.call(cbind, lapply(ctl, function(x) x[aaids]))
    rownames(ctl) <- aaids
    ctl[is.na(ctl)] <- 0

    ## remove all not in alphabet before taking col sum
    ctm <- matrix(0, nrow=length(alphabet), ncol=ncol(aa))
    rownames(ctm) <- alphabet
    ctm[alphabet[alphabet%in%rownames(ctl)],] <-
        ctl[alphabet[alphabet%in%rownames(ctl)],]
    ctm[is.na(ctm)] <- 0

    ## frequencies
    ctf <- t(t(ctm)/apply(ctm,2,sum))
    ctf
}

cols <- as.character(c(-10:10))
pfm1 <- getPFM(aam[ filt,cols,drop=FALSE])
pfm2 <- getPFM(aam[!filt,cols,drop=FALSE])

## TODO: adapt colors
ASN2 <- ASN
ASN2$cols <- aa.cols[ASN2$chars]

dfob <- createDiffLogoObject(pfm2, pfm1, alphabet=ASN2)


diffLogo(dfob)




