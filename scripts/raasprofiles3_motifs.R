
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
library(Biostrings)
AAS <- sort(unique(GENETIC_CODE))
AAT <- AAS[AAS!="*"]
diAAT <- sort(paste(AAT, rep(AAT,each=length(AAT)),sep=""))

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
mx <- 10
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
plotdev(file.path(mfig.path,paste0("AAS_position_encoded")),
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
plotdev(file.path(mfig.path,paste0("AAS_position_incorporated")),
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

cls <- clusterCluster(bdat$fromto, asite.bins, cl2.srt=asite.srt,
                      alternative="two.sided")
cls <- sortOverlaps(cls, p.min=1e-3)
plotdev(file.path(mfig.path,paste0("AAS_position_AAStype")),
        height=.2*nrow(cls$p.value)+1, width=mx/2+1, res=300)
par(mai=c(.5,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=p.min, p.txt=p.txt, ylab="AAS type", xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()
## sort by only enriched: TODO: implement in sortOverlaps
tmp <- cls
tmp$p.value[tmp$sign<0] <- 1
tmp <- sortOverlaps(tmp, p.min=1e-3, cut=TRUE)
cls <- sortOverlaps(cls, srt=rownames(tmp$p.value))

plotdev(file.path(mfig.path,paste0("AAS_position_AAStype_cut")),
        height=.2*nrow(cls$p.value)+1, width=mx/2+1, res=300)
par(mai=c(.5,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05)
plotOverlaps(cls, p.min=p.min, p.txt=p.txt, ylab="AAS type", xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.8)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(cls$p.value), labels=colnames(cls$p.value), las=2)
axis(2, at=nrow(cls$p.value):1, labels=rownames(cls$p.value), las=2)
mtext("position of AAS in peptide", 1, 1.5)
figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()


### AA MATRIX
aam <- do.call(rbind, strsplit(bdat$AA,""))
nc <- (ncol(aam)-1)/2
colnames(aam) <- -nc:nc

## AA MATRIX WITH AAS
aas <- aam
aas[,"0"] <- bdat$to

## HYPERGEO TESTS - tigther context
ovl <- aaProfile(aam[,as.character(-7:7)], abc=AAT)
plotdev(file.path(mfig.path,paste0("AAS_overlap")),
        height=5, width=5, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE)
dev.off()

plotdev(file.path(mfig.path,paste0("AAS_dotplot")),
        height=5, width=5, res=300)
par(mai=c(.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.05)
dotprofile(ovl, value="ratio", vcols=ttcols, xlab=NA,
           p.dot=p.min, lg2=TRUE, mxr=2,
           dot.sze=c(.3,2), ylab="amino acid", axis=1:2, show.total=TRUE)
dev.off()

## DEFINE CLASSES

classes <- cbind(
    kr=apply(aam[,as.character(-3:-1)], 1,
             function(x) any(x%in%c("K","R"))),
    Met=apply(aam[,as.character(c(-2,-1,1,2))], 1,
              function(x) sum(x%in%c("M"))) > 0,
    Trp=apply(aam[,as.character(c(-2,-1,1,2))], 1,
              function(x) any(x%in%c("W"))),
    Acidic= apply(aam[,as.character(c(-6:-1,1:6))], 1,
               function(x) sum(x%in%c("E","D")))>3,
    long=bdat$len > 1000  & bdat$median > -1,
    high=bdat$median > -1,
    Q=bdat$from=="Q",
    QA=bdat$fromto=="Q:A",
    QG=bdat$fromto=="Q:G",
    VQ=bdat$fromto=="V:Q",
    branched=bdat$from%in%c("I","L","V"),
    Nterm1=bdat$site%in%1,
    Nterm2=bdat$site%in%2)

fts <- unique(bdat$fromto)
ftcls <- matrix(NA, ncol=length(fts), nrow=nrow(aam))
colnames(ftcls) <-  paste0("fromto_",fts)
for ( i in seq_along(fts) ) 
    ftcls[,i] <- bdat$fromto==fts[i]
classes <- cbind(classes, ftcls)

frm <- unique(bdat$from)
ftcls <- matrix(NA, ncol=length(frm), nrow=nrow(aam))
colnames(ftcls) <- paste0("from_",frm)
for ( i in seq_along(frm) ) 
    ftcls[,i] <- bdat$from==frm[i]
#classes <- cbind(classes, ftcls)

frm <- unique(bdat$to)
ftcls <- matrix(NA, ncol=length(frm), nrow=nrow(aam))
colnames(ftcls) <- paste0("to_",frm)
for ( i in seq_along(frm) ) 
    ftcls[,i] <- bdat$to==frm[i]
#classes <- cbind(classes, ftcls)

cols <- as.character(c(-3:3))
rngs <- rep(list(cols), ncol(classes))
names(rngs) <- colnames(classes)

## custom ranges
rngs$QG <- as.character(-3:1)
rngs$QA <- as.character(-2:1)
rngs$Acidic <- rngs$Long <- rngs$High <-as.character(-10:10)
rngs$Met <- rngs$Trp <- as.character(-2:2)

psig <- 10^-c(3,5,10)

## use our internal AA colors
ASN2 <- ASN
ASN2$cols <- aa.cols[ASN2$chars]
ASN2$cols["V"] <- "#2eb774"


##pwm <- getPFM(aam)
##seqLogo(pwm, stackHeight=sumProbabilities, alphabet=ASN2)

for ( i in 1:ncol(classes) ) {

    id <- colnames(classes)[i]
    lb <- sub(".*_", "", id)
    filt <- classes[,id]
    cols <- rngs[[id]]
    axlab <- as.character(cols)
    axlab[axlab=="0"] <- ""
 
    if ( sum(filt)<2 ) next
    
    boxplot(bdat$median ~ filt)

    ## position weight matrices and diffLogo
    pfm1 <- getPFM(aam[ filt,cols,drop=FALSE])
    pfm2 <- getPFM(aam[!filt,cols,drop=FALSE])
    n1 <- sum(filt)
    n2 <- sum(!filt)
    dfob <- createDiffLogoObject(pfm1, pfm2, alphabet=ASN2)
    dfop <- enrichDiffLogoObjectWithPvalues(dfob, n1=n1, n2=n2)

    if ( any(dfop$pvals[cols!=0]<min(psig)) ) {

        mmai <- c(.5,.5,.1,.1)
        wd <- .4*length(cols) + .6
        ht <- 3
        mn <- dfob$ylim.negMax
        mx <- dfob$ylim.posMax
        plotdev(file.path(mfig.path,paste0("logos_", id)),
                height=ht, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=mn, ymax=mx)
        axis(1, at=1:length(cols), labels=axlab)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels="eAA", las=2)
        figlabel(paste0(lb,": ",sum(filt)), pos="topright", font=2)
        diffLogo_addPvals(dfop, ymin=mx)
        dev.off()
    }

    ## position weight matrices and diffLogo
    pfm1 <- getPFM(aas[ filt,cols,drop=FALSE])
    pfm2 <- getPFM(aas[!filt,cols,drop=FALSE])
    n1 <- sum(filt)
    n2 <- sum(!filt)
    dfob <- createDiffLogoObject(pfm1, pfm2, alphabet=ASN2)
    dfop <- enrichDiffLogoObjectWithPvalues(dfob, n1=n1, n2=n2)

    if ( any(dfop$pvals[cols!=0]<min(psig)) ) {

        mmai <- c(.5,.5,.1,.1)
        wd <- .4*length(cols) + .6
        ht <- 3
        mn <- dfob$ylim.negMax
        mx <- dfob$ylim.posMax
        plotdev(file.path(mfig.path,paste0("logos_", id,"_S")),
                height=ht, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=mn, ymax=mx)
        axis(1, at=1:length(cols), labels=axlab)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels="iAA", las=2)
        figlabel(paste0(lb,": ",sum(filt)), pos="topright", font=2)
        diffLogo_addPvals(dfop, ymin=mx)
        dev.off()
    }
    
}

plotdev(file.path(mfig.path,"protein_location_G"), ftyp,
        width=3, height=3, res=300)
hist(bdat$rpos[bdat$to=="G"])
dev.off()

## TODO:
## * RAAS profiles for EACH selection, by cancer types.

## * GENERAL PEPTIDE ANALYSIS

pp.min <- p.min #1e-20
pp.txt <- p.txt #1e-10

## from N-term
paa <- strsplit(bdat$BP, "")
paa <- do.call(rbind,lapply(paa, function(x) x[1:10]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT[!AAT%in%c("K","R","P","D","E")]))
plotdev(file.path(mfig.path,paste0("peptides_AA_Nterm")),
        height=5, width=5, res=300, type=ftyp)
par(mai=rep(.5,4), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=pp.min, p.txt=pp.txt, show.total=TRUE,
             xlab="Distance from peptide N-terminus", ylab="AA")
dev.off()

## from C-term
paa <- strsplit(bdat$BP, "")
paa <- do.call(rbind,lapply(paa, function(x) rev(x)[1:10]))
paa[is.na(paa)] <- "-"
ovl <- aaProfile(paa, abc=c(AAT[!AAT%in%c("K","R","P","D","E")]))
plotdev(file.path(mfig.path,paste0("peptides_AA_Cterm")),
        height=5, width=5, res=300, type=ftyp)
par(mai=rep(.5,4), mgp=c(1.3,.3,0), tcl=-.25)#, yaxs="i")
plotOverlaps(ovl, p.min=pp.min, p.txt=pp.txt, show.total=TRUE,
             xlab="Distance from peptide C-terminus", ylab="AA")
dev.off()
