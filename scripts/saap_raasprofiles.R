
source("~/work/mistrans/scripts/saap_utils.R")

library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

## TODO:
## * which AA are missing from mapped file and why? likely
## because BP didnt match, generate from SAAP/BP,
## * global or local background distribution?

####!!!
## TODO 20240215:
## * rows: plot only largest contributor! Q->G, T->V,
##   TODO: find largest RAAS effect size.
## * RAAS distributions: show all and scale/color by pval
## * formal approach to missing data: map measured peptides to the
##   full peptide space (main peptides or all possible).
## * plot by protein class, e.g only Albumins.
## * CHECK HANDLING OF DUPLICATE SAAP here and in tmt retrieval


## FROM->TO BY AA PROPERTY CLASSES
## DONT FILTER FOR CODONS, use hdat

## AA PROPERTY CLASSES
## https://en.wikipedia.org/wiki/Amino_acid#/media/File:ProteinogenicAminoAcids.svg
AAPROP <- list(charged=c("R","H","K","D","E"),
               polar=c("S","T","N","Q"),
               special=c("C","U","G","P"),
               hydrophobic=c("A","V","I","L","M","F","Y","W"))
tmp <- unlist(AAPROP)
AAPROP <- sub("[0-9]+","", names(tmp))
names(AAPROP) <- tmp
aaprop <- sub("hydrophobic", "hphobic", AAPROP)

## plot colors
docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))

gcols <- grey.colors(n=100, start=.9, end=0)

#### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")

in.file <- file.path(out.path,"saap_mapped.tsv")
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")

### PARAMETERES

p.min <- 1e-10
p.txt <- 1e-5

exclude.albumin <- FALSE #  TRUE # 
exclude.frequent <- FALSE # TRUE # 
frequent <- c("Q","W","H")

LAB <- "all"
fig.path <- file.path(proj.path,"figures","saap_raasprofiles")
if (  exclude.albumin ) {
    fig.path <- paste0(fig.path,"_noalb")
    LAB <- "no Alb./Hemog."
}
if ( exclude.frequent ) {
    tmp <- paste(frequent,collapse=",")
    fig.path <- paste0(fig.path,"_", gsub(",","",tmp))
    LAB <- paste("no", tmp)
}
dir.create(fig.path, showWarnings=FALSE)

### START

## PARSE & FILTER DATA
dat <- read.delim(in.file)
tmtf <- read.delim(tmt.file)


dat$Hemoglobin.Albumin <- as.logical(dat$Hemoglobin.Albumin)

#### TODO: missing values

## NOTE: tagged to removed in TMT, all but Healthy, but NOT in protein level
if ( FALSE ) {
    saap <- "LGMFNIQHGK" # flagged Keep=FALSE in tmt but not in hdat?
    saap <- "WYNLAIGSTGPWLK"
    ##saap <- "GVEAAGAMFLEAIPMSIPPEVK"
    tmtf[which(tmtf$SAAP==saap),c("Keep.SAAP","Dataset")]
    dat[which(dat$SAAP==saap),c("remove","Keep.SAAP","Dataset")]
}

## filter hard removes
hdat <- dat[which(!dat$remove),]


## get raw RAAS data TMT level
tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)
## remove excluded
cat(paste("removing", sum(!tmtf$Keep.SAAP),
          "tagged as false positive on TMT level\n"))
tmtf <- tmtf[tmtf$Keep.SAAP,]
## exclude NA or Inf
rm <- is.na(tmtf$RAAS) | is.infinite(tmtf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from TMT level\n"))
tmtf <- tmtf[!rm,]

## CLEAN PROTEIN LEVEL
tmtsaap <- unique(tmtf$SAAP)
rm <- !hdat$SAAP%in%tmtsaap

cat(paste("removing", sum(rm), "SAAP without values at TMT level\n"))
hdat <- hdat[!rm,]


## first find mutated AA pairs
## (column AASin input mixes I/L)
## and split mutated AA pairs into from/to
## TODO: get this from the columns in the mapped file, once clear
## why missing
saaps <- strsplit(hdat$SAAP,"")
bases <- strsplit(hdat$BP, "")
fromtol <- lapply(1:length(saaps), function(i) {
    pos <- which(saaps[[i]]!=bases[[i]])
    c(from=bases[[i]][pos], to=saaps[[i]][pos])
})
fromto <- do.call(rbind, fromtol)

## replace from/to columns
hdat$from <- fromto[,1]
hdat$to <- fromto[,2]

## fromto column
hdat$fromto <- apply(fromto, 1, paste0, collapse=":")

## generate columns by AA property
hdat$pfrom <- aaprop[fromto[,1]]
hdat$pto <- aaprop[fromto[,2]]
hdat$pfromto <- paste0(hdat$pfrom,":",hdat$pto)


#### FILTER

if (  exclude.albumin ) {
    ## TODO: also exclude from TMF level file
    hdat <- hdat[!hdat$Hemoglobin.Albumin,]
}

if ( exclude.frequent ) {

    ## TODO: MAKE SURE TO REMOVE THE SAME TYPES!!
    ## this is probably not required since we
    ## only look up SAAP that appear in hdat!

    rmsaap <- hdat$SAAP[hdat$from%in%frequent]
    tmtf <- tmtf[!tmtf$SAAP%in%rmsaap,]

    hdat <- hdat[!hdat$from%in%frequent,]
}

## find which are missing
##length(grep("NA",hdat$pfromto))
##adat <- hdat[grep("NA",hdat$pfromto, invert=TRUE),]

## sorting and labels
uds <- sort(unique(hdat$Dataset))
uds <- c("Healthy",uds[uds!="Healthy"])

srt <- c("charged","polar","hphobic","special")
srt <- paste(rep(srt,each=4), rep(srt,4), sep=":")
axex <- rep("",length(srt))
names(axex) <- srt
for ( i in seq_along(srt) ) {
    ft <- unlist(strsplit(srt[i],":"))
    axex[i] <- as.expression(bquote(.(ft[1]) %->% .(ft[2])))
}

fname <- file.path(fig.path,paste0("AAprop_wtests_"))
source("~/work/mistrans/scripts/saap_utils.R")
ovw <- raasProfile(x=hdat, id="SAAP", values=tmtf, 
                   rows="pfromto", cols="Dataset",
                   bg=TRUE, vid="RAAS", row.srt=srt, col.srt=uds,
                   use.test=w.test, do.plots=FALSE, xlab="TMT level RAAS",
                   fname=fname, verb=0)

## common histogram with color by sign and

## first, calculate densities
dns <- ovw$ALL
for ( i in 1:nrow(ovw$p.value) ) {
    for ( j in 1:ncol(ovw$p.value) ) {
        rw <- rownames(ovw$p.value)[i]
        cl <- colnames(ovw$p.value)[j]
        vls <- ovw$ALL[[rw]][[cl]]
        if (length(vls)>1 ) {
            dns[[rw]][[cl]] <- density(vls)
        }
    }
}
## get max
mxs <- rep(NA, length(dns))
for ( i in seq_along(dns) )
    mxs[i] <- max(unlist(lapply(dns[[i]],
                                function(z)
                                    ifelse("y"%in%names(z),
                                           max(unlist(z["y"])),-10))))

png(file.path(fig.path,paste0("AAprop_wtests_densities.png")),
    res=300, width=4, height=2, units="in")
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
hist(tmtf$RAAS, col="#77777755", border=NA, 
     xlim=c(-6,4),
     xlab="TMT level RAAS", ylab=NA, main=NA, axes=FALSE)
dns <- ovw$ALL
par(new=TRUE)
plot(1, col=NA, xlim=c(-6,4), ylim=c(0,.75),#max(mxs)),
     xlab=NA, ylab=NA, axes=FALSE)
for ( i in 1:nrow(ovw$p.value) ) {
    for ( j in 1:ncol(ovw$p.value) ) {
        rw <- rownames(ovw$p.value)[i]
        cl <- colnames(ovw$p.value)[j]
        vls <- ovw$ALL[[rw]][[cl]]
        if (length(vls)>1 ) {
            dns[[rw]][[cl]] <- density(vls)
            if ( ovw$p.value[i,j] <1e-3 )
                lines(dns[[rw]][[cl]],
                      col=ifelse(ovw$sign[i,j]==-1,"blue","red"),
                      lwd=-log10(ovw$p.value[i,j])/-log10(p.min))
            ## TODO: annotate
            
        }# else abline(v=vls)
    }
}
axis(2, labels=NA)
mtext("density",2,.5)
axis(1)
dev.off()



png(file.path(fig.path,paste0("AAprop_wtests_legend.png")),
    res=300, width=2, height=2, units="in")
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlapsLegend(p.min=p.min, p.txt=p.txt, type=2, col=ttcols)
dev.off()

png(file.path(fig.path,paste0("AAprop_wtests.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1, ylab=NA, xlab="", col=ttcols, show.total=TRUE)
axis(2, length(axex):1, labels=axex, las=2)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

## local raas colors
png(file.path(fig.path,"raas_profile_colors.png"),
    res=300, width=3, height=3, units="in")
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(c(ovw$median),
                          q=0, colf=viridis::viridis,
                          mn=-3, mx=0,
                          n=10,
                          xlim=c(-4,1),
                          plot=TRUE,
                          mai=c(.5,.5,.1,.1), xlab="All TMT level RAAS")
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(LAB, pos="bottomleft", cex=1)
dev.off()

png(file.path(fig.path,paste0("AAprop_wtests_raas_median.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
txt <- ovw$count
txt[txt==0] <- ""
txt.col <- ifelse(ovw$median< -2, "white","black")
image_matrix(ovw$median, col=lraas.col$col, cut=TRUE, breaks=lraas.col$breaks,
             axis=1, text=txt, text.col=txt.col, ylab=NA,
             xlab="", text.cex=.8)
axis(2, length(axex):1, labels=axex, las=2)
mtext("median TMT-level RAAS", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

png(file.path(fig.path,paste0("AAprop_wtests_count_unique.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
cnt <- ovw$unique
cnt[cnt==0] <- NA
txt <- ovw$unique
txt[txt==0] <- ""
txt.col <- ifelse(ovw$unique>quantile(c(ovw$unique),.95),"white","black")
image_matrix(cnt, col=gcols, axis=1:2,
             text=txt, text.col=txt.col, ylab="",
             xlab="", text.cex=.8)
##axis(4, nrow(ovw$p.value):1, labels=ACODONS[rownames(ovw$p.value)], las=2)
mtext("mistranslated", 2, 3.5, adj=-1.7, cex=1.2)
mtext("# unique SAAP", 3, .5, cex=1.2)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()

## by AA for reference - TODO: do this by dataset
## as previously

for ( ds in uds ) {
    ddat <- hdat[hdat$Dataset==ds,]
    dtmt <- tmtf[tmtf$Dataset==ds,]
    dtmt <- split(dtmt$RAAS,  dtmt$SAAP)

    ovw <- raasProfile(x=ddat, id="SAAP", values=dtmt, 
                       rows="to", cols="from",
                       bg=FALSE, vid="RAAS", 
                       use.test=w.test, do.plots=FALSE, 
                       verb=0)

    png(file.path(fig.path,paste0(ds,"_AA_wtests.png")),
        res=300, width=5, height=5, units="in")
    par(mai=c(.5,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
                 show.total=TRUE)
    figlabel(LAB, pos="bottomleft", cex=1.5)
    figlabel(pos="bottomright", text=ds, cex=2, font=2)
    dev.off()
}

### FIND MOST FREQUENT
## TODO: instead find strongest
## TODO: remove most frequent per Dataset, eg. N:M vs. T:V in PDAC

ptype <- unique(hdat$pfromto)
pmost <- sapply(ptype, function(x) {
    names(which.max(table(hdat$fromto[hdat$pfromto==x])))
})[srt]

tdat <- hdat[hdat$fromto%in%pmost,]

ovw <- raasProfile(x=tdat, id="SAAP", values=tmtf, 
                   rows="fromto", cols="Dataset",
                   bg=TRUE, vid="RAAS", row.srt=pmost, col.srt=uds,
                   use.test=w.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)

png(file.path(fig.path,paste0("AAprop_wtests_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
figlabel(LAB, pos="bottomleft", cex=1.5)
mtext("most frequent", 2, 3, cex=1.5)
dev.off()

## remove most frequent
tdat <- hdat[!hdat$fromto%in%pmost,]

ovw <- raasProfile(x=tdat, id="SAAP", values=tmtf, 
                   rows="pfromto", cols="Dataset",
                   bg=TRUE, vid="RAAS", row.srt=srt, col.srt=uds,
                   use.test=w.test, do.plots=FALSE, xlab="TMT level RAAS",
                   verb=0)

png(file.path(fig.path,paste0("AAprop_wtests_wo_best.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
figlabel(LAB, pos="bottomleft", cex=1.5)
mtext("w/o freq.", 1, 2, cex=1.5, adj=-1)
dev.off()
