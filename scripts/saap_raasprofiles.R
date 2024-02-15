
source("~/work/mistrans/scripts/saap_utils.R")

library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

## TODO:
## * which AA are missing from mapped file and why? likely
## because BP didnt match, generate from SAAP/BP,
## * global or local background distribution?

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

#### START

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")

in.file <- file.path(out.path,"saap_mapped.tsv")
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")

incl.albumin <- FALSE #  TRUE # 
exclude.frequent <- FALSE # TRUE # 
frequent <- c("Q")#,"T","S")

LAB <- "all"
fig.path <- file.path(proj.path,"figures","saap_raasprofiles")
if ( !incl.albumin ) {
    fig.path <- paste0(fig.path,"_noalb")
    LAB <- "no Alb./Hemog."
}
if ( exclude.frequent ) {
    fig.path <- paste0(fig.path,"_nofreq")
    LAB <- paste("no", paste(frequent,collapse=";"))
}
dir.create(fig.path, showWarnings=FALSE)

## PARSE & FILTER DATA
dat <- read.delim(in.file)
dat$Hemoglobin.Albumin <-
    as.logical(dat$Hemoglobin.Albumin)
hdat <- dat[!dat$remove,]



## get raw RAAS data TMT level
tmtf <- read.delim(tmt.file)
tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)
## remove excluded
cat(paste("removing", sum(!tmtf$Keep.SAAP),
          "tagged as false positive on TMT level\n"))
tmtf <- tmtf[tmtf$Keep.SAAP,]
## exclude NA or Inf
rm <- is.na(tmtf$RAAS) | is.infinite(tmtf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from TMT level\n"))
tmtf <- tmtf[!rm,]
tmt <- split(tmtf$RAAS, tmtf$SAAP)



## first find mutated AA pairs
## (column AASin input mixes I/L)
## and split mutated AA pairs into from/to
## TODO: get this from the columns in the mapped file, once clear
## why missing
saaps <- strsplit(hdat$SAAP,"")
bases <- strsplit(hdat$BP, "")
fromto <- lapply(1:length(saaps), function(i) {
    pos <- which(saaps[[i]]!=bases[[i]])
    c(from=bases[[i]][pos], to=saaps[[i]][pos])
})
fromto <- do.call(rbind, fromto)

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

if ( !incl.albumin ) {
    ## TODO: also exclude from TMF level file
    hdat <- hdat[!hdat$Hemoglobin.Albumin,]
}

if ( exclude.frequent ) {
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

fname <- file.path(fig.path,paste0("AAprop_all_wtests_"))
ovw <- raasProfile(x=hdat, id="SAAP", values=tmtf, 
                   rows="pfromto", cols="Dataset",
                   bg=TRUE, vid="RAAS", row.srt=srt, col.srt=uds,
                   use.test=w.test, do.plots=FALSE, xlab="TMT level RAAS",
                   fname=fname, verb=1)


png(file.path(fig.path,paste0("AAprop_all_wtests_legend.png")),
    res=300, width=2, height=2, units="in")
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlapsLegend(p.min=1e-10, p.txt=1e-5, type=2, col=ttcols)
dev.off()

png(file.path(fig.path,paste0("AAprop_all_wtests.png")),
    res=300, width=4, height=6, units="in")
par(mai=c(1,1.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=1e-10, p.txt=1e-5,
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

png(file.path(fig.path,paste0("AAprop_all_wtests_raas_median.png")),
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

png(file.path(fig.path,paste0("AAprop_all_wtests_count_unique.png")),
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

fname <- file.path(fig.path,paste0("AA_all_wtests_"))
ovw <- raasProfile(x=hdat, id="SAAP", values=tmtf, 
                   rows="fromto", cols="Dataset",
                   bg=TRUE, vid="RAAS", col.srt=uds,
                   use.test=w.test, do.plots=FALSE, xlab="TMT level RAAS",
                   fname=fname, verb=1)

png(file.path(fig.path,paste0("AA_all_wtests.png")),
    res=300, width=3, height=45, units="in")
par(mai=c(1,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovw, p.min=1e-10, p.txt=1e-5,
             text.cex=.8, axis=1:2, ylab=NA, xlab="", col=ttcols,
             show.total=TRUE)
figlabel(LAB, pos="bottomleft", cex=1.5)
dev.off()
