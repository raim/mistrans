
library(segmenTools)
source("~/work/mistrans/scripts/saap_utils.R")


##tmt.file <- file.path("~/data/mistrans/originalData/",
##                      "All_SAAP_TMTlevel_quant_df_withTonsil.xlsx")
##dat <- readxl::read_xlsx(tmta.file)
##dat <- as.data.frame(data)
## use pre-filtered data here
tmt.file <- file.path("~/data/mistrans/originalData/",
                      "All_SAAP_TMTlevel_quant_df.xlsx")
##                      "All_filtered_SAAP_TMTlevel_quant_df_withTonsil.xlsx")
dat <- readxl::read_xlsx(tmt.file)
dat <- as.data.frame(dat)

dat <- dat[dat[,"Keep SAAP"],]

p.min <- 1e-10
p.txt <- 1e-3

tpath <- "~/data/mistrans/figures/tonsil"
dir.create(tpath)
setwd(tpath)

raas <- as.numeric(dat$RAAS)
raas[!is.finite(raas)] <- NA

dtyp <- dat$"TMT/Tissue"
dtyp[grep("^S", dtyp)] <- "CANCER"
ddtyp <- dtyp
ddtyp[ddtyp=="tonsil" & dat$Digest == "Trypsin"] <- paste0("TONSIL_Trypsin")
ddtyp[ddtyp=="tonsil" & dat$Digest != "Trypsin"] <- paste0("TONSIL_Other")

digest <- dat$Digest
digest[dat$Digest=="Trypsin" & dtyp=="CANCER"] <- "Trypsin_cancer"
digest[dat$Digest=="Trypsin" & dtyp!="CANCER"] <- "Trypsin_tissues"
digest[dat$Digest=="Trypsin" & dtyp=="tonsil"] <- "Trypsin_tonsil"

ovl <- clusterCluster(ddtyp, dat$AAS, alternative="two.sided")
fromq <- grep("^Q", colnames(ovl$p.value), value=TRUE)
ovlq <- t(sortOverlaps2(t(ovl), axis=2, srt=fromq))
ovlc <- sortOverlaps(ovl, axis=1, p.min=p.min, cut=TRUE)

plotdev("aas_dataset_fromq", width=6, height=6)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlq, p.min=p.min, p.txt=p.txt, xlab=NA, ylab=NA, show.total=TRUE)
dev.off()

plotdev("aas_dataset_signif", width=.25*ncol(ovlc$p.value)+2, height=6)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlc, p.min=p.min, p.txt=p.txt, xlab=NA, ylab=NA,
             show.total=TRUE, text.cex=.9)
dev.off()


ovl <- clusterCluster(ddtyp[ddtyp!="CANCER"],
                      dat$AAS[ddtyp!="CANCER"], alternative="two.sided")
fromq <- grep("^Q", colnames(ovl$p.value), value=TRUE)
ovlq <- t(sortOverlaps2(t(ovl), axis=2, srt=fromq))
ovlc <- t(sortOverlaps2(t(ovl), axis=2, p.min=p.min, cut=TRUE))

plotdev("aas_dataset_fromq_tissues", width=6, height=6)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlq, p.min=p.min, p.txt=p.txt, xlab=NA, ylab=NA, show.total=TRUE)
figlabel("tissues", pos="bottomleft",font=2, cex=1.2)
dev.off()

plotdev("aas_dataset_signif_tissues", width=.25*ncol(ovlc$p.value)+2, height=6)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlc, p.min=p.min, p.txt=p.txt, xlab=NA, ylab=NA,
             show.total=TRUE, text.cex=.9)
figlabel("tissues", pos="bottomleft",font=2, cex=1.2)
dev.off()



ovl <- clusterCluster(digest, dat$AAS, alternative="two.sided")
ovlc <- t(sortOverlaps2(t(ovl), axis=2, p.min=p.min, cut=TRUE))
plotdev("aas_digest_signif", width=.25*ncol(ovlc$p.value)+2, height=2.5)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlc, p.min=p.min, p.txt=p.txt, xlab=NA, ylab=NA, show.total=TRUE, text.cex=.9)
dev.off()

ovl <- clusterCluster(digest[dtyp!="CANCER"],
                      dat$AAS[dtyp!="CANCER"], alternative="two.sided")
ovlc <- t(sortOverlaps2(t(ovl), axis=2, p.min=p.min, cut=TRUE))
plotdev("aas_digest_signif_tissues", width=.25*ncol(ovlc$p.value)+2, height=2.5)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlc, p.min=p.min, p.txt=p.txt, xlab=NA, ylab=NA, show.total=TRUE, text.cex=.9)
figlabel("tissues", pos="bottomleft",font=2, cex=1.2)
dev.off()

ovl <- clusterCluster(digest[dtyp=="tonsil"],
                      dat$AAS[dtyp=="tonsil"], alternative="two.sided")
ovlc <- t(sortOverlaps2(t(ovl), axis=2, p.min=1e-2, cut=TRUE))
plotdev("aas_digest_signif_tonsil", width=.25*ncol(ovlc$p.value)+2, height=2.5)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlc, p.min=p.min, p.txt=p.txt, xlab=NA, ylab=NA, show.total=TRUE, text.cex=.9)
figlabel("tonsil", pos="bottomleft",font=2, cex=1.2)
dev.off()

## RAAS

plotdev("tonsil_raas", width=6, height=4)
par(mai=c(1.2,.5,.1,.5), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(raas ~ dtyp, las=2, xlab=NA, ylab="RAAS")
abline(h=-1.35, col=2)
axis(4, at=-1.35, col.axis=2, col=2, las=2)
dev.off()

plotdev("digest_raas", width=3, height=4)
par(mai=c(1.2,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(raas ~ digest, las=2, xlab=NA, ylab="RAAS")
abline(h=-1.35, col=2)
dev.off()

plotdev("digest_raas_tonsil", width=3, height=4)
par(mai=c(1.2,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(raas[dtyp=="tonsil"] ~ dat$Digest[dtyp=="tonsil"],
        las=2, xlab=NA, ylab="RAAS")
figlabel("tonsil", pos="bottomleft",font=2, cex=1.2)
dev.off()


## inspect Q->G in tonsil
pases <- unique(dat$Digest)
for ( i in seq_along(pases) ) {
    dig <- pases[i]
    flt <- dat$AAS=="Q to G" &
        (dat$"TMT/Tissue"=="tonsil" & dat$Digest==dig)
    bps <- dat$BP[flt]
    rss <- round(raas[flt],1)
    sps <- dat$SAAP[flt]
    cat(paste0(">",dig,"\n",
               paste(paste(bps,",",sps,"\t",rss),collapse="\n"),"\n"))
}


## use pre-filtered data here
pat.file <- file.path("~/data/mistrans/originalData/",
                      "All_SAAP_patient_level_quant_df.txt")
dat <- read.delim(pat.file)


dat <- dat[dat[,"Keep.SAAP"],]

