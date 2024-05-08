
library(segmenTools)

tmt.file <- file.path("~/data/mistrans/originalData/",
                      "All_SAAP_TMTlevel_quant_df_withTonsil.xlsx")
tmt.file <- file.path("~/data/mistrans/originalData/",
                      "All_filtered_SAAP_TMTlevel_quant_df_withTonsil.xlsx")
dat <- readxl::read_xlsx(tmt.file)
dat <- as.data.frame(dat)

head(dat)


raas <- as.numeric(dat$RAAS)
raas[!is.finite(raas)] <- NA

dtyp <- dat$"TMT/Tissue"
dtyp[grep("^S", dtyp)] <- "cancer"
ddtyp <- dtyp
ddtyp[ddtyp=="tonsil" & dat$Digest == "Trypsin"] <- paste0("Trypsin_tonsil")

digest <- dat$Digest
digest[digest=="Trypsin" & dtyp=="tonsil"] <- "Trypsin_tonsil"

ovl <- clusterCluster(ddtyp, dat$AAS, alternative="two.sided")

fromq <- grep("^Q", colnames(ovl$p.value), value=TRUE)
ovlq <- sortOverlaps(ovl, axis=1, srt=fromq)
ovlc <- sortOverlaps(ovl, axis=1, p.min=1e-10, cut=TRUE)

plotdev("aas_dataset_fromq", width=6, height=6)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlq, p.min=1e-10, xlab=NA, ylab=NA, show.total=TRUE)
dev.off()

plotdev("aas_dataset_signif", width=.25*ncol(ovlc$p.value)+2, height=6)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlc, p.min=1e-10, xlab=NA, ylab=NA, show.total=TRUE, text.cex=.9)
dev.off()


ovl <- clusterCluster(digest, dat$AAS, alternative="two.sided")
ovlc <- sortOverlaps(ovl, axis=1, p.min=1e-5, cut=TRUE)
plotdev("aas_digest_signif", width=.25*ncol(ovlc$p.value)+2, height=2.5)
par(mai=c(.6,1.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovlc, p.min=1e-10, xlab=NA, ylab=NA, show.total=TRUE, text.cex=.9)
dev.off()


plotdev("tonsil_raas", width=6, height=4)
par(mai=c(1.2,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(raas ~ dtyp, las=2, xlab=NA, ylab="RAAS")
dev.off()

plotdev("digest_raas", width=3, height=4)
par(mai=c(1.2,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(raas ~ digest, las=2, xlab=NA, ylab="RAAS")
dev.off()

plotdev("digest_raas_tonsil", width=3, height=4)
par(mai=c(1.2,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
boxplot(raas[dtyp=="tonsil"] ~ dat$Digest[dtyp=="tonsil"],
        las=2, xlab=NA, ylab="RAAS")
dev.off()
