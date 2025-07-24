
file <- "/home/raim/data/mistrans/processedData/Supplemental_Data_7.SAAP_coordinates.tsv"

dat <- read.delim(file)

## only unique SAAP
dat <- dat[!duplicated(dat$SAAP),]

barplot(table(dat$site))

relpos <- dat$site/nchar(dat$SAAP)

hist(relpos)



segmenTools::dense2d(dat$site, dat$RAAS)
segmenTools::dense2d(relpos, dat$RAAS)

sum(dat$fromto=="Q:G")/nrow(dat)
sum(dat$site<4)/nrow(dat)
sum(dat$fromto=="Q:G" & dat$site<4)/nrow(dat)


tmt.file <- file.path("~/data/mistrans/originalData/",
                      "All_SAAP_TMTlevel_quant_df.xlsx")
##                      "All_filtered_SAAP_TMTlevel_quant_df_withTonsil.xlsx")
ton <- readxl::read_xlsx(tmt.file)
ton <- as.data.frame(ton)

bpsaap_data <- paste0(dat$BP, "_", dat$SAAP)
bpsaap_tonsil <- paste0(ton$BP, "_", ton$SAAP)

ton <- ton[match(bpsaap_data, bpsaap_tonsil),]

table(ton$Digest)
hist(relpos[ton$Digest=="Trypsin"], 
     breaks=seq(0,1,.05), xlab="relative position in peptide")
hist(relpos[ton$Digest!="Trypsin"], add=TRUE, col=2, breaks=seq(0,1,.05))

barplot(table(dat$site), xlab="abs. position in peptide (dist. from N-term)")
barplot(table(dat$site[ton$Digest!="Trypsin"]), add=TRUE, col=2)
