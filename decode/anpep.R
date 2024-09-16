library(readxl)

infile <- '~/data/mammary/originalData/wang19_table_EV1.xlsx'
fig.path <- '/home/raim/data/mistrans/figures/'

dat <- read_xlsx(infile, sheet=4)

colnames(dat)

nrm <- dat[,3:31]
nrm <- t(t(nrm)/apply(nrm,2,sum, na.rm=TRUE))

prts <- rbind(ANPEP=nrm[dat$"Gene name" == "ANPEP",],
              NPEPPS=nrm[dat$"Gene name" == "NPEPPS",])

png(file.path(fig.path,"aminopeptidase_N_bytissue.png"),
              units="in", width=10, height=5, res=300)
par(mai=c(1.6,1,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(prts, las=2, beside=TRUE, legend=TRUE)
mtext(expression(I[ANPEP]/I[total]), 2, 3.3)
#legend("topright", "ANPEP: aminopeptidase N")
dev.off()

##plot(log(prts[1,]), log(prts[2,]))
