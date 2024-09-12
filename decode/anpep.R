library(readxl)

infile <- '~/data/mammary/originalData/wang19_table_EV1.xlsx'
fig.path <- 'home/raim/data/mistrans/figures/'

dat <- read_xlsx(infile, sheet=4)

colnames(dat)

nrm <- dat[,3:31]
nrm <- t(t(nrm)/apply(nrm,2,sum, na.rm=TRUE))

gidx <- dat$"Gene name" == "ANPEP"

png(file.path(fig.path,"aminopeptidase_N_bytissue.png",
              units="in", width=10, height=5, res=300)
par(mai=c(1.6,1,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(nrm[gidx,], las=2)
mtext(expression(I[ANPEP]/I[total]), 2, 3.3)
legend("topright", "ANPEP: aminopeptidase N")
dev.off()

