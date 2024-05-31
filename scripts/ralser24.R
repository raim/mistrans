### FOR RALSER LAB PRESENTATION
## 1000+ yeasts - transcriptome, proteome and protein degradation
library(readxl)
library(segmenTools)
options(stringsAsFactors=FALSE)

data.path <- "~/data/"
mam.path <- file.path(data.path, "mammary")
yeast.path <- file.path(data.path, "yeast")

fig.path <- file.path(mam.path,"figures","ralser")
dir.create(fig.path)
ftyp <- "png"

## heatmap colors
ttcols <- colorRampPalette(c("#0000FF", "#FFFFFF","#FF0000"))(100)
gr.lab <- expression(growth~rate/h^-1)
dg.lab <- expression(degradation~rate/h^-1)
mdg.lab <- expression(median~degradation~rate/h^-1)
hl.lab <- expression("half-life"/h)

## data files
yfeature.file <- file.path(yeast.path,"feature_R64-1-1_20110208_allData.csv")
yhlf.file <- file.path(yeast.path,"processedData","christiano14.csv")
cls.file <- file.path(yeast.path,"processedData","redoxClusters.RData")
c24.file <- "~/data/yeast/originalData/caudal24_stables.xlsx"
m24.file <- "~/data/yeast/originalData/muenzner24_stables.xlsx"

## load yeast gene table
ygenes <- read.delim(yfeature.file)
ygenes <- ygenes[ygenes$type=="gene",]

### load half-lives
yhld <- read.delim(yhlf.file)
yhlf <- yhld[match(ygenes$ID, yhld$ENSG), "t1.2..hours."]
yhlf[which(yhlf==">= 100")] <- 120
yhlf <- as.numeric(yhlf)

## get 2012 clustering
load(cls.file)
cls.srt <- all$coarse$srt
cls.col <- all$coarse$col[cls.srt]
## ONLY MAJOR
cls.srt.major <- grep("NA",grep("[A-Z]",cls.srt,value=TRUE),
                      invert=TRUE,value=TRUE) 
CLM <- factor(ygenes$CL_rdx, levels=cls.srt.major)

## TRANSCRIPTOME, blindly using log2 fold change on sheet 7
## Caudal et al. 2024, Table S6
c24 <- as.data.frame(read_xlsx(c24.file, sheet=7))

plotdev(file.path(fig.path,"transcriptome_clusters"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.3,.05), mgp=c(1.3,.3,0), tcl=-.25)
idx <- match(ygenes$ID, c24$ORF)
bp <- boxplot(c24$log2FoldChange[idx] ~ CLM,
        col=cls.col[cls.srt.major], las=2, xlab=NA,
        ylab="log2FoldChange", axes=FALSE)
mtext("yeast co-expression cohorts", 1, 1.6)
axis(2)
Map(axis, 1, at=1:length(bp$n), labels=bp$names, las=2,
    col.axis=cls.col[bp$names])
Map(axis, 3, at=1:length(bp$n), labels=bp$n, las=2, cex.axis=.7)
box()
legend("bottomleft", "Caudal et al. 2024", bty="n", seg.len=0)
dev.off()

plotdev(file.path(fig.path,"transcriptome_clusters_arrow"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.3,.05), mgp=c(1.3,.3,0), tcl=-.25)
idx <- match(ygenes$ID, c24$ORF)
bp <- boxplot(c24$log2FoldChange[idx] ~ CLM,
        col=cls.col[cls.srt.major], las=2, xlab=NA,
        ylab="log2FoldChange", axes=FALSE)
axis(3, at=1:length(bp$n), labels=bp$n, las=2, cex.axis=.7)
axis(2, las=2)
##mtext("yeast co-expression cohorts", 1, 1.6)
arrows(x0=1.5, x1=length(bp$n)-1, y0=-9, xpd=TRUE, length=.05, lwd=3)
text(length(bp$n)/2, -9.5, labels="time", xpd=TRUE, font=2)
box()
legend("bottomleft", "Caudal et al. 2024", bty="n", seg.len=0)
dev.off()


## @Muenzner2024: SI Tables 17 (Dynamic SILAC) and 16 (growth rate) 
## strain wise half-lives
## kdp denotes the turnover rate in (1/h), HL the protein half-life
## (in h), and Intensity is an abundance estimate
m24 <- as.data.frame(read_xlsx(m24.file, sheet=17,skip=3))

## strain doubling times and growth rates
m24mu <- as.data.frame(read_xlsx(m24.file, sheet=16,skip=3))
mus <- log(2)/(m24mu[,2]/60) # growth rate/h-1
names(mus) <- m24mu[,1]

## TODO: growth rate vs. abundance
## for each strain, calculate Ribi/RP vs. S/C ratio
## and compare against growth rate.

## split half-live data by strain!
m24genes <- unique(m24$systematic_name)

m24strains <- split(m24, m24$strain)

## align all data and bring into matrix format

## intensities
m24int <- lapply(m24strains,
                 function(x) x[match(m24genes,x[,"systematic_name"]),
                               "Intensity"])
m24int <- do.call(cbind, m24int)
rownames(m24int) <- m24genes
## normalize by total intensity
## TODO: can the be done better? available total data?
m24int <- t(t(m24int)/apply(m24int, 2, sum, na.rm=TRUE))

## wrong half-lives w/o growth correction
m24hlfw <- lapply(m24strains,
                  function(x) x[match(m24genes,x[,"systematic_name"]),
                                "HL"])
m24hlfw <- do.call(cbind, m24hlfw)
rownames(m24hlfw) <- m24genes

## strain correlation of half-lives
com <- cor(m24hlfw, use="pairwise.complete")
image_matrix(com, breaks=seq(-1,1,length.out=length(ttcols)+1), col=ttcols)

## strain correlation of intensities
com <- cor(log10(m24int), use="pairwise.complete")
image_matrix(com, breaks=seq(.7,1,length.out=length(ttcols)+1), col=ttcols)

## CORRECTED FOR GROWTH RATE
## deg_apparent = deg + mu


m24degw <- log(2)/m24hlfw
## to test matrix addition: m24degw[] <- 0
m24deg <- t(t(m24degw) - mus)
m24hlf <- log(2)/m24deg ## CORRECTED HALF-LIVES

hist(c(m24hlfw), breaks=100, xlab="reported half-life/h")
hist(c(m24hlf), breaks=100, xlab="reported half-life/h")


brks <- seq(0,5,.05)
hist(mus, border=2, breaks=brks, freq=FALSE,
     xlab=expression(reported~degradation~rate/h^-1))#
hist(m24degw, breaks=brks, freq=FALSE, add=TRUE)

## corrected -> below 0
hist(m24deg, breaks=100)
hist(ash(c(m24hlf)), xlab="arcsinh corrected half-lives")


## NOTE: TOO SHORT HALF-LIVES, USING NON-CORRECTED DATA BELOW

## translation rate ~ int*(deg+mu)
m24trnsl <- m24int * log(2)/m24hlfw

### PLOT

## map to our main gene set
idx <- match(ygenes$ID, m24genes)
for ( i in 1:nrow(m24mu) ) {
    strain <- m24mu$strain[i]
    int <- m24int[,strain]
    plotdev(file.path(fig.path,paste0("proteome_clusters_",strain)),
            type=ftyp, height=3, width=3, res=200)
    par(mai=c(.5,.5,.3,.05), mgp=c(1.1,.3,0), tcl=-.25)
    bp <- boxplot(log10(int[idx]) ~ CLM,
                  ylim=c(-6, -1),
                  col=cls.col[cls.srt.major], las=2, xlab=NA,
                  ylab=expression(log[10](intensity/total)), axes=FALSE)
    axis(3, at=1:length(bp$n), labels=bp$n, las=2, cex.axis=.7)
    axis(2)
    mtext(paste("strain", strain), 1, 1, font=2, cex=1.3)
    dev.off()
}

idx <- match(ygenes$ID, m24genes)
for ( i in 1:nrow(m24mu) ) {
    strain <- m24mu$strain[i]
    hlf <- m24hlfw[,strain]
    plotdev(file.path(fig.path,paste0("halflives_clusters_",strain)),
            type=ftyp, height=3, width=3, res=200)
    par(mai=c(.5,.5,.3,.05), mgp=c(1,.3,0), tcl=-.25)
    bp <- boxplot(log10(hlf[idx]) ~ factor(ygenes$CL_rdx, levels=cls.srt.major),
                  col=cls.col[cls.srt.major],
                  ylim=log10(c(1,10)),
                  xlab=NA,
                  ylab=hl.lab, axes=FALSE)
    axis(3, at=1:length(bp$n), labels=bp$n, las=2, cex.axis=.7)
    axis(2, at=-1:10, labels=10^(-1:10))
    axis(2, at=log10(rep(1:10, 5) * 10^rep(-1:3, each=10)), tcl=-.125,
         labels=FALSE)
    mtext(paste("strain", strain), 1, 1, font=2, cex=1.3)
    dev.off()
}

idx <- match(ygenes$ID, m24genes)
for ( i in 1:nrow(m24mu) ) {
    strain <- m24mu$strain[i]
    hlf <- m24hlfw[,strain]
    plotdev(file.path(fig.path,paste0("degradation_cor_",strain)),
            type=ftyp, height=3, width=3, res=200)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(log10(log(2)/hlf[idx]), log10(log(2)/yhlf),
            xlab="log10 deg. rate, Muenzner 2024",
            ylab="log10 deg. rate, Christiano 2014")
    figlabel(pos="bottomleft", strain, font=2)
    abline(a=0, b=1, col=5)
    dev.off()
}

idx <- match(ygenes$ID, m24genes)
for ( i in 1:nrow(m24mu) ) {
    strain <- m24mu$strain[i]
    hlf <- m24hlfw[,strain]
    int <- m24int[,strain]
    trnsl <- m24trnsl[,strain]

    plotdev(file.path(fig.path,paste0("halflife_cor_",strain)),
            type=ftyp, height=3, width=3, res=200)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(log10(hlf[idx]), log10(yhlf),
            xlab="half-lives, Muenzner 2024/h",
            ylab="half-lives, Christiano 2014/h", axes=FALSE)
    abline(a=0, b=1, col=5)
    for ( ax in 1:2 ) {
        axis(ax, at=-1:10, labels=10^(-1:10))
        axis(ax, at=log10(rep(1:10, 5) * 10^rep(-1:3, each=10)), tcl=-.125,
             labels=FALSE)
    }
    figlabel(pos="bottomleft", strain, font=2)
    dev.off()
}

boxplot(log10(yhlf) ~ factor(ygenes$CL_rdx, levels=cls.srt.major),
        col=cls.col[cls.srt.major])


### COMPARISON TO GROWTH RATES
## TODO: reduce intensities to cluster means,
## TODO: account for strains, euploid vs. natural aneu vs. lab aneu,
## NOTE: anticorrelation of AB-D but no relation to growth rate,
##       BUT with broader deifnitions (A,AB) vs. (cd.n,D) there is no
##       anticorrelation but correlation of ratio to growth rate. 
## TODO: P = k*R/(dp+mu)
##       P*(dp+mu) = k*R,
##       P ~ mu or R ~ mu?

## strain-wise A/D ratio
idx <- match(ygenes$ID, m24genes)
int <- m24int[idx,]
rownames(int) <- ygenes$ID

cohorts <- split(ygenes$ID, ygenes$CL_rdx)
cohm <- matrix(0, ncol=length(cohorts), nrow=nrow(int))
colnames(cohm) <- names(cohorts)
rownames(cohm) <- rownames(int)
for ( j in seq_along(cohorts) ) 
    cohm[cohorts[[j]],j] <- 1        # boolean matrix
cohn <- t(t(cohm)/apply(cohm,2,sum)) # normalize by cohort size

counts <- t(int) # count matrix
counts[is.na(counts)] <- 0
## states = counts X cohorts
states <- as.data.frame(counts %*% cohn)
statel <- log10(states)

snames <- rownames(states)
shighlight <- snames
names(shighlight) <- snames
shighlight[!snames%in%c("BBV","CNL")] <- ""

plotdev(file.path(fig.path,"cohort_intensities_AD"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(states$A, states$D, density=FALSE, xlab=NA, ylab=NA)
mtext("cohort A",1,1.3, col=cls.col["A"], font=2)
mtext("cohort D",2,1.3, col=cls.col["D"], font=2)
shadowtext(states$A, states$D, labels=shighlight, pos=4, xpd=TRUE, col=1,cex=1)
mtext("cohort mean intensities", 3, 0.15)
dev.off()

plotdev(file.path(fig.path,"cohort_intensities_DC"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(states$D, states$C, density=FALSE, xlab=NA, ylab=NA)
mtext("cohort D",1,1.3, col=cls.col["D"], font=2)
mtext("cohort C",2,1.3, col=cls.col["C"], font=2)
shadowtext(states$D, states$C, labels=shighlight, pos=4, xpd=TRUE, col=1,cex=1)
mtext("cohort mean intensities", 3, 0.15)
dev.off()

plotdev(file.path(fig.path,"cohort_intensities_ABC"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(states$AB, states$C, density=FALSE, xlab=NA, ylab=NA)
mtext("cohort AB",1,1.3, col=cls.col["AB"], font=2)
mtext("cohort C",2,1.3, col=cls.col["C"], font=2)
shadowtext(states$AB, states$C, labels=shighlight, pos=4, xpd=TRUE, col=1,cex=1)
mtext("cohort mean intensities", 3, 0.15)
dev.off()

plotdev(file.path(fig.path,"cohort_intensities_ABD"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(states$AB, states$D, density=FALSE, xlab=NA, ylab=NA)
mtext("cohort AB",1,1.3, col=cls.col["AB"], font=2)
mtext("cohort D",2,1.3, col=cls.col["D"], font=2)
shadowtext(states$AB, states$D, labels=shighlight, pos=4, xpd=TRUE, col=1,cex=1)
mtext("cohort mean intensities", 3, 0.15)
dev.off()

## ratio vs. growth rate
rat <- log2(states$AB/states$D)
names(rat) <- rownames(states)

plotdev(file.path(fig.path,"cohort_growth_ABD"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(mus, rat[names(mus)], density=FALSE, xlab=gr.lab,
        ylab=expression(log[2](AB/D)),
        legpos="bottomright")
shadowtext(mus, rat[names(mus)], labels=shighlight[names(mus)], pos=4,
           xpd=TRUE, col=1, cex=1)
dev.off()

## instead use medians of broader classes 
rp <- apply(int, 2, function(x) median(x[ygenes$CL_rdx%in%c("A","AB")],
                                       na.rm=TRUE))
sc <- apply(int, 2, function(x) median(x[ygenes$CL_rdx==c("cd.n","D")],
                                       na.rm=TRUE))
plotdev(file.path(fig.path,"cohort_intensities_growth_stress"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(rp, sc, density=FALSE, xlab=gr.lab, ylab="stress/catabolism")
shadowtext(rp, sc, labels=shighlight, pos=4, xpd=TRUE, col=1,cex=1)
mtext("cohort mean intensities", 3, 0.15)
dev.off()

## ratio vs. growth rate
rat <- log2(rp/sc)
names(rat) <- rownames(states)

plotdev(file.path(fig.path,"cohort_growth_growth_stress"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(mus, rat[names(mus)], density=FALSE, xlab=expression(growth~rate/h^-1),
        ylab=expression(log[2](growth/stress)), legpos="bottomright")
shadowtext(mus, rat[names(mus)], labels=shighlight[names(mus)], pos=4,
           xpd=TRUE, col=1, cex=1)
dev.off()

## correlate median half life
mhlf <- apply(m24hlfw,2,median,na.rm=TRUE)
mdeg <- apply(m24degw,2,median,na.rm=TRUE)

plotdev(file.path(fig.path,"halflives_growthrate"),
        type=ftyp, height=3, width=3, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(mus, mdeg[names(mus)], density=FALSE,
        xlab=gr.lab, ylab=mdg.lab, legpos="bottomright")
abline(a=0,b=1, col=5)
shadowtext(mus, mdeg[names(mus)], labels=shighlight[names(mus)],
           pos=4, xpd=TRUE, col=1, cex=1)
dev.off()
