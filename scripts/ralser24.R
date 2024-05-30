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

## TRANSCRIPTOME, blindly using log2 fold change on sheet 7
c24 <- as.data.frame(read_xlsx(c24.file, sheet=7))

idx <- match(ygenes$ID, c24$ORF)
boxplot(c24$log2FoldChange[idx] ~ factor(ygenes$CL_rdx, levels=cls.srt.major),
        col=cls.col[cls.srt.major])



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
## deg = deg_apparent+mu
## to test matrix addition: m24hlfw[] <- log(2)

hist(c(m24hlfw), breaks=100, xlab="reported half-life/h")

m24degw <- log(2)/m24hlfw
m24deg <- t(t(m24degw) - mus)
m24hlf <- log(2)/m24deg ## CORRECTED HALF-LIVES

brks <- seq(0,5,.05)
hist(m24degw, breaks=brks, freq=FALSE,
     xlab=expression(reported~degradation~rate/h^-1))
hist(mus, border=2, breaks=brks, freq=FALSE)#, add=TRUE)

## corrected -> below 0
hist(m24deg, breaks=100)
hist(ash(c(m24hlf)), xlab="arcsinh corrected half-lives")

## NOTE: too short half-lives, using non-corrected data below


## translation rate ~ int*(deg+mu)
m24trnsl <- m24int * log(2)/m24hlfw

## map to our main gene set
idx <- match(ygenes$ID, m24genes)
for ( i in 1:nrow(m24mu) ) {
    strain <- m24mu$strain[i]
    hlf <- m24hlfw[,strain]
    int <- m24int[,strain]
    trnsl <- m24trnsl[,strain]

    boxplot(log10(hlf[idx]) ~ factor(ygenes$CL_rdx, levels=cls.srt.major),
            col=cls.col[cls.srt.major])
    
    plotCor(log10(log(2)/hlf[idx]), log10(log(2)/yhlf))
    abline(a=0, b=1)

    boxplot(log10(int[idx]) ~ factor(ygenes$CL_rdx, levels=cls.srt.major),
            col=cls.col[cls.srt.major])

    boxplot(log(trnsl[idx]) ~ factor(ygenes$CL_rdx, levels=cls.srt.major),
            col=cls.col[cls.srt.major])
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

plotCor(states$AB, states$D, density=FALSE, xlab="RP", ylab="S/C")

rat <- log2(states$AB/states$D)
names(rat) <- rownames(states)

plotCor(mus, rat[names(mus)], density=FALSE, xlab=expression(growth~rate/h^-1),
        ylab="log2 ratio AB/D")
shadowtext(mus, rat[names(mus)], labels=names(mus), pos=4, xpd=TRUE, col=1)

## instead use medians of broader classes 
rp <- apply(int, 2, function(x) median(x[ygenes$CL_rdx%in%c("A","AB")],
                                       na.rm=TRUE))
sc <- apply(int, 2, function(x) median(x[ygenes$CL_rdx==c("cd.n","D")],
                                       na.rm=TRUE))
plotCor(rp, sc, density=FALSE, xlab="growth", ylab="stress/catabolism")

rat <- log2(rp/sc)

plotCor(mus, rat[names(mus)], density=FALSE, xlab=expression(growth~rate/h^-1),
        ylab="log2 ratio growth/stress broad")
shadowtext(mus, rat[names(mus)], labels=names(mus), pos=4, xpd=TRUE, col=1)
