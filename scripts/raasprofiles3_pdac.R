
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)
library(gprofiler2)

### INVESTIGATING T->V in PDAC
## * counts and overlaps of BP/SAAP in other datasets,
## * TODO: load patient-level dataset to investigate
##   differences between tumour and adjacent tissues.


## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

## additional data
patient.file  <- file.path(dat.path,"All_SAAP_patient_level_quant_df.xlsx")
goslim.file  <- file.path(mam.path,"processedData","goslim.tsv")

pdfig.path <- file.path(fig.path,"PDAC")
dir.create(pdfig.path, showWarnings=FALSE)

corW <- corH <- 2.5
pmai <- c(.5,.5,.25,.25)
pmpg <- c(1.3,.3,0)

## LOCAL ENSEMBL GENE:PROTEIN MAPPING
g2p <- bdat$ensembl
names(g2p) <- bdat$gene

## ADD BP:SAAP COLUMN to select
tmtf$pair <- paste(tmtf$BP,tmtf$SAAP)
bdat$pair <- paste(bdat$BP,bdat$SAAP)

## PDAC ONLY
tmtp <- tmtf[tmtf$Dataset=="PDAC",]

## W/O PDAC
tmto <- tmtf[tmtf$Dataset!="PDAC",]

## T->V in all 
tmtv <- tmtf[tmtf$fromto=="T:V",]

## T->V in PDAC
tmtvp <- tmtf[tmtf$fromto=="T:V" & tmtf$Dataset=="PDAC",]

## T->V in all others
tmtvo <- tmtf[tmtf$fromto=="T:V" & tmtf$Dataset!="PDAC",]

## BP/SAAP LEVEL: T->V in PDAC
pdat <- bdat[bdat$pair%in%tmtp$pair,] # all BP/SAAP that appear in PDAC
tvdat <- pdat[pdat$fromto=="T:V",]

## COUNTS in PDAC vs OTHERS


## T->V/PDAC BP/SAAP in others: 1 in Healthy, 3 in LUAD
table(tmto$Dataset[which(tmto$pair%in%tvdat$pair)])

## proteins with multiple T->V
dup <- which(duplicated(tvdat$gene))[1]
tvdat$SAAP[which(tvdat$gene%in%tvdat$gene[dup])] # ACTG2 contains 3 distinct T->V

sort(table(tmtf$name[tmtf$fromto=="T:V"]))

sort(table(tmtf$Dataset[tmtf$fromto=="T:V" & tmtf$name=="ACTG1"]))

## ACTG1, unique T->V sites: 6 in PDAC, 1 in lung
sort(table(site$Dataset[site$fromto=="T:V" & site$name=="ACTG1"]))
sort(table(site$tissue[site$fromto=="T:V" & site$name=="ACTG1"]))

site[site$fromto=="T:V" & site$name=="ACTG1",]

### COUNTS for slides

sets <- c(all="tmtf",PDAC="tmtp","T:V"="tmtv", "T:V & PDAC"="tmtvp", "T:V & other"="tmtvo")

counts <- matrix(NA, ncol=length(sets), nrow=4)
colnames(counts) <- sets
rownames(counts) <- c("BP/SAAP", "AAS", "sites", "proteins") 
for ( set in sets) {

    tmp <- get(set)
    counts[1,set] <- length(unique(paste(tmp$BP, tmp$SAAP))) # unique BP/SAAP
    counts[2,set] <- length(unique(paste(tmp$unique.site, tmp$AAS))) # unique AAS
    counts[3,set] <- length(unique(tmp$unique.site)) # unique sites 
    counts[4,set] <- length(unique(tmp$ensembl)) # 1997 unique proteins
}
colnames(counts) <- names(sets)
knitr::kable(counts, format = "markdown")

## 25 unique  T->V protein sites in other tissues
length(unique(tmtf$unique.site[which(tmtf$Dataset!="PDAC" & tmtf$fromto=="T:V")]))
## 24 unique proteins in T->V/PDAC
length(unique(tmtf$ensembl[which(tmtf$Dataset!="PDAC" & tmtf$fromto=="T:V")]))

## T->V in other tissues
table(tmtf$BP[which(tmtf$Dataset!="PDAC" & tmtf$fromto=="T:V")],
      tmtf$TMT.Tissue[which(tmtf$Dataset!="PDAC" & tmtf$fromto=="T:V")])

## OVERLAP PLOT T->V in PDAC vs. OTHERS


ft <- unique(site$fromto)
aas.srt <- unique(c(grep("T:", ft, value=TRUE),grep(":V", ft, value=TRUE)))
aas.srt <- c("T:V", "Q:G")

dsets <- site$Dataset
##dsets[dsets=="Healthy"] <- site$Tissues[dsets=="Healthy"]
nsets <- stringr::str_count(dsets, ";")+1
##dsets[nsets > 1] <- ">1"
dsets[nsets > 2] <- ">2"
dsets[nsets > 3] <- ">3"
ds.srt <- unique(dsets)
ds.srt <- ds.srt[order(nchar(ds.srt))]

ovl <- clusterCluster(cl2=dsets, cl1=site$fromto,
                      cl1.srt=aas.srt, cl2.srt=ds.srt)
present <- colnames(ovl$overlap)[apply(ovl$overlap,2,sum)>0]
ovl <- sortOverlaps(ovl, axis=1, srt=present)

plotdev(file.path(pdfig.path,"TV_per_Dataset"),
        width=length(present)*.3+1, height=length(aas.srt)*.2+2)
par(mai=c(1.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovl, show.total=TRUE, xlab="", ylab="", p.min=p.min, p.txt=p.txt)
figlabel("count of unique site x AAS in Datasets", pos="bottom")
dev.off()

## same for tissues
dsets <- site$tissue
##dsets[dsets=="Healthy"] <- site$Tissues[dsets=="Healthy"]
nsets <- stringr::str_count(dsets, ";")+1
##dsets[nsets > 1] <- ">1"
##dsets[nsets > 2] <- ">2"
dsets[nsets > 3] <- ">3"
ds.srt <- unique(dsets)
ds.srt <- ds.srt[order(nchar(ds.srt))]

ovl <- clusterCluster(cl2=dsets, cl1=site$fromto,
                      cl1.srt=aas.srt, cl2.srt=ds.srt)
## NOTE: only T:V
present <- colnames(ovl$overlap)[ovl$overlap["T:V",]>0]
ovl <- sortOverlaps(ovl, axis=1, srt=present)

plotdev(file.path(pdfig.path,"TV_per_Tissue"),
        width=length(present)*.3+1, height=length(aas.srt)*.2+2.5)
par(mai=c(2,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovl, show.total=TRUE, xlab="", ylab="", p.min=p.min, p.txt=p.txt)
figlabel("count of unique site x AAS in Tissues", pos="bottom")
dev.off()


### COMPARE dataset sizes
## TODO: add T->V

fmat <- matrix(NA, ncol=2, nrow=nrow(tmtf))
##fmat[,3] <- tmtf$Dataset
fmat[,2] <- tmtf$Dataset
fmat[,1] <- tmtf$Dataset
##fmat[duplicated(tmtf$unique),2] <- "na"
fmat[duplicated(paste(tmtf$BP,tmtf$SAAP)),1] <- "na"

cnts <- apply(fmat, 2, table)
cnts <- do.call(rbind,lapply(cnts, function(x) x[uds]))
rownames(cnts) <- c("BP/SAAP", "#RAAS")
if ( interactive() ) {
    barplot(cnts, beside=TRUE, legend=TRUE)
    
}


## TODO:
## * Is there an enrichment for amino acids flanking the T-->V substitution ? 
## * Do we observe high frequency of T-->V in the "healthy" pancreas
##  from the label free dataset ?
## * Is the RAAS for T-->V different between the tumors and their
##  sounding tissues ?

### FUNCTIONS of T->V containing proteins in PDAC.


## get GOslim table for use with clusterAnnotation
got <- parseAnnotationList(genes[,c("ID","GOslim")]) 
## replace GO IDs by terms
terms <- read.delim(goslim.file)
trms <- terms[,2]
names(trms) <- terms[,1]
colnames(got) <- trms[colnames(got)]

if ( FALSE ) {
    ## TODO: raas profile?
    go.ovl <- raasProfile(x=tmtu, id="BP.SAAP", 
                          rows=got[tmtu$gene,], cols="Dataset",
                          bg=TRUE, value=value, 
                          col.srt=uds,
                          use.test=use.test, do.plots=FALSE,
                          xlab=xl.raas,
                          verb=0)
}

gois <- unique(tvdat$gene)

if ( FALSE )
    gores <- gprofiler2::gost(query = gois, organism = "human")

gotv <- got[gois,]
gotv <- gotv[, apply(gotv,2,sum)>0]

sort(apply(gotv, 2, sum))

pnms[g2p[names(which(gotv[,"cellular amino acid metabolic process"]))]]

## GLUD1: glutamate dehydrogenase 1,
## ALDH1A1: Aldehyde dehydrogenase 1 family, member A1,
## ACAT1: acetyl-CoA acetyltransferase 1, mitochondrial,
## VARS1: valyl-tRNA synthetase 1,
## GART: Trifunctional purine biosynthetic protein adenosine-3,
##       de novo purine biosynthesis.

## CPS1: Carbamoyl phosphate synthetase I (CPS I): transfers an
##       ammonia molecule to a molecule of bicarbonate that has been
##       phosphorylated by a molecule of ATP. The resulting carbamate
##       is then phosphorylated with another molecule of ATP. The
##       resulting molecule of carbamoyl phosphate leaves the enzyme.

## ETFA: Electron-transfer-flavoprotein, alpha subunit, also known as
##       ETF-α.[5] Together with Electron-transfer-flavoprotein, beta
##       subunit, encoded by the 'ETFB' gene, it forms the
##       heterodimeric electron transfer flavoprotein (ETF). The
##       native ETF protein contains one molecule of FAD and one
##       molecule of AMP, respectively.[6][7]

pats <- read_xlsx(patient.file)
pats <- as.data.frame(pats)
pats$sample <-
    sub("_[NT].*", "", pats[,"Sample name"])


pats <- pats[pats[,"Keep SAAP"],]

patv <- pats[pats$AAS=="T to V" & pats$Dataset=="PDAC",]

plotdev(file.path(pdfig.path,"TV_RAAS_per_tumor"),
        width=corW, height=3)
par(mai=c(.5,.5,.5,.05), mgp=c(1.3,.3,0), tcl=-.25)
bp <- boxplot(patv$RAAS ~ patv[,"Sample type"], ylab=xl.raas,
              xlab="T->V in PDAC; all patients")
axis(3, at=1:2, labels=bp$n, las=2)
dev.off()

## RAAS fold changes between tumor and normal pairs of BP/SAAP
patl <- split(patv, patv$sample)

x <- patl[[3]]

rp <- lapply(patl, function(x) {
    xl <- split(x, paste(x$BP, x$SAAP))
    
    rpaired <- lapply(xl, function(y) {
        if ( sum(y[,"Sample type"]=="Tumor") == sum(y[,"Sample type"]=="Normal") )
            data.frame(tumor= y[y[,"Sample type"]=="Tumor" ,"RAAS"],
                       normal=y[y[,"Sample type"]=="Normal","RAAS"],
                       sample=y[y[,"Sample type"]=="Tumor","sample"],
                       BP=    y[y[,"Sample type"]=="Tumor","BP"],
                       SAAP=  y[y[,"Sample type"]=="Tumor","SAAP"],
                       Dataset=  y[y[,"Sample type"]=="Tumor","Dataset"])
        else NULL
    })
    rpaired <- do.call(rbind, rpaired)
    if ( !is.null(rpaired) ) {
        if ( ncol(rpaired)!=6 )
            rpaired <- NULL
    }
    
    rpaired
})
rp <- do.call(rbind, rp)
rownames(rp) <- NULL

## SHOULD MOSTLY BE AVAILABLE IN PDAC
table(tmtf$Dataset[paste(tmtf$BP,tmtf$SAAP)%in%paste(rp$BP, rp$SAAP)])
## NOTE: BP but not SAAP present in many other Datasets
table(tmtf$Dataset[tmtf$BP%in%rp$BP])
table(tmtf$Dataset[tmtf$SAAP%in%rp$SAAP])

table(bdat$Dataset[bdat$BP%in%pdat$BP])
table(bdat$Dataset[bdat$SAAP%in%pdat$SAAP])

plotdev(file.path(pdfig.path,"TV_RAAS_tumor_normal"),
        width=3, height=3)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(rp[,"normal"], rp[,"tumor"], density=TRUE,
        xlab=expression(log[10](RAAS["normal"])),
        ylab=expression(log[10](RAAS["tumor"])),
        main="paired BP/SAAP per patient")
abline(a=0,b=1, col=5)
dev.off()

plotdev(file.path(pdfig.path,"TV_RAAS_tumor_normal_lg2fc"),
        width=3, height=3)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
hist(log2(10^rp[,"tumor"]/10^rp[,"normal"]), #density=TRUE,
        xlab=expression(log[2](RAAS["tumor"]/RAAS["normal"])),
        main="paired BP/SAAP per patient", breaks=100)
abline(v=0, col=2)
dev.off()

## 191 unique BP/SAAP
length(unique(paste(rp$BP, rp$SAAP)))
