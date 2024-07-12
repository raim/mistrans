
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

## PDAC ONLY
tmtp <- tmtf[tmtf$Dataset=="PDAC",]
pdat <- bdat[bdat$BP%in%tmtp$BP,]

## W/O PDAC
tmtwop <- tmtf[tmtf$Dataset!="PDAC",]

## all RAAS values of BP that appear in PDAC
tmtpd <- tmtf[tmtf$BP%in%pdat$BP,]

## T->V in PDAC
tvdat <- pdat[pdat$fromto=="T:V",]

## T->V in Precurser Level RAAS
tmtv <- tmtf[tmtf$fromto=="T:V",]

## COUNTS in PDAC vs OTHERS

## BP/SAAP unique to PDAC and T->V?

table(tmtv$Dataset) # T->V mostly in PDAC

table(tmtpd$Dataset) # BP in PDAC: many also present in other

## table
nrow(bdat) # 8402 unique BP/SAAP 
nrow(pdat) # 1650 unique BP/SAAP in PDAC
nrow(tvdat)#  192 unique BP/SAAP in T->V

length(unique(bdat$BP)) # 5303 unique BP
length(unique(pdat$BP)) #  621 unique BP in PDAC



sum( pdat$BP %in% tmtwop$BP) # 1269 BP also occur in other data sets
sum(!pdat$BP %in% tmtwop$BP) #  381 DO NOT occur in other data sets

length(unique(tvdat$name)) # 147 unique proteins with T-->V

dup <- which(duplicated(tvdat$gene))[1]
tvdat$SAAP[which(tvdat$gene==tvdat$gene[dup])] # ACTG2 contains 3 distinct T->V

## 7080 unique sites in all
length(unique(tmtf$unique.site))
## 786 unique sites in PDAC
length(unique(tmtf$unique.site[which(tmtf$Dataset=="PDAC")]))
## 188 unique T->V protein sites in PDAC
length(unique(tmtf$unique.site[which(tmtf$Dataset=="PDAC" & tmtf$fromto=="T:V")]))
## 25 unique  T->V protein sites in other tissues
length(unique(tmtf$unique.site[which(tmtf$Dataset!="PDAC" & tmtf$fromto=="T:V")]))

## T->V in other tissues
table(tmtf$BP[which(tmtf$Dataset!="PDAC" & tmtf$fromto=="T:V")],
      tmtf$TMT.Tissue[which(tmtf$Dataset!="PDAC" & tmtf$fromto=="T:V")])

## TODO: summary plot T->V in PDAC vs. OTHERS

bdat$fromto <- paste0(bdat$from,":",bdat$to)

ft <- unique(bdat$fromto)
aas.srt <- unique(c(grep("T:", ft, value=TRUE),grep(":V", ft, value=TRUE)))
aas.srt <- c("T:V", "Q:G")

dsets <- bdat$Datasets
##dsets[dsets=="Healthy"] <- bdat$Tissues[dsets=="Healthy"]
nsets <- stringr::str_count(dsets, ";")+1
##dsets[nsets > 1] <- ">1"
dsets[nsets > 2] <- ">2"
dsets[nsets > 3] <- ">3"
ds.srt <- unique(dsets)
ds.srt <- ds.srt[order(nchar(ds.srt))]

ovl <- clusterCluster(cl2=dsets, cl1=bdat$fromto,
                      cl1.srt=aas.srt, cl2.srt=ds.srt)
present <- colnames(ovl$overlap)[apply(ovl$overlap,2,sum)>0]
ovl <- sortOverlaps(ovl, axis=1, srt=present)

plotdev(file.path(pdfig.path,"TV_per_Dataset"),
        width=length(present)*.3+1, height=length(aas.srt)*.2+2)
par(mai=c(1.5,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovl, show.total=TRUE, xlab="", ylab="", p.min=p.min, p.txt=p.txt)
figlabel("count of unique BP/SAAP in Datasets", pos="bottom")
dev.off()

## same for tissues
dsets <- bdat$Tissues
##dsets[dsets=="Healthy"] <- bdat$Tissues[dsets=="Healthy"]
nsets <- stringr::str_count(dsets, ";")+1
##dsets[nsets > 1] <- ">1"
##dsets[nsets > 2] <- ">2"
dsets[nsets > 3] <- ">3"
ds.srt <- unique(dsets)
ds.srt <- ds.srt[order(nchar(ds.srt))]

ovl <- clusterCluster(cl2=dsets, cl1=bdat$fromto,
                      cl1.srt=aas.srt, cl2.srt=ds.srt)
## NOTE: only T:V
present <- colnames(ovl$overlap)[ovl$overlap["T:V",]>0]
ovl <- sortOverlaps(ovl, axis=1, srt=present)

plotdev(file.path(pdfig.path,"TV_per_Tissue"),
        width=length(present)*.3+1, height=length(aas.srt)*.2+2.5)
par(mai=c(2,.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovl, show.total=TRUE, xlab="", ylab="", p.min=p.min, p.txt=p.txt)
figlabel("count of unique BP/SAAP in Tissues", pos="bottom")
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

pats <- pats[pats[,"Keep SAAP"],]

patv <- pats[pats$AAS=="T to V" & pats$Dataset=="PDAC",]

plotdev(file.path(pdfig.path,"TV_RAAS_per_tumor"),
        width=corW, height=3)
par(mai=c(.5,.5,.5,.05), mgp=c(1.3,.3,0), tcl=-.25)
bp <- boxplot(patv$RAAS ~ patv[,"Sample type"], ylab=xl.raas,
              xlab="T->V in PDAC; all patients")
axis(3, at=1:2, labels=bp$n, las=2)
dev.off()
