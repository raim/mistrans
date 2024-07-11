
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)
library(gprofiler2)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

## additional data
goslim.file  <- file.path(mam.path,"processedData","goslim.tsv")

pfig.path <- file.path(fig.path,"PDAC")
dir.create(pfig.path, showWarnings=FALSE)

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
tmtv <- tmtp[tmtp$BP%in%tvdat$BP,]

## COUNTS in PDAC vs OTHERS

## BP/SAAP unique to PDAC and T->V?
table(tmtv$Dataset) # T->V in PDAC: not in other Datasets!
table(tmtpd$Dataset) # BP in PDAC: many also present in other

## table
nrow(pdat) # 1650 unique BP/SAAP in PDAC
nrow(tvdat) # 192 unique BP/SAAP in T->V


sum( pdat$BP %in% tmtwop$BP) #  1269 BP also occur in other data sets
sum(!pdat$BP %in% tmtwop$BP) # 381 DO NOT occur in other data sets

length(unique(tvdat$name)) # 147 unique proteins with T-->V

dup <- which(duplicated(tvdat$gene))[1]
tvdat$SAAP[which(tvdat$gene==tvdat$gene[dup])] # ACTG2 contains 3 distinct T->V



### FUNCTIONS of T->V containing proteins in PDAC.

sort(table(tvdat$name))


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

