
library(Biostrings) # for blosum62
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")

pfig.path <- file.path(fig.path,"rfam")
dir.create(pfig.path, showWarnings=FALSE)

corW <- corH <- 2.5
pmai <- c(.5,.5,.25,.25)
pmpg <- c(1.3,.3,0)

## additional data
secis.path <- file.path(out.path,"secis_elements.out")
rfam.path <- file.path(mam.path,"processedData",
                       "Homo_sapiens.GRCh38.cdna.large.rfam_cut_ga.out")

library(cmchainer)
source("~/programs/gIScanner/R/cmchainer.R") # TODO: mv to segmenTools


## SECIS ELEMENTS
rfs <- parseHomHits(files=secis.path, type="rfam11")

rfs$transcript <- sub("\\.","",rfs$target)

## NO OVERLAP OF SECIS TARGETS WITH AAS PROTEINS
rfs$transcript%in%bdat$transcript

## ALL RFAM 
rfs <- parseHomHits(files=rfam.path, type="cmscan")
rfs$transcript <- sub("\\.","",rfs$query)

## NO OVERLAP OF ANY RFAM WITH AAS PROTEINS
sum(rfs$transcript%in%bdat$transcript)
