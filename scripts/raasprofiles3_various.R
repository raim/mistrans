
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

vfig.path <- file.path(fig.path,"various")
dir.create(vfig.path, showWarnings=FALSE)

## TODO: raas profiles for alternative row/classes
## * motifs vs. disordered,
## * 
