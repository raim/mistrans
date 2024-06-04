
library(Biostrings) # for blosum62
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
source("~/work/mistrans/scripts/raasprofiles3_init.R")

dfig.path <- file.path(fig.path,"domains")
dir.create(dfig.path, showWarnings=FALSE)


tmtf$clan.ebi[tmtf$clan.ebi==""] <- "none" 
tmtf$pfam[tmtf$pfam==""] <- "none" 

## TODO: replace by clan/pfam names

ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="clan", cols="Dataset",
                   bg=TRUE, value="RAAS", ##row.srt=mid,
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0)

ovc <- sortOverlaps(ovw, p.min=1e-3, cut=TRUE, sign=1)

plotProfiles(ovc,
             fname=file.path(dfig.path,paste0("pfams_",SETID,"_all")),
             mai=c(.8,1,.6,.6), ttcols=ttcols, value="median",
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols, ffam="monospace")
