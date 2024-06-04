
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

##  additional data
## pfam hits
pfam <- file.path(mam.path,"processedData",
                  "Homo_sapiens.GRCh38.pep.large_annotations.csv")
## pfam download
pfdlf <- file.path(mam.path,"originalData", "9606.tsv.gz")
## pfam clans
clans <- file.path(mam.path,"originalData", "pfam", "Pfam-A.clans.tsv.gz")


## pfam
## fixed width format, where \s was replaced with ; by sed
## NOTE: three sets of from/to: hmm, ali, env;
## env may overlap between multiple - could be fused!
pfm <- read.csv(file=pfam, sep=";", fill=FALSE, header=FALSE,comment.char="#")
pfmh <- c(
    "target", "accession" , "tlen",            
    "query"   , "accession" , "qlen",                 
    "E-value"      , "score"     , "bias",                 
    "#"            , "of"        , "c-Evalue",             
    "i-Evalue"     , "score"     , "bias",            
    "from"         , "to"        , "from",                 
    "to"           , "FROM"      , "TO",                   
    "acc")#          , "description")
colnames(pfm) <- pfmh

## pfam clans: used to collapse domains below!
pclan <- read.delim(clans, row.names=1 , header=FALSE)
pfm$clan <- pclan[sub("\\..*","", pfm$accession),"V3"]
pfm$clan[pfm$clan==""] <- pfm$target[pfm$clan==""]

use.pclan <- TRUE# FALSE
if ( !use.pclan )
    pfm$clan <- pfm$target

pfl <- split(pfm, pfm$query)
names(pfl) <- sub("\\.[0-9]+", "", names(pfl))


### START ANALYSIS


#### PROTEIN DOMAINS

## MEDIAN SITE AND PROTEIN RAAS

### median raas per unique mane protein site

sitl <- split(tmtf$RAAS, paste(tmtf$mane, tmtf$pos))
site <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## mod. column names and add protein and site info
colnames(site) <- paste0("RAAS.", colnames(site))
site$mane <- sub(" .*", "", rownames(site))
site$pos <- as.numeric(sub(".* ", "", rownames(site)))

## add gene names
site$name <- ens2nam [site$mane]

## add uniprot id
site$uniprot <- unlist(lapply(ens2u[site$mane], paste, collapse=";"))

## RAAS COLOR - NOTE: broad codon definition of RAAS colors
site$RAAS.color <- num2col(site$RAAS.median,
                           limits=c(RAAS.MIN, RAAS.MAX), colf=arno)


### order PFAM clans in proteins
pfl <- pfl[site$mane]


## for each pfam domain collect all RAAS and do dot profile
pfml <- split(pfm, pfm$clan)

dom <- pfml[[1]]

## TODO: merge overlapping domains?

prts <- unique(dom$query)
 
pfpr <- split(dom, sub("\\..*", "", dom$query))

## all positions that belong to a given domain
pfdm <- lapply(pfpr,
               function(x) 
                   c(unlist(apply(x, 1, function(y) y["FROM"]:y["TO"]))))
pfdm <- lapply(pfdm, unique)

## get all RAAS for these proteins
tmtp <- tmtf[tmtf$ensembl%in%names(pfdm),]
tmtl <- split(tmtp, tmtp$ensembl)
