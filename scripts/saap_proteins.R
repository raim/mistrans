
## read protein fasta, iupred3, s4pred, and AAS,
## and generate per protein plots


source("~/work/mistrans/scripts/saap_utils.R")

library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)


#### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")

mam.path <- "~/data/mammary/"

in.file <- file.path(out.path,"saap_mapped3.tsv")
## protein fasta
pfasta <- file.path(dat.path,"all_proteins.fa")
## coding region fasta
tfasta <- file.path(mam.path,"processedData","coding.fa")

## protein-transcript map
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
## CDS structure
cdsmap <- file.path(mam.path,"processedData","protein_cds_structure.dat")
cdspos <- file.path(mam.path,"processedData","protein_cds_coordinates.tsv")


## structure predictions
## s4pred
s4pred <- file.path(mam.path,"processedData",
                    "Homo_sapiens.GRCh38.pep.large_s4pred.fas.gz")
## iupred3/anchor2
iupred <- file.path(mam.path,"processedData","iupred3")

## pfam hits
pfam <- file.path(mam.path,"processedData",
                  "Homo_sapiens.GRCh38.pep.large_annotations.csv")

## genome feature file
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

## PARSE & FILTER DATA
dat <- read.delim(in.file)
## remove SAAP/BP w/o protein
rm <- dat$ensembl==""
dat <- dat[!rm,]
aasl <- split(dat, dat$ensembl) 

## TMT Level RAAS Values
tmtf <- read.delim(tmt.file)

## s4pred
s4p <- readFASTA(s4pred, grepID=TRUE)
## skip version number - TODO: use version numbers!
names(s4p) <- sub("\\.[0-9]+", "", names(s4p))

## iupred3
iufiles <- list.files(pattern=paste0(".*iupred3.tsv.gz"), path=iupred)
names(iufiles) <- sub("\\..*","", iufiles)


## pfam
## fixed width format, where \s was replaced with ; by sed
## TODO: load pfam description
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

pfl <- split(pfm, pfm$query)
names(pfl) <- sub("\\.[0-9]+", "", names(pfl))

## NOTE: three sets of from/to: hmm, ali, env;
## env may overlap between multiple - could be fused!
