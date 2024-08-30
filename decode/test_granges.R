
library(rtracklayer)

## attempt to get the correct path for input data;
## these should be downloaded or generated as outlined
## in README.md
proj.path <- file.path(Sys.getenv("DECODE")) 
if ( proj.path=="" ) # author's local path
    proj.path <- "/home/raim/data/mistrans"
SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")

add.path <- file.path(proj.path,"additionalData")
mam.path <- file.path("/home/raim/data/mammary")

bw.file <- file.path(mam.path, "originalData",
                       "gwipsvizRiboseq.bw")
gff3.file <- file.path(mam.path, "originalData",
                       "Homo_sapiens.GRCh38.110.gff3.gz")


gff <- import(gff3.file, format = "gff3")
cds <- gff[gff$type=="CDS",] # only CDS



bw <- rtracklayer::import.bw(bw.file) 
seqlevelsStyle(bw) <- "NCBI"

ov <- findOverlaps(bw, gff)


##cds.val <- mcols(bw)[subjectHits(ov), ]
##
##cds_with_bed <- cbind(as.data.frame(cds_gr[queryHits(overlaps), ]),
##                      as.data.frame(cds_bed_values))
