
## TEST GENOMIC AND TRANSCRIPT COORDINATES OF AAS from map_peptides.R

dat.path <- "/home/raim/data/"
mam.path <- file.path(dat.path, "mammary")
mís.path <- file.path(dat.path, "mistrans")


## AAS mapping by map_peptides.R
in.file <- file.path(mis.path, "processedData", "saap_mapped.tsv")

## protein fasta
pfas <- file.path(mis.path, "processedData", "all_proteins.fa")
## transcript fasta
tfas <- file.path(mam.path, "originalData", "Homo_sapiens.GRCh38.cdna.all.fa.gz")
## chromosome fasta
gfas <- file.path(mam.path, "originalData", "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
