
## all packages required for the amino acid substitution analysis

## from CRAN
install.packages("viridis")
install.packages("readxl")

## from github
## https://github.com/raim/segmenTools (release RAAS_preprint).
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("raim/segmenTools@RAAS_preprint")

## for motifs.R 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DiffLogo")

## for kraq.R
BiocManager::install("Biostrings")
BiocManager::install("seqinr")
BiocManager::install("MatrixGenerics")
BiocManager::install("universalmotif")

## for rna.R
## https://josephcrispell.github.io/projects/basicplotter (2024-08-24),
remotes::install_github("JosephCrispell/basicPlotteR")


## for model.R: random forest model
install.packages("dplyr")
install.packages("ggplot2")
install.packages("xgboost")
