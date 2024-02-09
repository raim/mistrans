#!/bin/bash

## SCRIPTS to analyze location and function of mistranslation
## events, derived by Shiri Tsour from the slavovlab.

## INPUT GENOME DATA IS GENERATED
## by from genomeBrowser/mammary/setup.sh

## map each peptide to position in protein and transcript
R --vanilla < scripts/map_peptides2.R > log/map.txt

## copy required files to a new directory
## much smaller for easier transfer
if [ false ]; then
    cut -f 12 ~/data/mistrans/processedData/saap_mapped.tsv | sed '/^$/d' | tail -n +2 > ~/data/mistrans/processedData/mapped_ensembl_proteins.dat
    while read p; do
	echo "$p" 
	find ~/data/mammary/processedData/iupred3 -name "${p}*_iupred3.tsv.gz" \
	     >>  trasnfer_files.txt
    done < ~/data/mistrans/processedData/mapped_ensembl_proteins.dat
    cat trasnfer_files.txt |sort |uniq | xargs -I % cp % iupred3_selected/
fi

## GENERATE SUBSETS OF PROTEIN/TRANSCRIPT FASTA AND
## RETRIEVE GENOME LEVEL DATA FROM BIGWIG FILES
 

 
R --vanilla < scripts/saap_analysis.R 
## functional enrichment of SAAP-harboring proteins
R --vanilla < scripts/saap_function.R



cp -a ~/Documents/sejour23_fig1.jpg .

pandoc -t beamer tsour23pre.md --filter pandoc-citeproc -o tsour23pre.pdf

## TODO: instead, generate a blastdb of all transcripts and proteins,
## blast base peptides within those proteins and record locations
