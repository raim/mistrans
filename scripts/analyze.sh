#!/bin/bash

## SCRIPTS to analyze location and function of mistranslation
## events, derived by Shiri Tsour from the slavovlab.

## INPUT GENOME DATA IS GENERATED
## by from genomeBrowser/mammary/setup.sh

## map each peptide to position in protein and transcript
R --vanilla < scripts/map_peptides2.R > log/map.txt
R --vanilla < scripts/saap_analysis.R 
## functional enrichment of SAAP-harboring proteins
R --vanilla < scripts/saap_function.R



cp -a ~/Documents/sejour23_fig1.jpg .

pandoc -t beamer tsour23pre.md --filter pandoc-citeproc -o tsour23pre.pdf

## TODO: instead, generate a blastdb of all transcripts and proteins,
## blast base peptides within those proteins and record locations
