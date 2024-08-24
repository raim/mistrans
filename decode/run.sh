#!/bin/bash

## OUTPUT DATA PATHS
## SET THIS PATH TO WHERE YOU WANT OUTPUT
export DECDATA=${MYDATA}/decode_results

## generate output paths
mkdir -p $DECDATA/log
mkdir -p $DECDATA/figures
mkdir -p $DECDATA/originalData
mkdir -p $DECDATA/processedData

## initialize RAAS TABLES
## NOTE that this script is called redundantly
R --vanilla <  raas_init.R &> ${DECDATA}/log/raas_init.txt  

### CALCULATE RAAS PROFILES
R --vanilla <  codons.R &> ${DECDATA}/log/codons.txt  # FIGURE 2
R --vanilla <  aminoacids.R &> ${DECDATA}/log/aminoacids.txt  # FIGURE 3
R --vanilla <  motifs.R &> ${DECDATA}/log/motifs.txt  # FIGURE 4
R --vanilla <  kraq.R &> ${DECDATA}/log/kraq.txt  # EXTENDED DATA FIGURE 8: motif selection
R --vanilla <  structure.R &> ${DECDATA}/log/structure.txt 
R --vanilla <  function.R &> ${DECDATA}/log/function.txt 
R --vanilla <  proteins.R &> ${DECDATA}/log/proteins.txt 
R --vanilla <  model.R &> ${DECDATA}/log/model.txt 
R --vanilla <  rna.R &> ${DECDATA}/log/rna.txt
