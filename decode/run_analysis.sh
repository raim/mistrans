#!/bin/bash


# Ensure R is installed
if ! command -v R &> /dev/null; then
    echo "R needs to be isntalled  and in your PATH."
    exit 1
fi

## Note, that all scripts also load functions defined in
## raas_utils.R

## Initialize BP/SAAP and RAAS TABLES
## NOTE, that this script is called redundantly from each script below.
## It loads all required libraries, and defines all input files and output
## data paths. This script also generates the
## "Supplemental_Data_7.SAAP_coordinates.tsv" as well as tables for
## unique protein sites used in the random forest model by Andrew Leduc.
$myR --vanilla <  raas_init.R &> $DECDATA/log/raas_init.txt  

### CALCULATE RAAS PROFILES and GENERATE PUBLICATION FIGURES
R --vanilla < codons.R &> $DECDATA/log/codons.txt  # FIGURE 2
R --vanilla < aminoacids.R &> $DECDATA/log/aminoacids.txt  # FIGURE 3
R --vanilla < motifs.R &> $DECDATA/log/motifs.txt  # FIGURE 4, FIGURE 5
R --vanilla < function.R &> $DECDATA/log/function.txt # FIGURE 4
R --vanilla < kraq.R &> $DECDATA/log/kraq.txt  # Ext. Data Figure 8: motifs
R --vanilla < proteins.R &> $DECDATA/log/proteins.txt # FIGURE 5
R --vanilla < structure.R &> $DECDATA/log/structure.txt  # Ext. Data Figure 10
R --vanilla < rna.R &> $DECDATA/log/rna.txt # Suppl. Figure 2

## TODO: add 1d protein plots and chimeraX generating code
