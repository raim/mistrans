#!/bin/bash


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
$myR --vanilla < codons.R &> $DECDATA/log/codons.txt  # FIGURE 2
$myR --vanilla < aminoacids.R &> $DECDATA/log/aminoacids.txt  # FIGURE 3
$myR --vanilla < motifs.R &> $DECDATA/log/motifs.txt  # FIGURE 4, FIGURE 5
$myR --vanilla < function.R &> $DECDATA/log/function.txt # FIGURE 4
$myR --vanilla < kraq.R &> $DECDATA/log/kraq.txt  # Ext. Data Figure 8: motifs
$myR --vanilla < proteins.R &> $DECDATA/log/proteins.txt # FIGURE 5
$myR --vanilla < structure.R &> $DECDATA/log/structure.txt  # Ext. Data Figure 10
$myR --vanilla < rna.R &> $DECDATA/log/rna.txt # Suppl. Figure 2

## TODO: add 1d protein plots and chimeraX generating code
