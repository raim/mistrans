#!/bin/bash

## ANALYSIS OF AMINO ACID SUBSTITUTION SITES
##  as provided by Shiri Tsour.

## OUTPUT DATA PATHS
## SET THIS PATH TO WHERE YOU WANT TO SAVE
## INPUT (supplemental data) and OUTPUT (figures, tables, log files)
## TODO: parse DECDATA in the R scripts
export DECDATA=/home/raim/data/decode

## generate output paths
mkdir -p $DECDATA/log
mkdir -p $DECDATA/figures
mkdir -p $DECDATA/processedData
mkdir -p $DECDATA/originalData

## BLAST 2.15.0+ is REQUIRED, set the path to executables here
blastdir=${HOME}/programs/ncbi-blast-2.15.0+/bin

### DATA COLLECTION

## MANUALLY DOWNLOAD SUPPLEMENTAL DATA TABLES 2-4 AND STORE IN
## $DECDATA/originalData

## MANUALLY DOWNLOAD analysis_pipeline.zip AND STORE IN
## $DECDATA/processedData

## 1) CONVERT SAAP/BP INPUT DATA (by Shiri Tsour) to text files
##Supplemental_Data_2.SAAP_proteins.xlsx == All_SAAP_protein_filter_df.txt,
##Supplemental_Data_3.SAAP_precursor_quant.xlsx == All_SAAP_TMTlevel_quant_df.txt,
##Supplemental_Data_4.SAAP_reporter_quant.xlsx == All_SAAP_patient_level_quant_df.txt
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "locale=C format=automatic separator='	' eol=unix" $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.xlsx $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.txt 
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "locale=C format=automatic separator='	' eol=unix" $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.xlsx $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.txt 
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "locale=C format=automatic separator='	' eol=unix" $DECDATA/originalData/Supplemental_Data_2.SAAP_proteins.xlsx $DECDATA/originalData/Supplemental_Data_2.SAAP_proteins.txt
## gzip to save disk space
gzip $DECDATA/originalData/Supplemental_Data_2.SAAP_proteins.txt $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.txt $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.txt

## 2) COLLECT BP SEQUENCES
gunzip -c $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.txt.gz | cut -f 4 | \
    sort | uniq | grep -v "^BP$" | grep -v -e '^$' > $DECDATA/processedData/tmp.txt
gunzip -c $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.txt.gz | cut -f 4 | \
    sort | uniq | grep -v "^BP$" | grep -v -e '^$' > $DECDATA/processedData/tmp2.txt
cat $DECDATA/processedData/tmp.txt $DECDATA/processedData/tmp2.txt | \
    sort |uniq | awk '{print ">" $0 ORS $0}' - \
		     > $DECDATA/processedData/unique_bp.fas
rm -f $DECDATA/processedData/tmp.txt $DECDATA/processedData/tmp2.txt

## 3) GET Ensembl PROTEIN SEQUENCES
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz $DECDATA/originalData/

## 4) ADD PATIENT-SPECIFIC SINGLE AMINO REPLACEMENTS 
#All_SAAP_protein_filter_df.txt
R --vanilla < get_mutated_proteins.R  &> ${DECDATA}/log/mutated_proteins.txt 

## 5) BLAST BP IN Ensembl+MUTATION PROTEINS

$blastdir/makeblastdb -in ${DECDATA}/processedData/all_proteins.fa -parse_seqids -title "ensembl hg38 proteins" -dbtype prot  &> ${DECDATA}/log/blast_database.txt
## NOTE: blast hits are pre-filtered for query length=alignment length and >75% identity
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${DECDATA}/processedData/unique_bp.fas -db ${DECDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore" 2> ${DECDATA}/log/unique_bp_blast.txt | awk '{if($5==$6 && $3>75) print}'  |grep -v "^#" > ${DECDATA}/processedData/unique_bp_blast.tsv 


## 6) SELECT BEST-MATCHING PROTEIN HIT FOR EACH BP
R --vanilla < get_protein_match.R &> ${DECDATA}/log/protein_match.txt

### DATA ANALYSIS AND PLOTS

## Note, that all scripts also load functions defined in
## raas_utils.R

## Initialize BP/SAAP and RAAS TABLES
## NOTE, that this script is called redundantly from each script below.
## It loads all required libraries, and defines all input files and output
## data paths. This script also generates the
## "Supplemental_Data_7.SAAP_coordinates.tsv" as well as tables for
## unique protein sites used in the random forest model by Andrew Leduc.
R --vanilla <  raas_init.R &> ${DECDATA}/log/raas_init.txt  

### CALCULATE RAAS PROFILES and GENERATE PUBLICATION FIGURES
R --vanilla < codons.R &> ${DECDATA}/log/codons.txt  # FIGURE 2
R --vanilla < aminoacids.R &> ${DECDATA}/log/aminoacids.txt  # FIGURE 3
R --vanilla < motifs.R &> ${DECDATA}/log/motifs.txt  # FIGURE 4, FIGURE 5
R --vanilla < function.R &> ${DECDATA}/log/function.txt # FIGURE 4
R --vanilla < kraq.R &> ${DECDATA}/log/kraq.txt  # Ext. Data Figure 8: motifs
R --vanilla < proteins.R &> ${DECDATA}/log/proteins.txt # FIGURE 5
R --vanilla < structure.R &> ${DECDATA}/log/structure.txt  # Ext. Data Figure 10
R --vanilla < rna.R &> ${DECDATA}/log/rna.txt # Suppl. Figure 2

## TODO: add 1d protein plots and chimeraX generating code
