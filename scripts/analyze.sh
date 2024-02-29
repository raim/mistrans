#!/bin/bash


export MISDATA=${MYDATA}/mistrans
export MAMDATA=${MYDATA}/mammary
SRC=$GENBRO/data/mammary
THIS=${HOME}/work/mistrans

## SCRIPTS to analyze location and function of mistranslation
## events, derived by Shiri Tsour from the slavovlab.

## INPUT GENOME DATA IS GENERATED
## by from genomeBrowser/mammary/setup.sh

## analyze data structure, different number of replicates per unique SAAP
R --vanilla < ${THIS}/scripts/saap_means.R > log/means.txt


## ANNOTATE ALL BP/SAAP

## 1) collect all BP and SAAP/BP pairs

## BP as fasta for blast
cut -f 5 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 5 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort |uniq | awk '{print ">" $0 ORS $0}' - > ${MISDATA}/processedData/unique_bp.fas

## BP/SAAP as simple table, basis for search in proteins
cut -f 4,5 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 4,5 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort | uniq > ${MISDATA}/processedData/unique_saap.tsv

## 2) collect all proteins tagged with mutations and add these to protein DB;
##    generates ${MISDATA}/processedData/all_proteins.fa 
R --vanilla < ${THIS}/scripts/get_mutated_proteins.R

## 3) blast all BP against ensembl proteins + mutations
blastdir=${HOME}/programs/ncbi-blast-2.15.0+/bin
## generate local blastdb 
$blastdir/makeblastdb -in ${MISDATA}/processedData/all_proteins.fa -parse_seqids -title "ensembl hg38 proteins" -dbtype prot
## blast - filter full length hit alignment length=query length,
## and at least 75% identity with awk.
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/unique_bp.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  | awk '{if($5==$6 && $3>75) print}'  |grep -v "^#" > ${MISDATA}/processedData/unique_bp.tsv

## 4) find best matching protein
R --vanilla < ${THIS}/scripts/get_protein_match.R > log/match.txt


## map each peptide to position in protein and transcript
R --vanilla < ${THIS}/scripts/map_peptides3.R > log/map3.txt


## GENERATE SUBSETS OF PROTEIN/TRANSCRIPT FASTA, and
## RETRIEVE GENOME LEVEL DATA FROM BIGWIG FILES

### CALCULATE RAAS PROFILES
## TODO: log files
## with Albumin
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=FALSE/' scripts/raasprofiles.R | R --vanilla &
## without Albumin
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=FALSE/' scripts/raasprofiles.R | R --vanilla &
## unique SAAP without Albumin 
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=TRUE/' scripts/raasprofiles.R | R --vanilla &
## unique SAAP with Albumin 
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=TRUE/' scripts/raasprofiles.R | R --vanilla &

## same for healthy tissues
## with Albumin
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=FALSE/' scripts/raasprofiles.R | R --vanilla &
## without Albumin
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=FALSE/' scripts/raasprofiles.R | R --vanilla &
## unique SAAP without Albumin 
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=TRUE/' scripts/raasprofiles.R | R --vanilla &
## unique SAAP with Albumin 
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=TRUE/' scripts/raasprofiles.R | R --vanilla



R --vanilla < ${THIS}/scripts/saap_analysis.R  > log/analysis.txt
## functional enrichment of SAAP-harboring proteins
R --vanilla < ${THIS}/scripts/saap_function.R


### GENERATE RESULTS SLIDES AND PAPER FIGURES


cp -a ~/Documents/sejour23_fig1.jpg .

pandoc -t beamer tsour23pre.md --filter pandoc-citeproc -o tsour23pre.pdf

## TODO: instead, generate a blastdb of all transcripts and proteins,
## blast base peptides within those proteins and record locations
