#!/bin/bash


export MISDATA=${MYDATA}/mistrans
export MAMDATA=${MYDATA}/mammary
SRC=$GENBRO/data/mammary
THIS=${HOME}/work/mistrans

## SCRIPTS to analyze location and function of mistranslation
## events, derived by Shiri Tsour from the slavovlab.

## INPUT GENOME DATA IS GENERATED
## by from genomeBrowser/mammary/setup.sh

## additional data
wget https://degronopedia.com/degronopedia/download/data/DEGRONOPEDIA_degron_dataset.xlsx -P $MISDATA/originalData/

## analyze data structure, different number of replicates per unique SAAP
R --vanilla < ${THIS}/scripts/saap_means.R > log/means.txt


## ANNOTATE ALL BP/SAAP

## 1) collect all BP and SAAP/BP pairs

## BP as fasta for blast - TODO: remove header!
cut -f 5 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 5 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort |uniq | awk '{print ">" $0 ORS $0}' - > ${MISDATA}/processedData/unique_bp.fas

## how many? 7991 unique BP
grep -n ">"  ${MISDATA}/processedData/unique_bp.fas |wc -l

## SAAP as fasta for blast - TODO: remove header!
cut -f 4 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 4 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort |uniq | awk '{print ">" $0 ORS $0}' - > ${MISDATA}/processedData/unique_saap.fas

## how many? 15061 unique SAAP
grep -n ">"  ${MISDATA}/processedData/unique_saap.fas |wc -l

## BP/SAAP as simple table, basis for search in proteins
cut -f 4,5 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 4,5 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort | uniq > ${MISDATA}/processedData/unique_saap.tsv

## how many? 15669 unique SAAP/BP, 15061 unique SAAP
wc -l  ${MISDATA}/processedData/unique_saap.tsv
cut -f 1  ${MISDATA}/processedData/unique_saap.tsv |sort|uniq|wc -l

## 2) collect all proteins tagged with mutations and add these to protein DB;
##    generates ${MISDATA}/processedData/all_proteins.fa 
R --vanilla < ${THIS}/scripts/get_mutated_proteins.R

## 3) blast all BP against ensembl proteins + mutations
blastdir=${HOME}/programs/ncbi-blast-2.15.0+/bin
## generate local blastdb 
$blastdir/makeblastdb -in ${MISDATA}/processedData/all_proteins.fa -parse_seqids -title "ensembl hg38 proteins" -dbtype prot
## blast - filter full length hit alignment length=query length,
## and at least 75% identity with awk.
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/unique_bp.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  | awk '{if($5==$6 && $3>75) print}'  |grep -v "^#" > ${MISDATA}/processedData/unique_bp_blast.tsv

## 3.A) blast SP against ensembl+mutations
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/unique_saap.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  | awk '{if($5==$6 && $3>75) print}'  |grep -v "^#" > ${MISDATA}/processedData/unique_saap_blast.tsv

## some statistics on blast
## TODO: expand this QC analysis a bit,
## do SAAP have the expected mismatches?
R --vanilla < ${THIS}/scripts/saap_blast_stats.R

### 4) COLLECT DATA FOR ALL BP/SAAP

## 4.A) find best matching protein
##    GENERATES ${MISDATA}/processedData/bp_mapped.tsv
R --vanilla < ${THIS}/scripts/get_protein_match.R > ${THIS}/log/match.txt

## 4.B) map each peptide to it's positions in proteins and transcripts, and
## add a variety of collected information on protein structure, e.g.
## iupred3, anchor2, s4pred, codon, ...
##    GENERATES ${MISDATA}/processedData/saap_mapped3.tsv, and
##    QC figures in ${MISDATA}/figures/saap_mapping3/
R --vanilla < ${THIS}/scripts/map_peptides3.R > ${THIS}/log/map3.txt

### 5) ANALYSIS

## PROTEINS: protein level analysis
## TODO: revisit previous saap_function.R and saap_analysis.R
R --vanilla < ${THIS}/scripts/saap_proteins.R > ${THIS}/log/proteins.txt

## export sequence context of AAS
R --vanilla < ${THIS}/scripts/extract_aas_context.R > ${THIS}/log/context.txt

## TODO: run deep learning script here! currently done interactively
## from above script, but instead sequence input variation should be done
## in python script.

## MOTIFS: get and analyze sequences surrounding the ASS
## TODO: split this script, use output from above extract_aas_context.R
R --vanilla < ${THIS}/scripts/get_aas_context.R > ${THIS}/log/context.txt

## kplogo
## TODO: run over all from:to classes and find a way to collect and plot
## results by RAAS!
~/programs/kpLogo/bin/kpLogo seqcontext_all_.fa  -alphabet protein -o kplogo/all
~/programs/kpLogo/bin/kpLogo seqcontext_fromto_Q:G.fa -alphabet protein -o kplogo/QG
~/programs/kpLogo/bin/kpLogo seqcontext_fromto_T:V.fa -alphabet protein -o kplogo/TV
~/programs/kpLogo/bin/kpLogo seqcontext_methionine_TRUE.fa -alphabet protein -o kplogo/M
~/programs/kpLogo/bin/kpLogo seqcontext_tryptophane_TRUE.fa -alphabet protein -o kplogo/W


## GENERATE SUBSETS OF PROTEIN/TRANSCRIPT FASTA, and
## RETRIEVE GENOME LEVEL DATA FROM BIGWIG FILES

### CALCULATE RAAS PROFILES
## TODO: log files
## with Albumin
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3.R | R --vanilla &
## without Albumin
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3.R | R --vanilla &
## unique SAAP without Albumin 
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3.R | R --vanilla &
## unique SAAP with Albumin 
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3.R | R --vanilla &

## same for healthy tissues
## with Albumin
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3.R | R --vanilla &
## without Albumin
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3.R | R --vanilla &
## unique SAAP without Albumin 
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3.R | R --vanilla &
## unique SAAP with Albumin 
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3.R | R --vanilla


## OLD and OBSOLETE - TODO: redo functional analysis
R --vanilla < ${THIS}/scripts/saap_analysis.R  > log/analysis.txt
## functional enrichment of SAAP-harboring proteins
R --vanilla < ${THIS}/scripts/saap_function.R


### GENERATE RESULTS SLIDES AND PAPER FIGURES


cp -a ~/Documents/sejour23_fig1.jpg .

pandoc -t beamer tsour23pre.md --filter pandoc-citeproc -o tsour23pre.pdf

## TODO: instead, generate a blastdb of all transcripts and proteins,
## blast base peptides within those proteins and record locations
