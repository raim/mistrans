#!/bin/bash


export MISDATA=${MYDATA}/mistrans
export MAMDATA=${MYDATA}/mammary
SRC=$GENBRO/data/mammary
THIS=${HOME}/work/mistrans

## SCRIPTS to analyze location and function of mistranslation
## events, derived by Shiri Tsour from the slavovlab.

## INPUT GENOME DATA IS GENERATED
## by from genomeBrowser/mammary/setup.sh

### ADDITIONAL DATA

## DEGRONS
wget https://degronopedia.com/degronopedia/download/data/DEGRONOPEDIA_degron_dataset.xlsx -P $MISDATA/originalData/

## CODON DATA
## TODO: add download url for dana14_codons.ods/csv

## @HernandezAlias2023: Using protein-per-mRNA differences among human
## tissues in codon optimization
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-02868-2/MediaObjects/13059_2023_2868_MOESM3_ESM.xlsx -P $MISDATA/originalData/ -O hernandez-alias23_file3.xlsx

#W @Wu2019
## We calculated the codon stability coefficient (CSC) as the Pearson
## correlation coefficient between mRNA stability and codon
## occurrence. [...] The CSC scores do not present strong correlation
## with codon usage (Figure 1—figure supplement 1C).
wget https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDUzOTYvZWxpZmUtNDUzOTYtZmlnMS1kYXRhMi12Mi5jc3Y-/elife-45396-fig1-data2-v2.csv?_hash=oV0Fjo95uQOzu5LreFXU9sbiAG2ub8ZzXLyP%2B0iTk98%3D -P $MISDATA/originalData/ -O elife-45396-fig1-data2-v2.csv



### SAAP/RAAS ANALYSIS

## analyze data structure, different number of replicates per unique SAAP
R --vanilla < ${THIS}/scripts/saap_means.R > log/means.txt


## ANNOTATE ALL BP/SAAP

## 1) collect all BP and SAAP/BP pairs

## BP as fasta for blast - TODO: remove header!
cut -f 5 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 5 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort |uniq | awk '{print ">" $0 ORS $0}' - > ${MISDATA}/processedData/unique_bp_2.fas

## how many? 7991 unique BP
grep -n ">"  ${MISDATA}/processedData/unique_bp.fas |wc -l
## with tonsil in TMT dataset
grep -n ">"  ${MISDATA}/processedData/unique_bp_2.fas |wc -l

## SAAP as fasta for blast - TODO: remove header!
cut -f 4 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 4 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort |uniq | awk '{print ">" $0 ORS $0}' - > ${MISDATA}/processedData/unique_saap_2.fas

## how many? 15061 unique SAAP
grep -n ">"  ${MISDATA}/processedData/unique_saap.fas |wc -l
## with tonsil in TMT dataset
grep -n ">"  ${MISDATA}/processedData/unique_saap_2.fas |wc -l

## BP/SAAP as simple table, basis for search in proteins
cut -f 4,5 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 4,5 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort | uniq > ${MISDATA}/processedData/unique_saap_2.tsv

## how many? 15669 unique SAAP/BP, 15061 unique SAAP
wc -l  ${MISDATA}/processedData/unique_saap.tsv
cut -f 1  ${MISDATA}/processedData/unique_saap.tsv |sort|uniq|wc -l
## with tonsil in TMT dataset
wc -l  ${MISDATA}/processedData/unique_saap_2.tsv
cut -f 1  ${MISDATA}/processedData/unique_saap_2.tsv |sort|uniq|wc -l

## 2) collect all proteins tagged with mutations and add these to protein DB;
##    generates ${MISDATA}/processedData/all_proteins.fa 
R --vanilla < ${THIS}/scripts/get_mutated_proteins.R

## 3) blast all BP against ensembl proteins + mutations
blastdir=${HOME}/programs/ncbi-blast-2.15.0+/bin
## generate local blastdb 
$blastdir/makeblastdb -in ${MISDATA}/processedData/all_proteins.fa -parse_seqids -title "ensembl hg38 proteins" -dbtype prot
## blast - filter full length hit alignment length=query length,
## and at least 75% identity with awk.
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/unique_bp_2.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  | awk '{if($5==$6 && $3>75) print}'  |grep -v "^#" > ${MISDATA}/processedData/unique_bp_blast_2.tsv

## 3.B) blast SP against ensembl+mutations
## NOTE: using BP as fasta title: gives warning of >50 valid AA in title
## this should not be a problem, but TODO: check whether these are correctly
## reflected (ini full length) in blast output, or cut at 50.
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/unique_saap_2.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  | awk '{if($5==$6 && $3>75) print}'  |grep -v "^#" > ${MISDATA}/processedData/unique_saap_blast_2.tsv

## some statistics on blast
## TODO: expand this QC analysis a bit,
## do SAAP have the expected mismatches?
R --vanilla < ${THIS}/scripts/saap_blast_stats.R

## 3.C) blast all main peptides against 20k core proteins
cat ${MISDATA}/originalData/All_main_tryptic_peptide_list.txt | sort |uniq | awk '{print ">" $0 ORS $0}' > ${MISDATA}/processedData/main_peptides.fas
## call blast, only report full length hits with 100% identity
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/main_peptides.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  | awk '{if($5==$6 && $3==100) print}'  |grep -v "^#" > ${MISDATA}/processedData/main_peptides_blast.tsv

### 4) COLLECT DATA FOR ALL BP/SAAP

## 4.A) find best matching protein
##    GENERATES ${MISDATA}/processedData/bp_mapped_3.tsv
R --vanilla < ${THIS}/scripts/get_protein_match.R > ${THIS}/log/match_3.txt

## collect ONLY required iupred3 data for transfer to intron
if [ false ]; then
    cd ~/data/mammary/processedData
    all=$(cut -f 2 ~/data/mistrans/processedData/bp_mapped_3.tsv |sort|uniq|grep ENSP|grep -v "_")
    for i in $all
    do
	echo "$i"
	# or do whatever with individual element of the array
	find iupred3/ -name "${i}*.tsv.gz" -exec cp -a {} iupred3_selected \;
    done
    zip -r iupred3_selected iupred3_selected
    cd -
fi  


## 4.B) map each peptide to it's positions in proteins and transcripts, and
## add a variety of collected information on protein structure, e.g.
## iupred3, anchor2, s4pred, codon, ...
##    GENERATES ${MISDATA}/processedData/saap_mapped_5.tsv, and
##    QC figures in ${MISDATA}/figures/saap_mapping5/
R --vanilla < ${THIS}/scripts/map_peptides3.R > ${THIS}/log/map_5.txt

### 5) ANALYSIS


### CALCULATE RAAS PROFILES
R --vanilla <  ${THIS}/scripts/raasprofiles3_codons.R # FIGURE 2
R --vanilla <  ${THIS}/scripts/raasprofiles3_aminoacids.R # FIGURE 3
R --vanilla <  ${THIS}/scripts/raasprofiles3_motifs.R # FIGURE 4
R --vanilla <  ${THIS}/scripts/raasprofiles3_kraq.R ## TODO: reproduce diAA
R --vanilla <  ${THIS}/scripts/raasprofiles3_structure.R
R --vanilla <  ${THIS}/scripts/raasprofiles3_domains.R
R --vanilla <  ${THIS}/scripts/raasprofiles3_proteins.R

### TODO: move those two scripts to genomeBrowser
R --vanilla < ${THIS}/scripts/halflives.R
R --vanilla < ${THIS}/scripts/ralser24.R

## all protein profiles
R --vanilla <  ${THIS}/scripts/saap_proteins.R



### OUTDATED OLD CODE:

## MOTIFS: get and analyze sequences surrounding the ASS

## export sequence context of AAS
## generates file ${MISDATA}/processedData/saap_context.tsv
R --vanilla < ${THIS}/scripts/export_aas_context.R > ${THIS}/log/context_export.txt

## extract sequence sets from saap_context.tsv for motif analysis

## binomial distribution (hypergeo) tests of AA around AAS
## including miscleavage etc.
R --vanilla < ${THIS}/scripts/get_aas_context.R > ${THIS}/log/context_analysis.txt

## TODO: fuse two context scripts; and move randomization etc. to deep learning
## python script; run kplogo; select motifs and sequences for RAAS analysis.

## kplogo
## TODO: run over all from:to classes and find a way to collect and plot
## results by RAAS!
cd ~/data/mistrans/processedData/motifs
~/programs/kpLogo/bin/kpLogo seqcontext_all_.fa  -alphabet protein -o kplogo/all
~/programs/kpLogo/bin/kpLogo seqcontext_fromto_Q:G.fa -alphabet protein -o kplogo/QG
~/programs/kpLogo/bin/kpLogo seqcontext_fromto_T:V.fa -alphabet protein -o kplogo/TV
~/programs/kpLogo/bin/kpLogo seqcontext_methionine_TRUE.fa -alphabet protein -o kplogo/M
~/programs/kpLogo/bin/kpLogo seqcontext_tryptophane_TRUE.fa -alphabet protein -o kplogo/W

## TODO: log files
## with Albumin
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## without Albumin
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## unique SAAP without Albumin 
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## unique SAAP with Albumin 
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &

## same for healthy tissues
## with Albumin
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## without Albumin
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## unique SAAP without Albumin 
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## unique SAAP with Albumin 
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla


## OLD and OBSOLETE - TODO: redo functional analysis
R --vanilla < ${THIS}/scripts/saap_analysis.R  > log/analysis.txt
## functional enrichment of SAAP-harboring proteins
R --vanilla < ${THIS}/scripts/saap_function.R


### GENERATE RESULTS SLIDES AND PAPER FIGURES


cp -a ~/Documents/sejour23_fig1.jpg .

pandoc -t beamer tsour23pre.md --filter pandoc-citeproc -o tsour23pre.pdf

## TODO: instead, generate a blastdb of all transcripts and proteins,
## blast base peptides within those proteins and record locations
