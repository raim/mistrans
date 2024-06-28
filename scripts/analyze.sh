#!/bin/bash


export MISDATA=${MYDATA}/mistrans
export MAMDATA=${MYDATA}/mammary
SRC=$GENBRO/data/mammary
THIS=${HOME}/work/mistrans

mkdir $MISDATA/log
mkdir $MISDATA/figures
mkdir $MISDATA/originalData
mkdir $MISDATA/processedData

## SCRIPTS to analyze location and function of mistranslation
## events, derived by Shiri Tsour from the slavovlab.

### NOTE : INPUT GENOME DATA IS GENERATED
### by genomeBrowser/mammary/setup.sh in $MAMDATA

### ADDITIONAL DATA

## DEGRONS
## NOT USED
##wget https://degronopedia.com/degronopedia/download/data/DEGRONOPEDIA_degron_dataset.xlsx -P $MISDATA/originalData/

## CODON DATA

## @HernandezAlias2023: Using protein-per-mRNA differences among human
## tissues in codon optimization
## NOT USED
##wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-02868-2/MediaObjects/13059_2023_2868_MOESM3_ESM.xlsx -P $MISDATA/originalData/ -O hernandez-alias23_file3.xlsx

## @Dana2014
## TODO: add source and download url for dana14_codons.ods/csv,
## optionally (if present) used in raasprofiles3_codons.R

## @Wu2019
## We calculated the codon stability coefficient (CSC) as the Pearson
## correlation coefficient between mRNA stability and codon
## occurrence. [...] The CSC scores do not present strong correlation
## with codon usage (Figure 1—figure supplement 1C).
cd $MISDATA/originalData/
wget https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDUzOTYvZWxpZmUtNDUzOTYtZmlnMS1kYXRhMi12Mi5jc3Y-/elife-45396-fig1-data2-v2.csv?_hash=oV0Fjo95uQOzu5LreFXU9sbiAG2ub8ZzXLyP%2B0iTk98%3D  -O elife-45396-fig1-data2-v2.csv
cd -

### SAAP/RAAS ANALYSIS


### TODO:
## * move all old data and scripts to a backup folder,
## * consider which parts of genomeBrowser/data/mammary should/could be
##   moved here,
## * rename file versions (scripts and output data) to avoid confusion,
## * freshly download and record Shiri's data,
## * generate .txt files from xlsx, and rerun the whole pipeline.

## Shiri Tsour: main input files from google drive
## downloaded from shared google drive on 20240627
## * All_SAAP_patient_level_quant_df.xlsx
## * All_SAAP_TMTlevel_quant_df.xlsx
## * All_SAAP_protein_filter_df.xlsx

## Sent by Shiri Tsour via slack:
## * main_peptide_quant_df.xlsx
## * tonsil_main_peptide_quant_df.xlsx
## * All_main_nontryptic_peptide_list.txt
## * All_main_tryptic_peptide_list.txt
## * All_main_nontryptic_noKR_peptide_list.txt
## * All_main_ArgC_LysC_peptide_list.txt

## convert xlsx files to text
## NOTE: when copy-pasting this line to bash you need to manually
## add the tab charcter in `separator='<tab>' by typing ctrl-v <TAB>
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "locale=C format=automatic separator='	' eol=unix" $MISDATA/originalData/All_SAAP_TMTlevel_quant_df.xlsx $MISDATA/originalData/All_SAAP_TMTlevel_quant_df.txt 
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "locale=C format=automatic separator='	' eol=unix" $MISDATA/originalData/All_SAAP_patient_level_quant_df.xlsx $MISDATA/originalData/All_SAAP_patient_level_quant_df.txt 
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "locale=C format=automatic separator='	' eol=unix" $MISDATA/originalData/All_SAAP_protein_filter_df.xlsx $MISDATA/originalData/All_SAAP_protein_filter_df.txt



## ANALYZE DATA STRUCTURE:
## different number of measurements per unique SAAP
R --vanilla < ${THIS}/scripts/saap_means.R > log/saap_means.txt


## ANNOTATE ALL BP/SAAP

## 1) collect all BP and SAAP/BP pairs

## BP as fasta for blast - TODO: remove header!
cut -f 5 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 5 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort |uniq | awk '{print ">" $0 ORS $0}' - > ${MISDATA}/processedData/unique_bp.fas

## how many? 8720 unique BP
grep -n ">"  ${MISDATA}/processedData/unique_bp.fas |wc -l

## SAAP as fasta for blast - TODO: remove header!
cut -f 4 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 4 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort |uniq | awk '{print ">" $0 ORS $0}' - > ${MISDATA}/processedData/unique_saap.fas

## how many? 15910 unique SAAP
grep -n ">"  ${MISDATA}/processedData/unique_saap.fas |wc -l


## BP/SAAP as simple table, basis for search in proteins
cut -f 4,5 ${MISDATA}/originalData/All_SAAP_TMTlevel_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp.tsv
cut -f 4,5 ${MISDATA}/originalData/All_SAAP_patient_level_quant_df.txt | sort | uniq > ${MISDATA}/processedData/tmp2.tsv
cat ${MISDATA}/processedData/tmp.tsv ${MISDATA}/processedData/tmp2.tsv | sort | uniq > ${MISDATA}/processedData/unique_saap.tsv

## how many? 16526 unique SAAP/BP, 15910 unique SAAP
wc -l  ${MISDATA}/processedData/unique_saap.tsv
cut -f 1  ${MISDATA}/processedData/unique_saap.tsv |sort|uniq|wc -l

## 2) collect all proteins tagged with mutations and add these to protein DB;
##    generates ${MISDATA}/processedData/all_proteins.fa 
R --vanilla < ${THIS}/scripts/get_mutated_proteins.R  > ${MISDATA}/log/mutated_proteins.txt 2>&1

## how many? 131328 proteins!
grep -n ">"  ${MISDATA}/processedData/all_proteins.fa|wc -l

## 3) blast all BP against ensembl proteins + mutations
blastdir=${HOME}/programs/ncbi-blast-2.15.0+/bin
## generate local blastdb 
$blastdir/makeblastdb -in ${MISDATA}/processedData/all_proteins.fa -parse_seqids -title "ensembl hg38 proteins" -dbtype prot > ${MISDATA}/log/blast_database.txt 2>&1
## blast - filter full length hit alignment length=query length,
## and at least 75% identity with awk.
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/unique_bp.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore" 2> ${MISDATA}/log/unique_bp_blast.txt | awk '{if($5==$6 && $3>75) print}'  |grep -v "^#" > ${MISDATA}/processedData/unique_bp_blast.tsv 

## 3.B) blast SP against ensembl+mutations
## NOTE: using BP as fasta title: gives warning of >50 valid AA in title
## this should not be a problem, but TODO: check whether these are correctly
## reflected (ini full length) in blast output, or cut at 50.
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/unique_saap.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  2> ${MISDATA}/log/unique_saap_blast.txt | awk '{if($5==$6 && $3>75) print}'  |grep -v "^#" > ${MISDATA}/processedData/unique_saap_blast.tsv 

## some statistics on blast
## TODO: expand this QC analysis a bit,
## do SAAP have the expected mismatches?
R --vanilla < ${THIS}/scripts/saap_blast_stats.R > ${MISDATA}/log/saap_blast_stats.txt 2>&1

## 3.C) blast all main peptides against 20k core proteins
cat ${MISDATA}/originalData/All_main_tryptic_peptide_list.txt | sort |uniq | awk '{print ">" $0 ORS $0}' > ${MISDATA}/processedData/main_peptides.fas
## call blast, only report full length hits with 100% identity
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/main_peptides.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  2> ${MISDATA}/log/main_peptides_blast.txt | awk '{if($5==$6 && $3==100) print}'  |grep -v "^#" > ${MISDATA}/processedData/main_peptides_blast.tsv

### 4) COLLECT DATA FOR ALL BP/SAAP

## 4.A) find best matching protein
##    GENERATES ${MISDATA}/processedData/bp_mapped_3.tsv
R --vanilla < ${THIS}/scripts/get_protein_match.R > ${MISDATA}/log/protein_match.txt 2>&1

## collect ONLY required iupred3 data for transfer to intron
if [ false ]; then
    cd ~/data/mammary/processedData
    all=$(cut -f 2 ~/data/mistrans/processedData/bp_mapped.tsv |sort|uniq|grep ENSP|grep -v "_")
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
R --vanilla < ${THIS}/scripts/map_peptides.R > ${MISDATA}/log/map_peptides.txt 2>&1

### 5) ANALYSIS


### CALCULATE RAAS PROFILES
R --vanilla <  ${THIS}/scripts/raasprofiles3_codons.R > ${MISDATA}/log/codons.txt 2>&1 # FIGURE 2
R --vanilla <  ${THIS}/scripts/raasprofiles3_aminoacids.R > ${MISDATA}/log/aminoacids.txt 2>&1 # FIGURE 3
R --vanilla <  ${THIS}/scripts/raasprofiles3_motifs.R > ${MISDATA}/log/motifs.txt 2>&1 # FIGURE 4
R --vanilla <  ${THIS}/scripts/raasprofiles3_kraq.R > ${MISDATA}/log/kraq.txt 2>&1 ## TODO: reproduce diAA
R --vanilla <  ${THIS}/scripts/raasprofiles3_structure.R > ${MISDATA}/log/structure.txt 2>&1
R --vanilla <  ${THIS}/scripts/raasprofiles3_function.R > ${MISDATA}/log/function.txt 2>&1
R --vanilla <  ${THIS}/scripts/raasprofiles3_proteins.R > ${MISDATA}/log/proteins.txt 2>&1

## all protein profiles
R --vanilla <  ${THIS}/scripts/saap_proteins.R &> ${MISDATA}/log/protein_plots.txt 

## all protein chimeraX codes for pdb
R --vanilla <  ${THIS}/scripts/raasprofiles3_pdbscan.R > ${MISDATA}/log/protein_plots.txt 2>&1

### TODO: move those two scripts to genomeBrowser
R --vanilla < ${THIS}/scripts/halflives.R
R --vanilla < ${THIS}/scripts/ralser24.R

## TODO: collect publication figures here!

