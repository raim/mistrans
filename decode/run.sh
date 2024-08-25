#!/bin/bash

## ANALYSIS OF AMINO ACID SUBSTITUTION SITES
##  as provided by Shiri Tsour.

## SOFTWARE DEPENDENCIES

## additionally to standard bash tools (cat, grep, cut, awk, gunzip, rsync),
## this script additionally requires gnumeric's ssconvert, NCBI's blast,
## and a recent R version. Please define the paths to these tools here!

blastdir=${HOME}/programs/ncbi-blast-2.15.0+/bin
ssconvert=/usr/bin/ssconvert
myR=/usr/bin/R

## NUMBER OF CPUs used for blast
nthreads=7

## OUTPUT DATA PATHS
## SET THIS PATH TO WHERE YOU WANT TO SAVE INPUT and OUTPUT
export DECDATA=/home/raim/data/decode

## generate output paths
mkdir -p $DECDATA/log
mkdir -p $DECDATA/figures
mkdir -p $DECDATA/processedData
mkdir -p $DECDATA/originalData

### DATA COLLECTION

## 1) REQUIRED MANUAL DOWNLOADS

## download additionalData.zip which provides human genome data files
## generated in by the genomeBrowser project AND unzip in $DECDATA/

## cp -a ~/data/mistrans/additionalData.zip $DECDATA/
## cd $DECDATA; unzip additionalData.zip

echo "testing existence of required input data"

required=${DECDATA}/additionalData/features_GRCh38.110.tsv.gz
if [ ! -f "$required" ]; then
    echo "MANUAL STEPS REQIRED: Please download the additionalData.zip"
    echo "and unzip as instructed."
    exit 1
fi



## 2) COLLECT BP/SAAP SEQUENCES
## produces $DECDATA/processedData/unique_bp.fas and
##          $DECDATA/processedData/unique_saap.tsv
echo "extracting BP and SAAP sequences"

## write BP into a fasta file for use with blast
gunzip -c $DECDATA/additionalData/All_SAAP_TMTlevel_quant_df.txt.gz \
    | cut -f 5 | sort | uniq | grep -v "^BP$" | \
    grep -v -e '^$' > $DECDATA/processedData/tmp.txt
gunzip -c $DECDATA/additionalData/All_SAAP_patient_level_quant_df.txt \
    | cut -f 5 | sort | uniq | grep -v "^BP$" | \
    grep -v -e '^$' > $DECDATA/processedData/tmp2.txt
cat $DECDATA/processedData/tmp.txt $DECDATA/processedData/tmp2.txt | \
    sort |uniq | awk '{print ">" $0 ORS $0}' - \
		     > $DECDATA/processedData/unique_bp.fas
rm -f $DECDATA/processedData/tmp.txt $DECDATA/processedData/tmp2.txt


## BP/SAAP as simple table, basis for search in proteins
gunzip -c $DECDATA/additionalData/All_SAAP_TMTlevel_quant_df.txt.gz \
    | cut -f 4,5 | sort | uniq | grep -v -P "^SAAP\tBP" | \
    grep -v -e '^$' > $DECDATA/processedData/tmp.tsv
gunzip -c $DECDATA/additionalData/All_SAAP_patient_level_quant_df.txt \
    | cut -f 4,5 | sort | uniq | grep -v -P "^SAAP\tBP" | \
    grep -v -e '^$' > $DECDATA/processedData/tmp2.tsv
cat $DECDATA/processedData/tmp.tsv $DECDATA/processedData/tmp2.tsv | \
    sort | uniq > $DECDATA/processedData/unique_saap.tsv
rm -f $DECDATA/processedData/tmp.txt $DECDATA/processedData/tmp2.txt


## 3) GET Ensembl PROTEIN SEQUENCES
echo "downloading Ensembl protein fasta"
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz $DECDATA/originalData/

## 4) ADD PATIENT-SPECIFIC SINGLE AMINO REPLACEMENTS
## produces $DECDATA/processedData/all_proteins.fa
echo "adding mutations to fasta"
$myR --vanilla < get_mutated_proteins.R  &> $DECDATA/log/mutated_proteins.txt 


## 5) BLAST BP IN Ensembl+MUTATION PROTEINS
## produces $DECDATA/processedData/unique_bp_blast.tsv
echo "RUNNING BLAST SEARCH OF BP"

$blastdir/makeblastdb -in $DECDATA/processedData/all_proteins.fa -parse_seqids \
		      -title "ensembl hg38 proteins" -dbtype prot  \
    &> $DECDATA/log/blast_database.txt
## NOTE: blast hits are pre-filtered for query length=alignment length
##       and >75% identity
format="6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"
${blastdir}/blastp  -num_threads ${nthreads} -task blastp-short \
	   -query  $DECDATA/processedData/unique_bp.fas \
	   -db $DECDATA/processedData/all_proteins.fa  \
	   -outfmt "$format" 2> $DECDATA/log/unique_bp_blast.txt \
    | awk '{if($5==$6 && $3>75) print}' \
    | grep -v "^#" > $DECDATA/processedData/unique_bp_blast.tsv 


## 6) SELECT BEST-MATCHING PROTEIN HIT FOR EACH BP
## produces $DECDATA/processedData/bp_mapped.tsv
echo "getting best blast hit for each BP"

$myR --vanilla < get_protein_match.R &> $DECDATA/log/protein_match.txt

## 7) COLLECT ALL DATA FOR UNIQUE BP/SAAP PAIRS
## produces $DECDATA/processedData/saap_mapped.tsv

## NOTE: this script requires large additional data produced
## by genomeBrowser project and partially
## requiring high performance computing. The script can be run
## after doing these calculations and providing the correct path
## to the resulting data directory; see
## https://gitlab.com/raim/genomeBrowser/-/blob/master/data/mammary/setup.sh

## Here, we provide the produced file saap_mapped.tsv in the
## downloaded additionalData folder instead, to allow running
## the R analysis scripts that produce the final published figures, below.
echo "collecting BP protein, transcript and genome coordinates" \
     "and various structural data"

$myR --vanilla < map_peptides.R &> $DECDATA/log/map_peptides.txt


### DATA ANALYSIS AND PLOTS
echo "RUNNING R ANALYSIS PIPELINE, producing all figures and data tables" 

run_analysis.sh

