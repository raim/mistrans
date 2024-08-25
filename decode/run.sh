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

## download SUPPLEMENTAL DATA TABLES 2-4 AND COPY TO
## $DECDATA/originalData

## download additionalData.zip which provides human genome data files
## generated in by the genomeBrowser project AND unzip in $DECDATA/

required1=${DECDATA}/originalData/Supplemental_Data_3.SAAP_precursor_quant.xlsx
if [ ! -f "$required1" ]; then
    echo "MANUAL STEPS REQIRED: Please download the supplemental data tables"
    echo "and copy to required folder as instructed."
    exit 1
fi
required2=${DECDATA}/additionalData/features_GRCh38.110.tsv.gz
if [ ! -f "$required2" ]; then
    echo "MANUAL STEPS REQIRED: Please download the additionalData.zip"
    echo "and unzip as instructed."
    exit 1
fi


## 2) CONVERT SAAP/BP INPUT DATA (by Shiri Tsour) to text files
$ssconvert --export-type=Gnumeric_stf:stf_assistant \
	   -O "locale=C format=automatic separator='	' eol=unix" \
	   $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.xlsx \
	   $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.tsv 
$ssconvert --export-type=Gnumeric_stf:stf_assistant \
	   -O "locale=C format=automatic separator='	' eol=unix" \
	   $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.xlsx \
	   $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.tsv 
$ssconvert --export-type=Gnumeric_stf:stf_assistant \
	   -O "locale=C format=automatic separator='	' eol=unix" \
	   $DECDATA/originalData/Supplemental_Data_2.SAAP_proteins.xlsx \
	   $DECDATA/originalData/Supplemental_Data_2.SAAP_proteins.tsv
## gzip to save disk space
gzip $DECDATA/originalData/Supplemental_Data_2.SAAP_proteins.tsv \
     $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.tsv \
     $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.tsv

exit 1

## 3) COLLECT BP/SAAP SEQUENCES
## produces $DECDATA/processedData/unique_bp.fas and
##          $DECDATA/processedData/unique_saap.tsv

## write BP into a fasta file for use with blast
gunzip -c $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.tsv.gz \
    | cut -f 4 | sort | uniq | grep -v "^BP$" | \
    grep -v -e '^$' > $DECDATA/processedData/tmp.txt
gunzip -c $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.tsv.gz \
    | cut -f 4 | sort | uniq | grep -v "^BP$" | \
    grep -v -e '^$' > $DECDATA/processedData/tmp2.txt
cat $DECDATA/processedData/tmp.txt $DECDATA/processedData/tmp2.txt | \
    sort |uniq | awk '{print ">" $0 ORS $0}' - \
		     > $DECDATA/processedData/unique_bp.fas
rm -f $DECDATA/processedData/tmp.txt $DECDATA/processedData/tmp2.txt


## BP/SAAP as simple table, basis for search in proteins
gunzip -c $DECDATA/originalData/Supplemental_Data_3.SAAP_precursor_quant.tsv.gz \
    | cut -f 3,4 | sort | uniq | grep -v -P "^SAAP\tBP" | \
    grep -v -e '^$' > $DECDATA/processedData/tmp.tsv
gunzip -c $DECDATA/originalData/Supplemental_Data_4.SAAP_reporter_quant.tsv.gz \
    | cut -f 3,4 | sort | uniq | grep -v -P "^SAAP\tBP" | \
    grep -v -e '^$' > $DECDATA/processedData/tmp2.tsv
cat $DECDATA/processedData/tmp.tsv $DECDATA/processedData/tmp2.tsv | \
    sort | uniq > $DECDATA/processedData/unique_saap.tsv
rm -f $DECDATA/processedData/tmp.txt $DECDATA/processedData/tmp2.txt


## 4) GET Ensembl PROTEIN SEQUENCES
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz $DECDATA/originalData/

## 5) ADD PATIENT-SPECIFIC SINGLE AMINO REPLACEMENTS
## produces $DECDATA/processedData/all_proteins.fa
$myR --vanilla < get_mutated_proteins.R  &> $DECDATA/log/mutated_proteins.txt 


## 6) BLAST BP IN Ensembl+MUTATION PROTEINS
## produces $DECDATA/processedData/unique_bp_blast.tsv

$blastdir/makeblastdb -in $DECDATA/processedData/all_proteins.fa -parse_seqids \
		      -title "ensembl hg38 proteins" -dbtype prot  \
    &> $DECDATA/log/blast_database.txt
## NOTE: blast hits are pre-filtered for query length=alignment length
##       and >75% identity
format="6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"
${blastdir}/blastp  -num_threads 7 -task blastp-short \
	   -query  $DECDATA/processedData/unique_bp.fas \
	   -db $DECDATA/processedData/all_proteins.fa  \
	   -outfmt $format 2> $DECDATA/log/unique_bp_blast.txt \
    | awk '{if($5==$6 && $3>75) print}' \
    | grep -v "^#" > $DECDATA/processedData/unique_bp_blast.tsv 


## 6) SELECT BEST-MATCHING PROTEIN HIT FOR EACH BP
## produces $DECDATA/processedData/bp_mapped.tsv
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

if false; then
    $myR --vanilla < map_peptides.R &> $DECDATA/log/map_peptides.txt
fi



### DATA ANALYSIS AND PLOTS

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
