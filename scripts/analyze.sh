#!/bin/bash

## SCRIPTS to analyze location and function of mistranslation
## events, as provided by Shiri Tsour and Nikolai Slavov.


### TODO:
## * consider which parts of genomeBrowser/data/mammary should/could be
##   moved here,

### MAIN INPUT DATA

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


## program paths
blastdir=${HOME}/programs/ncbi-blast-2.15.0+/bin

## data paths
export MISDATA=${MYDATA}/mistrans
export MAMDATA=${MYDATA}/mammary
SRC=$GENBRO/data/mammary
THIS=${HOME}/work/mistrans

## generate output paths
mkdir -p $MISDATA/log
mkdir -p $MISDATA/figures
mkdir -p $MISDATA/originalData
mkdir -p $MISDATA/processedData


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
R --vanilla < ${THIS}/scripts/get_mutated_proteins.R  &> ${MISDATA}/log/mutated_proteins.txt 

## how many? 131328 proteins!
grep -n ">"  ${MISDATA}/processedData/all_proteins.fa|wc -l

## 3) blast all BP against ensembl proteins + mutations
## generate local blastdb 
$blastdir/makeblastdb -in ${MISDATA}/processedData/all_proteins.fa -parse_seqids -title "ensembl hg38 proteins" -dbtype prot  &> ${MISDATA}/log/blast_database.txt
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
R --vanilla < ${THIS}/scripts/saap_blast_stats.R &> ${MISDATA}/log/saap_blast_stats.txt

## SEARCH SECIS RFAM in ensembl transcripts
## SECIS1-5: RF00031, RF01988, RF01989, RF01990, RF01991
myrfam=${MYDATA}/rfam/
declare -a rfams=("RF00031" "RF01988" "RF01989" "RF01990" "RF01991")

## now loop through the above array
for i in "${rfams[@]}"; do
    echo $i
    cmsearch --cpu 7 --tblout  ${MISDATA}/processedData/${i}.out --notextw ${myrfam}/${i}.cm $MAMDATA/originalData/Homo_sapiens.GRCh38.cdna.all.fa.gz &>> ${MISDATA}/processedData/${i}.log
done
## NOTE: the file secis_elements.out can be parsed by cmchainer::parseHomHits
grep -h -v "^#" ${MISDATA}/processedData/RF*.out > ${MISDATA}/processedData/secis_elements.out

### 4) COLLECT DATA FOR ALL BP/SAAP

## 4.A) find best matching protein
##    GENERATES ${MISDATA}/processedData/bp_mapped.tsv
R --vanilla < ${THIS}/scripts/get_protein_match.R &> ${MISDATA}/log/protein_match.txt

## collect ONLY required iupred3 data for transfer to laptop (intron)
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
##    GENERATES ${MISDATA}/processedData/saap_mapped.tsv, and
##    QC figures in ${MISDATA}/figures/saap_mapping/
R --vanilla < ${THIS}/scripts/map_peptides.R &> ${MISDATA}/log/map_peptides.txt

## sanity check derived coordinates: load protein, transcript and genome fasta files
## to check correct nucleotides/codons at the selected positions.
R --vanilla < ${THIS}/scripts/check_coordinates.R &> ${MISDATA}/log/check_coordinates.txt

### 5) ANALYSIS

## initializiation for all scripts below; NOTE, that this script
## is (re-)called directly from the scripts below; TODO: call init
## only once and use in scripts.
## PRODUCES site_raas.txv,  a site-specific raas values with
## protein annotation, for use in random forest model.
R --vanilla <  ${THIS}/scripts/raasprofiles3_init.R &> ${MISDATA}/log/init.txt 

### CALCULATE RAAS PROFILES
R --vanilla <  ${THIS}/scripts/raasprofiles3_codons.R &> ${MISDATA}/log/codons.txt  # FIGURE 2
R --vanilla <  ${THIS}/scripts/raasprofiles3_aminoacids.R &> ${MISDATA}/log/aminoacids.txt  # FIGURE 3
R --vanilla <  ${THIS}/scripts/raasprofiles3_motifs.R &> ${MISDATA}/log/motifs.txt  # FIGURE 4
R --vanilla <  ${THIS}/scripts/raasprofiles3_kraq.R &> ${MISDATA}/log/kraq.txt  # EXTENDED DATA FIGURE 8: motif selection
R --vanilla <  ${THIS}/scripts/raasprofiles3_structure.R &> ${MISDATA}/log/structure.txt 
R --vanilla <  ${THIS}/scripts/raasprofiles3_function.R &> ${MISDATA}/log/function.txt 
R --vanilla <  ${THIS}/scripts/raasprofiles3_proteins.R &> ${MISDATA}/log/proteins.txt 
R --vanilla <  ${THIS}/scripts/raasprofiles3_model.R &> ${MISDATA}/log/model.txt 
R --vanilla <  ${THIS}/scripts/raasprofiles3_rna.R &> ${MISDATA}/log/rna.txt

## COMPARE MULTIPLE BP per SAAP
R --vanilla <  ${THIS}/scripts/raasprofiles3_multipleBP.R &> ${MISDATA}/log/multipleBP.txt 


## PROTEIN 3D: chimeraX codes for pdb
R --vanilla <  ${THIS}/scripts/raasprofiles3_pdbscan.R &> ${MISDATA}/log/protein_pdb.txt 


## PROTEIN 1D: plots of AAS along protein 1D

## first blast all main peptides against 20k core proteins
## NOTE: runs a few hours!
cat ${MISDATA}/originalData/All_main_tryptic_peptide_list.txt | sort |uniq | awk '{print ">" $0 ORS $0}' > ${MISDATA}/processedData/main_peptides.fas
## call blast, only report full length hits with 100% identity
${blastdir}/blastp  -num_threads 7 -task blastp-short -query  ${MISDATA}/processedData/main_peptides.fas -db ${MISDATA}/processedData/all_proteins.fa   -outfmt "6 qseqid sacc pident mismatch length qlen slen sstart send  evalue bitscore"  2> ${MISDATA}/log/main_peptides_blast.txt | awk '{if($5==$6 && $3==100) print}'  |grep -v "^#" > ${MISDATA}/processedData/main_peptides_blast.tsv

## plot 1D for all proteins
R --vanilla <  ${THIS}/scripts/saap_proteins.R &> ${MISDATA}/log/protein_plots.txt 

### COLLECT PUBLICATION FIGURES HERE!



## Figure 2:
## * codon frequency as table file,
## * full codon dotprofiles,

ftyp=pdf #png

results=${MISDATA}/results_${ftyp}/
mkdir $results

cp -a ${THIS}/scripts/README.md $results/

cp -a ${MISDATA}/processedData/saap_mapped.tsv $results/
cp -a ${MISDATA}/processedData/sites_raas.tsv $results/

cp -a ${MISDATA}/figures/raasprofiles3/legend_dotplot_acols_slim.${ftyp} $results/


mkdir $results/codons
cp -a ${MISDATA}/figures/raasprofiles3/codons/codon_frequencies.tsv $results/codons/
cp -a ${MISDATA}/figures/raasprofiles3/codons/codons_raas_fbg.${ftyp} $results/codons/
cp -a ${MISDATA}/figures/raasprofiles3/codons/codons_raas_wu19.${ftyp} $results/codons/
## NOTE: markdown includes figures from _codons.R 
pandoc ${THIS}/scripts/results_codons.md -t beamer -o $results/codons/codon_plot_A4.pdf
pdfcrop --clip $results/codons/codon_plot_A4.pdf  $results/codons/Extended_Data_Figure_4.pdf

## Figure 3:
## 3e: encoded/incorporated selection,
## 3b: AAS type sorted/colored by AA properties,
## 3c: AA property selection.
## Extended 5a: incorporated by encoded over all Datasets,
## Extended 5b: AA property,
## Extended 5f: all encoded/incorporated.

mkdir $results/aminoacids
cp -a  ${MISDATA}/figures/raasprofiles3/aminoacids/AA_cancer_cut_dotplot_manual.${ftyp} $results/aminoacids/
cp -a  ${MISDATA}/figures/raasprofiles3/aminoacids/AA_cancer_all_dotplot.${ftyp} $results/aminoacids/

cp -a  ${MISDATA}/figures/raasprofiles3/aminoacids/AAprop_cancer_cut_dotplot_manual.${ftyp} $results/aminoacids/
cp -a  ${MISDATA}/figures/raasprofiles3/aminoacids/AAprop_cancer_dotplot_manual.${ftyp} $results/aminoacids/

cp -a  ${MISDATA}/figures/raasprofiles3/aminoacids/fromAA_cancer_cut_dotplot_manual_rotated.${ftyp} $results/aminoacids/
cp -a  ${MISDATA}/figures/raasprofiles3/aminoacids/fromAA_cancer_dotplot.${ftyp} $results/aminoacids/

cp -a  ${MISDATA}/figures/raasprofiles3/aminoacids/toAA_cancer_cut_dotplot_manual_rotated.${ftyp} $results/aminoacids/
cp -a  ${MISDATA}/figures/raasprofiles3/aminoacids/toAA_cancer_dotplot.${ftyp} $results/aminoacids/

mkdir $results/motifs
cp -a ${MISDATA}/figures/raasprofiles3/motifs/selected/AA_logos_KRAQ.${ftyp} $results/motifs/
cp -a ${MISDATA}/figures/raasprofiles3/motifs/selected/AA_logos_CCxCC.${ftyp} $results/motifs/
cp -a ${MISDATA}/figures/raasprofiles3/motifs/selected/AA_logos_MMxMM.${ftyp} $results/motifs/
cp -a ${MISDATA}/figures/raasprofiles3/motifs/selected/AA_logos_WWxWW.${ftyp} $results/motifs/
cp -a ${MISDATA}/figures/raasprofiles3/motifs/motifs_cancer_dotplot.${ftyp} $results/motifs/
cp -a ${MISDATA}/figures/raasprofiles3/motifs/classes_conservation_disorder_raas.${ftyp} $results/motifs/

##cp -a ${MISDATA}/figures/raasprofiles3/motifs/selected/logos_fromto_*_encoded.${ftyp} $results/motifs/
## NOTE: markdown includes figures from _kraq.R and _motifs.R
pandoc ${THIS}/scripts/results_motifs.md  -o $results/motifs/Extended_Data_Figure_8_A4.pdf
pdfcrop --clip $results/motifs/Extended_Data_Figure_8_A4.pdf $results/motifs/Extended_Data_Figure_8.pdf

mkdir $results/structure
cp -a ${MISDATA}/figures/raasprofiles3/proteins/protein_intensities_all.${ftyp} $results/structure/
cp -a ${MISDATA}/figures/raasprofiles3/proteins/protein_halflives_all.${ftyp} $results/structure/
cp -a ${MISDATA}/figures/raasprofiles3/proteins/protein_lengths_all.${ftyp} $results/structure/
cp -a ${MISDATA}/figures/raasprofiles3/proteins/protein_Tmelt_all.${ftyp} $results/structure/
cp -a ${MISDATA}/figures/raasprofiles3/proteins/proteins_raas.tsv $results/structure/

cp -a ${MISDATA}/figures/raasprofiles3/structure/structure_cor_iupred3_RAAS.${ftyp} $results/structure/
cp -a ${MISDATA}/figures/raasprofiles3/structure/structure_cor_MMSeq2_RAAS.${ftyp} $results/structure/

mkdir $results/function
cp -a ${MISDATA}/figures/raasprofiles3/function/type_go_cancer_ptgt_high_dotplot.${ftyp} $results/function
## NOTE: markdown includes figures from _function.R 
pandoc ${THIS}/scripts/go_dissection.md -t beamer -o $results/function/go_dissection.pdf

cd $MISDATA
zip -r results_${ftyp} results_${ftyp}

##zip file of all protein plots,
##zip file of all chimeraX commands.
