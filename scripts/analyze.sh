#!/bin/bash

## download ensembl human proteins
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

## map each peptide to position
R --vanilla < map_peptides.R

cp -a ~/Documents/sejour23_fig1.jpg .

pandoc -t beamer tsour23pre.md --filter pandoc-citeproc -o tsour23pre.pdf

## TODO: instead, generate a blastdb of all transcripts and proteins,
## blast base peptides within those proteins and record locations
