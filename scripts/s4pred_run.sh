#!/bin/sh


i=${1}
seq_per_call=${2}
ncpu=${3}

fasta=/home/${USER}/data/mammary/processedData/Homo_sapiens.GRCh38.pep.large.fa
spred=/home/${USER}/programs/s4pred/


## use awk to get range of sequences
from=`echo $i  $seq_per_call | awk '{print $1 * $2}' -`
to=`echo $from $seq_per_call | awk '{print $1 + $2 -1}' -`
echo $i from $from to $to

## retrieve sequences from fasta file
srun awk "/^>/ {n++} (n>=$from && n<=$to) {print}" $fasta > s4pred_${i}.fas

## call python
srun python ${spred}/run_model.py s4pred_${i}.fas --threads $ncpu -z -s -o /home/${USER}/data/mammary/processedData/s4pred/
