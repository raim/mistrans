#!/bin/bash


module load python/3.6.6

ncpu=20
fasta=/home/${USER}/data/mammary/processedData/Homo_sapiens.GRCh38.pep.large.fa
spred=/home/${USER}/programs/s4pred/

## 10min for 100 sequences in tests with $ncpu=20
## 600 per hour -> do 600 per run to split into ca. 200 jobs

total=`grep ">" $fasta |wc -l`
seq_per_call=2000
start=0
end=`echo $total $seq_per_call | awk '{print int($1/$2)}' -`

for (( i=$start; i<=$end; i++ )); do

    ## use awk to get range of sequences
    from=`echo $i  $seq_per_call | awk '{print $1 * $2}' -`
    to=`echo $from $seq_per_call | awk '{print $1 + $2 -1}' -`
    echo $i from $from to $to
    ## retrieve sequences from fasta file
    awk "/^>/ {n++} (n>=$from && n<=$to) {print}" $fasta > s4pred_${i}.fas

    ## run script
    sbatch -a $i --job-name=s4pred_${i} --output=s4pred_${i}.out -p short --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB python ${spred}/run_model.py s4pred_${i}.fas --threads $ncpu -z -s -o /home/${USER}/data/mammary/processedData/s4pred/
done

