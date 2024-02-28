#!/bin/bash



## copy only files required for SAAP analysis to a new directory
## much smaller for easier transfer
if [ false ]; then
    cut -f 12 ~/data/mistrans/processedData/saap_mapped.tsv | sed '/^$/d' | tail -n +2 > ~/data/mistrans/processedData/mapped_ensembl_proteins.dat
    while read p; do
	echo "$p" 
	find ~/data/mammary/processedData/iupred3 -name "${p}*_iupred3.tsv.gz" \
	     >>  trasnfer_files.txt
    done < ~/data/mistrans/processedData/mapped_ensembl_proteins.dat
    cat trasnfer_files.txt |sort |uniq | xargs -I % cp % iupred3_selected/
fi
