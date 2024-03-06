#!/bin/bash



## copy only files required for SAAP analysis to a new directory
## much smaller for easier transfer

cut -f 13 ~/data/mistrans/processedData/saap_mapped3.tsv | sed '/^$/d' | tail -n +2 |sort|uniq > ~/data/mistrans/processedData/mapped_ensembl_proteins3.dat
while read p; do
    echo "$p" 
    find ~/data/mammary/processedData/iupred3 -name "${p}*_iupred3.tsv.gz" \
	 >>  transfer_files3.txt
done < ~/data/mistrans/processedData/mapped_ensembl_proteins3.dat
cat transfer_files3.txt |sort|uniq| xargs -I % cp % iupred3_selected/

zip -r iupred3_selected iupred3_selected
