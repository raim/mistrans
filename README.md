Protein mutation data by Shiri Tsour, 2023.

Ideas

# 5'/3' bias

* map peptides and mutations via ensembl/refseq ID to full length proteins,
* analyze 5'/3' bias, and compare to 5'/3' bias of codon usage.

simplest approache: 

* get all unique ensembl proteins,
* for each peptide, get its coordinate with the protein.

# secondary structure

* map to secondary structure prediction, ideally with
  local classification which AA has which role,
* alpha v beta may have different translation speeds and different mutation rates,
* alpha can switch to beta, and beta is involved in protein aggregation;
  such switching may positively or negatively affect protein stability.
