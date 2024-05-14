Protein mutation data by Shiri Tsour, 2023.

Ideas

## TODO

# ANDREW WAS HERE

* hypothesis tree:
    - synthesis: codons, ribosome density, translation kinetics,
	- degradation: degrons,
	- dominance: catalytic sites, protein interactions, 

* unstructured domains:
    - iupred3,  PONDR® VSL2B,
	- 3500 proteomes from viruses and the three domains of life: https://doi.org/10.1080/07391102.2012.675145,
	- https://disprot.org/
* catalytic sites:
    - https://www.ebi.ac.uk/thornton-srv/m-csa/ 

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
