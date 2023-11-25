---
title: | 
    | Functional Analysis of MTPs
date: \today
author: Rainer Machn&eacute;, Boston
output: beamer_presentation
classoption: t
bibliography: /home/raim/ref/tata.bib
citation_package: natbib
header-includes:
    - \graphicspath{{/home/raim/data/mistrans/figures/}} 
---


# Codons

![](mtp_per_protein.png){width=60%}

# Codons

\tiny


* number_annot_WP: glycolysis,
* number_annot_GO:CC: 
    - \tiny **extracellular vesicles**,
	- cytoskeleton,
    - mitochondrion,
	- lysosome, vacuole,
	- proteasome,
	- translation and a bit of ribosome.
	
	
## likely false-positives with genomic mutations

(number_annot_GO:CC/MF/BP),

* top hit MTP/protein: `ENSP00000418649`, no cross-annotation via
  mapping pipeline (likely the `freeBayes`-based hits?,
* Many Immunglobulin V regions! 
  
## non-unique mapping

* eg. Top_leading_protein `ENSP00000216962` (w/o mutation!) has
different Blast_ref_protein entries: `NP_005600.1| glycogen
phosphorylase, muscle form isoform 1` vs `NP_002853.2| glycogen
phosphorylase, brain form`.

## TODO

* GOslim: get a better handle to reduce significant hits
  for the most comprehensive annotation,
* map to transcripts to get codons,
* map to genome (splicing, accurate TSS),\newline e.g. control:
  genomic mutations and MTP vs nucleosome data.
