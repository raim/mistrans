---
title: | 
    | Relative Location of Translation Errors
date: \today
author: Rainer Machn&eacute;, Boston
output: beamer_presentation
classoption: t
bibliography: /home/raim/ref/tata.bib
citation_package: natbib
header-includes:
    - \graphicspath{{/home/raim/data/mistrans/figures/}} 
---

# Relative Location of Translation Errors?

#### Previous Data:

\scriptsize

* **Codon optimization and translation speed** were reported
  to have a 5' $\rightarrow$ 3' gradient, such that 5' ends
  are encoded by rare codons and translated slowly. 
* @Tuller2010 suggested this avoid ribosome stalling in highly transcribed
  genes, the **translation ramp hypothesis**,
* @Sejour2023elife show that this translation ramp affects mostly lowly
  expressed genes and conclude that it is likely NOT a selected property.

\vspace{3ex}
  
\phantom{noch viel mittiger}![](doc/sejour23_fig1.jpg){width=50%}\newline 
\phantom{noch viel mittiger}\tiny RRT: ribosome residence time [@Gardin2014]

\vspace{-3ex}

#### Hypothesis:

mistranslation frequency may depend on
tRNA frequency/codon optimality, and thus could
also show a N $\rightarrow$ C gradient.


# Protein Stability Changes: $\Delta\Delta G$ 

\vspace{3ex}

<!--
#### $\alpha$-helix vs $\beta$-sheet vs coil vs DNA structure

* Ancient protein folds, likely older than the genetic code,
* $\alpha$-helix ~3.6 AA per turn, ie. ~10.8 bp, ie. ~1 turn of
  the DNA helix, **coincidence or correlation**?
* $\alpha \rightarrow \beta \rightarrow$ protein aggregation and *stability*?
-->

Are we looking for stabilizing or de-stabilizing mistranslations?

* **Specific** effects, **stabilizing**:
    - Destroying degron signals,
	- ubiquitinylation sites, 
	- ...
* **Unspecific** effects, **destabilizing**:
    - Leading to aggregation,
    - $\Delta\Delta G > 0$,
	- $\alpha \rightarrow \beta$ switch.

# Relative Location of Translation Errors: Q&D Scan

\footnotesize

* `All_MTP_BP_sharedPep_quant.xlsx` $\Rightarrow$
  `Homo_sapiens.GRCh38.pep.all.fa` via Ensembl ID (`^ENSP...`).
<!--	- \scriptsize column: `Top_leading_protein` contains ENSEMBL protein IDs.-->
* Located **2.3k unique MTPs** within the Ensembl protein sequences via 
 `gregexpr` of **BP**:


![](relative_position_hist.png){width=80%}


# Relative Location of Translation Errors: Q&D Scan

\footnotesize

* `All_MTP_BP_sharedPep_quant.xlsx` $\Rightarrow$
  `Homo_sapiens.GRCh38.pep.all.fa` via Ensembl ID (`^ENSP...`).
<!--	- \scriptsize column: `Top_leading_protein` contains ENSEMBL protein IDs.-->
* Located **2.3k unique MTPs** within the Ensembl protein sequences via 
 `gregexpr` of **BP**:


![](relative_position_hist_main_peptides.png){width=80%}


# Relative Location of Translation Errors: Q&D Scan

\footnotesize

* `All_MTP_BP_sharedPep_quant.xlsx` $\Rightarrow$
  `Homo_sapiens.GRCh38.pep.all.fa` via Ensembl ID (`^ENSP...`).
<!--	- \scriptsize column: `Top_leading_protein` contains ENSEMBL protein IDs.-->
* Located **3.3k unique MTPs** within the Ensembl protein sequences via 
 `gregexpr` of **BP**:


![](absolute_position_length.png){width=80%}

# Relative Location of Translation Errors: Q&D Scan

\footnotesize

* `All_MTP_BP_sharedPep_quant.xlsx` $\Rightarrow$
  `Homo_sapiens.GRCh38.pep.all.fa` via Ensembl ID (`^ENSP...`).
<!--	- \scriptsize column: `Top_leading_protein` contains ENSEMBL protein IDs.-->
* Located **2.3k unique MTPs** within the Ensembl protein sequences via 
 `gregexpr` of **BP**:


![](relative_position_hist_small.png){width=80%}

# Relative Location of Translation Errors: Q&D Scan

\footnotesize

* `All_MTP_BP_sharedPep_quant.xlsx` $\Rightarrow$
  `Homo_sapiens.GRCh38.pep.all.fa` via Ensembl ID (`^ENSP...`).
<!--	- \scriptsize column: `Top_leading_protein` contains ENSEMBL protein IDs.-->
* Located **2.3k unique MTPs** within the Ensembl protein sequences via 
 `gregexpr` of **BP**:


![](relative_position_hist_large.png){width=80%}


# Relative Location of Translation Errors: Q&D Scan

\footnotesize

* `All_MTP_BP_sharedPep_quant.xlsx` $\Rightarrow$
  `Homo_sapiens.GRCh38.pep.all.fa` via Ensembl ID (`^ENSP...`).
<!--	- \scriptsize column: `Top_leading_protein` contains ENSEMBL protein IDs.-->
* Located **2.3k unique MTPs** within the Ensembl protein sequences via 
 `gregexpr` of **BP**:


![](relative_position_RAAS.png){width=80%}


# Relative Location of Translation Errors: Q&D Scan

\scriptsize

\vspace{3ex}

#### Q&D Method:

* File `All_MTP_BP_sharedPep_quant.xlsx`:\scriptsize
    - \scriptsize 79,604 MTP/BP pairs, mapped to blast hits,
	- column: `Top_leading_protein` contains ENSEMBL protein IDs.
* Downloaded ensembl's `Homo_sapiens.GRCh38.pep.all.fa.gz`, and
    - \scriptsize Filter 61,987 rows where `Top_leading_protein`
      appears in fasta,
    - **Note**: includes mutation suffixes `ENSP..._Q13R` are implemented
	(`segmenTools` new function `mutationPositions`).
* Take each mutated peptide only once,
    - \scriptsize
	**2,344 rows with a unique `MTP sequence`**.
* Use R's `gregexpr` to locate base peptide `BP` in protein sequence,\newline
  record position and total length of the protein.


# Relative Location of Translation Errors: Q&D Scan

\scriptsize

\vspace{3ex}

#### TODO:

* Use main peptide distribution to get a local-specific enrichment
  score and p-value,
    - \scriptsize Why are a lot of "Main peptides" not found in the ENSEMBL sequence?
* Get transcripts and map all proteins and mutations to codons.

#### Not yet tested:

* Codon context matters, e.g. *codon pairs with rare arginine codons
  and successive proline codons were among the slowest codon pairs
  translated in vivo* [@Chevance2014],
* Translation speed is higher for $\alpha$-helices than for
  $\beta$-sheets [@Thanaraj1996],
* $\alpha\rightarrow\beta$ **fold switches** [@Dalal1997] may lead to
  aggregates [@Pan1993] and change stability.
    - \scriptsize can be predicted via `JPred4` [@Mishra2021]
* IUPred3 for disordered regions [].


# Codon Optimality and Ribosome Dwell Times

\scriptsize

* Human codon usage tables:
    - \scriptsize Ensembl-based: http://genomes.urv.es/CAIcal/CU_human_nature
	via E-CAI server [@Puigbo2008],
	- NCBI-based at KEGG: http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606m
	- CAI vs tAI vs gtAI [@Anwar2023],
	- cTE vs nTE [@Pechmann2013].
	

* Spatially resolved codon optimality:
	- \scriptsize $\alpha$ v $\beta$ translation speed [@Thanaraj1996], 
	and stability; **fold switches** [@Pan1993; @Dalal1997; @Mishra2021], 
    - Supply vs. demand model, and differences **within**
	alpha-helices, related to co-folding [@Pechmann2013].

SIMPLE: $\Rightarrow$ **CAI/tAI of MTP vs. Main peptides?** $\Leftarrow$


* Ribosome profiling data - dwell times ~ codon opt.:
    - \scriptsize
	non-trivial normalization required [@OConnor2016],
	- meta-analysis yeast to mammals, and new mouse data
	[@Gobet2020] (Naef again),
	- meta-analysis regarding cycloheximide bias [@Sharma2021],
	- DeepShape and a Codon Residence Index [@Cui2019].
	
* Translation ramp:	
	- \scriptsize
	Is the translation ramp [@Tuller2010] really required for highly
	expressed genes are a side-effect of fast 5' evolution @Sejour2023elife?


# What are Optimal Codons?

![](/home/raim/Documents/pechmann13_fig1b.png){width=25%}![](/home/raim/Documents/pechmann13_fig2f.png){width=25%}![](/home/raim/Documents/pechmann13_fig4c_left.png){width=25%}![](/home/raim/Documents/pechmann13_fig5d.png){width=25%}

\scriptsize

* **Translation ramp**: initial short dip, just ~10 AA,
  required for highly expressed genes/to avoid stalling [@Tuller2010], or
  just a side-effect of fast 5' evolution [@Sejour2023elife],
* @Pechmann2013: tRNA vs. codons supply and demand model,
  and heterogeneity along CDS, even **within** $\alpha$-helices,
* @Torrent2018: **tRNA availability changes during stress**,
* Measured ribosome residence or dwell times per codon/tRNA:

![](/home/raim/Documents/gobet20_fig1c_A.png){width=32%}![](/home/raim/Documents/gobet20_fig5a.png){width=30%}

# Previous: Position Bias along Proteins

![](relative_position_hist.png){width=45%}




\footnotesize

* Position bias towards C-terminal end.


# Previous: Position Bias along Proteins

![](relative_position_hist_small.png){width=45%}![](relative_position_hist_large.png){width=45%}

![](absolute_position_length.png){width=45%}![](relative_position_RAAS.png){width=45%}


\footnotesize

* Position bias towards C-terminal end,
* Stronger in small proteins.


# MTP Functional Analysis

... requires unique set of proteins: analyze MTP per protein.

![](mtp_length.png){width=33%}![](mtp_raas_number.png){width=33%}![](mtp_raas_density.png){width=33%}

![](mtp_per_protein.png){width=50%}![](mtp_raas_density_discrete.png){width=50%}

\footnotesize 


# MTP Functions: Immune Response


(number_annot_GO:CC/MF/BP),

* top hit MTP/protein: `ENSP00000418649`, no cross-annotation via
  mapping pipeline (likely the `freeBayes`-based hits?,
* Many Immunglobulin V regions! 

GO:MF antigen binding,
GO:CC immunoglobulin complex
Go:BP immune response

# MTP Functions: Extracellular

\tiny

GO: external extracellular vesicle junction adhesion secretory

# MTP Functions: Cytoskeleton

# MTP Functions: Metabolism

\tiny

GO: glycoly* gluconeogene* mitochon* organell* vacuol* lysos* proteaso* translation ribosom*


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

* Use a consistent annotation set for all functional analyses,
* Hands-on, **add to gene table**: complete column with gene names,
  add columns with GO terms,
* GOslim: get a better handle to reduce significant hits
  for the most comprehensive annotation,
* map to transcripts to get codons,
* map to genome (splicing, accurate TSS),\newline e.g. control:
  genomic mutations and MTP vs nucleosome data.


# References {.allowframebreaks} 

\tiny

<!--
pandoc -t beamer tsour23pre.md --filter pandoc-citeproc -o tsour23pre.pdf
-->
