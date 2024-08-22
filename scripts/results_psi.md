---
fontfamily: helvet
header-includes: 
  - \renewcommand{\familydefault}{\sfdefault}
  - \graphicspath{{/home/raim/data/mistrans/figures/raasprofiles3}}
  - \usepackage[left=2cm,right=2cm,bottom=2cm,top=2cm]{geometry}
  - \pagenumbering{gobble}
bibliography: /home/raim/ref/tata.bib
---

![](rna/psi_hypergeotest.pdf){height=35%}![](rna/psi_raas_codons.pdf){height=35%}

\vspace{-48ex}\Large\textbf{a}\hspace{31ex}\textbf{b} 

\vspace{33ex}

![](rna/psi_raas_codon1_A549.pdf){height=35%}![](rna/psi_raas_codon2_NTERA.pdf){height=35%}

\vspace{-33ex}\Large\textbf{c}\hspace{31ex}\textbf{d}

<!--
\normalsize

Caption: (a) 180 sites in AAS codons significantly overlap with sites of
modified U residues in transcripts, as defined by nanopore sequencing
[@McCormick2024pre]. To further evaulate the p-value of a cumulative
hypergeometic distribution test (see Methods for details) 
calculated p-value, we plotted the p value distribution over different
counts in panel and indicate the real count and its associated p-value
in red (a)

Methods:

180 sites in AAS codons overlap with sites of modified U residues in
transcripts, as defined by nanopore sequencing [@McCormick2024pre]. To
test whether this overlap is significant, we mapped all modified
U sites from 6 cell lines to the coding sequences of ENSEMBL MANE
transcripts, yielding 39,723 unique sites in 6,771 transcripts. We
also reduced the set of AAS from 7,069 to 6,967 sites that map to MANE
transcripts. We then counted all U in the union of all 7,471
transcripts with AAS and/or $\psi$ sites as the total set, and all
4,250 U in codons at AAS sites as the test set, and asked whether it
is significant to find at least 180 U by a cumulative hypergeometric
distribution test ($p[X>179]$, in R: \texttt{phyper(q=180-1, m=39723,
n=2879635-39723, k=4250, lower.tail=FALSE}). To further evaulate the
calculated p-value, we plotted the p value distribution over different
counts in panel and indicate the real count and its associated p-value
in red (a). Notably, we find no global relation of the measured
fraction of modified U (column \texttt{mm.DirectMINUSmm.IVT} in the
original data set) to the median RAAS of these sites (b). Scanning for
correlations over all six cell lines and three codon positions, we
find only two subsets that show a significant positive (c) or negative
(d) correlation, but both at very low fraction of modified sites, and
the authors of the original data set would not consider these sites as
true positive pseudouridylation sites.
-->
