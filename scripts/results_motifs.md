---
header-includes: 
  - \graphicspath{{/home/raim/data/mistrans/figures/raasprofiles3}}
bibliography: /home/raim/ref/tata.bib
---


# <!--  Sequence Context: logos by AAS type -->



![](motifs/AAS_overlap.pdf){height=50%}![](kraq/peptides_AAS_AAStype_tight.pdf){height=50%}

\vspace{-3ex}

\phantom{das muss in die mitte}![](legend_pvals_horizontal.pdf){width=20%}

\tiny

Extended Data Figure Motifs: \textbf{Amino Acid and Substitution Site
Enrichment Profiles.} (a) Lysine (K) and Arginine (R), the trypsin
cleavage sites are enriched directly upstream of the substitution
sites. Tryptophane (W), methionine (M), and at lower signficance,
glycine (G) and cytosine (C) are enriched directly adjacent to
substitution sites. (b) Substitutions by glycine (G) or alanine (A)
are enriched within the 3 N-temrinal amino acids of base peptides (N1
to N3), i.e., directly after the trypsin cleavage sites (K or
R). Various substitutions involving arginine (N), methionine (M) or
glutamate (E) as either the encoded or the incorporated amino acid are
enriched distant from the N- and C-termini (>9). Only substitution
types with significant enrichments ($p\le10^{-10}$) are shown in (b).
\textcolor{red}{These observations lead to the definition of sequences
for Figure 4e, where sequence difference logos were calculated for all
substitutions by G or A within the 3 N-terminal sites of base peptides
(motif \texttt{KRAQ}), and for all substitutions that had at least one
W (\texttt{WWxWW}), M (\texttt{MMxMM}), C (\texttt{CCxCC}) or G
(\texttt{GGxGG}) within 2 positions of the substitution site. The
selection \texttt{GGxGG} showed no specific enrichments at the
subsitution site and is not shown.} TODO: omit AAS from filtering in this and logos?


# 


![](motifs/logos/logos_fromto_Q:G_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_Q:A_encoded.pdf){width=19%}
<!-- ![](motifs/logos/logos_fromto_I:G_encoded.pdf){width=19%} prev sig on intron -->
![](motifs/logos/logos_fromto_E:N_encoded.pdf){width=19%}
<!-- ![](motifs/logos/logos_fromto_S:N_encoded.pdf){width=19%} prev sig -->
![](motifs/logos/logos_fromto_T:V_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_N:M_encoded.pdf){width=19%}
<!--  ![](motifs/logos/logos_fromto_T:Q_encoded.pdf){width=19%} prev sig on exon -->

![](motifs/logos/logos_fromto_M:G_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_L:G_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_N:G_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_C:G_encoded.pdf){width=19%}

![](motifs/logos/logos_fromto_E:C_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_I:Q_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_L:Q_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_P:T_encoded.pdf){width=19%}


![](motifs/logos/logos_fromto_E:M_encoded.pdf){width=19%}
<!-- ![](motifs/logos/logos_fromto_M:D_encoded.pdf){width=19%} prev sig -->
![](motifs/logos/logos_fromto_W:D_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_H:D_encoded.pdf){width=19%}
![](motifs/logos/logos_fromto_T:D_encoded.pdf){width=19%}

\tiny

Extended Data Figure Motifs: \textbf{Sequence Difference Logos for AAS
Types.} Sequence difference logos were calculated for all unique
sequences surrounding substitution sites, subset for all observed
substitution types (encoded$\to$incorporated amino acids) using the R
package \texttt{DiffLogo} (version 2.26.0) [@Nettling2015]. Plots were
only generated if any of the positions -3 to +3 around a substitution
site showed a significant enrichment with $p<10^{-10}$, and all these
logos are shown. The logos were grouped by common patterns: cntnd. on
nextpage

#

\tiny

cntd. from previous page. (i) Substitutions by glycine or alanine
(Q$\to$A, Q$\to$G, M$\to$G, L$\to$G) are enriched directly upstream
with lysine (K) or arginine (R), i.e. they are preferentially observed
at the N-terminus of base peptides, next to the trypsin cleavage sites
(K or R). (ii) Substitutions of glutamine (Q$\to$A) or ariginine
(N$\to$G) are flanked by cysteine (C) enrichments. (iii) substitutions
E$\to$N, T$\to$V and N$\to$M are flanked by methionine (M)
enrichments. (iv) Substitutions EC, IQ, LQ and PT are flanked by
tryptophane (W) enrichments. \textcolor{red}{To generalize these
grouped observations, we defined four sequence classes used for
the difference logos in Figure 4e, where sequence difference
logos (each against all other sequences in our set) were calculated
for all substitutions by G or A within the 3 N-terminal sites of base
peptides (motif \texttt{KRAQ}), and for all substitutions that had at
least one W (\texttt{WWxWW}), M (\texttt{MMxMM}), C (\texttt{CCxCC})
or G (\texttt{GGxGG}) within 2 positions of the substitution site. The
selection \texttt{GGxGG} showed no specific enrichments at the
subsitution site and is not shown.}

# Methods

\tiny

\textbf{Sorted Enrichment Profiles.}  Amino acid enrichments around
unique substitution sites (ExtFig Motifs(a)) and the location of
substitution sites along base peptides (ExtFig Motifs(b)) were
evaluated using functions of the R package \texttt{segmenTools}
(release tag at github `AAS_preprint`), described in detail by
@Behle2022. Shortly: enrichments were calculated by cumulative
hypergeometric distribution tests (function
\texttt{clusterCluster}). Significantly enriched ($p\le10^{-5}$) amino
acids or sites were sorted along the x-axis from top to bottom
(function \texttt{sortOverlaps}). The field color scales with
$-log_2(p)$, cut at $p\le10^{-10}$ and the white text indicates
$p\le10^{-5}$ (function \texttt{plotOverlaps}). ~~The red line
indicates the cut-off below which no significant enrichments were
observed (all $p>10^{-10}$).~~

\textbf{RAAS Distribution Profiles.}
The global distribution of $\log_{10}(\text{RAAS})$ is approximately
log-normal.  Therefore, the distributions for various subsets of the
observed substitutions were compared to the distribution of all other
values from the same distribution (globally for Figures w,x, or within
cancer types for Figures y,z, by two-sided Student's t-tests (R's
\textit{t.test} function with the Welch approximation, called from our
custom function \texttt{raasProfile}). The p-value~~, sidedness (test
statistic $t<0$ or $t>0$)~~ and $\log_{10}$ of the median RAAS value
were recorded, and used to generate the RAAS distribution profiles
(function \textit{dotprofile}), where the dot size scales with
$-\log_{10}(p)$ up to a cut-off, and dot colors scale with the median
RAAS value, as indicated by individual Figure legends.

\textbf{Sequence Context Analysis.} To analyzse the sequence context
of amino acid substitution sites, we generated a custom protein
database of all proteins defined in the Ensembl human genome release
\texttt{GRCh38.110}, supplemented with all proteins with
patient-specific mutations from SHIRI's PIPELINE. All unique base
peptides were searched against this database using NCBI BLAST (version
2.15.0+), and only blast hits with $100 \%$ sequence identity and no
mismatches were considered. For these hits, we obtained the amino acid
context from the original Ensembl protein and the codon from the
matched Ensembl transcript. If the substitution itself covered a
patient-specific mutation the codon was not considered.

# Methods

\tiny


\textbf{Additional Data Sources.}  \texttt{IUPred3} [@Erdos2021] was
used to calculated a "disorder" score for each protein in our
database, with options \texttt{-a} to add the \texttt{ANCHOR2}
prediction of disordered protein interaction sites, and \texttt{-s} to
use the "medium" smoothing type. PFAM domains for each protein were
predicted with \texttt{HMMER3} (version 3.4 (Aug. 2023)), using a
reporting cutoff \texttt{-E 0.01} for a comprehensive
output. Alternative scores for disordered (\texttt{flDPnn}, @Hu2021)
and disordered protein interaction sites (\texttt{disordRDPbind},
@Peng2015), as well as a \texttt{MMSeq2}-based sequence conservation
score were obtained \textit{via} the \texttt{DescribePROT} database
[@Zhao2021].


\textbf{Sequence Difference Logos} were generated with the R package
\texttt{DiffLogo} (version 2.26.0) [@Nettling2015] with a customized
\texttt{DiffLogo} function that allows to set the y-axis limits of the
plots, and p-values were calculated with the ... omitted (set to 1) for
positions affected by our sequence selection.


# References{.allowframebreaks}

\scriptsize
