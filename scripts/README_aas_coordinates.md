

The file `aas_coordinates.tsv` provides summary results for peptides
with identified amino acid substitutions. The encoded "base peptides"
(BP) have been mapped to the human genome version GCRh38, ensembl
release 110, by blast, and only BP with 100% identity matches to
ensembl proteins are reported. All coordinates (position) are 1-based.

# Columns:

* `BP`: "base peptide" as encoded by the genome or by patient-specific
  transcripts, these were blasted against the Ensembl proteome fasta
  file (GRCh38, release 110), supplemented with mutated proteins
  detected by one of two variant calling pipelines.
* `SAAP`: "substituted amino acid peptide", same as BP with one AA difference.
* `site`: position of the substituted amino acid in the BP and the SAAP,
* `fromto`: the encoded and substituted amino acid (difference between
  BP and SAAP) in the format `encoded:substituted`.
* `codon`: the codon for the encoded amino acid (`from`) in 
   the matching ensembl `transcript`.
* `RAAS`: the log10 of the abundance ratio of the SAAP and BP
  peptides, this is a median if the same BP and SAAP were measured in
  multiple samples (`RAAS.n`>1).
* `RAAS.n`: number of measurement used for the median.
* `Datasets`, `Tissues`: all cancer and tissue types were the BP and
  SAAP were measured.
* `name`: human-readable gene name (`Name=` tag in the gff3 description column.
* `protein`: ensembl protein ID of the main selected blast hit of the the BP.
* `protein.position`: position of the substituted amino acid in the protein
  indicated in column `protein`.
* `transcript`: the ensembl ID for the transcript encoding for
  the protein in column ensembl.
* `transcript.position`: the position of the `codon` (2nd position) in
  the ensembl transcript sequence as provided in the fasta file:
  `Homo_sapiens.GRCh38.cdna.all.fa.gz`.
* `gene`: the ensembl ID of the gene that encodes for the transcript and protein.
* `chr`, `coor`, `strand`: the chromosome, chromosome position and
  strand of the `codon` (2nd position).
