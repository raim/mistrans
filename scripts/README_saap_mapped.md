

The file `saap_mapped.tsv` provides the blast-based re-mapping of base
peptides (BP) to the ensembl human genome version GCRh38, release 110,
and including proteins with patient-specific single amino acid
substitutions detected by one of two variant calling pipelines.

Most interesting for external use are the columns `pos`, `tpos` and
`chr`/`coor`/`strand` which provide the exact positions of the
detected amino acid substitutions in the ensembl protein, transcript
and genome sequences. See below for a detailed description of all columns.

# Columns:

* `BP`: "base peptide" as encoded by the genome or by patient-specific
  transcripts, these were blasted against the Ensembl proteome (GRCh38, release
  110) fasta, supplemented with mutated proteins detected by one of two
  variant calling pipelines.
* `SAAP`: "substituted amino acid peptide", same as BP with one AA difference.
* `protein`: ensembl protein ID of the main selected blast hit, which includes
  proteins with patient-specific substitutions, indicated by suffixes (`_Q23G`);
  Note, that these names were shortened to 47 chars but are unique, full original
  names (from the variant calling pipeline) are available if required.
* `identity`, ..., `bitscore`: blast output.
* `ensembl`: the ensembl protein ID without suffixes.
* `gene`: the ensembl ID of the gene that encodes for the transcript and protein.
* `transcript`: the ensembl ID for the transcript encoding for
  the protein in column ensembl.
* `MANE`: the ensembl transcript ID for the MANE subset of transcripts
  associated with the same gene, mostly but not always identical with
  the `transcript` column. If it differs the blast hit was NOT for the
  MANE isoform of the same gene.
* `MANE.protein`: the ensembl ID for the protein encoded by the MANE transcript.
* `numGO`: number of GO terms associated with the gene.
* `match`: classification of the blast hit quality, `good` where
  full-length matches with 100% identity and 0 mismatches, only those
  should be used for analysis of protein, transcript or genome
  coordinates. Note, that we were not able to account for mutations
  detected by the second variant calling pipeline. Doing this should
  remove the `bad` (<100 % identity, >0 mismatches, or not a full
  length hit) and `wrong` (>2 mismatches) blast hits.
* `nogene`: no ensembl gene could be mapped.
* `IG`: the blast hit is tagged as an immunoglobulin (pattern
   `biotype=IG_.*` in the ensembl gff3 file); since these are expected
   to be mutated in individual immune cells but not detected by
   patient-specific genome or transcript sequencing these ARE EXCLUDED
   from all analyses.
* `albumin`: contains "albumin" in the gff3 file's description column.
* `globin`:  contains "globin" in the gff3 file's description column.
* `extracellular`: the gene is is annotated with the GO term
   `GO:0005576` for "extracellular region".
* `exclude`: outdated tag to exclude the BP/SAAP from analyses, if
  either `IG`, `albumin` or `globin` is TRUE.
* `name`: human-readable gene name (`Name=` tag in the gff3 description column.
* `site`: position of the substituted amino acid in the BP and the SAAP,
* `pos`: position of the substituted amino acid in the protein
  indicated in column `protein` (blast hit).
* `len`: length of the protein.
* `rpos`: `pos`/`len`, relative position of the substituted amino acid
  in the protein.
* `from`: the encoded amino acid (in the BP).
* `to`: the substituted amino acid (in the SAAP).
* `codon`: the codon coding for the encoded amino acid (`from`) in 
   the matching ensembl `transcript`.
* `tpos`: the position of the `codon` (2nd position) in the ensembl
  transcript sequence as provided in the fasta file:
  `Homo_sapiens.GRCh38.cdna.all.fa.gz`.
* `chr`, `coor`, `strand`: the chromosome, chromosome position and
  strand of the `codon` (2nd position).
* `s4pred`: protein secondary structure prediction by S4PRED at the
  substituted amino acid; C: coil, E: beta-sheet, H: alpha helix.
* `C.protein`: total number of S4pred-predicted C in the `protein`.
* `E.protein`: total number of S4pred-predicted E in the `protein`.
* `H.protein`: total number of S4pred-predicted H in the `protein`.
* `iupred3`: iupred3-based disordered score at the substituted amino acid.
* `iupred3.protein`: mean iupred3 score of the whole protein.
* `anchor2`: anchor2 (via iupred3) score of disordered protein
  interaction at the substituted amino acid.
* `anchor2.protein`: mean anchor2 score of the whole protein.
* `MMSeq2`: MMSeq2-based sequence conservation score, via the
  describePROT database.
* `ASAquick`: surface accessibility score, via the describePROT database.
* `DisoRDPbind`: disordered protein binding score (similar to
  `anchor2`), via the describePROT database.
* `SCRIBER`:  protein-binding prediction score, via the describePROT database.
* `flDPnn`: disordered score (similar to `iupred3`), via the
  describePROT database.
* `AA`: -25 to +25 amino acid sequence around the substituted amino
  acid, `-` fills up positions beyond the protein ends, such that the
  central position still corresponds to the substituted amino acid.
* `NT`: the transcript sequence coding for the amino acids in `AA`.
* `pfam`: custom `hmmer3`-based PFAM domains that overlap with
  substituted amino acid.
* `clan`: PFAM clans of the domains in column `pfam`.
* `pfam.ebi`: PFAM domains downloaded from InterPro that overlap with
  the substituted amino acid.
* `clan.ebi`: PFAM clans of the domains in column `pfam.ebi`.
