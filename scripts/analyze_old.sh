

### OUTDATED OLD CODE:

## MOTIFS: get and analyze sequences surrounding the ASS

## export sequence context of AAS
## generates file ${MISDATA}/processedData/saap_context.tsv
R --vanilla < ${THIS}/scripts/export_aas_context.R > ${MISDATA}/log/context_export.txt

## extract sequence sets from saap_context.tsv for motif analysis

## binomial distribution (hypergeo) tests of AA around AAS
## including miscleavage etc.
R --vanilla < ${THIS}/scripts/get_aas_context.R > ${MISDATA}/log/context_analysis.txt

## TODO: fuse two context scripts; and move randomization etc. to deep learning
## python script; run kplogo; select motifs and sequences for RAAS analysis.

## kplogo
## TODO: run over all from:to classes and find a way to collect and plot
## results by RAAS!
cd ~/data/mistrans/processedData/motifs
~/programs/kpLogo/bin/kpLogo seqcontext_all_.fa  -alphabet protein -o kplogo/all
~/programs/kpLogo/bin/kpLogo seqcontext_fromto_Q:G.fa -alphabet protein -o kplogo/QG
~/programs/kpLogo/bin/kpLogo seqcontext_fromto_T:V.fa -alphabet protein -o kplogo/TV
~/programs/kpLogo/bin/kpLogo seqcontext_methionine_TRUE.fa -alphabet protein -o kplogo/M
~/programs/kpLogo/bin/kpLogo seqcontext_tryptophane_TRUE.fa -alphabet protein -o kplogo/W

## TODO: log files
## with Albumin
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## without Albumin
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## unique SAAP without Albumin 
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## unique SAAP with Albumin 
sed 's/^healthy.*/healthy=FALSE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &

## same for healthy tissues
## with Albumin
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## without Albumin
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=FALSE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## unique SAAP without Albumin 
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=TRUE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla &
## unique SAAP with Albumin 
sed 's/^healthy.*/healthy=TRUE/;s/^exclude.albumin.*/exclude.albumin=FALSE/;s/^only.unique.*/only.unique=TRUE/' ${THIS}/scripts/raasprofiles3_codons.R | R --vanilla


## OLD and OBSOLETE - TODO: redo functional analysis
R --vanilla < ${THIS}/scripts/saap_analysis.R  > log/analysis.txt
## functional enrichment of SAAP-harboring proteins
R --vanilla < ${THIS}/scripts/saap_function.R


### GENERATE RESULTS SLIDES AND PAPER FIGURES


cp -a ~/Documents/sejour23_fig1.jpg .

pandoc -t beamer tsour23pre.md --filter pandoc-citeproc -o tsour23pre.pdf

## TODO: instead, generate a blastdb of all transcripts and proteins,
## blast base peptides within those proteins and record locations
