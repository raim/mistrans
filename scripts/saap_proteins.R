
## read protein fasta, iupred3, s4pred, and AAS,
## and generate per protein plots

## TODO:
## * BP+SAAP intensities, sum over positions.
## * 

## TODO: why is iupred3 for ENSP00000352639 wrong?
## sequence seems to stem from a differen short protein,
## multiple annotated proteins, e.g. ENSP00000464724.1
## TODO: why is iupred3 for ENSP00000364986 wrong?
## sequence seems to stem from ENSP00000460206.1
## TODO: ENSP00000374778 et al. stop codon missing?

source("~/work/mistrans/scripts/saap_utils.R")
source("~/programs/genomeBrowser/src/genomeBrowser_utils.R")

library(readxl)
library(viridis)
library(segmenTools)
library(seqinr)
library(Biostrings) # genetic code
##library(Rpdb)
options(stringsAsFactors=FALSE)


docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))
ttcolf <- colorRampPalette(ttcols)

#### PATHS AND FILES

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")
fig.path <- file.path(proj.path,"figures","proteins")
dir.create(fig.path)
dir.create(file.path(fig.path, "selected"))
dir.create(file.path(fig.path, "tight"))
dir.create(file.path(fig.path, "zoom"))

mam.path <- "~/data/mammary/"

## AAS 
in.file <- file.path(out.path,"saap_mapped.tsv")
## RAAS values
tmt.file <- file.path(proj.path,"originalData",
                      "All_SAAP_TMTlevel_quant_df.txt")
## AAS windows
window.file <- file.path(out.path,"aas_windows.tsv")

## protein fasta
pfasta <- file.path(out.path,"all_proteins.fa")
## coding region fasta
tfasta <- file.path(mam.path,"processedData","coding.fa")
## codon counts
codon.file <- file.path(mam.path,"processedData","coding_codons.tsv")

## protein-transcript map
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
## CDS structure
cdsmap <- file.path(mam.path,"processedData","protein_cds_structure.dat")
cdspos <- file.path(mam.path,"processedData","protein_cds_coordinates.tsv")

## protein complexes
humap.file <- file.path(mam.path,"originalData","humap2_complexes_20200809.txt")

## uniprot/ensembl mapping
uni.file <- file.path(mam.path,"originalData","uniprot_ensembl.dat")

## STRUCTURE PREDICTIONS

## s4pred
s4pred <- file.path(mam.path,"processedData",
                    "Homo_sapiens.GRCh38.pep.large_s4pred.fas.gz")
## iupred3/anchor2
iupred <- file.path(mam.path,"processedData","iupred3")

## phastcons data (retrieved by Andrew Leduc)
phastcons.file <- file.path(mam.path,"processedData","first3k_proteins.RData")

## describe prot
descrp <- file.path(mam.path,"processedData","describePROT")

## pfam hits
pfam <- file.path(mam.path,"processedData",
                  "Homo_sapiens.GRCh38.pep.large_annotations.csv")
## pfam download
pfdlf <- file.path(mam.path,"originalData", "9606.tsv.gz")
## pfam clans
clans <- file.path(mam.path,"originalData", "pfam", "Pfam-A.clans.tsv.gz")

## genome feature file
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

## PROTEOMIC PEPTIDES
mainp.file <- file.path(out.path,"main_peptides_blast.tsv")
maini.file <- file.path(dat.path,"tonsil_main_peptide_quant_df.xlsx") # tonsil
maina.file <- file.path(paste0(dat.path,"_20240627"),
                        "main_peptide_quant_df.xlsx") # other tissues

## ALL BP



## PARAMETERS

## plot
ftyp <- "png"

## only plot subset of substitutions with high RAAS
min.raas <- -1

## raas color scheme
RAAS.MIN <- -4
RAAS.MAX <-  1
colors <- "arno" # "inferno" #"rocket" # "viridis" # 

COLF <- get(colors)

## shapely/rasmol via 
shapely.cols <- character()
shapely.cols[c("D","E")]   <- "#E60A0A"  # acidic
shapely.cols[c("K","R")]   <- "#145AFF"  # basic
shapely.cols[c("S","T")]   <- "#FA9600"  # S/T - P-targets
shapely.cols[c("C","M")]   <- "#E6E600"  # sulfur-containing
shapely.cols[c("P")]     <- "#DC9682"    # stiff backbone 
shapely.cols[c("F","Y")]   <- "#3232AA"  # aromatic
shapely.cols[c("W")]     <- "#B45AB4"    # large aromatic
shapely.cols[c("H")]     <- "#8282D2"    # aromatic, nitrogen
shapely.cols[c("N","Q")]   <- "#00DCDC"  # nitrogen-containing
shapely.cols[c("G")]     <- "#EBEBEB"    # just backbone
shapely.cols[c("L","V","I")] <- "#0F820F"# branched chain
shapely.cols[c("A")]     <- "#C8C8C8"    # minimal residue

## darker grays for white background
shapely.cols[c("G")]     <- "#A0A0A0"   
shapely.cols[c("A")]     <- "#909090"    # minimal residue

## PARSE & FILTER DATA

## proteins
genes <- read.delim(feature.file)

### TODO: move raasprofiles3_init.R here as well to
### align with other analyses

## AAS
dat <- read.delim(in.file)

## remove SAAP/BP w/o protein
rm <- dat$ensembl==""
dat <- dat[!rm,]

## TMT Level RAAS Values
tmtf <- read.delim(tmt.file)


## exclude NA or Inf
## NOTE: since we don't to tests against the global
## distribution and only calculate protein-specific
## RAAS means, we don't need to filter RAAS data here.
rm <- is.na(tmtf$RAAS) | is.infinite(tmtf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from TMT level\n"))
tmtf <- tmtf[!rm,]

### FILTER: filter only tmtf for global distribution tests

## convert to logical
tmtf$Keep.SAAP <-  as.logical(tmtf$Keep.SAAP)
dat$Keep.SAAP <- !dat$IG

## TAG BAD BLAST MATCHES
dat$Keep.SAAP[dat$match!="good"] <- FALSE



## NOTE/TODO: KRAS removed via TMT leve file filter,
## since its tagged as "Potential.uncaptured.transcript".
## Here, for the protein plotter, BUT NOT in RAAS profile
## analysis, I remove this filter to also plot KRAS.
##pid="ENSP00000308495" # KRAS
##bp="LVVVGAGGVGK"

## OVERRULE TMT LEVEL KEEP SAAP COLUMN
##tmtf$Keep.SAAP <- TRUE

alls <- rbind(dat[,c("Keep.SAAP","SAAP")],
              tmtf[,c("Keep.SAAP","SAAP")])
alls <- split(alls$Keep.SAAP, alls$SAAP)

## remove if tagged so for any dataset

keep <- unlist(lapply(alls, all))
dat$keep <- keep[dat$SAAP]
tmtf$keep <- keep[tmtf$SAAP]


## remove excluded
cat(paste("removing", sum(!dat$keep),
          "tagged as false positive on protein level\n"))
dat <- dat[which(dat$keep),]

## get raw RAAS data TMT level
## remove excluded
cat(paste("removing", sum(!tmtf$keep),
          "tagged as false positive on TMT level\n"))
tmtf <- tmtf[tmtf$keep,]

## CROSS-INDEX BP/SAAP PAIRS
dat$ID <- paste(dat$BP, dat$SAAP)
tmtf$ID <- paste(tmtf$BP, tmtf$SAAP)

tmtf$IDX <- match(tmtf$ID, dat$ID)
dat$IDX <- 1:nrow(dat)

### ADD RAAS STATS to protein mapping table

## split RAAS by BP/SAAP
## TODO: use site-specific split as in raasprofiles3_proteins.R
tmtl <- split(tmtf$RAAS, paste(tmtf$BP, tmtf$SAAP))

## match to main data
tmt2dat <- match(paste(dat$BP, dat$SAAP), names(tmtl))
tmtl <- tmtl[tmt2dat]

## RAAS statistics
## TODO: add w/t test p-value
tmts <- lapply(tmtl, function(x) {
    x <- 10^x
    mn <- mean(x)
    md <- median(x)
    sd <- sd(x)
    cv <- sd/mn
    c(n=length(x),
      mean=log10(mn),
      median=log10(md),
      sd=log10(sd),
      cv=cv)
})
tmts <- as.data.frame(do.call(rbind, tmts))

## RAAS COLORS
plotdev(file.path(fig.path,paste0("legend_raas_proteins")),
        res=300, width=4, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,.3,0), tcl=-.25)
lraas.col <- selectColors(tmtf$RAAS,
                          mn=RAAS.MIN, mx=RAAS.MAX,colf=COLF,
                          n=50, plot=TRUE,
                          mai=c(.5,.5,.1,.1),
                          xlab=expression(TMT~level~log[10]*RAAS))
axis(1, at=seq(-4,4,.5), labels=FALSE)
figlabel(colors, pos="bottomleft", cex=1)
##figlabel(LAB, pos="bottomright", cex=1)
dev.off()

## add to main data
dat <- cbind(dat, tmts)

## remove SAAP without RAAS
dat <- dat[dat$n>0,]

## add RAAS coloring
intv <- findInterval(dat$median, lraas.col$breaks)
## raise to min/max
intv[intv==0] <- 1
intv[intv>length(lraas.col$col)] <- length(lraas.col$col)
dat$color <- lraas.col$col[intv]

## RAAS bins for order in plots
raas <- dat$median
raas[raas< -4] <- -3.99999999
raas[raas> 1] <- 1
raas.bins <- cut(raas, seq(-4,1,.1))
dat$RAAS.bins <- as.character(raas.bins)
raas.srt <- levels(raas.bins)
dat$raas <- raas

## LIST OF AAS BY PROTEIN
aasl <- split(dat, dat$ensembl)

## AAS windows
windows <- read.delim(window.file, header=FALSE)

## s4pred
s4p <- segmenTools::readFASTA(s4pred, grepID=TRUE)
## skip version number - TODO: use version numbers!
names(s4p) <- sub("\\.[0-9]+", "", names(s4p))

## iupred3
iufiles <- list.files(pattern=paste0(".*iupred3.tsv.gz"), path=iupred)
names(iufiles) <- sub("\\..*","", iufiles)

## phastcons : list with protein IDs
## NOTE: not complete and not used below
phastcons <- readRDS(phastcons.file)
names(phastcons) <- sub("\\..*", "", names(phastcons))

## describeProt
## NOTE: only MMSeq2 result in second column is used
dpfiles <- list.files(pattern=paste0(".*.tsv.gz"), path=descrp)
names(dpfiles) <- sub("\\.tsv.gz","",dpfiles)

## load transcript and protein fasta
## GET ENSEMBL PROTEINS - from project mammary
pfas <- segmenTools::readFASTA(pfasta, grepID=TRUE)

## get matching transcripts
tfas <- segmenTools::readFASTA(tfasta, grepID=TRUE)

## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
## reverse map transcript-protein
pamrt <- matrix(rownames(trmap), ncol=1)
rownames(pamrt) <- trmap[,1]

## rename by protein names via trmap, for quick access
## of protein-specific transcript
names(tfas) <- pamrt[names(tfas),1]

## rename phastcons by protein names
names(phastcons) <- pamrt[names(phastcons),1]


## load global gene-wise codon counts
codons <- read.delim(codon.file, row.names=1)
rownames(codons) <- pamrt[rownames(codons),1]

## codon counts in our AAS protein set
cod  <- codons[names(aasl),]

## rm NA proteins: encoded on scaffold or mtDNA
ina <- is.na(cod[,1]) 
cat(paste(sum(ina), "proteins not present in codons\n"))
cod <- cod[!ina,]

## calculate codon frequencies
codt <- apply(cod, 2, sum)
## transcript codon frequencies per AA
codl <- split(codt, GENETIC_CODE[sub(".*\\.","",names(codt))])
codl <- lapply(codl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
##codl <- codl[names(aaprop)[names(aaprop)%in%names(codl)]] ## SORT BY AA PROP
codf <- lapply(codl, function(x) x/sum(x)) # codon frequency
cods <- unlist(codf)
names(cods) <- sub("\\.","-",names(cods))
## ommit AA- prefix to use this directly with codons 
names(cods) <- sub(".*-","",names(cods))

## exons/CDS
cds <- read.delim(cdsmap, header=FALSE, row.names=1)
cdl <- lapply(strsplit(cds[,1],";"), as.numeric)
names(cdl) <- rownames(cds)

## pfam
## fixed width format, where \s was replaced with ; by sed
## NOTE: three sets of from/to: hmm, ali, env;
## env may overlap between multiple - could be fused!
pfm <- read.csv(file=pfam, sep=";", fill=FALSE, header=FALSE,comment.char="#")
pfmh <- c(
    "target", "accession" , "tlen",            
    "query"   , "accession" , "qlen",                 
    "E-value"      , "score"     , "bias",                 
    "#"            , "of"        , "c-Evalue",             
    "i-Evalue"     , "score"     , "bias",            
    "from"         , "to"        , "from",                 
    "to"           , "FROM"      , "TO",                   
    "acc")#          , "description")
colnames(pfm) <- pfmh

## pfam clans: used to collapse domains below!
pclan <- read.delim(clans, row.names=1 , header=FALSE)
pfm$clan <- pclan[sub("\\..*","", pfm$accession),"V3"]
pfm$clan[pfm$clan==""] <- pfm$target[pfm$clan==""]

use.pclan <- TRUE# FALSE
if ( !use.pclan )
    pfm$clan <- pfm$target
## split by protein
pfl <- split(pfm, pfm$query)
names(pfl) <- sub("\\.[0-9]+", "", names(pfl))


## UNIPROT <-> ENSEMBL MAPPING
## NOTE: ~90 duplicated ensembl IDs
uni2ens <- read.delim(uni.file, header=FALSE)
uni2ens[,2] <- sub("\\..*", "", uni2ens[,2]) # remove ensembl version tag
uni2ens[,1] <- sub("-.*", "", uni2ens[,1]) # remove uniprot version tag
## remove duplicate ensembl IDs
## TODO: is the list sorted? best uniprot hit?
uni2ens <- uni2ens[!duplicated(uni2ens[,2]),]


## PFAM 37.0 ANNOTATION from EBI
pfd <- read.delim(pfdlf, skip=3, header=FALSE)
colnames(pfd) <- c("seq id", "alignment start",
                   "alignment end",
                   "FROM",#"envelope start",
                   "TO",#"envelope end",
                   "hmm acc", "hmm name", "type",
                   "hmm start", "hmm end", "hmm length",
                   "bit score", "E-value", "clan")
## add ensembl ID
pfd$ensembl <- uni2ens[match(pfd[,1], uni2ens[,1]),2]
## remove missing ensembl IDs
pfd <- pfd[!is.na(pfd$ensembl),]
## split to list per ensembl ID
pfdl <- split(pfd, pfd$ensembl)

### PROTEOMICS
## main peptides detected in tonsil data
mainp <- read.delim(mainp.file, header=FALSE)
colnames(mainp) <-  c("MP","protein","identity", "mismatches",
                      "alen", "qlen", "slen", "sstart", "send", "e", "bitscore")

## INSPECT PEPTIDE INTENSITIES
## TODO: do this elsewhere and generate plots
if ( Sys.info() ["nodename"]=="exon" & interactive() ) {
    ## intensities in tissues w/o tonsil
    maina <- as.data.frame(read_xlsx(maina.file))
        
    colnames(maina)[1] <- "MP"
    
    mc <- apply(maina[,2:ncol(maina)], 1, function(x) sum(x>0))
    
    ## TODO: plot to file
    hist(mc, xlab="# tissues", ylab="# peptides")

    ## TODO: global analyses as previously
    ## * intensities per protein, compare peptides starting with A
    ##   with all peptides for that protein.
    
    ## add intensities to mainp
    sum(!mainp$MP%in%rownames(maina))
    sum(!rownames(maina)%in%mainp$MP)

    ## TODO: find merge params
    mains <- merge(mainp, maina, by="MP")

    ## tag peptides that start with A
    
    ## grep("^A",mains$MP)

    ## TODO: compare with protein intensities between ^A and other peptides
    ## ?do this only on tryptic peptides in maini??
    
}

## intensities in tonsil
maini <- as.data.frame(read_xlsx(maini.file))
rownames(maini) <- maini[,1]
maini <- maini[,2:ncol(maini)]

## filter ONLY trypsin peptides
maini <- maini[maini$Trypsin>0,]

## TODO: plot to file

plotdev(file.path(fig.path,paste0("legend_peptide_intensities")),
        res=300, width=4, height=4)
par(mai=c(.25,.5,.10,.15), mgp=c(1.4,.3,0), tcl=-.25)
mcc <- log10(maini$Trypsin)
mcc[mcc==-Inf] <- 0
mcc <- selectColors(mcc, mn=6, mx=11, reverse=TRUE,
                    xlab=expression(log[10](intensity)))
dev.off()


mc <- mcc$x.col
mc[maini$Trypsin==0] <- NA
names(mc) <- rownames(maini)

## add color to blast results
mainp$color <- mc[mainp$MP]

## ma


## AAS per protein
app <- unlist(lapply(aasl, nrow))
appc <- app
appc   [app>30 ] <- 30

if ( interactive() )
    hist(appc, breaks=1:30, xlab="AAS per protein")



## proteins of interest:
## TDP-43,
## proteasome: PSM*,
## TODO: glycolysis, AA metabolism,
## glycolysis: ENO1/2/3, GAPDH,
## histones:
## globin: HBA1/2, HBB,
## Ras/Rap: RAB*, RAP1A, RAP2A,
## collagen,
## CCR4-NOT: CNOT10,
## calmodulin: CALM1, CALML3
## TODO: RAAS profiles for protein complexes, 

POI <- c(pamrt[genes[genes$name=="KRAS","canonical"],],
         pamrt[genes[genes$name=="TARDBP","canonical"],],
         ## proteasome 20S subunit alpha 1 
         pamrt[genes[genes$name=="PSMA1","canonical"],],
         ## actins
         pamrt[genes[genes$name=="ACTA1","canonical"],],
         pamrt[genes[genes$name=="ACTB","canonical"],],
         pamrt[genes[genes$name=="ACTC1","canonical"],],
         pamrt[genes[genes$name=="ACTG1","canonical"],],
         pamrt[genes[genes$name=="ACTG2","canonical"],])

shiri.selection <- c("IFI30",
                     "PSMC5",
                     "TST",
                     "UGDH",
                     "POTEF",
                     "ENSP00000478289",
                     "EPHX1",
                     "RPSA",
                     "EEF1G",
                     "ARHGDIB",
                     "MYL12A",
                     "CKM",
                     "CORO1B",
                     "PTMA",
                     "ALDH6A1",
                     "GFAP",
                     "CLIC4",
                     "CTRB1",
                     "GPD1",
                     "AMY2A",
                     "EEF1A1",
                     "RBP2",
                     "RPL8",
                     "GSTA1",
                     "NARS1",
                     "MZB1",
                     "ASS1",
                     "DCXR",
                     "ENSP00000252485",
                     "TGFBI",
                     "ENSP00000438588",
                     "PGM1",
                     "PALLD",
                     "MYL9",
                     "ENSP00000446637",
                     "CKB",
                     "FBLN1",
                     "HLA-C")


## map gene names and uniprot
gidx <-  match(trmap[names(aasl),1], genes$MANE)
pnms <- genes$name[gidx]
pnms[is.na(pnms)] <- names(aasl)[is.na(pnms)]
names(pnms)  <- names(aasl)

puni <- genes$swissprot[gidx]
puni[is.na(puni)] <- names(aasl)[is.na(puni)]
names(puni)  <- names(aasl)

## INVESTIGATE SOME CASES
pid=names(which(pnms=="KRAS")) # several adjacent T:S from different BP
## TODO: fuse BP
pid=names(which(pnms=="KRAS")) # several adjacent T:S from different BP

pid=names(which(pnms=="S100A8"))
pids=names(pnms)[grep("S100",pnms)]

## TODO: track PSMB5 - why missing from tmtf in RAAS profile scripts?
pid=names(which(pnms=="PSMB5"))

## AHNAK: righ RAAS protein
pid=names(which(pnms=="AHNAK"))
## ACTG2: many AAS protein
pid=names(which(pnms=="ACTG2"))

## just shiri's and my collection
pids <- c(POI, names(pnms)[pnms%in%shiri.selection])


### ZOOMed PLOTs
pid=names(which(pnms=="GSTA1"))
pid=names(which(pnms=="MZB1"))
pid=names(which(pnms=="ALDOB"))
pid=names(which(pnms=="ANXA5"))
zooms <- list("ENSP00000335620"=c(1,90),#GSTA1
              "ENSP00000296511"=c(125,175),#ANXA5
              "ENSP00000497767"=c(135,250),#ALDOB
              "ENSP00000303920"=c(120,155),#MZB1
              "ENSP00000295137"=c(200,300)) #ACTG2
pids <- names(zooms)

## plot all proteins INCL. QC
pids <- names(aasl)#POI #

do.tight.plot <- FALSE
do.full.plot <- TRUE
for ( pid in pids ) {

    ffile <- file.path(fig.path, pnms[pid])
    tfile <- file.path(fig.path, "tight", paste0(pnms[pid], "_tight"))
    zfile <- file.path(fig.path, "zoom", paste0(pnms[pid], "_zoom"))
    sfile <- file.path(fig.path, "selected", pnms[pid])
    stfile <- file.path(fig.path, "selected", paste0(pnms[pid], "_tight"))

    cat(paste("getting data for", pnms[pid], pid, "\n"))

    ##if ( file.exists(ffile)l) next

    if ( !pid%in%names(aasl) ) {
        cat(paste("NOT FOUND\n"))
        next
    }
    
    ## collect all data for protein
    aas <- aasl[[pid]]
    pf <- pfl[[pid]]
    s4 <- gsub("H","0",gsub("E","/",gsub("C","-",s4p[[pid]]$seq)))
    iu <- NULL
    if ( pid %in% names(iufiles) )
        iu <- read.delim(file.path(iupred,iufiles[pid]), header=FALSE)
    dp <- NULL
    if ( puni[pid] %in% names(dpfiles) ) 
        dp <- read.delim(file.path(descrp,dpfiles[puni[pid]]), header=FALSE)


    ## secondary structure prediction as feature blocks
    s4b <- strsplit(s4p[[pid]]$seq,"")[[1]]
    s4coors <- NULL
    for ( tp in c("H","E") ) {
        tpidx <- which(s4b==tp)
        if ( length(tpidx) ) {
            tpends <- c(tpidx[which(diff(tpidx)>1)], tail(tpidx,1))
            tpstarts <- c(head(tpidx,1),tpidx[which(diff(tpidx)>1)+1])
            tpcoor <- cbind.data.frame(chr=1, start=tpstarts, end=tpends,
                                       strand=ifelse(tp=="H","+","-"),
                                       type=tp, color="#AAAAAA")
            s4coors <- rbind.data.frame(s4coors, tpcoor)
        }
    }

    ## sequences
    psq <- pfas[[pid]]$seq
    tsq <- tfas[[pid]]$seq
    phc <- NULL
    if ( pid %in% names(phastcons) )
        phc <- phastcons[[pid]]

    ## check lengths
    tlen <- nchar(tsq)/3 -1
    plen <- nchar(psq)
    slen <- nchar(s4)
    ilen <- nrow(iu)
    ## QC
    ## check lengths
    if ( length(unique(c(tlen,plen,slen,ilen)))>1 )
        cat(paste("WARNING:", pid, "lengths differ",
                  paste(unique(round(c(tlen,plen,slen,ilen),1)),
                        collapse=";"), "\n"))
    ## check sequence via translate
    if ( !is.null(iu) )
        if ( paste0(iu[,2],collapse="")!=psq )
            cat(paste("WARNING:", pid, "wrong iupred3 seq\n"))
    if ( !is.null(dp) )
        if ( paste0(dp[,1],collapse="")!=psq )
            cat(paste("WARNING:", pid, "wrong describePROT seq\n"))
    if ( is.null(tsq) ) {
        cat(paste("WARNING:", pid, "no transcript\n"))
    } else {
        tpsq <- sub("\\*$","",
                    paste0(seqinr::translate(strsplit(tsq,"")[[1]]),
                           collapse=""))
        if ( tpsq!=psq )
            cat(paste("WARNING:", pid, "wrong translation\n"))
    }
    ## check phastcons sequence
    if ( !is.null(phc) ) {
        if ( psq!=paste(phc[,"AA"], collapse="") )
            cat(paste("WARNING:", pid, "wrong phastcons sequences\n"))
    }        
    ## calculate vector of codon frequencies
    cmt <- tsq
    if ( !is.null(tsq) ) {
        cmt <- matrix(strsplit(tsq, "")[[1]],ncol=3, byrow=TRUE)
        cmt <- cods[apply(cmt,1,paste, collapse="")]
    }
    
    ## get domains
    pf <- pfl[[pid]]

    ## convert all to genomeBrowser style tables
    ## and use genomeBrowser plot functions.

    ## dummy plot coordinates for protein length
    coors <- c(chr=1,start=0, end=plen+1)
   
    ## fuse domains of same type
    ## using segmenTools segmentMerge via bedtools
    if ( !is.null(pf) ) {
        pf <- pf[,c("clan","FROM","TO","E-value")]
        colnames(pf) <- c("type","start","end","score")
        pf <- pf[order(pf$start),]
        pf <- cbind(pf,chr=1,strand=".")
        pf <- segmentMerge(pf, type="type", verb=0)
        pf$strand <- "."
        pf$name <- pf$type
        pf$color <- "#000000"
    }
    
    ## ADD DOWNLOADED PFAM/CLAN
    if ( pid%in%names(pfdl) ) {
        npf <- cbind(chr=1,
                     start=pfdl[[pid]]$FROM,
                     end=pfdl[[pid]]$TO,
                     strand=".",
                     type="official",
                     name=pfdl[[pid]]$"hmm name",
                     color="#ff0000")
        pf <- rbind(pf, npf)                    
    }

    ## ADD all "Main peptides"
    mps  <- mainp[mainp$protein==pid,, drop=FALSE]
    mpd <- NULL
    if ( nrow(mps)>0 )
        mpd <- data.frame(type="tonsil, tryptic",
                          start=mps$sstart,
                          end=mps$send,
                          chr=1, strand=".",
                          color=mps$color)

    ## AAS windows
    wins <- windows[windows[,1]==pid,]
    wind <- NULL
    if ( nrow(wins)>0 )
        wind  <- data.frame(type="AAS window",
                            start=wins[,2],
                            end=wins[,3],
                            chr=1, strand=".", color="#000000")

    mpd <- rbind(mpd, wind)
        
    ## PREPARE AAS PLOT
    ## order by RAAS such that higher RAAS will be on top of overlapping
    ##aas <- aas[order(aas$median),]
    aaco <- paste0(aas$from,":",aas$to)
    ## TODO: type 0,1,2,3 for closeby, not just identical
    aatp <- as.numeric(duplicated(round(aas$pos/10)))
    tagp <- tagDuplicates(aas$pos)
    aatp[aatp>0] <- tagp[aatp>0]
    ## NOTE: cex ~ n of measurements
    aad <- data.frame(type=round(aas$median,1), #aas$RAAS.bins, #aaco,#
                      name=aaco, #aas$to, #
                      start=aas$pos,
                      end=aas$pos,
                      chr=1, strand=".",
                      cex=log10(aas$n)+1,
                      color=aas$color,
                      codon=aas$codon,
                      raas=aas$median)
    ##aad <- aad[order(aas$median),,drop=FALSE]
    ##aad <- aad[aas$median >= min.raas, ]

    ## MOVING AVERAGE
    ## * 1: num of observations
    ## * 2: mean RAAS
    
    ## sort AAS by position
    aas <- aas[order(aas$pos),]

    mar <- man <- map <- mas <- rep(NA, plen)
    window <- 11
    use.test <- t.test
    for ( j in 1:plen ) {
        rng <- (j-floor(window/2)):(j+floor(window/2))
        rng <- rng[rng>0 & rng<= plen]
        ## get all BP/SAAP whose AAS is in this range
        idx <- which(aas$pos %in% rng)
        ## get all TMT LEVEL RAAS
        tidx <- which(tmtf$IDX%in% aas$IDX[idx])
        ## number, median RAAS and t-test p.value
        man[j] <- length(tidx)
        mar[j] <- log10(median(10^tmtf$RAAS[tidx]))
        if ( length(tidx)>2 ) {
            rt <- use.test(tmtf$RAAS[tidx], tmtf$RAAS)
            map[j] <- rt$p.value
            mas[j] <- sign(rt$statistic)
        }
    }



    ## skip plot if no RAAS is higher than minimum
    doit <- pid%in%POI | pnms[pid]%in%shiri.selection
    if ( sum(aas$median >= min.raas)==0 & !doit ) next

    cat(paste("PLOTTING", pnms[pid], pid, "\n"))

    wscale <- 1/30
    if ( plen<90 ) wscale <- 1/10
    if ( plen>2000 ) wscale <- 1/100
    mmai <- c(.05,1,.05,.1)
    
    heights <- c(.1,  # title
                 .25,  # PFAM/CLAN
                 .200,# iupred/anchor2
                 #.075,# conservation
                 .1,  # main peptides
                 .15, # secondary structure and cleavage
                 .5, # AAS type
                 .1,  # AAS context
                 .075,# AA
                 .075,# codon
                 .075 # CDS structure
                 )
    if ( do.full.plot ) {
        plotdev(ffile, width=min(c(100,plen*wscale)),
                height=2*sum(heights), type=ftyp, res=200)
        layout(mat=t(t(1:length(heights))), heights=heights)
        par(mai=mmai, xaxs="i", xpd=TRUE)
        
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        text(x=plen/2, y=1,
             labels=paste(pnms[pid],"/",puni[pid]), cex=2, xpd=TRUE)
        
        ## domains as arrows
        if ( !is.null(pf) ) {
            plotFeatures(pf, coors=coors, tcx=1.5, names=TRUE,
                         typord=TRUE, axis2=FALSE)
        } else plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        mtext("Pfam\nclan", 2, 2)
        
        ## domains as blocks
        ##plotFeatureBlocks(pf, coors=coors)
        
        

        ## iupred3/anchor2 as heatmap
        if ( !is.null(iu) ) {
            
            ## standardize MMSeq2
            cons <- dp[,2]
            cons <- (cons - 0)/(4.5-0)
            ## QC: mmseq2 must have samelength as iupred3
            if ( length(cons) != nrow(iu) )
                cons <- rep(NA, nrow(iu))
            iud <- cbind(chr=1,
                         coor=iu[,"V1"],
                         MMSeq2=cons,
                         iupred3=iu[,c("V3")],
                         anchor2=iu[,c("V4")])
            brks <- 0:100/100
            ## TODO: saver y-axis labelling in plotHeat!
            tm <- plotHeat(iud, coors=coors, breaks=brks,
                           colors=c(ttcolf(length(brks)-1)), axis2=TRUE)
        } else {
            plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        }

        ## main peptides
        mmaiaa <- mmai
        mmaiaa[1] <- 0.01
        par(mai=mmaiaa)
        if ( !is.null(mpd) ) {
            plotFeatures(mpd, coors=coors, names=FALSE, arrows=TRUE, 
                         typord=TRUE, axis2=TRUE, arrow=list(code=3, pch=NA))
        } else {
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        }
        mmaiaa <- mmai
        mmaiaa[3] <- 0.01
        par(mai=mmaiaa)
        
        ## secondary structure, K|R and AAS
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        text(1:plen, y=1, labels=strsplit(s4,"")[[1]], cex=1)
        axis(2, at=1, labels="PSSP", las=2)
        ## indicate K|R cleavage sites - TODO: omit Keil rules
        arg <- which(strsplit(psq,"")[[1]]%in%c("R"))
        if ( length(arg) )
            arrows(x0=arg, y0=1.5, y1=1, length=.05, lwd=1)
        axis(2, at=1.4, labels="K|R", las=2)
        lys <- which(strsplit(psq,"")[[1]]%in%c("K"))
        if ( length(lys) )
            arrows(x0=lys, y0=1.5, y1=1, length=.05, lwd=1.5, col=2)
        
        ## indicate ALL AAS
        arrows(x0=aas$pos, y0=0, y1=1, length=.05, lwd=3)
        arrows(x0=aas$pos, y0=0, y1=1, length=.05, lwd=2.5, col=aas$color)
        axis(2, at=.6, labels="AAS", las=2)
        
        par(mai=mmai, xpd=TRUE)
        
        ##abline(v=aas$pos, col=aas$color) ## TODO: arrows, color by RAAS etc.
        ## DETAILS for high RAAS AAS
        if ( FALSE ) { # custom with n~size
            plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
            for ( tp in unique(aad$type) )
                shadowtext(x=aad$start[aad$type==tp],
                           y=1,#-rep(tp,sum(aad$type==tp))*.25,
                           labels=aaco[aad$type==tp],
                           cex=aad$size[aad$type==tp],
                           col=aad$color[aad$type==tp])
        }
        plotFeatures(aad, coors=coors, names=TRUE, arrows=FALSE, tcx=1.5,
                     plotorder=order(aad$type), 
                     types=as.character(rev(seq(-10,10,.1))), cuttypes=TRUE,
                     typord=TRUE, axis2=TRUE)
        mtext("AAS\ntype", 2, 2)
        ## indicate ALL AAS
        ##arrows(x0=aas$pos, y0=-.1, y1=.15, length=.05, lwd=3, xpd=TRUE)
        ##arrows(x0=aas$pos, y0=-.1, y1=.15, length=.05,lwd=2.5, col=aas$color,
        ##       xpd=TRUE)
        

        ## AA context
        ## add AA sequences around AAS
        ## TODO: fuse overlapping
        pseq <- strsplit(psq,"")[[1]]
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        if ( nrow(aad)>0 )
            for ( j in 1:nrow(aad) ) {
                if ( aad$raas[j] < .1 ) next
                x <- aad$start[j]
                prng <- (x-1):(x+1)
                prng <- prng[prng>0 & prng <= plen]
                text(x=x, y=1, labels=paste(pseq[prng], collapse=""),
                     family="monospace")
            }
        mtext("context", 2, las=2, cex=.8)
        
        ## AA by amino acid color code
        mmaiaa <- mmai
        mmaiaa[1] <- 0.01
        aam <- cbind(chr=1,coor=1:length(pseq),
                     vals=match(pseq, names(shapely.cols)))
        par(mai=mmaiaa)
        plotHeat(aam, coors=coors, colors=shapely.cols)
        axis(3, at=aas$pos, labels=FALSE)
        mtext("amino acid", 2, las=2)
        ## codons as heatmap
        mmaiaa <- mmai
        mmaiaa[3] <- 0.01
        if ( !is.null(cmt) ) {
            cmm <- cbind(chr=1,coor=1:length(cmt),vals=cmt)
            brks <- 0:100/100
            par(mai=mmaiaa)
            plotHeat(cmm, coors=coors, breaks=brks,
                     colors=c(viridis(length(brks)-1)))
        ##axis(3, at=aas$pos, labels=FALSE)
        } else {
            plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        }
        mtext("codon frq.", 2, las=2)
        ## EXONS
        par(mai=mmai)
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        axis(3, at=c(1,1+cdl[[pid]]/3), tcl=1, label=FALSE)
        axis(1, at=c(1,1+cdl[[pid]]/3), tcl=1, label=FALSE)
        mtext("CDS", 2, las=2)
        
        dev.off()

        if ( doit | file.exists(paste0(sfile,".",ftyp)) )
            file.copy(paste0(ffile,".",ftyp),
                      paste0(sfile,".",ftyp), overwrite = TRUE)
    }
    
### TIGHT PLOT
    if ( do.tight.plot ) {
        cat(paste("TIGHT PLOT", pnms[pid], pid, "\n"))

        ## plot width
        zlen <- plen
        midp <- plen/2
        
        wscale <- 1/30
        if ( zlen<90 ) wscale <- 1/10
        if ( zlen>2000 ) wscale <- 1/100
        mmai <- c(.05,1,.05,.1)
        
        heights <- c(.1,  # title
                     .15,  # PFAM/CLAN
                     .2,  # RAAS CURVES
                     .2,  # main peptides
                     .1, # secondary structure and cleavage
                     .5, # AAS type
                     .075, # CDS structure
                     .1 # coordinates
                     )
        plotdev(tfile, width=min(c(100,zlen*wscale)),
                height=2*sum(heights), type=ftyp, res=200)
        layout(mat=t(t(1:length(heights))), heights=heights)
        par(mai=mmai, xaxs="i", xpd=TRUE)
        
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        text(x=midp, y=1,
             labels=paste(pnms[pid],"/",puni[pid]), cex=2, xpd=TRUE)
        
        ## domains as arrows
        if ( !is.null(pf) ) {
            plotFeatures(pf, coors=coors, tcx=1.5, names=TRUE,
                         typord=TRUE, axis2=FALSE)
        } else plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        mtext("Pfam\nclan", 2, 2)
        
        ## MOVING AVERAGE
        ## TODO: cex proportional to p-value
        ## 
        plot(1:plen, man, type="l",
             xlim=c(coors[2:3]), xlab=NA, ylab=NA, axes=FALSE)
        axis(2, las=2)
        abline(h=0, col="gray", xpd=FALSE)
        mtext("# RAAS", 2, 2.5, cex=.8, las=2)
        
        ## main peptides
        mmaiaa <- mmai
        mmaiaa[1] <- 0.01
        par(mai=mmaiaa)
        if ( !is.null(mpd) ) {
            plotFeatures(mpd, coors=coors, names=FALSE, arrows=TRUE,
                         types=c("tonsil, tryptic","AAS window"),
                     typord=TRUE, axis2=TRUE, arrow=list(code=3, pch=NA))
        } else {
            plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        }
        mmaiaa <- mmai
        mmaiaa[3] <- mmaiaa[1] <- 0.01
        par(mai=mmaiaa)
        
        ## secondary structure, K|R and AAS
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        text(1:plen, y=1, labels=strsplit(s4,"")[[1]], cex=1)
        axis(2, at=1, labels="PSSP", las=2)
        ## indicate K|R cleavage sites - TODO: omit Keil rules
        arg <- which(strsplit(psq,"")[[1]]%in%c("R"))
        if ( length(arg) )
        arrows(x0=arg, y0=1.5, y1=1, length=.05, lwd=1)
        axis(2, at=1.4, labels="K|R", las=2)
        lys <- which(strsplit(psq,"")[[1]]%in%c("K"))
        if ( length(lys) )
            arrows(x0=lys, y0=1.5, y1=1, length=.05, lwd=1.5, col=2)
        
        ## indicate ALL AAS
        ##arrows(x0=aas$pos, y0=0, y1=1, length=.05, lwd=1.5)
        ##arrows(x0=aas$pos, y0=0, y1=1, length=.05, lwd=1, col=aas$color)
        ##axis(2, at=.6, labels="AAS", las=2)
        
        mmaiaa <- mmai
        mmaiaa[3] <- mmaiaa[1] <- 0.01
        par(mai=mmaiaa, xpd=TRUE)
        
        plot(aad$start, aad$raas, pch=19, col=NA,
             xlim=c(coors[2:3]), xlab=NA, ylab=NA, axes=FALSE)
        axis(2)
        lines(1:plen, mar, lwd=2)
        points(aad$start, aad$raas, pch=19, col=aad$color, cex=3.5)
        points(aad$start, aad$raas, pch=1, col="white", cex=3.5)
        mtext(expression(log[10](RAAS)), 2, 2, cex=.8)
        if ( FALSE ) {
            plotFeatures(aad, coors=coors, names=TRUE, arrows=FALSE, tcx=1.5,
                         plotorder=order(aad$type), 
                         types=rev(seq(-10,10,.1)), cuttypes=TRUE,
                         typord=TRUE, axis2=FALSE)
            mtext("AAS\ntype", 2, 2)
        }
    
        ## EXONS
        par(mai=mmai)
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        axis(3, at=c(1,1+cdl[[pid]]/3), tcl=1, label=FALSE)
        axis(1, at=c(1,1+cdl[[pid]]/3), tcl=1, label=FALSE)
        mtext("CDS", 2, las=2)
        
        ## coordinates
        amai <- mmai
        par(mai=amai, xpd=FALSE)
        plot(1, xlim=c(coors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
        xat <- pretty(coors[2:3])
        if ( 0%in%xat )
            xat[xat==0] <- 1
        axis(1, tcl=.25, at=xat, mgp=c(0,-1.25,0))
        mtext("position", 2, .5, las=2, cex=.7)

        dev.off()
    }

    
### ZOOM PLOT

    cat(paste("ZOOM PLOT", pnms[pid], pid, "\n"))

    zlen <- plen
    midp <- plen/2
    zoors <- coors
    if ( pid %in% names(zooms) ) {
        zoors <- c(chr=1,
                   start=zooms[[pid]][1], end=zooms[[pid]][2]+1)
        zlen <- zoors[3] - zoors[2]
        midp <- mean(zoors[2:3])

    }

    ## remove pfam clans
    pf <- pf[pf$color=="#000000",]
   
    wscale <- 1/25
    if ( zlen<90 ) wscale <- 1/10
    if ( zlen>2000 ) wscale <- 1/100
    mmai <- c(.05,.5,.05,.1)
    
    heights <- c(.1,  # title
                 .15,  # PFAM/CLAN
                 .15, # secondary structure 
                 .3,  # RAAS CURVES
                 .4, # RAAS values
                 .1 # position labels
                 )
    plotdev(zfile, width=min(c(100,zlen*wscale)),
            height=2*sum(heights), type=ftyp, res=200)
    layout(mat=t(t(1:length(heights))), heights=heights)
    par(mai=mmai, xaxs="i", xpd=FALSE)
    
    plot(1, xlim=c(zoors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
    text(x=midp, y=1,
         labels=paste(pnms[pid],"/",puni[pid]), cex=2, xpd=TRUE)
    
    ## domains as arrows
    if ( !is.null(pf) ) {
        plotFeatures(pf, coors=zoors, tcx=1.5, names=TRUE,
                     typord=TRUE, axis2=FALSE)
    } else plot(1, xlim=c(zoors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
    mtext("Pfam", 2, 1.5)
        
    ## secondary structure as blocks
    par(mai=mmai)
    if ( !is.null(s4coors) ) {
        plotFeatureBlocks(data=s4coors, coors=zoors, axis2=FALSE, abline=TRUE)
        axis(2, at=c(.75,-.75), labels=c(expression(alpha),
                                       expression(beta)), las=2, cex.axis=1.2)
        ##mtext("structure", 2, 2, cex=.7, las=2)
    } else plot(1, xlim=c(zoors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)

    ## RAAS MEASUREMENT FREQUENCY
    plot(1:plen, man, type="l",
         xlim=c(zoors[2:3]), xlab=NA, ylab=NA, axes=FALSE)
    axis(2, las=2)
    abline(h=0, col="gray", xpd=FALSE)
    mtext("# RAAS", 2, 2, cex=.8)#, las=2)

    ## RAAS VALUES
    mmaiaa <- mmai
    mmaiaa[3] <- mmaiaa[1] <- 0.01
    par(mai=mmaiaa, xpd=TRUE)
    plot(aad$start, aad$raas, pch=19, col=NA,
         xlim=c(zoors[2:3]), xlab=NA, ylab=NA, axes=FALSE)
    axis(2)
    lines(1:plen, mar, lwd=2)
    points(aad$start, aad$raas, pch=19, col=aad$color, cex=2.5)
    points(aad$start, aad$raas, pch=1, col="white", cex=2.5)
    mtext(expression(log[10](RAAS)), 2, 2, cex=.8)

    ## coordinate axis
    amai <- mmai
    par(mai=amai, xpd=FALSE)
    plot(1, xlim=c(zoors[2:3]), col=NA, axes=FALSE, xlab=NA, ylab=NA)
    xat <- pretty(zoors[2:3])
    if ( 0%in%xat )
        xat[xat==0] <- 1
    axis(1, tcl=.25, at=xat, mgp=c(0,-1.25,0))
    mtext("position", 2, .5, las=2, cex=.7)
    dev.off()
    
    ## note that the file.exists part allows to update interesting proteins
    ## that were manually copied to the selected folder
    if ( doit | file.exists(paste0(stfile,".",ftyp)) )
        file.copy(paste0(tfile,".",ftyp),
                  paste0(stfile,".",ftyp), overwrite = TRUE)

}

## HOTSPOT SCAN



##library(rtracklayer)
##bw = import.bw('~/data/mammary/originalData/hg38.phastCons7way.bw')
##library(Rpdb)
##pdbf <- file.path(mam.path, "originalData", "afold",
##                  "AF-P25786-F1-model_v4.pdb")
##x = read.pdb(pdbf)
##visualize(x, type="p")

