
library(Biostrings) # for blosum62
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
source("~/work/mistrans/scripts/raasprofiles3_init.R")

pfig.path <- file.path(fig.path,"proteins")
dir.create(pfig.path, showWarnings=FALSE)
          
### START ANALYSIS


#### PROTEINS and COMPLEXES

## TODO: align use of raasProfile vs. listProfile

## MEDIAN SITE AND PROTEIN RAAS

### median raas per unique mane protein site

sitl <- split(tmtf$RAAS, paste(tmtf$mane, tmtf$pos))
site <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## mod. column names and add protein and site info
colnames(site) <- paste0("RAAS.", colnames(site))
site$mane <- sub(" .*", "", rownames(site))
site$pos <- as.numeric(sub(".* ", "", rownames(site)))

## add gene names
site$name <- ens2nam [site$mane]

## add uniprot id
site$uniprot <- unlist(lapply(ens2u[site$mane], paste, collapse=";"))

## RAAS COLOR
site$RAAS.color <- num2col(site$RAAS.median,
                           limits=c(RAAS.MIN, RAAS.MAX), colf=arno)

## protein median of site median RAAS
sitl <- split(site$RAAS.median, site$mane)
pbstat <- listProfile(sitl, y=tmtf$RAAS, use.test=use.test, min=3)

## protein median raas per protein w/o site-specific median first
ptl <- split(tmtf$RAAS, tmtf$mane)
ptstat <- listProfile(ptl, y=tmtf$RAAS, use.test=use.test, min=3)

if ( interactive() ) {
    ## test alternative measures of consistent RAAS
    ## CV:
    dense2d(site$RAAS.cv, -log10(site$RAAS.p.value),
            xlab="CV", ylab=expression(-log10(p)))
    dense2d(site$RAAS.cv, site$RAAS.sd, xlab="CV", ylab="SD")
    hist(site$RAAS.sd)
}

## ORDER SITES

## order proteins by RAAS
pbstat$rank <- rank(pbstat$median)

## make sure its ordered
## TODO: ORDER PROTEINS BY SIZE OR NUMBER OF AAS
site <-site[order(site$mane, site$pos),]

## order by protein RAAS rank (for hotspot plot)
site$rank <- pbstat[site$mane,"rank"]
site <-site[order(site$rank, site$pos),]

##
nsites <- as.numeric(sub(".*\\.","",tagDuplicates(site$mane)))
nsites[is.na(nsites)] <- 1
site$n <- nsites

### WRITE OUT SITE FILE
## TODO: do this upstream and make site file central protein level file
write.table(site, file=file.path(out.path, "aas_sites.tsv"), sep="\t",
            na="", row.names=FALSE, quote=FALSE)

if ( interactive() ) { # inspect some proteins
    site[grep("PSMA1", site$name),]
}


## list sites per protein

## ONLY PROTEINS>1
multip <- split(site$pos, site$mane)
multip <- names(lengths(multip)[lengths(multip)>1])
sites <- site[site$mane%in%multip,]
cpos <- cumsum(sites$pos)

## AAS and RAAS along all concatenated proteins

intv <- findInterval(sites$RAAS.median, vec=vbrks)
## raise to min/max
intv[intv==0] <- 1
intv[intv>length(vcols)] <- length(vcols)
raas.col <- vcols[intv]
plotdev(file.path(pfig.path,"hotspots_problem"), type=ftyp,
        width=10, height=3, res=200)
layout(t(t(1:3)), heights=c(.1,.4,.6))
par(mai=c(0.1,.5,.05,.1),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
## TODO: number of RAAS
par(mai=c(0,.5,.05,.1))
dense2d(cpos, rep(1, length(cpos)), axes=FALSE, ylab=NA, xlab=NA)
par(mai=c(0.1,.5,.01,.1))
plot(x=cpos, y=sites$RAAS.n, col=raas.col, pch=19, cex=.5,
     ylab="#RAAS per site",
     axes=FALSE, log="y", xpd=TRUE)
axis(2)
axis(3, at=cpos[sites$n==1], col=2, col.axis=2, labels=FALSE, tcl=.25)
par(mai=c(.35,.5,0.1,.1))
plot(cpos, sites$RAAS.median, type="h", col=raas.col, xpd=TRUE,
     xlab="position in concatenated proteins with >1 AAS, sorted by protein RAAS", ylab=xl.site)
points(cpos, sites$RAAS.median, col=raas.col, pch=19, cex=.2, xpd=TRUE)
## TODO: protein lines
axis(3, at=cpos[sites$n==1], col=2, col.axis=2, labels=FALSE, tcl=.5)
dev.off()

## MEDIAN RAAS FOR EACH UNIQUE PROTEIN POSITION



## median RAAS vs. number of sites
if ( interactive() )
    dense2d(pbstat$median, log10(pbstat$n))

## volcano
plotdev(file.path(pfig.path,paste0("proteins_volcano_sites")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
res <- volcano(pbstat, value="median",
               p.txt=-log10(0.01), v.txt=c(-2,-2), cut=5, mid=-2,
               ids=pnms, xlab=xl.prots)
mtext("protein RAAS, median of site medians", 3,0)
dev.off()


## PROTEIN MEDIAN RAAS per protein w/o site-specific median first

plotdev(file.path(pfig.path,paste0("proteins_volcano_all")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
res <- volcano(ptstat, value="median",
               p.txt=12, v.txt=c(-2,-1), cut=50, mid=0,
               ids=pnms, xlab=xl.prota)
mtext("protein RAAS, median of all RAAS", 3,0)
dev.off()

## compare median of all RAAS vs. median of site median RAAS
if ( interactive() ) {
    plotCor(pbstat$median, ptstat$median,
            xlab="median of site median RAAS", ylab="median of all RAAS")
}


## PROTEIN HALF-LIVES, @Mathieson2018
hlvd <- readxl::read_xlsx(math18.file)

## mean half-live over all replicates and cell types
## TODO: consider distribution
cidx <- grep("half_life", colnames(hlvd), value=TRUE)
hlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
names(hlv) <- unlist(hlvd[,1])


xl.hlf <- expression(protein~"half-life"/h)

## halflives site
plotdev(file.path(pfig.path,paste0("protein_halflives_site")),
        type=ftyp, res=300, width=3.5,height=3.5)
idx <- match(pnms[rownames(pbstat)], names(hlv))
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
plotCor(pbstat$median, log10(hlv[idx]),
        xlab=xl.prots, ylab=xl.hlf, axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
## TODO: log tick marks
##axis(2, at=rep(1:10, 4) * 10^rep(0:3, each=10), labels=FALSE)
dev.off()

## halflives all
plotdev(file.path(pfig.path,paste0("protein_halflives_all")),
        type=ftyp, res=300, width=3.5,height=3.5)
idx <- match(pnms[rownames(ptstat)], names(hlv))
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
plotCor(ptstat$median, log10(hlv[idx]), signif = 1,
        xlab=xl.prota, ylab=xl.hlf, axes=FALSE)
axis(1)
axis(2, at=1:10, labels=10^(1:10))
dev.off()

if ( interactive() ) {
    ## TODO: compare iupred3 and half-lives for all proteins, not just RAAS
    idx <- match(dat$name,names(hlv))
    plotCor(dat$iupred3.protein, log10(hlv[idx]))
}

for ( ctype in c("Bcells", "NK", "hepatocytes", "monocytes", "neurons") ) {
    cidx <- grep("half_life", colnames(hlvd), value=TRUE)
    cidx <- cidx[grep(ctype, cidx, ignore.case=TRUE)]
    hlv <- apply(hlvd[,cidx], 1, mean, na.rm=TRUE)
    names(hlv) <- unlist(hlvd[,1])
    
    plotdev(file.path(pfig.path,paste0("protein_halflives_sites_",ctype)),
            type=ftyp, res=300, width=3.5,height=3.5)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
    idx <- match(pnms[rownames(pbstat)], names(hlv))
    plotCor(pbstat$median, log10(hlv[idx]),
            xlab=xl.prots, ylab=xl.hlf, axes=FALSE)
    axis(1)
    axis(2, at=1:10, labels=10^(1:10))
    figlabel(ctype, pos="bottomleft")
    dev.off()
    
    plotdev(file.path(pfig.path,paste0("protein_halflives_all_",ctype)),
            type=ftyp, res=300, width=3.5,height=3.5)
    idx <- match(pnms[rownames(ptstat)], names(hlv))
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
    plotCor(ptstat$median, log10(hlv[idx]), signif = 1,
        xlab=xl.prota, ylab=xl.hlf, axes=FALSE)
    axis(1)
    axis(2, at=1:10, labels=10^(1:10))
    figlabel(ctype, pos="bottomleft")
    dev.off()
}

### PROTEIN RAAS vs. 20S Targets - @Pepelnjak2024
## TODO: p20 data is on peptide level with many duplicated prots
    
p20 <- data.frame(readxl::read_xlsx(pepe24.file, sheet=2))

xl.20s <- expression(median~log[2]("20S"/control))

## MEDIAN LG2FC PER PROTEIN
p20s <- p20[!is.na(p20$Log2.FC.),]
p20l <- split(p20s[,c("Log2.FC.")],
              p20s[,"UniprotID"])
p20p <- unlist(lapply(p20l, median, na.rm=TRUE))


## get uniprot mapping
## TODO: akugb with uni2e and ens2u  in init script.
p2ens <- uni2ens
p2ens <- p2ens[p2ens[,2]%in%tmtf$mane,]
p2ens <- p2ens[p2ens[,1]%in%names(p20p),]

cat(paste(sum(duplicated(p2ens[,1])),"uniprot with multiple ensembl\n"))
cat(paste(sum(duplicated(p2ens[,2])),"uniprot with multiple ensembl\n"))
p2el <- unlist(split(p2ens[,2], p2ens[,1]))

e2u <- names(p2el)
names(e2u) <- p2el

if ( length(unique(lengths(p2el)))>1 )
    stop("non-unique uniprot 2 ensembl mapping")

##
lg2fc <- p20p
names(lg2fc) <- p2el[names(lg2fc)]
lg2fc <- lg2fc[rownames(pbstat)]
plotdev(file.path(pfig.path,paste0("p20_protein_lg2fc_site")),
        type=ftyp, res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
plotCor(pbstat$median, lg2fc, xlab=xl.prots,
        ylab=xl.20s, signif=2, round=2)
dev.off()
lg2fc <- p20p
names(lg2fc) <- p2el[names(lg2fc)]
lg2fc <- lg2fc[rownames(ptstat)]
plotdev(file.path(pfig.path,paste0("p20_protein_lg2fc_all")),
        type=ftyp, res=300, width=3.5,height=3.5)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
plotCor(ptstat$median, lg2fc, xlab=xl.prota,
        ylab=xl.20s, signif=2, round=2)
dev.off()

## brute force
if ( FALSE ) {
    ## get uniprot mapping
    p2ens <- uni2ens
    p2ens <- p2ens[p2ens[,2]%in%tmtf$mane,]
    p2ens <- p2ens[p2ens[,1]%in%p20[,1],]
    
    cat(paste(sum(duplicated(p2ens[,1])),"uniprot with multiple ensembl\n"))
    cat(paste(sum(duplicated(p2ens[,2])),"uniprot with multiple ensembl\n"))
    p2el <- unlist(split(p2ens[,2], p2ens[,1]))
    p20$ensembl <- p2ens[match(p20[,1],p2ens[,1]),2]
    
    
    log2fc <- p20$Log2.FC.[match(rownames(pbstat), p20$ensembl)]
    logp <- -log10(p20$adj.pvalue[match(rownames(pbstat), p20$ensembl)])
    slogp <- logp*sign(log2fc)
    plotdev(file.path(pfig.path,paste0("p20_peptide_volcano")),
            type=ftyp, res=300, width=3.5,height=3.5)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
    dense2d(log2fc, logp, xlab=xl.20s, ylab="-log10(p)")
    dev.off()
    
    plotdev(file.path(pfig.path,paste0("p20_peptide_protein_raas_pval_all")),
            type=ftyp, res=300, width=3.5,height=3.5)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
    plotCor(pbstat$median, slogp, xlab=xl.prots,
            ylab="p20/tryptic: sign(lg2fc)*-log10(p)")
    dev.off()
    
    plotdev(file.path(pfig.path,paste0("p20_peptide_protein_raas_lg2fc_all")),
            type=ftyp, res=300, width=3.5,height=3.5)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
    plotCor(pbstat$median, log2fc, xlab=xl.prots,
            ylab=xl.20s)
    dev.off()
}


### PROTEIN RAAS vs. IUPRED

## get whole protein mean iupred3 score
iu3 <- split(dat$iupred3.protein, dat$MANE.protein)
## QC: all protein level
table(unlist(lapply(iu3, function(x) length(unique))))
iu3 <- unlist(lapply(iu3, unique))

## NOTE: no correlation to protein-wide iupred3
plotdev(file.path(pfig.path,paste0("protein_iupred3_all")),
        type=ftyp, res=300, width=3.5,height=3.5)
idx <- match(rownames(ptstat), names(iu3))
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.25,0), tcl=-.25)
plotCor(ptstat$median, log(hlv[idx]), signif = 1,
        xlab=xl.prota,
        ylab="iupred3")
dev.off()

### PROTEIN WINDOWS - 50 AA WINDOWS


## add position windows to raas table
## see cbind(1:200,round(seq(1:200)/50),
##           round((25+seq(1:200))/50))
## for window assignment of positions
tmtf$win1 <- round(tmtf$pos/50) # first is 1:25
tmtf$win2 <- round((tmtf$pos+25)/50) # first is 1:50

## search windows with high RAAS
windows <- append(split(tmtf$RAAS, paste(tmtf$mane, tmtf$win1)),
                  split(tmtf$RAAS, paste(tmtf$mane, tmtf$win2)))
hist(lengths(windows), breaks=100)

## RAAS statistics for windows
winstat <- listProfile(windows, y=tmtf$RAAS, min=3)
winstat$mane <- sub("\\..*","",rownames(winstat))
ids <- pnms[winstat$mane]
names(ids) <- rownames(winstat)

plotdev(file.path(pfig.path,paste0("windows_volcano_all")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
volcano(winstat, cut=50, v.txt=c(-2,0), p.txt=10,
        ids=ids,
        xlab=xl.all)
mtext("50 AA windows, median of all RAAS", 3,0)
dev.off()

## protein window stat based on site-specific RAAS medians

site$win1 <- round(site$pos/50) # first is 1:25
site$win2 <- round((site$pos+25)/50) # first is 1:50

## search windows with high RAAS
windows <- append(split(site$RAAS.median,
                        paste0(site$mane,"-w1-",site$win1)),
                  split(site$RAAS.median,
                        paste0(site$mane,"-w2-",site$win2)))
hist(lengths(windows), breaks=100)

winstat <- listProfile(windows, y=tmtf$RAAS, min=2)
winstat$mane <- sub("-w.*","",rownames(winstat))

ids <- strsplit(rownames(winstat),"-")
ids <- unlist(lapply(ids, function(x) pnms[x[1]]))
names(ids) <- rownames(winstat)

plotdev(file.path(pfig.path,paste0("windows_volcano_sites")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
volcano(winstat, value="median", cut=20, v.txt=c(-2,-2),
        p.txt=6, ids=unlist(ids),
        xlab=xl.site)
        
mtext("50 AA windows, median of site medians", 3,0)
show <- (winstat$median>0 & winstat$p<0.001) | winstat$median>1.6
shadowtext(winstat$median[show],
           -log10(winstat$p[show]),
           labels=ids[show],
           pos=4, xpd=TRUE, col=1, cex=.8)
show <- pnms[winstat$mane] %in% c("PGM1","PSMA1")
points(winstat$median[show],
       -log10(winstat$p[show]), col=3)
shadowtext(winstat$median[show],
           -log10(winstat$p[show]),
           labels=ids[show],
           pos=4, xpd=TRUE, col=3, font=2)
dev.off()

## PROTEIN COMPLEXES

## TODO: first boil down RAAS to median per site to avoid
## over-estimation by single sites, eg. Q->G in PSMA1.

## rename names by ensembl in position-wise RAAS medians

## distribution in protein complexes

## hu.MAP v2
if ( FALSE ) {
    humap <- read.csv(humap.file)
    hulst <- strsplit(humap[,3]," ")
    huids <- humap[,1]
}

## CORUM v4.1
humap <- read.delim(corum.file)
hulst <- strsplit(humap$subunits.UniProt.IDs.,";")
huids <- sub(" complex$","",humap[,2])

names(hulst)<- huids

## replace complex uniprot genes by ALL ensembl proteins that map to it
huens <- lapply(hulst, function(x) {
    unique(uni2ens[uni2ens[,1]%in%x,2])
})


## COMPLEXES - MEDIAN SITE RAAS
## replace complex ensembl IDs with rows in RAAS table
huraas <- lapply(huens, function(x) c(na.omit(pbstat[x,"median"])))

## RAAS statistics
hustat <- listProfile(huraas, y=tmtf$RAAS, min=3)

plotdev(file.path(pfig.path,paste0("complex_volcano_sites")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
res <- volcano(hustat, value="median",
               p.txt=3, v.txt=c(-Inf,-1.5), cut=50, mid=-1,
               xlab=xl.site)
mtext("protein complexes, median of site medians", 3,0)
dev.off()


## COMPLEXES - ALL RAAS
## replace complex ensembl IDs with all RAAS in TMT level table
huraas <- lapply(huens, function(x) tmtf$RAAS[tmtf$mane%in%x])

## RAAS statistics
hustat <- listProfile(huraas, y=tmtf$RAAS, min=3)

plotdev(file.path(pfig.path,paste0("complex_volcano_all")),
        type=ftyp, res=300, width=4.4,height=4)
par(mai=c(.5,.5,.25,.5), mgp=c(1.3,.25,0), tcl=-.25)
res <- volcano(hustat, value="median",
               p.txt=20, v.txt=c(-Inf,-1), cut=50, mid=0,
               xlab=xl.all)
mtext("protein complexes, median of all RAAS", 3,0)
dev.off()

## NOTE: EZR, MSN, 

## investigate complexes
if ( FALSE ) {
    humap[grep("Spliceosome, E", humap$ComplexName),c("ComplexName",
                                                      "Complex.comment",
                                                      "subunits.Gene.name.")]
    humap[grep("Drosha", humap$ComplexName),c("ComplexName",
                                              "Complex.comment",
                                              "subunits.Gene.name.")]
}

### WRITE CHIMERAX ATTRIBUTE FILES

## LOAD PDB2ENSEMBL
## and write RAAS color display for selected/all structures
pdb2ens <- read.csv("~/data/mistrans/originalData/pdb_chain_ensembl.csv",
                    skip=1)

## filter all structures where we have AAS
pdb2ens <- pdb2ens[pdb2ens$TRANSLATION_ID%in%site$mane,]


pdb.path <- file.path(out.path,"pdb_attributes")
dir.create(pdb.path)

pdbid <- "6kwy"

pdb <- pdb2ens[pdb2ens[,1]==pdbid,]

show.str <- "show "

## TODO: write RAAS coloring string

cat(paste0("attribute: raas\nrecipient: residues\n\n"),
    file=file.path(pdb.path, paste0(pdbid, "_raas.defattr")))
for ( chain in unique(pdb$CHAIN) ) {
    ens <- unique(pdb$TRANSLATION_ID[pdb$CHAIN==chain])
    sts <- site[site$mane==ens,]
    for ( i in 1:nrow(sts) ) {
        sid <-
            paste0("/",chain, ":",sts$pos[i])
        show.str <- paste0(show.str, " ", sid)
        cat(paste0("\t", sid,"\t",sts$RAAS.median[i],"\n"),
            file=file.path(pdb.path, paste0(pdbid, "_raas.defattr")),
            append=TRUE)
    }
}

## commandline string to copy paste into chimeraX
cat(paste0("open ",pdbid,"; set bgColor black; color all #E6E6FA; hide all atoms; show all cartoons;",show.str,"; open /home/raim/data/mistrans/processedData/pdb_attributes/",pdbid,"_raas.defattr; color byattribute raas palette ^RdYlBu\n\n"), file=file.path(pdb.path, paste0(pdbid, "_raas.txt")))

## TODO: can we set RAAS colors palette in chimeraX?
paste(arno(5),collapse=";")

### CLUSTER BY REGIONS
