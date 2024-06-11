
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

hfig.path <- file.path(fig.path,"hotspots")
dir.create(hfig.path, showWarnings=FALSE)

corW <- corH <- 2.5

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

## 

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

plotdev(file.path(hfig.path,paste0("windows_volcano_all")),
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

plotdev(file.path(hfig.path,paste0("windows_volcano_sites")),
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

plotdev(file.path(hfig.path,paste0("complex_volcano_sites")),
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

plotdev(file.path(hfig.path,paste0("complex_volcano_all")),
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

### AAS along proteins
head(tmtf$unique.site)

tmtl <- split(tmtf$RAAS, tmtf$unique.site)
tmtl <- listProfile(tmtl, y=tmtf$RAAS, use.test=use.test, min=3)
tmts <- tmtf[!duplicated(tmtf$unique.site),]

## TODO: calculate a moving average of RAAS in sliding window with 1 bp
## 


### 3D STRUCTURE - WRITE CHIMERAX ATTRIBUTE FILES

## LOAD PDB2ENSEMBL
## and write RAAS color display for selected/all structures
pdb2ens <- read.csv("~/data/mistrans/originalData/pdb_chain_ensembl.csv",
                    skip=1)

## filter all structures where we have AAS
pdb2ens <- pdb2ens[pdb2ens$TRANSLATION_ID%in%site$mane,]


pdb.path <- file.path(out.path,"pdb_attributes")
dir.create(pdb.path)

pdbid <- "6kwy" # proteasome cryoEM
pdbid <- "6qzp" # ribosome cryoEM
pdbid <- "1xk4" # EF hand example, S100A8,  AAS hotspot in heterodimer
                # with S100A9

### loop through ALL pdbs and write chimeraX command strings
### for PDBs containing proteins with AA substitutions.
for ( pdbid in unique(pdb2ens[,1]) ) {

    pdb <- pdb2ens[pdb2ens[,1]==pdbid,]

    show.str <- ""
    col.str <- ""
    sel.str <- ""
    
    defattr.file <- file.path(pdb.path, paste0(pdbid, "_raas.defattr"))
    
    if ( FALSE ) cat(paste0("attribute: raas\nrecipient: residues\n\n"),
                     file=defattr.file)
    for ( chain in unique(pdb$CHAIN) ) {
        ens <- na.omit(unique(pdb$TRANSLATION_ID[pdb$CHAIN==chain]))
        ## get AAS for this pdb chain!
        sts <- site[site$mane==ens,]
        for ( i in 1:nrow(sts) ) {
            sid <-
                paste0("/",chain, ":",sts$pos[i])
            sel.str <- paste0(sel.str," ",sid)
            show.str <- paste(show.str, "show", sid, "atoms;")
            col.str <- paste(col.str, "color",sid, sts$RAAS.color[i],";")
            if ( FALSE) cat(paste0("\t", sid,"\t",sts$RAAS.median[i],"\n"),
                            file=defattr.file,
                            append=TRUE)
        }
    }
    
    ## commandline string to copy paste into chimeraX
    if ( FALSE) cat(paste0("open ",pdbid,"; set bgColor black; color all #E6E6FA; hide all atoms; show all cartoons;",show.str,"; open /home/raim/data/mistrans/processedData/pdb_attributes/",pdbid,"_raas.defattr; color byattribute raas palette ^RdYlBu\n\n"), file=file.path(pdb.path, paste0(pdbid, "_raas_attribute.txt")))
    
    
    ## TODO: can we set RAAS colors palette in chimeraX?
    cat(paste0("open ",pdbid,
               "; set bgColor black",
               "; color all #E6E6FA; hide all atoms; show all cartoons;",
               col.str,
               show.str,
               "; select ",sel.str,"\n\n"),
        file=file.path(pdb.path, paste0(pdbid, "_raas.txt")))
}



if ( !interactive() ) quit("no")

###
library(iPAC)
library(SpacePAC)


## PSMA1

uid <- "P25786"
pdbid <- "6kwy"
chain <- "E"

## RAP1A
uid <- "P62834"
pdbid <- "1gua" 
chain <- "A"

pdbid <- "1xk4" # EF hand example, S100A8,  AAS hotspot in heterodimer
                # with S100A9

## TODO: use bio3d package instead of iPAC,
## TODO: use dssp to extract secondary structures from PDB,

## get ensembl IDs for this structure
pdb <- pdb2ens[pdb2ens[,1]==pdbid,]
eids <- unique(pdb$TRANSLATION_ID)

## get structure and fasta
## TODO: loop over all fasta that map to a given structure
## and fuse results?
CIF <- paste0("https://files.rcsb.org/view/",pdbid,".cif")
Fasta <- paste0("https://www.uniprot.org/uniprot/", uid, ".fasta")
Positions <- iPAC::get.AlignedPositions(CIF, Fasta, chain)

## AAS 
usites <- site[which(site$uniprot==uid),]

## TODO: add RAAS to position
raas3d <- Positions$Position
colnames(raas3d) <- c("AA","pos","chain","x","y","z")

## get RAAS
raas3d$RAAS <- rep(NA, nrow(raas3d))
raas3d$RAAS <- median(sites$RAAS.median)
midx <- match(usites$pos, raas3d$Can.Count)
raas3d$RAAS[midx[!is.na(midx)]] <- usites$RAAS.median[!is.na(midx)]

dst <- as.matrix(dist(raas3d[,c("x","y","z")],
                      diag = TRUE, upper = FALSE))

## NOTE: from SpacePAC
source("~/work/mistrans/scripts/saap_utils.R")
plotRaas3d(raas3d, center=usites$pos,
           color=usites$RAAS.color, radius=usites$n-1,
           name=paste(unique(usites$name),collapse=";"))


## generate mutation matrix
seq <- readFASTA(Fasta)
Mutations <- matrix(0,
                    ncol=nchar(seq[[1]]$seq),
                    nrow=nrow(usites))
for ( i in 1:(nrow(usites)) )
    Mutations[i,usites$pos[i]] <- 1
## OBLIGATORY COLUMN NAMES
## strange errors without
colnames(Mutations) <- paste0("V",1:ncol(Mutations))


icl <- ClusterFind(mutation.data=Mutations, 
                   position.data=Positions$Positions,
                   create.map = "Y",Show.Graph = "Y")

## load local copy of spacepac to learn and debug
if ( FALSE ) {
    allrs <- paste0(file.path("~/programs/SpacePAC/R/",
                              list.files(path="~/programs/SpacePAC/R/",
                                         pattern="*.R")))
    for ( file in allrs ) source(file)

    ## residue distance matrix
    dst <- create.Distance.Matrix(Positions$Positions)
}


mycs <- SpaceClust(Mutations, Positions$Positions,
                   numsims =1000, simMaxSpheres = 3,
                   radii.vector = 1:7, method = "SimMax")

mycs$optimal.sphere[,8]

mycp <- SpaceClust(Mutations, Positions$Positions,
                   radii.vector = c(1:10), alpha = .99,
                   method = "Poisson", multcomp="none")
mycp$result.poisson[,"P.Value"]


### spacepac example

PIK3CA.CIF <- "https://files.rcsb.org/view/2ENQ.cif"
PIK3CA.Fasta <- "https://www.uniprot.org/uniprot/P42336.fasta"
PIK3CA.Positions <- get.AlignedPositions(PIK3CA.CIF, PIK3CA.Fasta, "A")

## MUTATION MATRIX: COLUMNS are residues in FASTA,
## rows contain 1 for individual mutations

data(KRAS.Mutations) 
data(PIK3CA.Mutations) 
my.clusters <- SpaceClust(PIK3CA.Mutations, PIK3CA.Positions$Positions,
                          numsims =1000, simMaxSpheres = 3,
                          radii.vector = c(1,2,3,4), method = "SimMax")
my.clusters <- SpaceClust(PIK3CA.Mutations, PIK3CA.Positions$Positions,
                          radii.vector = c(1,2,3,4),
                          alpha = .05, method = "Poisson")
### CLUSTER BY REGIONS
