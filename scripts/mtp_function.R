
library(readxl)
library(segmenTools)
options(stringsAsFactors=FALSE)
library(gprofiler2)

proj.path <- "/home/raim/data/mistrans"
dat.path <- file.path(proj.path,"originalData")
out.path <- file.path(proj.path,"processedData")
fig.path <- file.path(proj.path,"figures")

##dat <- read_xlsx("All_MTP_BP_sharedPep_quant.xlsx")
##colnames(dat) <- gsub(" ","_", colnames(dat))

##dat <- read_xlsx(file.path(dat.path,"All_MTP_BP_sharedPep_quant_03Oct23.xlsx"))
##colnames(dat) <- gsub(" ","_", colnames(dat))

## READ MTP TABLE
## output from map_peptides.R
dat <- read.delim(file.path(out.path,"mtp_mapped.tsv"))


### REDUCE TABLE TO UNIQUE HUMAN PROTEINS

mids <- dat$proteinID
mids.sze <- table(mids) # sort(table(mids), decreasing=TRUE)

## NOTE: many mutations in expected, eg. ENSP00000474524
## "V region of the variable domain of immunoglobulin heavy chains
## that participates in the antigen recognition"

## analyze RAAS distribution of multiply mutated/MTPd
raas <- as.numeric(dat$RAAS)
raas[raas==-Inf] <- min(raas[is.finite(raas)]) -1
raas[raas== Inf] <- max(raas[is.finite(raas)]) +1

hist(raas)

raas.lst <- split(raas, f=mids)
raas <- cbind(##n=mids.sze,
              n=unlist(lapply(raas.lst, length)),#sum(raas$n!=raas$n2)==0
              mean=unlist(lapply(raas.lst, mean)),
              sd=unlist(lapply(raas.lst, sd)))
raas <- as.data.frame(raas)

head(raas$n)

### UNIQUE PROTEIN MATRIX

## ADD GENE DESCRIPTIONS
desc <- sub("\\[Homo sapiens\\]","",sub(".*\\|","",dat$Blast_ref_protein))
desc.lst <- split(desc, f=mids)
desc.lst <- lapply(desc.lst, function(x) trimws(x[x!=""]))
desc.lst <- lapply(desc.lst, unique)

## check unique? NOTE: several with non-unique ensembl<->refseq mapping
## i.e. Blast_ref_protein vs. Top_leading_protein
table(unlist(lapply(desc.lst, length)))
which(unlist(lapply(desc.lst, length))==2)

if ( any(as.numeric(names(table(unlist(lapply(desc.lst,length)))))>1) )
    warning("non-unique mapping between ensembl and refseq")
desc <- unlist(lapply(desc.lst, function(x) x[1]))

## TODO: collect more info here,
## * protein length

## ADD PROTEIN LENGTH
len <- dat$len
len.lst <- split(len, f=mids)
len.lst <- lapply(len.lst, function(x) trimws(x[x!=""]))
len.lst <- lapply(len.lst, unique)

## check unique? NOTE: several with non-unique ensembl<->refseq mapping
## i.e. Blast_ref_protein vs. Top_leading_protein
table(unlist(lapply(len.lst, length)))
which(unlist(lapply(len.lst, length))==2)

if ( any(as.numeric(names(table(unlist(lapply(len.lst,length)))))>1) )
    warning("non-unique mapping between ensembl and refseq")

len <- as.numeric(unlist(len.lst))

## MTP density
dns <- raas$n/len

raas <- cbind.data.frame(raas, description=desc, length=len, density=dns)

## sort by MTP/protein
raas <- raas[order(unlist(raas$density), decreasing=TRUE),]

## plot MTP/protein distribution
rn <- unlist(raas$n)
rn[rn>10] <- 11
## counts
rcnt <- table(rn)
## percent of total MTPs
rprc <- rcnt/nrow(dat)
## percent of total proteins
pprc <- rcnt/nrow(raas)

png(file.path(fig.path,"mtp_per_protein.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(rcnt, xlab="unique MTP/protein",
        ylab="# proteins")
legend("topright",c(paste0("unique MTP: ", nrow(dat)),
                    paste0("unique proteins: ", nrow(raas))))
dev.off()

png(file.path(fig.path,"mtp_per_protein_percent_mtp.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(rprc, xlab="unique MTP/protein",
        ylab="% of total unique MTP")
legend("topright",c(paste0("unique MTP: ", nrow(dat)),
                    paste0("unique proteins: ", nrow(raas))))
dev.off()
png(file.path(fig.path,"mtp_per_protein_percent_protein.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(pprc, xlab="unique MTP/protein",
        ylab="% of total unique proteins")
legend("topright",c(paste0("unique MTP: ", nrow(dat)),
                    paste0("unique proteins: ", nrow(raas))))
dev.off()

## plot RAAS per protein variance
## TODO: understand RAAS per protein distribution or get good value
## test whether correlating
plotdev(file.path(fig.path,"mtp_raas_mean_sd"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(raas$mean, raas$sd, xlab="mean RAAS", ylab="SD RAAS")
dev.off()

mx <- max(abs(c(raas$sd,raas$mean)), na.rm=TRUE)
brks <- seq(-mx-1, mx+1, 1)
hist(raas$mean, breaks=brks)
hist(raas$mean[raas$n>1], breaks=brks, col=2, add=TRUE)
hist(raas$sd[raas$n>1], xlab="standard deviation of RAAS", breaks=brks,
     border=4, add=TRUE)

hist(log10(raas$density), breaks=100,
     xlab=expression(log[10](MTP/AA)))

## RAAS vs. MTP DENSITY
plotdev(file.path(fig.path,"mtp_raas_density"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(raas$mean, log10(raas$density), xlab="mean RAAS",
        ylab=expression(log[10](MTP/AA)))
dev.off()

plotdev(file.path(fig.path,"mtp_raas_number"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(raas$mean, log10(raas$n), xlab="mean RAAS",
        ylab=expression(log[10](MTP/protein)))
dev.off()

plotdev(file.path(fig.path,"mtp_length"),
        height=3.5, width=3.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
dense2d(raas$len, raas$n, log="xy", xlab="protein length/AA",
        ylab="MTP/protein")
dev.off()


## ENRICHMENT TESTS
go.path <- file.path(fig.path,"annotation")
dir.create(go.path, showWarnings=FALSE)

## categorizations

rn <- raas$n
rn[rn>3] <- 4
cls.mat <- cbind(number=rn)
                
rs <- raas$mean
rs[abs(rs)>10] <- 10 * sign(rs[abs(rs)>10])
brks <- seq(-10,10,4)
bins <- cut(raas$mean, breaks=brks)
bn <- as.character(bins)
bn[is.na(bn)] <- "na"

cls.mat <- cbind(number=rn,
                 raas=bn)
cls.srts <- list(number=1:4,
                 raas=levels(bins))
cls.labs <- c(number="MTP/protein",
              raas="mean RAAS")
rownames(cls.mat) <- rownames(raas)

ovl <- clusterCluster(cls.mat[,1], cls.mat[,2])
## calculate optimal figure height: result fields + figure margins (mai)
nh <- nrow(ovl$p.value) *.3 + 1
nw <- ncol(ovl$p.value) *.4 + 1

## RAAS vs. MTP densitiy
plotdev(file.path(fig.path,"mtp_raas_density_discrete"),
        height=nh, width=nw, res=200)
par(mai=c(.75,.5,.4,.4), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovl, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=cls.labs[1])
mtext(cls.labs[2], 1, 2.5)
dev.off()
       
## TODO: look for manually defined functions over all annotations
my.functions <-
    list(immune=c("antigen binding","immunoglobulin","immune response"),
         extracellular=c("external", "extracellular",
                         "vesicle", "junction", "adhesion", "secretory"),
         cytoskeleton=c("actin[^g]","myosin","tubul","adhesion","fiber",
                        "polymer","filament"),
         metabolism=c("glycoly", "gluconeogene", "proteolysis", "catabolism",
                      "vacuol", "lysos", "proteaso", "translation", "ribosom"),
         organelle=c("mitochon", "organell","vacuol", "lysos"))

for ( ct in 1:ncol(cls.mat) ) {


    ## get clustering
    cid <- colnames(cls.mat)[ct]
    cls <- cls.mat[,ct]

    cl.srt <- cls.srts[[ct]]

    cat(paste("CALCULATING ANNOTATION ENRICHMENTS for", cid, "\n"))
    

    ## ENRICHMENT OVER ALL CATEGORIES in gprofiler2
    ovll <- runGost(cls, organism="hsapiens", cls.srt=cl.srt) 
    ## plot enrichments
    for ( ctgy in names(ovll) ) {
    
        ovl <- ovll[[ctgy]]

        ## filter large annotations such as GO
        p.filt <- 1e-1
        cut <- TRUE
        if ( length(grep("^GO:",ctgy))>0 ) {

            ## 
            trm.srt <- ovl$num.query[,1]
            keep <- trm.srt > 10 & trm.srt <5e3
            trm.srt <- names(trm.srt)[keep]
            ovl <- sortOverlaps(ovl, srt=trm.srt)
            
            p.filt <- 1e-3
            cut <- TRUE
        } else if ( ctgy=="TF" ) {
            p.filt <- 1e-3
            cut <- TRUE
        }

        
        ## sort along clusters
        ovc <- sortOverlaps(ovl, p.min=p.filt, cut=cut)
        
        ## calculate optimal figure height: result fields + figure margins (mai)
        nh <- nrow(ovc$p.value) *.2 + 1
        nw <- ncol(ovc$p.value) *.4 + 5
        
        ## plot figure
        plotdev(file.path(go.path,paste0(cid,"_annot_",ctgy)),
                height=nh, width=nw, res=200)
        par(mai=c(.75,4.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        plotOverlaps(ovc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
                     xlab=NA, ylab=NA)
        mtext(cls.labs[ct], 1, 1, adj=-1.5, cex=1.2, font=2)
        axis(2, at=nrow(ovc$p.value)+1.5, label=ctgy,
             cex.axis=1.2, font=2, xpd=TRUE, las=2)
        dev.off()
    }
    
    ## ENRICHMENT OVER PRE-SELECTED TERMS in gprofiler2
    ovll <- runGost(cls, organism="hsapiens", cls.srt=cl.srt,
                    terms=my.functions) 
    ## plot enrichments
    for ( ctgy in names(ovll) ) {
    
        ovl <- ovll[[ctgy]]

        ## filter large annotations such as GO
        p.filt <- 1e-2
        cut <- TRUE
        
        ## sort along clusters
        ovc <- sortOverlaps(ovl, p.min=p.filt, cut=cut)
        
        ## calculate optimal figure height: result fields + figure margins (mai)
        nh <- nrow(ovc$p.value) *.2 + 1
        nw <- ncol(ovc$p.value) *.4 + 5
        
        ## plot figure
        plotdev(file.path(go.path,paste0(cid,"_annot_",ctgy)),
                height=nh, width=nw, res=200)
        par(mai=c(.75,4.5,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
        plotOverlaps(ovc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
                     xlab=NA, ylab=NA)
        mtext(cls.labs[ct], 1, 1, adj=-1.5, cex=1.2, font=2)
        axis(2, at=nrow(ovc$p.value)+1.5, label=ctgy,
             cex.axis=1.2, font=2, xpd=TRUE, las=2)
        dev.off()
    }
}
