
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

## add gene descriptions
desc <- sub("\\[Homo sapiens\\]","",sub(".*\\|","",dat$Blast_ref_protein))
desc.lst <- split(desc, f=mids)
desc.lst <- lapply(desc.lst, function(x) trimws(x[x!=""]))
desc.lst <- lapply(desc.lst, unique)

## check unique? NOTE: several with non-unique ensembl<->refseq mapping
## i.e. Blast_ref_protein vs. Top_leading_protein
table(unlist(lapply(desc.lst, length)))
which(unlist(lapply(desc.lst, length))==2)

if ( any(table(unlist(lapply(desc.lst,length)))>1) )
    warning("non-unique mapping between ensembl and refseq")
desc <- unlist(lapply(desc.lst, function(x) x[1]))

raas <- cbind.data.frame(raas, description=desc)

## sort by MTP/protein
raas <- raas[order(unlist(raas$n), decreasing=TRUE),]

## plot MTP/protein distribution
rn <- unlist(raas$n)
rn[rn>10] <- 11
## counts
rcnt <- table(rn)
## percent of total MTPs
rprc <- rcnt/nrow(dat)

png(file.path(fig.path,"mtp_per_protein.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(rcnt, xlab="unique MTPs/protein",
        ylab="# proteins")
legend("topright",c(paste0("unique MTP: ", nrow(dat)),
                    paste0("unique proteins: ", nrow(raas))))
dev.off()

png(file.path(fig.path,"mtp_per_protein_percent.png"),
    res=300, width=5, height=3.5, units="in")
par(mai=c(.5,.5,.15,.1), mgp=c(1.3,.3,0), tcl=-.25)
barplot(rprc, xlab="unique MTPs/protein",
        ylab="% of total unique MTPs")
legend("topright",c(paste0("unique MTP: ", nrow(dat)),
                    paste0("unique proteins: ", nrow(raas))))
dev.off()

## plot RAAS per protein variance
## TODO: understand RAAS per protein distribution or get good value
## test whether correlating
plot(raas$mean, raas$sd)

mx <- max(abs(c(raas$sd,raas$mean)), na.rm=TRUE)
brks <- seq(-mx-1, mx+1, 1)
hist(raas$mean, breaks=brks)
hist(raas$mean[raas$n>1], breaks=brks, col=2, add=TRUE)
hist(raas$sd[raas$n>1], xlab="standard deviation of RAAS", breaks=brks,
     border=4, add=TRUE)


## ENRICHMENT TESTS
go.path <- file.path(fig.path,"annotation")
dir.create(go.path)

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
cls.labs <- c(number="MTPs/protein",
              raas="mean RAAS")

ovl <- clusterCluster(cls.mat[,1], cls.mat[,2])
plotOverlaps(ovl, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=cls.labs[2], ylab=cls.labs[1])

for ( ct in 1:ncol(cls.mat) ) {

    ##cls <- get(paste0("cls",ctypes[ct]))
    ##cls.srt <- get(paste0("cls",ctypes[ct],".srt"))

    ## by MTP/protein bins
    cid <- colnames(cls.mat)[ct]
    cls <- cls.mat[,ct]

    cl.srt <- cls.srts[[ct]]
    cl.lst <- split(rownames(raas), f=as.character(cls))[cl.srt]
    cl.sze <- unlist(lapply(cl.lst, length))[cl.srt]
    
    gores <- gost(query=cl.lst, organism = "hsapiens", significant=FALSE)
    
    ## PLOT ENRICHMENTS
    for ( ctgy in rev(unique(gores$result$source)) ) { #c("CC","BP","MF") ) {
        
        go <- gores$result
        go <- go[go$source==ctgy,]

        cat(paste("analyzing category", ctgy, "with",
                  length(unique(go$term_id)), "entries\n"))
    
    
        terms <- go[,c("term_id","term_name","term_size")]
        terms <- terms[!duplicated(terms$term_id),]
        rownames(terms) <- terms$term_id
        
        ## filter large annotations such as GO
        p.filt <- 1e-1
        cut <- TRUE
        if ( length(grep("^GO:",ctgy))>0 ) {
            terms <- terms[terms$term_size >10 & terms$term_size<5e3,]
            p.filt <- 1e-3
            cut <- TRUE
        } else if ( ctgy=="TF" ) {
            p.filt <- 1e-3
            cut <- TRUE
        } 
        
        ## construct overlap matrix
        govl <- matrix(NA, nrow=nrow(terms), ncol=length(cl.lst))
        colnames(govl) <- names(cl.lst)
        rownames(govl) <- terms$term_name
        
        gopvl <- gocnt <- govl
        
        ## generate overlap structure:
        ## fill matrices with overlap p-values and counts, num.query, num.target
        ## TODO: do this more efficiently, one loop and vectors.
        for ( i in 1:nrow(govl) ) 
            for ( j in 1:ncol(govl) ) {
                idx <- which(go$query==colnames(govl)[j] &
                             go$term_id==terms$term_id[i])
                if ( length(idx)>1 ) {
                    stop("too many fields; shouldn't happen")
                    break
                } else if ( length(idx)==0 ) {
                    gopvl[i,j] <- 1
                    gocnt[i,j] <- 0
                } else {
                    gopvl[i,j] <- go$p_value[idx]
                    gocnt[i,j] <- go$intersection_size[idx]
                }
            }
        
        ## construct overlap class
        ovl <- list()
        ovl$p.value <- gopvl
        ovl$count <- gocnt
        ovl$num.query <- as.matrix(terms[,"term_size",drop=FALSE])
        ovl$num.target <- t(as.matrix(cl.sze))
        class(ovl) <- "clusterOverlaps"
        
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
        dev.off()
    }
}
