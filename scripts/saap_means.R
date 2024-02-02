
### TODO:
## * mean of mean of means: better use median?
## * 

library(readxl)
library(segmenTools)
options(stringsAsFactors=FALSE)

mam.path <- "/home/raim/data/mammary"
proj.path <- "/home/raim/data/mistrans"

gen.path <- file.path(mam.path, "originalData")
dat.path <- file.path(proj.path,"originalData")
fig.path <- file.path(proj.path,"figures","saap_means")
out.path <- file.path(proj.path,"processedData")


## figures
dir.create(fig.path, showWarnings=FALSE)

tmt.file <- file.path(dat.path, "All_filtered_SAAP_TMTlevel_quant_df.xlsx")
pat.file <- file.path(dat.path, "All_filtered_SAAP_patient_level_quant_df.xlsx")


## ANALYZING MEAN TMT RAAS AT DIFFERENT LEVELS

tmtf <- read_xlsx(tmt.file)
tmtf$RAAS <- as.numeric(tmtf$RAAS)
## exclude NA or Inf
rm <- is.na(tmtf$RAAS) | is.infinite(tmtf$RAAS)
cat(paste("removing", sum(rm), "with NA or Inf RAAS from TMT level\n"))
tmtf <- tmtf[!rm,]

summary.types <- cbind(s4=tmtf$SAAP,
                       s3=paste0(tmtf$SAAP,"_",tmtf$BP),
                       s2=paste0(tmtf$Dataset,"_",tmtf$SAAP),
                       s1=paste0(tmtf$Dataset,"_",tmtf$SAAP,"_",tmtf$BP))

for ( i in 1:ncol(summary.types) ) {

    id <- colnames(summary.types)[i]
    CL <- summary.types[,i]

    ## split by summary type
    tmt <- split(tmtf$RAAS, CL)

    ## delog for CV
    tmtl <- lapply(tmt, function(x) x^10) # delog
    tlen <- unlist(lapply(tmtl, length))
    tmds <- unlist(lapply(tmtl, median, na.rm=FALSE))
    tmns <- unlist(lapply(tmtl, mean, na.rm=FALSE))
    tdev <- unlist(lapply(tmtl, sd, na.rm=FALSE))
    tcvs <- tdev/tmns

    
    
    ## calculate mean RAAS from raw data for each of the SAAP
    ## don't delog first, to reproduce shiri's values
    tmnr <- unlist(lapply(tmt, mean, na.rm=FALSE))
    tmdr <- lapply(tmt, median, na.rm=TRUE)
    tmdr <- unlist(lapply(tmdr, function(x) ifelse(length(x)==0, NA, x)))
    tsdr <- unlist(lapply(tmt, sd, na.rm=TRUE))

    ## NOTE: CV for log-normal data, see @Canchola2017
    ## cv = sqrt(10^(log(100)*var) -1)
    cvl <- function(x) sqrt(10^(log(10)*var(x,na.rm=TRUE)) -1)
    tcvr <- unlist(lapply(tmt, cvl))
    tcvr[tcvr>10] <- 10 # remove huge outlier
     
    ## collect number of times this SAAP/BP was found per cancer type
    tmtf$count <- tlen[CL]


    ## CV on de-logged level,
    png(file.path(fig.path,paste0(id,"_raas_means_delogged_CV.png")),
         res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(tcvs, log10(tlen), axes=FALSE,
            xlab=expression(CV==sdev/mean), cex=.6,
            ylab="count of unique TMT rows")
    axis(1)
    axis(2, at=log10(c(1:10,1:10*10,1:10*100)), labels=NA, tcl=par("tcl")/2)
    axis(2, at=log10(c(1,10,50,100,200)), labels=c(1,10,50,100,200))
    dev.off()

    ## delogged means deviate strongly
    png(file.path(fig.path,paste0(id,"_raas_means_delogged.png")),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(log10(tmns), tmnr, xlab=expression(log[10](mean(10^RAAS))),
            ylab="mean RAAS from TMT level file")
    abline(a=0, b=1, col=1)
    dev.off()

    ## regression to mean problem
    png(file.path(fig.path,paste0(id,"_raas_means_datasets.png")),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    df <- data.frame(mean=tmnr, n=tlen)
    df <- df[df$n>0,]
    dense2d(df$mean, log10(df$n),
            ylab="count of unique TMT rows", cex=.6,
            xlab="mean TMT level RAAS", axes=FALSE, xlim=c(-6,4))
    axis(1)
    axis(2, at=log10(c(1:10,1:10*10,1:10*100)), labels=NA, tcl=par("tcl")/2)
    axis(2, at=log10(c(1,10,50,100,200)), labels=c(1,10,50,100,200))
    dev.off()
    
    ## regression to mean 
    png(file.path(fig.path,paste0(id,"_raas_means_datasets_all.png")),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(tmtf$RAAS, log10(tmtf$count), cex=.2,
            xlab="TMT level RAAS", ylab="count of unique TMT rows", axes=FALSE,
            xlim=c(-6,4))
    axis(1)
    axis(2, at=log10(c(1:10,1:10*10,1:10*100)), labels=NA, tcl=par("tcl")/2)
    axis(2, at=log10(c(1,10,50,100,200)), labels=c(1,10,50,100,200))
    dev.off()

    ## total counts at each level
    png(file.path(fig.path,paste0(id,"_raas_means_datasets_counts.png")),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    lgcnt <- table(log10(tmtf$count))
    plot(names(lgcnt), lgcnt, type="h", axes=FALSE,
         xlab="count of unique TMT rows", ylab=NA)
    corners = par("usr") 
    text(x = corners[1]-mean(corners[1:2])/4, y = mean(corners[3:4]),
         "total count", srt=-90, xpd=TRUE)
    axis(2)
    axis(1, at=log10(c(1:10,1:10*10,1:10*100)), labels=NA, tcl=par("tcl")/2)
    axis(1, at=log10(c(1,10,50,100,200)), labels=c(1,10,50,100,200))
    dev.off()
    
    
    ## SD
    png(file.path(fig.path,paste0(id,"_raas_means_SD.png")),
         res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(tsdr, log10(tlen), axes=FALSE,
            xlab=expression(standard~deviation), cex=.6,
            ylab="count of unique TMT rows")
    axis(1)
    axis(2, at=log10(c(1:10,1:10*10,1:10*100)), labels=NA, tcl=par("tcl")/2)
    axis(2, at=log10(c(1,10,50,100,200)), labels=c(1,10,50,100,200))
    dev.off()
    
    ## CV for log-normal data
    png(file.path(fig.path,paste0(id,"_raas_means_CV.png")),
         res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dense2d(tcvr, log10(tlen), axes=FALSE,
            xlab=NA,
            cex=.6, ylab="count of unique TMT rows")
    mtext(expression(CV==sqrt(10^(ln(10)*var(x)) -1)), 1, 1.5)
    axis(1)
    axis(2, at=log10(c(1:10,1:10*10,1:10*100)), labels=NA, tcl=par("tcl")/2)
    axis(2, at=log10(c(1,10,50,100,200)), labels=c(1,10,50,100,200))
    dev.off()

    ## regression to the MEDIAN 
    png(file.path(fig.path,paste0(id,"_raas_medians_datasets.png")),
        res=300, width=3.5, height=3.5, units="in")
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    df <- data.frame(mean=tmdr, n=tlen)
    df <- df[df$n>0,]
    dense2d(df$mean, log10(df$n),
            ylab="count of unique TMT rows", cex=.6,
            xlab="median TMT level RAAS", axes=FALSE, xlim=c(-6,4))
    axis(1)
    axis(2, at=log10(c(1:10,1:10*10,1:10*100)), labels=NA, tcl=par("tcl")/2)
    axis(2, at=log10(c(1,10,50,100,200)), labels=c(1,10,50,100,200))
    dev.off()
}
###
