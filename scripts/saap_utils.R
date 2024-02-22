
# from->to as expression, e.g. for axis labels
ftlabels <- function(srt) {
    axex <- rep("",length(srt))
    names(axex) <- srt
    for ( i in seq_along(srt) ) {
        ft <- unlist(strsplit(srt[i],":"))
        axex[i] <- as.expression(bquote(.(ft[1]) %->% .(ft[2])))
    }
    axex
}

## wilcox test with normalized U-statistic, to
## be used for statistical profile plots
w.test <- function(x,y) {
    res <- wilcox.test(x,y)
    ## normalized U-statistic
    tt <- res$statistic/(sum(!is.na(x))*sum(!is.na(y))) -0.5
    rt <- list()
    rt$statistic <- unlist(tt)
    rt$p.value <- unlist(res$p.value)
    rt
}


## plot p.values, RAAS, and count.
plotovl <- function(ovl, text="count", cut=TRUE, value="median",
                    txt.cut=-2, ...) {

    ## TODO: auto-select txt.cut, in image_matrix based on brightness,
    ## as in 
    txt <- ovl[[text]]
    txt[txt==0] <- ""
    txt.col <- ifelse(ovl[[value]]< txt.cut, "white","black")
    image_matrix(ovl[[value]], cut=cut,
                 text=txt, text.col=txt.col, ...)
}

## volcano plot function for overlap class
## produced by raasProfile
volcano <- function(ovl, cut=15, value="mean", p.txt=6, v.txt, mid,
                    density=TRUE, ...) {

    
    tmv <- gdata::unmatrix(t(ovl[[value]]))
    ##tmv <- c(ovl[[value]])
    tpv <- -log10(gdata::unmatrix(t(ovl$p.value)))
    tpv[tpv>cut] <- cut

    na <- is.na(tmv) | is.na(tpv)
    tpv <- tpv[!na]
    tmv <- tmv[!na]

    if ( density ) {
        dense2d(tmv, tpv, ylab=expression(-log[10](p)), ...)
    } else {
        plot(tmv, tpv, ylab=expression(-log[10](p)), ...)
    }
    axis(4)

    ## indicate high
    show <- tpv>p.txt
    if ( !missing(v.txt) )
        show <- show & (tmv<v.txt[1] | tmv>v.txt[2])

    show <- which(show)
    if ( length(show) ){
        if ( missing(mid) )
            mid <- mean(ovl[[value]],na.rm=TRUE)
        shadowtext(tmv[show], tpv[show], labels=names(tpv)[show],
                   pos=ifelse(tmv[show]<mid, 2,4), col="red",
                   font=2, cex=.8, xpd=TRUE)
    }
}

## calculate statistical profiles, similar to segmenTools::clusterProfiler,
## but working on lists of unequal lengths instead of a matrix
raasProfile.old <- function(x=cdat, id="SAAP", values=tmt,
                        bg=FALSE, vid,
                        rows="to", cols="aacodon",
                        row.srt, col.srt,
                        use.test=use.test,
                        do.plots=interactive(),
                        xlab="value", fname="profile_",
                        verb=FALSE) {

    ## TODO: use full genetic code for consistent columsn and rows!

    if ( !missing(row.srt) ) aas <- row.srt
    else aas <- sort(unique(x[,rows]))
    if ( !missing(col.srt) ) acod <- col.srt
    else acod <- sort(unique(x[,cols]))

    ## FILTER SRTS
    aas <- aas[aas%in%x[,rows]]
    acod <- acod[acod%in%x[,cols]]
    
    codt <- list()
    
    tt <- matrix(0, ncol=length(acod), nrow=length(aas))
    colnames(tt) <- acod
    rownames(tt) <- aas
    tc <- tp <- tm <- td <- tt
    tp[] <- 1
    tm[] <- td[] <- NA


    min.obs=2
    allvals <- list()
    for ( i in seq_along(aas) ) {

        aa <- aas[i]
        allvals[[aa]] <- list()

        for ( j in seq_along(acod) ) {

            cod <- acod[j]
        
            ## get all SAAP for this codon
            idx <- x[,cols]==cod & x[,rows]==aa

            ## NOTE: unique SAAP!!
            csap <- x[idx,id]
            csap <- unique(csap)

            if ( length(idx)==0 ) next
            if ( sum(idx,na.rm=TRUE)==0 ) next
            
            if (verb>0) cat(paste("testing",cod,"to", aa, "\n"))

            ## if a local background is chosen, tmt
            ## should be a table with the bg column
            ## and we filtering TMT only for this subset
            if ( bg ) {

                ## TODO: CHECK THIS: split by SAAP, get all RAAS
                ## duplicate SAAP effect?
                if ( verb>0 )
                    cat(paste("getting subset of values:",
                              cols, cod, vid, id,"\n"))
                vals <- values[values[,cols]==cod,] ## FILTER DATASET
                vals <- split(vals[,vid], vals[,id])## split RAAS/vid by SAAP/id
            } else vals <- values
                
            
            ## get all RAAS values
            rm <- !csap%in%names(vals)
            if ( sum(rm) ) 
                cat(paste("WARNING:", aa, cod, ":", sum(rm), "of", length(csap),
                          "IDs not found in values:",
                          paste(csap[rm],collapse=","),"\n"))
            csap <- csap[!rm]
            
            ## get all RAAS values
            y <- unlist(vals[csap])
            X <- unlist(vals)

            allvals[[aa]][[cod]] <- y
            
            tc[i,j] <- length(y)

            if ( length(y)==0 ) next
            
            tm[i,j] <- mean(y)
            td[i,j] <- median(y)
           
            if ( sum(!is.na(y)) >= min.obs ) {
                tts <- use.test(y, X)
                codt[[cod]] <- tts
                tt[i,j] <- tts$statistic
                tp[i,j] <- tts$p.value
                
                if ( do.plots ) {
                    outfile <- paste0(fname,cod,"_",aa,".png")
                    hcol <- ifelse(tts$statistic<0, "#0000ff","#ff0000")
                    if ( verb>0 ) cat(paste("plotting", outfile, "\n"))
                    png(outfile, width=3, height=2,
                        units="in", res=100)
                    par(mai=c(.5,.15,.25,.15), tcl=-.25, mgp=c(1.3,.3,0))
                    brks <- seq(min(X), max(X), length.out=20)
                    hist(y, freq=FALSE, border=NA,
                         col=paste0(hcol,77),
                         xlim=range(X), main=NA,
                         xlab=xlab, breaks=brks, axes=FALSE)
                    par(new=TRUE)
                    hist(X, freq=FALSE, col="#77777777", border=NA,
                         xlim=range(X), main=NA, xlab=NA, ylab=NA, axes=FALSE,
                         breaks=brks)
                    abline(v=median(X), col=1, lwd=2)
                    abline(v=td[i,j], col=hcol, lwd=2)
                    axis(1)
                    mtext(as.expression(bquote(.(cod) %->% .(aa))),
                          3,-.1, font=2, cex=1.5)
                    legend(ifelse(tts$statistic<0, "topright","topleft"),
                           c(paste0("n=", length(y)),
                             paste0("u/t=", round(tts$statistic,1)),
                             paste0("p=", signif(tts$p.value,1))),
                           bty="n", seg.len=0, x.intersp=0)
                    dev.off()
                }
            }
        }
    }
    if ( verb>0 ) cat(paste("DONE\n"))
    
    ## construct overlap object
    ova <- list()
    ova$p.value <- tp
    ova$statistic <- tt
    ova$count <- tc
    ova$mean <- tm
    ova$median <- td

    ## unique SAAP from main data table
    ova$unique <- table(x[,rows], x[,cols])[rownames(tp),colnames(tp)]
    
    sg <- sign(tt)
    sg[is.na(sg)] <- 1
    sg[sg==0] <- 1
    ova$sign <- sg
    
    ##ova$statistic[is.na(ova$statistic)] <- ""
    
    ## TODO: add counts
    ova$num.target <-
        t(as.matrix(table(x[,cols])[colnames(tp)]))
    ova$num.query <-
        as.matrix(table(x[,rows]))[rownames(tp),,drop=FALSE]


    ## attach all raw values
    ova$ALL <- allvals
    
    class(ova) <- "clusterOverlaps"
    ova
}

##ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
##                    delog=TRUE, rows="pfromto", cols="Dataset",
##                    bg=TRUE, row.srt=srt, col.srt=uds,
##                    use.test=t.test, do.plots=TRUE, xlab="TMT level RAAS",
##                    fname=fname, verb=1)

## calculate statistical profiles, similar to segmenTools::clusterProfiler,
## but working on lists of unequal lengths instead of a matrix
raasProfile <- function(x=tmtf, id="SAAP", 
                        value="RAAS", delog=TRUE, replace=TRUE,
                        bg=FALSE, na.rm=FALSE,
                        rows="to", cols="aacodon",
                        row.srt, col.srt,
                        use.test=use.test,
                        do.plots=interactive(),
                        xlab="value", fname="profile_",
                        verb=FALSE) {

    ## check presence
    if ( any(!c(value, rows, cols)%in%colnames(tmtf)) )
        stop("on of the requested values ",
             paste(value, rows, cols, collapse=";"),
             " is not present in the data")
    
    ## sorting of row and column classes
    if ( !missing(row.srt) ) aas <- row.srt
    else aas <- sort(unique(x[,rows]))
    if ( !missing(col.srt) ) acod <- col.srt
    else acod <- sort(unique(x[,cols]))

    ## remove rows with missing values
    if ( na.rm ) {
        ina <- is.na(x[,value])
        if ( sum(ina)>0 ) {
            wrn <- paste0(sum(ina), " rows with NA values removed")
            cat(paste("WARNING: ", wrn, "\n"))
            warning(wrn)
            x <- x[!ina,]
        }
    }
    
    ## filter for available classes
    aas <- aas[aas%in%x[,rows]]
    acod <- acod[acod%in%x[,cols]]

    ## result structures
    tt <- matrix(0, ncol=length(acod), nrow=length(aas))
    colnames(tt) <- acod
    rownames(tt) <- aas
    tc <- tp <- tm <- td <- uq <- tt
    tp[] <- 1
    tm[] <- td[] <- NA

    min.obs=2
    allvals <- list()
    for ( i in seq_along(aas) ) {

        aa <- aas[i]
        allvals[[aa]] <- list()

        for ( j in seq_along(acod) ) {

            cod <- acod[j]

            ## get all rows for the current intersection
            idx <- x[,cols]==cod & x[,rows]==aa
            
            if ( any(is.na(idx)) ) {
                wrn <- paste0(sum(is.na(idx)), " NA in ", aa, "/", cod)
                cat(paste0("WARNING: ", wrn,"\n"))
                warning(wrn)
                idx[is.na(idx)] <- FALSE
            }

            ## local background: column type!
            bgidx <- !idx
            if ( bg ) 
                bgidx <- x[,cols]==cod
            if ( !replace )
                bgidx <- bgidx & !idx

            if ( any(is.na(bgidx)) ) {
                wrn <- paste0(sum(is.na(bgidx)), " NA in ",
                              aa, "/", cod, "background")
                cat(paste0("WARNING: ", wrn,"\n"))
                warning(wrn)
                idx[is.na(idx)] <- FALSE
            }

            ## collect unique counter by ID column
            if ( !missing(id) )                 
                uq[i,j] <- length(unique(x[idx,id]))
            
            if (verb>0) cat(paste("testing",cod,"to", aa, "\n"))

            ## get all RAAS values
            y <- unlist(x[  idx, value])
            X <- unlist(x[bgidx, value])

            ## store all values for this class
            allvals[[aa]][[cod]] <- y

            ## value count
            tc[i,j] <- length(y)

            if ( length(y)==0 ) {
                if ( verb>0 )
                    cat(paste0("WARNING: no values for ",aa, "/", cod,"\n"))
                next
            }

            ## MEAN AND MEDIAN
            ## delog for mean and median
            ly <- y
            if ( delog ) 
                ly <- 10^y
            ymn <- mean(ly)
            ymd <- median(ly)
            if ( delog ) {
                ymn <- log10(ymn)
                ymd <- log10(ymd)
            }
            tm[i,j] <- ymn
            td[i,j] <- ymd

            ## T-TEST or RANK SUM TEST
            if ( sum(!is.na(y)) >= min.obs ) {
                
                tts <- use.test(y, X)

                tt[i,j] <- tts$statistic
                tp[i,j] <- tts$p.value
                
                if ( do.plots ) {
                    outfile <- paste0(fname,cod,"_",aa,".png")
                    hcol <- ifelse(tts$statistic<0, "#0000ff","#ff0000")
                    if ( verb>0 ) cat(paste("plotting", outfile, "\n"))
                    png(outfile, width=3, height=2,
                        units="in", res=100)
                    par(mai=c(.5,.15,.25,.15), tcl=-.25, mgp=c(1.3,.3,0))
                    brks <- seq(min(X), max(X), length.out=20)
                    hist(y, freq=FALSE, border=NA,
                         col=paste0(hcol,77),
                         xlim=range(X), main=NA,
                         xlab=xlab, breaks=brks, axes=FALSE)
                    par(new=TRUE)
                    hist(X, freq=FALSE, col="#77777777", border=NA,
                         xlim=range(X), main=NA, xlab=NA, ylab=NA, axes=FALSE,
                         breaks=brks)
                    abline(v=median(X), col=1, lwd=2)
                    abline(v=td[i,j], col=hcol, lwd=2)
                    axis(1)
                    mtext(as.expression(bquote(.(cod) %->% .(aa))),
                          3,-.1, font=2, cex=1.5)
                    legend(ifelse(tts$statistic<0, "topright","topleft"),
                           c(paste0("n=", length(y)),
                             paste0("u/t=", round(tts$statistic,1)),
                             paste0("p=", signif(tts$p.value,1))),
                           bty="n", seg.len=0, x.intersp=0)
                    dev.off()
                }
            }
        }
    }
    if ( verb>0 ) cat(paste("DONE\n"))
    
    ## construct overlap object
    ova <- list()
    ova$p.value <- tp
    ova$statistic <- tt
    ova$count <- tc
    ova$mean <- tm
    ova$median <- td
    ova$unique <- uq

    ova$delog <- delog

    ## sign of change, used in plotOverlaps
    sg <- sign(tt)
    sg[is.na(sg)] <- 1
    sg[sg==0] <- 1
    ova$sign <- sg
    
     ## add counts
    ova$num.target <-
        t(as.matrix(table(x[,cols])[colnames(tp)]))
    ova$num.query <-
        as.matrix(table(x[,rows]))[rownames(tp),,drop=FALSE]


    ## attach all raw values
    ova$ALL <- allvals
    
    class(ova) <- "clusterOverlaps"
    ova
}
