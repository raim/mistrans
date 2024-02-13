w.test <- function(x,y) {
    res <- wilcox.test(x,y)
    ## normalized U-statistic
    tt <- res$statistic/(sum(!is.na(x))*sum(!is.na(y))) -0.5
    rt <- list()
    rt$statistic <- unlist(tt)
    rt$p.value <- unlist(res$p.value)
    rt
}



## volcano plot function for overlap class
## produced by raasProfile
volcano <- function(ovl, cut=15, value="mean", p.txt=6, v.txt, mid=0, ...) {

    
    tmv <- gdata::unmatrix(t(ovl[[value]]))
    ##tmv <- c(ovl[[value]])
    tpv <- -log10(gdata::unmatrix(t(ovl$p.value)))
    tpv[tpv>cut] <- cut

    na <- is.na(tmv) | is.na(tpv)
    tpv <- tpv[!na]
    tmv <- tmv[!na]
    
    dense2d(tmv, tpv, ylab=expression(-log[10](p)), ...)
    axis(4)

    show <- tpv>p.txt
    if ( !missing(v.txt) )
        show <- show & (tmv<v.txt[1] | tmv>v.txt[2])

    show <- which(show)
    if ( length(show) )
        shadowtext(tmv[show], tpv[show], labels=names(tpv)[show],
                   pos=ifelse(tmv[show]<mid, 2,4), col="red",
                   font=2, cex=.8, xpd=TRUE)
}

raasProfile <- function(x=cdat, id="SAAP", values=tmt,
                        rows="to", cols="aacodon",
                        use.test=use.test,
                        do.plots=interactive(), verb=FALSE) {

    ## TODO: use full genetic code for consistent columsn and rows!
    aas <- sort(unique(x[,rows]))
    acod <- sort(unique(x[,cols]))
    codt <- list()
    
    tt <- matrix(0, ncol=length(acod), nrow=length(aas))
    colnames(tt) <- acod
    rownames(tt) <- aas
    tc <- tp <- tm <- td <- tt
    tp[] <- 1
    tm[] <- td[] <- NA


    min.obs=2
    for ( i in seq_along(aas) ) {
        for ( j in seq_along(acod) ) {

            aa <- aas[i]
            cod <- acod[j]
        
            ## get all SAAP for this codon
            idx <- x[,cols]==cod & x[,rows]==aa
            
            csap <- unique(x[idx,id])

            if ( sum(idx)==0 ) next
            
            ##cat(paste("testing",cod,"to", aa, "\n"))
            
            ## get all RAAS values
            rm <- !csap%in%names(values)
            if ( sum(rm) & verb ) 
                cat(paste0(aa, cod, ": ", sum(rm), " ID not found in values\n"))

            csap <- csap[!rm]
            y <- unlist(values[csap])
            X <- unlist(values) 
            
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
                    hist(y, freq=FALSE, border=2, xlim=range(X), ylim=c(0,.6))
                    hist(X, freq=FALSE, add=TRUE, xlim=range(X))
                }
            }
        }
    }
    
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
    
    class(ova) <- "clusterOverlaps"
    ova
}
