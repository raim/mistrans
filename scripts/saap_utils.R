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

            ## if a local background is chosen, tmt
            ## should be a table with the bg column
            ## and we filtering TMT only for this subset
            if ( bg ) {
                cat(paste("getting subset of values",cod,"\n"))
                vals <- values[values[,cols]==cod,]
                vals <- split(vals[,vid], vals[,id])
            } else vals <- values
                
            
            ## get all RAAS values
            rm <- !csap%in%names(vals)
            if ( sum(rm) & verb ) 
                cat(paste0(aa, cod, ": ", sum(rm), " ID not found in values\n"))

            csap <- csap[!rm]
            y <- unlist(vals[csap])
            X <- unlist(vals) 
            
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
