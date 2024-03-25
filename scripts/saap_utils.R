
require(segmenTools)

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

dotprofile <- function(x, value, vcols=viridis::viridis(100),
                       vbrks, p.dot=1e-10, dot.sze=c(.3,2),
                       lg2=FALSE, mxr,
                       show.total=FALSE,
                       test=FALSE, ...) {

    ## values for coloring
    vals <- x[[value]]
    if ( lg2 ) {
        vals <- log2(vals)
        if ( missing(mxr) ) mxr <- max(abs(vals))
        vals[vals >  mxr] <-  mxr
        vals[vals < -mxr] <- -mxr
    }

    ## breaks
    if ( missing(vbrks) )
        vbrks <- seq(min(vals,na.rm=TRUE), max(vals, na.rm=TRUE), length.out=length(vcols)+1)

    navals <- vals
    navals[] <- NA
    image_matrix(navals, breaks=vbrks, col=vcols, ...)
    
    p <- -log10(x$p.value)
    p[p>-log10(p.dot)] <- -log10(p.dot)
    z <- p/-log10(p.dot)


    
    ## intersect data to colors
    cols <- vcols[findInterval(t(vals),
                               seq(min(vbrks), max(vbrks),
                                   length.out=length(vbrks)),
                               all.inside = TRUE)]
    d.sze <- dot.sze[1]+dot.sze[2]*c(t(z))
    if ( test )
        points(x = rep(1:ncol(z), nrow(z)),
               y = rep(nrow(z):1, each= ncol(z)), cex=d.sze)
    points(x = rep(1:ncol(z), nrow(z)),
           y = rep(nrow(z):1, each= ncol(z)),
           cex=d.sze, pch=19,
           col=cols)
    
    toty <- totx <- FALSE
    if (is.logical(show.total)) {
        if (show.total) 
            toty <- totx <- TRUE
    }
    else if (is.character(show.total)) {
        if (show.total == "x") 
            totx <- TRUE
        if (show.total == "y") 
            toty <- TRUE
        if (show.total %in% c("xy", "yx")) 
            toty <- totx <- TRUE
    }
    if (toty) 
        if ("num.query" %in% names(x)) 
            axis(4, at = length(x$num.query):1, labels = x$num.query, 
                las = 2, lwd = 0, lwd.ticks = 1)
    if (totx) 
        if ("num.target" %in% names(x)) 
            axis(3, at = 1:length(x$num.target), labels = x$num.target, 
                las = 2, lwd = 0, lwd.ticks = 1)
    if (toty | totx) 
        figlabel("total", region = "figure", pos = "topright", 
            cex = par("cex"))
}

## testing to use enrichment of overlap tables
## as dot size -> too big for low numbers
plotOverlaps2 <- function(ovw,
                          col,
                          max.val=3, 
                          p.min=1e-10, p.txt=1e-5,
                          dot.sze=c(.2,3), bg=TRUE,
                          ...) {
    
    ## p-value scaling
    p <- -log10(ovw$p.value)
    p[p>-log10(p.min)] <- -log10(p.min)
    p <- p/-log10(p.min)

    p <- p*ovw$sign
    
    ## enrichment ratio scaling
    r <- abs(log2(ovw$ratio))
    r[r>max.val] <- max.val
    r <- r/max.val
    

    ## point size
    d.cex <- dot.sze[1]+dot.sze[2]*c(t(r))


    ## p-value colors
    if ( missing(col) ) {
        docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
        upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
        col <- unique(c(rev(docols), upcols))
    }
    brks <- seq(-1, 1, length.out = length(col)+1)
    cols <- col[findInterval(t(p), brks, all.inside = TRUE)]

    navals <- ovw[["ratio"]]
    navals[] <- NA
    image_matrix(navals,  col=col, breaks=brks, ...)

    if ( bg )
        points(x = rep(1:ncol(r), nrow(r)),
               y = rep(nrow(r):1, each= ncol(r)), cex=d.cex+.5, pch=1)

    points(x = rep(1:ncol(r), nrow(r)),
           y = rep(nrow(r):1, each= ncol(r)),
           cex=d.cex, pch=19,
           col=cols)
    return(cbind(c(p),cols))
}

plotProfiles <- function(ovw, mai=c(.6,.5,.5,.5),
                         fw=.25, fh=.2,
                         ttcols=ttcols, p.min=p.min, p.txt=p.txt,
                         dot.sze=dot.sze, p.dot=p.dot, ## TODO: legend
                         value="median", vcols, vbrks,
                         count="unique", gcols=gcols,
                         fname="profile", mtxt, mtxt.line=1.3, llab, rlab,
                         col.lines, # column classes - vertical lines
                         verb=0) {

    ## replace row labels with arrows
    rows <- rownames(ovw$p.value)
    axex <- rep("",length(rows))
    names(axex) <- rows
    for ( i in seq_along(rows) ) {
        if ( length(grep(":",rows[i])) ) {
            ft <- unlist(strsplit(rows[i],":"))
            axex[i] <- as.expression(bquote(.(ft[1]) %->% .(ft[2])))
        } else axex[i] <- as.expression(bquote(.(rows[i])))
    }
  

    
    ## calculate optimal figure height: result fields + figure margins (mai)
    nh <- nrow(ovw$p.value) *fh + mai[1] + mai[3]
    nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]


    ## p-value profiles
    outf <- paste0(fname,"_wtests")
    if (verb>0) cat(paste("plotting test profile",outf,p.min,p.txt,"\n"))
    segmenTools::plotdev(outf,
                         height=nh, width=nw, res=300)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovw, p.min=p.min, p.txt=p.txt,
                 text.cex=.8, axis=1, ylab=NA, xlab=NA,
                 col=ttcols, show.total=TRUE)
    axis(2, length(axex):1, labels=axex, las=2)
    if ( !missing(mtxt) ) mtext(mtxt, 2, mtxt.line)
    if ( !missing(llab) ) figlabel(llab, pos="bottomleft", cex=1.2)
    if ( !missing(rlab) ) figlabel(rlab, pos="bottomright", cex=.8)
    if ( !missing(col.lines) )
        abline(v=.5+col.lines, lwd=1, xpd=FALSE)
    dev.off()

    if ( missing(vbrks) ) {
        vals <- ovw[[value]]
        vals <- vals[!is.na(vals)]
        vbrks <- seq(min(vals), max(vals), length.out=101)
        vcols <- viridis::viridis(100)
    }

    ## combined effect size and p-value plot
    if ( verb ) cat(paste("plotting dotplot\n"))
    segmenTools::plotdev(paste0(fname,"_dotplot"),
                         height=nh, width=nw, res=300)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
    dotprofile(x=ovw, value=value, vbrks=vbrks,
               vcols=vcols, p.dot=p.dot, dot.sze=dot.sze,
               axis=1, xlab=NA, ylab=NA)
    box()
    axis(2, length(axex):1, labels=axex, las=2)
    axis(3, at=1:ncol(ovw$num.target), labels=ovw$num.target[1,],las=2)
    axis(4, at=nrow(ovw$num.query):1, labels=ovw$num.query[,1],las=2)
    if ( !missing(mtxt) ) mtext(mtxt, 2, mtxt.line)
    if ( !missing(llab) ) figlabel(llab, pos="bottomleft", cex=1.2)
    if ( !missing(rlab) ) figlabel(rlab, pos="bottomright", cex=.8)
    if ( !missing(col.lines) )
        abline(v=.5+col.lines, lwd=1, xpd=FALSE)
    dev.off()
    ##return(p)

    if (verb>0) cat(paste("plotting",value,"RAAS\n"))
    
    segmenTools::plotdev(paste0(fname,"_raas_",value),
                         height=nh, width=nw, res=300)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
    txt <- ovw[["count"]]
    txt[txt==0] <- ""
    q1 <- quantile(c(ovw[[value]]), probs=.4, na.rm=TRUE)
    txt.col <- ifelse(ovw[[value]]<q1, "white","black")
    image_matrix(ovw[[value]], col=vcols, breaks=vbrks, axis=1,
                 text=txt, text.col=txt.col, xlab=NA, ylab=NA, text.cex=.8)
    axis(2, length(axex):1, labels=axex, las=2)
    if ( !missing(mtxt) ) mtext(mtxt, 2, mtxt.line)
    if ( !missing(llab) ) figlabel(llab, pos="bottomleft", cex=1.2)
    if ( !missing(rlab) ) figlabel(rlab, pos="bottomright", cex=.8)
    if ( !missing(col.lines) )
        abline(v=.5+col.lines, lwd=1, xpd=FALSE, col="white")
    dev.off()

    if (verb>0) cat(paste("plotting unique number\n"))
   
    segmenTools::plotdev(paste0(fname,"_count_unique"),
                         height=nh, width=nw, res=300)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
    cnt <- ovw$unique
    cnt[cnt==0] <- NA
    txt <- ovw$unique
    txt[txt==0] <- ""
    txt.col <- ifelse(ovw$unique>quantile(c(ovw$unique),.95),"white","black")
    image_matrix(cnt, col=gcols, axis=1:2,
                 text=txt, text.col=txt.col, ylab=NA, xlab=NA, text.cex=.8)
    ##axis(2, length(axex):1, labels=axex, las=2)
    if ( !missing(mtxt) ) mtext(mtxt, 2, mtxt.line)
    if ( !missing(llab) ) figlabel(llab, pos="bottomleft", cex=1.2)
    if ( !missing(rlab) ) figlabel(rlab, pos="bottomright", cex=.8)
    if ( !missing(col.lines) )
        abline(v=.5+col.lines,
               lwd=1, xpd=FALSE)
    dev.off()

    if (verb>0) cat(paste("plotting volcano\n"))

    segmenTools::plotdev(paste0(fname,"_volcano"),
                         height=3, width=4, res=300)
    par(mai=c(.5,.75,.1,.75), mgp=c(1.3,.3,0), tcl=-.25)
    volcano(ovw, cut=100, p.txt=-log10(p.txt),
            v.txt=c(-Inf,-1), density=TRUE,
            xlab=paste(value,"TMT RAAS"), value=value)
    abline(v=get(value, mode="function")(ovw[[value]],na.rm=TRUE))
    if ( !missing(llab) ) figlabel(llab, pos="bottomleft", cex=1.2)
    if ( !missing(rlab) ) figlabel(rlab, pos="bottomright", cex=.8)
    dev.off()
    
    if (verb>0) cat(paste("plotting histograms\n"))

    ## common histogram with color by sign and
    ## first, calculate densities
    dns <- ovw$ALL
    for ( i in 1:nrow(ovw$p.value) ) {
        for ( j in 1:ncol(ovw$p.value) ) {
            rw <- rownames(ovw$p.value)[i]
            cl <- colnames(ovw$p.value)[j]
            vls <- ovw$ALL[[rw]][[cl]]
            if (length(vls)>1 ) {
            dns[[rw]][[cl]] <- density(vls)
            }
        }
    }
    ## get max
    mxs <- rep(NA, length(dns))
    for ( i in seq_along(dns) )
        mxs[i] <- max(unlist(lapply(dns[[i]],
                                    function(z)
                                        ifelse("y"%in%names(z),
                                               max(unlist(z["y"])),-10))))
    
    segmenTools::plotdev(paste0(fname,"_densities"),
                         res=300, width=4, height=2)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    dns <- ovw$ALL
    hist(unlist(dns), col="#77777755", border=NA, 
         xlim=c(-6,4),
         xlab=expression(TMT~level~log[10]*RAAS), ylab=NA, main=NA, axes=FALSE)
    par(new=TRUE)
    plot(1, col=NA, xlim=c(-6,4), ylim=c(0,max(mxs)),
         xlab=NA, ylab=NA, axes=FALSE)
    for ( i in 1:nrow(ovw$p.value) ) {
        for ( j in 1:ncol(ovw$p.value) ) {
            rw <- rownames(ovw$p.value)[i]
            cl <- colnames(ovw$p.value)[j]
            vls <- ovw$ALL[[rw]][[cl]]
            if (length(vls)>1 ) {
                dns[[rw]][[cl]] <- density(vls)
                if ( ovw$p.value[i,j] <1e-3 )
                    lines(dns[[rw]][[cl]],
                          col=ifelse(ovw$sign[i,j]==-1,"blue","red"),
                          lwd=-log10(ovw$p.value[i,j])/-log10(p.min))
                ## TODO: annotate
            }
        }
    }
    axis(2)#, labels=NA)
    mtext("density",2,1.3)
    axis(1)
    dev.off()
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
                    lg2=FALSE, mxr, density=TRUE, xlab, adjust, ...) {

    ## values for coloring
    vals <- ovl[[value]]
    if ( lg2 ) {
        vals <- log2(vals)
        if ( missing(mxr) ) mxr <- max(abs(vals))
        vals[vals >  mxr] <-  mxr
        vals[vals < -mxr] <- -mxr
    }

    
    tmv <- gdata::unmatrix(t(vals))
    ##tmv <- c(ovl[[value]])

    pvals <- ovl$p.value
    if ( !missing(adjust) )
        pvals[] <- p.adjust(c(pvals), method=adjust)
    
    
    tpv <- -log10(gdata::unmatrix(t(pvals)))
    tpv[tpv>cut] <- cut

    na <- is.na(tmv) | is.na(tpv)
    tpv <- tpv[!na]
    tmv <- tmv[!na]

    if ( length(tmv)<2 ) 
        density <- FALSE

    if ( missing(xlab) ) xlab <- value
    
    if ( density ) {
        dense2d(tmv, tpv, ylab=expression(-log[10](p)), xlab=xlab, ...)
    } else {
        plot(tmv, tpv, ylab=expression(-log[10](p)), xlab=xlab, ...)
    }
    axis(4)

    ## indicate high
    show <- tpv>p.txt
    if ( !missing(v.txt) )
        show <- show & (tmv<v.txt[1] | tmv>v.txt[2])

    show <- which(show)
    if ( length(show) ){
        if ( missing(mid) )
            mid <- mean(vals,na.rm=TRUE)
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
                        bg=FALSE, bg.dir="col", na.rm=FALSE,
                        rows="to", cols="aacodon",
                        row.srt, col.srt, filter=TRUE,
                        use.test=use.test,
                        do.plots=FALSE,
                        xlab="value", fname="profile_",
                        verb=FALSE) {

    ## check presence
    if ( any(!c(value, rows, cols)%in%colnames(x)) )
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
    if ( filter ) {
        aas <- aas[aas%in%x[,rows]]
        acod <- acod[acod%in%x[,cols]]
    }
    
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
            if ( bg ) {
                if ( bg.dir=="col" )
                    bgidx <- x[,cols]==cod
                else if ( bg.dir=="row" )
                    bgidx <- x[,rows]==aa
            }
            if ( !replace ) ## do not include foreground in background!
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

## calculated two-sided hypergeo tests for an amino acid
## table, for each column
#' @param x a matrix of character values, e.g. amino acids, where
#' rows are protein sequences and columns relative positions
#' @param abc list of characters to consider, e.g. all amino acids
aaProfile <- function(x, abc, k=1, p.min, verb=0) {

    if ( missing(abc) )
        abc <- sort(unique(c(x)))

    aam <- matrix(1, ncol=ncol(x), nrow=length(abc))
    rownames(aam) <- abc
    colnames(aam) <- colnames(x)
    aac <- aap <- aar <- aam
    aac[] <- 0
    
    for ( i in 1:nrow(aam)   ) {
        if ( verb>0 ) cat(paste("testing", rownames(aam)[i], "\n"))
        for ( j in  1:ncol(aam) ) {

            if ( verb>1 ) cat(paste("\t", colnames(aam)[j], "\n"))
            
            aa <- rownames(aam)[i]
            pos <- colnames(aam)[j]
            
            m <- sum(c(x)==aa) # white balls
            n <- sum(x%in%abc)-m # black balls
            N <- m+n # total number of balls

            q <- sum(x[,j]==aa)    # white balls drawn - AT CURRENT POSITION
            k <- sum(x[,j]%in%abc) # number of balls drawn
            
            aac[i,j] <- q
            ## enriched
            aam[i,j] <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)
            ## deprived
            aap[i,j] <- phyper(q=q, m=m, n=n, k=k, lower.tail=TRUE)
            ## frequency ratio
            aar[i,j] <- q/k * N/m
        }
    }
    ## FILTER SIGNIFICANT
    if ( !missing(p.min) ) {
        pval <- aam
        pval[aap<aam] <- aap[aap<aam]
        sig <- apply(pval, 1, function(x) any(x<p.min))
        aam <- aam[sig,]
        aap <- aap[sig,]
        aac <- aac[sig,]
        aar <- aar[sig,]
    }

    ## construct overlap class
    ovl <- list()
    ovl$count <- aac
    ovl$ratio <- aar
    ovl$p.value <- aam
    ovl$p.value[aap<aam] <- aap[aap<aam]
    ovl$sign <- ifelse(aap<aam,-1,1)
    ovl$num.query <-
        t(t(table(c(x))[rownames(aam)]))
    ovl$num.target <-
        t(apply(x, 2, function(x) sum(x%in%rownames(aam))))

    class(ovl) <- "clusterOverlaps"
    
    ovl
}
