
## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
source("~/work/mistrans/scripts/raasprofiles3_init.R")

afig.path <- file.path(fig.path,"aminoacids")
dir.create(afig.path, showWarnings=FALSE)

### START ANALYSIS

xl.raas <- expression(log[10](RAAS)) # *bar(RAAS))
xl.raaa <- expression(log[10](RAAS))
xl.raau <- expression(log[10]*bar(RAAS[unique]))


## RAAS profiles by AAS site in peptides
## TODO: median RAAS of unique BP/SAAP vs. rel.position,
##       similar to iupred3
for ( ds in auds ) {

    tmtd <- tmtf
    if ( ds!="all" ) {
        if ( ds!="cancer" )
            tmtd <- tmtf[tmtf$Dataset==ds,]
        else
            tmtd <- tmtf[tmtf$Dataset!="Healthy",]
    }
    plotdev(file.path(afig.path,paste0("AASsite_",SETID,"_cor_",ds)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tmtd$RAAS.median, tmtd$rsite,
            xlab=NA,
            ylab="relative AAS position in peptide")
    mtext(xl.raau, 1, 1.6)
    figlabel(ds, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    dev.off()
}


fname <- file.path(dpath,paste0("AASsite_",SETID,"_wtests_"))
ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
                    delog=TRUE, rows="rsite.bins", cols="Dataset",
                    bg=TRUE, row.srt=rev(rsrt), col.srt=uds,
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    fname=fname, verb=0)
plotProfiles(ovw, fname=file.path(afig.path,paste0("AASsite_",SETID)),
             mai=c(.8,1,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             ttcols=ttcols, value="median",
             mtxt="relative AAS position", mtxt.line=3.5,
             rlab=LAB, llab="", ftyp=ftyp,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols, verb=1)
  
## RAAS profiles by AA properties
fname <- file.path(dpath,paste0("AAprop_",SETID,"_wtests_"))
ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
                    delog=TRUE, rows="pfromto", cols="Dataset",
                    bg=TRUE, row.srt=srt, col.srt=uds,
                    use.test=use.test, do.plots=TRUE,
                    xlab=xl.raas,
                    fname=fname, verb=0)
ovwp <- sortOverlaps(ovw, axis=2, p.min=p.min, cut=TRUE)
ovwp <- sortOverlaps(ovwp, axis=2, srt=sort(rownames(ovwp$p.value)))



## plot all
plotProfiles(ovw, fname=file.path(afig.path,paste0("AAprop_",SETID)),
             mai=c(.7,1.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="", ftyp=ftyp,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols, verb=1)
if ( nrow(ovwp$p.value)>0 )
    plotProfiles(ovwp,
                 fname=file.path(afig.path,paste0("AAprop_",SETID,"_cut")),
                 mai=c(.7,1.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="", ftyp=ftyp,
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols)

## TODO: median vs. p-value as a measure of effect size
if ( interactive() ) {
    plot(ovw$count, ovw$median, xlab="count", ylab="RAAS median")
    plot(ovw$count, -log(ovw$p.value), xlab="count", ylab=expression(-log(p)))
    plot(ovw$median, -log(ovw$p.value))
}


## by "from" amino acid
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="from", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0)
plotProfiles(ovw, fname=file.path(afig.path,paste0("fromAA_",SETID)),
             mtxt="encoded AA",
             mai=c(.8,.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="", ftyp=ftyp,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## by "to" amino acid
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="to", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0)
plotProfiles(ovw, fname=file.path(afig.path,paste0("toAA_",SETID)),
             mai=c(.8,.5,.5,.5),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             mtxt="substituted AA", ttcols=ttcols, value="median",
             rlab=LAB, llab="", ftyp=ftyp,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)


## plot 16 plots for pfrom-to combos, for each to all AA
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="fromto", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=TRUE,
                   fname=file.path(dpath,paste0("AA_",SETID,"_")),
                   xlab=xl.raas,
                   verb=0)

## by AA property transition
for ( ptype in unique(tmtf$pfromto) ) {

    ## sort by to
    qsrt <- sort(unique(tmtf$fromto[tmtf$pfromto==ptype]))
    ft <- unlist(lapply(strsplit(qsrt,":"), function(x) x[2]))
    qsrt <- qsrt[order(ft)]
                 
    ovwp <- sortOverlaps(ovw, srt=qsrt)

    plotProfiles(ovwp, fname=file.path(afig.path,paste0(ptype,"_",SETID)),
                 mai=c(.8,1.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="", ftyp=ftyp,
                 mtxt=ptype, mtxt.line=3,
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols)
    ## cut
    ovwp <- sortOverlaps(ovwp, axis=2, p.min=p.txt, cut=TRUE)
    if ( nrow(ovwp$p.value)>1 ) #TODO: fix error with nrow=1 in 
        plotProfiles(ovwp,
                     fname=file.path(afig.path,paste0(ptype,"_",SETID,"_cut")),
                     mai=c(.8,1.5,.5,.5), ttcols=ttcols, value="median",
                     p.min=p.min, p.txt=p.txt,
                     dot.sze=dot.sze, p.dot=p.dot,
                     rlab=LAB, llab="", ftyp=ftyp,
                     mtxt=ptype, mtxt.line=3,
                     vcols=vcols, vbrks=vbrks,
                     gcols=gcols)
    
}

## plot only significantly high RAAS
ovwp <-ovw

## only significant
ovwp <- sortOverlaps2(ovwp, axis=2, p.min=p.min, cut=TRUE)

if ( nrow(ovwp$p.value)>0 ) {

    ## sort plot
    ## by from AA and internally by to AA
    oldsrt <- rownames(ovwp$p.value)
    taa <- sub(".*:","",oldsrt)
    faa <- sub(":.*","",oldsrt)
    laa <- split(taa, faa)
    laa <- lapply(laa, function(x) x[order(match(x, aa.srt))])
    laa <- laa[order(match(names(laa),aa.srt))]
    laa <- unlist(laa)
    newsrt <- gsub("[0-9]","",paste0(names(laa),":",laa))
    ovwp <- sortOverlaps2(ovwp, axis=2, srt=newsrt)
    
    
    
    faa <- sub(":.*","",rownames(ovwp$p.value))
    ft.cols <- aa.cols[faa]
    names(ft.cols) <- rownames(ovwp$p.value)
    
    plotProfiles(ovwp,
                 fname=file.path(afig.path,paste0("AA_",SETID,"_cut")),
                 mai=c(.8,.6,.6,.6), fw=.2, fh=.2,
                 ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="",  ftyp=ftyp,
                 ffam="monospace",
                 axis2.col=ft.cols,
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols)

    ## calculate optimal figure height: result fields + figure margins (mai)
    nr <- nrow(ovwp$p.value)
    nc <- ncol(ovwp$p.value)
    mai <- c(.8,.6,.6,.6)
    nh <- nr *fh + mai[1] + mai[3]
    nw <- nc *fw + mai[2] + mai[4]
    ffam <- "monospace"

    rsrt <- rownames(ovwp$p.value)
    taa <- sub(".*:","", rsrt)
    faa <- sub(":.*","", rsrt)

    fname <- file.path(afig.path,paste0("AA_",SETID,"_cut"))
    ## combined effect size and p-value plot
    plotdev(paste0(fname,"_dotplot_manual"),
            height=nh, width=nw, res=300, type=ftyp)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
    dotprofile(x=ovwp, value="median", vbrks=vbrks,
               vcols=vcols, p.dot=p.dot,
               dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
    box()
    axis(1,  1:ncol(ovwp$p.value),
         labels=colnames(ovwp$p.value), las=2, family=ffam)
    ##axis(2, length(axex):1, labels=axex, las=2, family=ffam)
    shadowtext(rep(-2, nr), nr:1, faa, col=aa.cols[faa], xpd=TRUE,
               font=2, r=.1)
    text(rep(-1,nr), nr:1, expression(""%->%""), xpd=TRUE)
    shadowtext(rep( 0, nr), nr:1, taa, col=aa.cols[taa], xpd=TRUE,
               font=2, r=.1)
    axis(3, at=1:ncol(ovwp$num.target),
         labels=format(ovwp$num.target[1,], big.mark=",", trim=TRUE),las=2)
    axis(4, at=nrow(ovwp$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
         labels=format(ovwp$num.query[,1], big.mark=",", trim=TRUE),las=2)

    dev.off()
 }


## TODO: plot on p-value correction
plot(ovw$p.value, qvalue::qvalue(c(ovw$p.value))$qvalues)


## by AA->AA
for ( ds in auds ) {

    tmtd <- tmtf
    if ( ds!="all" ) {
        if ( ds!="cancer" )
            tmtd <- tmtf[tmtf$Dataset==ds,]
        else
            tmtd <- tmtf[tmtf$Dataset!="Healthy",]
    }

    ## NOTE: bg=FALSE, no column-wise background!
    ovw <- raasProfile(x=tmtd, id="SAAP", 
                       rows="to", cols="from",
                       bg=FALSE, value="RAAS", 
                       use.test=use.test, do.plots=FALSE, 
                       verb=0)

    plotProfiles(ovw, fname=file.path(afig.path,paste0("AA_",SETID,"_",ds)),
                 mai=c(.8,.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=ds, ftyp=ftyp,
                 vcols=vcols, vbrks=vbrks,
                 gcols=gcols)
}


### BY STRUCTURAL FEATURES
for ( ds in auds ) {

    tmtd <- tmtf
    dsl <- ""
    if ( ds!="all" ) {
        if ( ds!="cancer" )
            tmtd <- tmtf[tmtf$Dataset==ds,]
        else
            tmtd <- tmtf[tmtf$Dataset!="Healthy",]
        dsl <- ds
    }
    ovl <- clusterCluster(tmtd$iupred3.bins, paste0(tmtd$raas.bins), 
                          cl1.srt=c(rev(levels(iupred3.bins)),"na"),
                          cl2.srt=raas.srt)
    plotdev(file.path(afig.path,paste0("structure_iupred3_",SETID,"_ovl_",ds)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.9,.9,.5,.5), mgp=c(3.3,.3,0), tcl=-.25)
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt,
                 xlab=NA, ylab=NA,
                 show.total=TRUE)
    mtext("disordered score", 2, 3.5)
    mtext(xl.raau, 1, 3.5)
    figlabel(dsl, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    dev.off()
    
    plotdev(file.path(afig.path,paste0("structure_iupred3_",SETID,"_cor_",ds)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tmtd$RAAS.median, tmtd$iupred3,
            xlab=NA, ylab="disordered score")
    mtext(xl.raau, 1, 1.6)
    figlabel(dsl, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    dev.off()
}

## TODO: IUPRED3 vs. RAAS bins
##ovw <- raasProfile(x=tmtf, id="SAAP", 
##                   rows="iupred3.bins", cols="Dataset",
##                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
##                   col.srt=uds,
##                   use.test=use.test, do.plots=FALSE, xlab="TMT level RAAS",
##                   verb=0)
##

## IUPRED3
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="iupred3.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("iupred3_",SETID,"_")))

plotProfiles(ovw, fname=file.path(afig.path,paste0("structure_iupred3_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="disordered score", mtxt.line=3.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## ANCHOR2
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="anchor2.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(rsrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, fname=file.path(afig.path,"anchor2_"))

plotProfiles(ovw, fname=file.path(afig.path,paste0("structure_anchor2_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, ftyp=ftyp,
             mtxt="ANCHOR2", mtxt.line=3.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)

## S4 PRED SECONDARY STRUCTURE
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="sstruc", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(ssrt), col.srt=uds,
                   use.test=use.test, do.plots=TRUE,
                   xlab=xl.raas,
                   verb=0, fname=file.path(dpath,"structure_s4pred_"))

plotProfiles(ovw, fname=file.path(afig.path,paste0("structure_s4pred_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, ftyp=ftyp,
             mtxt="SPRED4", mtxt.line=3.3,
             vcols=vcols, vbrks=vbrks,
             gcols=gcols)



if ( FALSE ) {

    ## developing dot plots for combined effect/p plot
plotOverlaps(ovw, p.min=p.min, p.txt=p.txt, txt.col=NA,
                 text.cex=.8, axis=1:2, ylab=NA, xlab=NA,
                 col=ttcols, show.total=TRUE)
    mai=c(.8,.9,.5,.5)
    value <- "median"
    vcols <- vcols
    vbrks <- vbrks
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25)
    navals <- ovw[[value]]
    navals[] <- NA
    image_matrix(navals, breaks=vbrks,
                 col=vcols, axis=1, xlab=NA, ylab=NA)
    p <- -log10(ovw$p.value)
    p[p>-log10(p.min)] <- -log10(p.min)
    z <- p/-log10(p.min)
    points(x = rep(1:ncol(z), nrow(z)),
           y = rep(nrow(z):1, each= ncol(z)),
           cex=1.5*c(t(z)), pch=19,
           col=num2col(t(ovw[[value]]),limits=range(vbrks),
                       colf=viridis::viridis, n=length(vcols)))
    axis(2, length(axex):1, labels=axex, las=2)
    axis(4, at=nrow(ovw$num.query):1, labels=ovw$num.query[,1],las=2)
    axis(3, at=1:ncol(ovw$num.target), labels=ovw$num.target[1,],las=2)
    if ( !missing(mtxt) ) mtext(mtxt, 2, mtxt.line)
    if ( !missing(llab) ) figlabel(llab, pos="bottomleft", cex=1.2)
    if ( !missing(rlab) ) figlabel(rlab, pos="bottomright", cex=.8)

    
}
