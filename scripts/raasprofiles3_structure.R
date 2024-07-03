
### ANALYZE AMINO ACID AND PROTEIN STRUCTURAL PROPERTIES at AAS

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")
    
sfig.path <- file.path(fig.path,"structure")
dir.create(sfig.path, showWarnings=FALSE)

### START ANALYSIS

xl.raas <- expression(log[10](RAAS)) # *bar(RAAS))
xl.raaa <- expression(log[10](RAAS))
xl.raau <- expression(log[10]*bar(RAAS[unique]))
 
### BY STRUCTURAL FEATURES


## loop over bins
nms <- sapply(rank.vals, as.expression)
names(nms) <- rank.vals
nms["iupred3"] <- "disordered"
nms["DisoRDPbind"] <- "binding"
nms["MMSeq2"] <- "conserved"
nms["loglen"] <- expression(log[10](length))
for ( val in rank.vals ) {
    rows <- paste0(val,".bins")
    bsrt <- c("na",levels(get(rows)))

    ovw <- raasProfile(x=tmtf, id="SAAP", 
                       rows=rows, cols="Dataset",
                       bg=TRUE, value="RAAS", row.srt=rev(bsrt),
                       col.srt=uds,
                       use.test=use.test, do.plots=FALSE,
                       xlab=xl.raas,
                       verb=0)

    ina <- "na" %in% rownames(ovw$p.value)

    mmai <- c(.05,.5,.05,.7)
    nw <- ncol(ovw$p.value)*.2 + mmai[2] + mmai[4]
    nh <- nrow(ovw$p.value)*.2 + mmai[1] + mmai[3]
    fname <- file.path(sfig.path,paste0(val,"_bin_",SETID))
    plotdev(paste0(fname,"_dotplot_manual"),
            height=nh, width=nw, res=300, type=ftyp)
    par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25)#, family="monospace")
    dotprofile(x=ovw, value="median", vbrks=abrks,
               vcols=acols, p.dot=p.dot,
               dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
    box()
    mtext(nms[val], 2, 1.3)
    axis(4, at=nrow(ovw$p.value):1,
         labels=format(ovw$num.query[,1], big.mark=",", trim=TRUE),
         las=2, family="monospace")
    polygon(y=c(1+as.numeric(ina), nrow(ovw$p.value), nrow(ovw$p.value)),
            x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
    if ( ina )
        axis(2, at=1, label="na", las=2)
    dev.off()

}


## loop over ranks
for ( val in rank.vals ) {
    rows <- paste0(val,".rank")
 
    ovw <- raasProfile(x=tmtf, id="SAAP", 
                       rows=rows, cols="Dataset",
                       bg=TRUE, value="RAAS", row.srt=rev(ranksrt),
                       col.srt=uds,
                       use.test=use.test, do.plots=FALSE,
                       xlab=xl.raas,
                       verb=0, 
                       fname=file.path(dpath,paste0("iupred3_",SETID,"_")))

    ina <- "na" %in% rownames(ovw$p.value)
    
    mmai <- c(.05,.5,.05,.7)
    nw <- ncol(ovw$p.value)*.2 + mmai[2] + mmai[4]
    nh <- nrow(ovw$p.value)*.2 + mmai[1] + mmai[3]
    fname <- file.path(sfig.path,paste0(val,"_rank_",SETID))
    plotdev(paste0(fname,"_dotplot_manual"),
            height=nh, width=nw, res=300, type=ftyp)
    par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25)#, family="monospace")
    dotprofile(x=ovw, value="median", vbrks=abrks,
               vcols=acols, p.dot=p.dot,
               dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
    box()
    mtext(nms[val], 2, 1.3)
    axis(4, at=nrow(ovw$p.value):1,
         labels=format(ovw$num.query[,1], big.mark=",", trim=TRUE),
         las=2, family="monospace")
    polygon(y=c(1+as.numeric(ina), nrow(ovw$p.value), nrow(ovw$p.value)),
            x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
    if ( ina )
        axis(2, at=1, label="na", las=2)
    dev.off()
}




## TODO: align this with protein profiles "site" matrix,
## use unique site raas instead of unique Dataset/BP/SAAP Raas.
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

    ## calculate RAAS BINS
    ## NOTE: this approach leads to the same protein site being
    ## included multiple times!
    
    ovl <- clusterCluster(tmtd$iupred3.bins, paste0(tmtd$raas.bins), 
                          cl1.srt=c(rev(levels(iupred3.bins)),"na"),
                          cl2.srt=raas.srt)
    plotdev(file.path(sfig.path,paste0("structure_iupred3_",SETID,"_ovl_",ds)),
            type=ftyp, res=300, width=3, height=3)
    par(mai=c(.9,.9,.5,.5), mgp=c(3.3,.3,0), tcl=-.25)
    plotOverlaps(ovl, p.min=p.min, p.txt=p.txt,
                 xlab=NA, ylab=NA, show.total=TRUE)
    mtext("disordered score", 2, 3.5)
    mtext(xl.raau, 1, 3.5)
    figlabel(dsl, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    dev.off()
    
    plotdev(file.path(sfig.path,paste0("structure_iupred3_",SETID,"_cor_",ds)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tmtd$RAAS.median, tmtd$iupred3,
            xlab=NA, ylab="disorder, IUpred3")
    mtext(xl.raau, 1, 1.6)
    figlabel(dsl, pos="bottomleft", font=2, cex=1.2)
    figlabel(LAB, pos="bottomright", cex=.7)
    dev.off()

    plotdev(file.path(sfig.path,paste0("structure_flDPnn_",SETID,"_cor_",ds)),
            type=ftyp, res=300, width=3,height=3)
    par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
    plotCor(tmtd$RAAS.median, tmtd$flDPnn,
            xlab=NA, ylab="disorder, flDPnn")
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

## ASAquick
iusrt <- levels(ASAquick.bins)
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="ASAquick.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(iusrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("ASAquick_",SETID,"_")))

plotProfiles(ovw,
             fname=file.path(sfig.path,paste0("structure_ASAquick_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="ASAquick", mtxt.line=3.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)
## SCRIBER
iusrt <- levels(SCRIBER.bins)
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="SCRIBER.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(iusrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("SCRIBER_",SETID,"_")))

plotProfiles(ovw, fname=file.path(sfig.path,paste0("structure_SCRIBER_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="SCRIBER", mtxt.line=3.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)

## DisoRDPbind
iusrt <- levels(DisoRDPbind.bins)
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="DisoRDPbind.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(iusrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("DisoRDPbind_",SETID,"_")))

plotProfiles(ovw,
             fname=file.path(sfig.path,paste0("structure_DisoRDPbind_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="DisoRDPbind", mtxt.line=3.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)

## MMSeq2
iusrt <- c("na",levels(MMSeq2.bins))
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="MMSeq2.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(iusrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("MMSeq2_",SETID,"_")))

plotProfiles(ovw, fname=file.path(sfig.path,paste0("structure_MMSeq2_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="conservation,\nMMSeq2", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)

## flDPnn
iusrt <- levels(flDPnn.bins)
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="flDPnn.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(iusrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, 
                   fname=file.path(dpath,paste0("flDPnn_",SETID,"_")))

plotProfiles(ovw, fname=file.path(sfig.path,paste0("structure_flDPnn_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="disorder, flDPnn", mtxt.line=3.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)

## IUPRED3
iusrt <- levels(iupred3.bins)
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="iupred3.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(iusrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0)

plotProfiles(ovw, fname=file.path(sfig.path,paste0("iupred3_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="disorder, IU3", mtxt.line=3.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)


## ANCHOR2
ansrt <- levels(anchor2.bins)
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="anchor2.bins", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(ansrt),
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0, fname=file.path(sfig.path,"anchor2_"))

plotProfiles(ovw, fname=file.path(sfig.path,paste0("structure_anchor2_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, ftyp=ftyp,
             mtxt="ANCHOR2", mtxt.line=3.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)

## S4 PRED SECONDARY STRUCTURE
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="sstruc", cols="Dataset",
                   bg=TRUE, value="RAAS", row.srt=rev(ssrt), col.srt=uds,
                   use.test=use.test, do.plots=TRUE,
                   xlab=xl.raas,
                   verb=0, fname=file.path(dpath,"structure_s4pred_"))

plotProfiles(ovw, fname=file.path(sfig.path,paste0("structure_s4pred_",SETID)),
             mai=c(.8,.9,.5,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, ftyp=ftyp,
             mtxt="SPRED4", mtxt.line=3.3,
             vcols=acols, vbrks=abrks,
             gcols=gcols)


### COMPARE STRUCTURAL MEASURES at AAS
## iupred3 vs.flDPnn
## anchor2 vs. DisoRDPbind

fname <- file.path(sfig.path,paste0("structure_cor_iupred3_flDPnn"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(bdat$iupred3, bdat$flDPnn, xlab="iupred3", ylab="flDPnn")
dev.off()

fname <- file.path(sfig.path,paste0("structure_cor_anchor2_DisoRDPbind"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(bdat$anchor2, bdat$DisoRDPbind, xlab="anchor2", ylab="DisoRDPbind",
        legpos="bottomright")
dev.off()

fname <- file.path(sfig.path,paste0("structure_cor_MMSeq2_RAAS"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(bdat$MMSeq2, bdat$RAAS,
        xlab="conservation, MMSeq2",
        ylab=xl.raas, legpos="topright")
dev.off()

fname <- file.path(sfig.path,paste0("structure_cor_MMSeq2_iupred3"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(bdat$MMSeq2, bdat$iupred3,
        xlab="conservation, MMSeq2",
        ylab="disordered score, iupred3")
dev.off()

fname <- file.path(sfig.path,paste0("structure_cor_flDPnn_RAAS"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(bdat$flDPnn, bdat$RAAS,
        xlab="disordered score, flDPnn",
        ylab=xl.raas, legpos="bottomright")
dev.off()
fname <- file.path(sfig.path,paste0("structure_cor_iupred3_RAAS"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(bdat$iupred3, bdat$RAAS,
        xlab="disordered score, iupred3",
        ylab=xl.raas, legpos="bottomleft")
dev.off()

## TODO: better analyze alpha-helical assocation
## get actual prediction score into the bp_mapped.tsv!
plotCor(bdat$RAAS, bdat$C.protein)
plotCor(bdat$RAAS, bdat$H.protein)
plotCor(bdat$RAAS, bdat$E.protein)



