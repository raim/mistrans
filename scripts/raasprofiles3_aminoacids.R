
### ANALYZE AMINO ACID AND PROTEIN STRUCTURAL PROPERTIES at AAS

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
 



## RAAS profiles by relative position of AAS in peptide
rssrt <- levels(tmtf$rsite.bins)
fname <- file.path(dpath,paste0("AASsite_",SETID,"_wtests_"))
ovw  <- raasProfile(x=tmtf, id="SAAP", value="RAAS",
                    delog=TRUE, rows="rsite.bins", cols="Dataset",
                    bg=TRUE, row.srt=rev(rssrt),
                    col.srt=uds,
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
             vcols=acols, vbrks=abrks,
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
             vcols=acols, vbrks=abrks,
             gcols=gcols, verb=1)
source("~/work/mistrans/scripts/saap_utils.R")
if ( nrow(ovwp$p.value)>0 )
    plotProfiles(ovwp,
                 fname=file.path(afig.path,paste0("AAprop_",SETID,"_cut")),
                 mai=c(.7,1.5,.5,.5), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="", ftyp=ftyp,
                 vcols=acols, vbrks=abrks,
                 plot.legend=TRUE,
                 gcols=gcols)


### FULL MANUAL DOTPLOT BY AAprop for MANUSCRIP MAIN
## calculate optimal figure height: result fields + figure margins (mai)
nr <- nrow(ovw$p.value)
nc <- ncol(ovw$p.value)
mai <- c(.8,1.75,.7,.6)
nh <- nr *fh + mai[1] + mai[3]
nw <- nc *fw + mai[2] + mai[4]
ffam <- "monospace"#"sans"

if ( SETID=="tissues" ) ffam <- "sans"

rsrt <- rownames(ovw$p.value)
tap <- sub(".*:","", rsrt)
fap <- sub(":.*","", rsrt)

### NOTE: using tighter color range
fname <- file.path(afig.path,paste0("AAprop_",SETID))
## combined effect size and p-value plot
plotdev(paste0(fname,"_dotplot_manual"),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
dotprofile(x=ovw, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
box()
axis(1,  1:ncol(ovw$p.value),
     labels=colnames(ovw$p.value), las=2, family=ffam)
##axis(2, length(axex):1, labels=axex, las=2, family=ffam)
shadowtext(rep(-6, nr), nr:1, fap, col=aaprop.cols[fap], xpd=TRUE,
           font=2, r=.1)
text(rep(-3.75,nr), nr:1, expression(""%->%""), xpd=TRUE)
shadowtext(rep(-1.5, nr), nr:1, tap, col=aaprop.cols[tap], xpd=TRUE,
           font=2, r=.1)
axis(3, at=1:ncol(ovw$num.target),
     labels=format(ovw$num.target[1,], big.mark=",", trim=TRUE),las=2)
axis(4, at=nrow(ovw$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
     labels=format(ovw$num.query[,1], big.mark=",", trim=TRUE),las=2)
dev.off()

## TIGHT MANUAL DOTPLOT BY AAprop for MANUSCRIP MAIN
# calculate optimal figure height: result fields + figure margins (mai)
nr <- nrow(ovwp$p.value)
nc <- ncol(ovwp$p.value)
mai <- c(.8,1.75,.1,.6)
nh <- nr *fh + mai[1] + mai[3]
nw <- nc *fw + mai[2] + mai[4]

rsrt <- rownames(ovwp$p.value)
tap <- sub(".*:","", rsrt)
fap <- sub(":.*","", rsrt)

### NOTE: using tighter color range
fname <- file.path(afig.path,paste0("AAprop_",SETID))
## combined effect size and p-value plot
plotdev(paste0(fname,"_cut_dotplot_manual"),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
dotprofile(x=ovwp, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
box()
axis(1,  1:ncol(ovwp$p.value),
     labels=colnames(ovwp$p.value), las=2, family=ffam)
##axis(2, length(axex):1, labels=axex, las=2, family=ffam)
shadowtext(rep(-6, nr), nr:1, fap, col=aaprop.cols[fap], xpd=TRUE,
           font=2, r=.1)
text(rep(-3.75,nr), nr:1, expression(""%->%""), xpd=TRUE)
shadowtext(rep(-1.5, nr), nr:1, tap, col=aaprop.cols[tap], xpd=TRUE,
           font=2, r=.1)
##axis(3, at=1:ncol(ovwp$num.target),
##     labels=format(ovwp$num.target[1,], big.mark=",", trim=TRUE),las=2)
axis(4, at=nrow(ovwp$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
     labels=format(ovwp$num.query[,1], big.mark=",", trim=TRUE),las=2)
dev.off()


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
## SORT ROWS BY MEDIAN RAAS
nsrt <- names(sort(apply(ovw$median, 1, median)))
ovw <- sortOverlaps(ovw, axis=2, srt=nsrt)
plotProfiles(ovw, fname=file.path(afig.path,paste0("fromAA_",SETID)),
             mtxt="Encoded AA",
             mai=c(.8,.5,.05,.6),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             ttcols=ttcols, value="median",
             rlab=LAB, llab="", ftyp=ftyp,
             ##axis2.col=aap.cols,
             vcols=acols, vbrks=abrks,
             gcols=gcols, plot.all=TRUE, ffam=ffam)

## from AA - ONLY SIGNIFICANT

ovwp <- sortOverlaps(ovw, axis=2, p.min=p.min, cut=TRUE)
nsrt <- names(sort(apply(ovwp$median, 1, median)))
ovwp <- sortOverlaps(ovwp, axis=2, srt=nsrt)

nr <- nrow(ovwp$p.value)
nc <- ncol(ovwp$p.value)
mai <- c(.8,.5,.1,.6)
nh <- nr *fh + mai[1] + mai[3]
nw <- nc *fw + mai[2] + mai[4]
fap <- rownames(ovwp$p.value)

fname <- file.path(afig.path,paste0("fromAA_",SETID))
plotdev(paste0(fname,"_cut_dotplot_manual"),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
dotprofile(x=ovwp, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
box()
axis(1,  1:ncol(ovwp$p.value),
     labels=colnames(ovwp$p.value), las=2, family=ffam)
shadowtext(rep(-.1, nr), nr:1, fap, col=aap.cols[fap], xpd=TRUE,
           font=2, r=.1)
mtext("Encoded AA", 2, 1.3, family=ffam)
##axis(3, at=1:ncol(ovwp$num.target),
##     labels=format(ovwp$num.target[1,], big.mark=",", trim=TRUE),las=2)
axis(4, at=nrow(ovwp$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
     labels=format(ovwp$num.query[,1], big.mark=",", trim=TRUE),las=2)
dev.off()

## ROTATE
ovwr <- t(ovwp)
nr <- nrow(ovwr$p.value)
nc <- ncol(ovwr$p.value)
mai <- c(.4,.8,.6,.1)
nh <- nr *fh + mai[1] + mai[3]
nw <- nc *fw + mai[2] + mai[4]
fap <- colnames(ovwr$p.value)

fname <- file.path(afig.path,paste0("fromAA_",SETID))
plotdev(paste0(fname,"_cut_dotplot_manual_rotated"),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
dotprofile(x=ovwr, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
box()
axis(2, nrow(ovwr$p.value):1,
     labels=rownames(ovwr$p.value), las=2, family=ffam)
shadowtext(1:nc, rep(0, nc), fap, col=aap.cols[fap], xpd=TRUE,
           font=2, r=.1)
mtext("Encoded AA", 1, 1, family=ffam)
axis(3, at=1:ncol(ovwr$num.target),
     labels=format(ovwr$num.target[1,], big.mark=",", trim=TRUE),las=2)
##axis(4, at=nrow(ovwr$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
##     labels=format(ovwr$num.query[,1], big.mark=",", trim=TRUE),las=2)
dev.off()


## BY "TO" AMINO ACID
ovw <- raasProfile(x=tmtf, id="SAAP", 
                   rows="to", cols="Dataset",
                   bg=TRUE, value="RAAS", col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0)


## SORT ROWS BY MEDIAN RAAS
nsrt <- names(sort(apply(ovw$median, 1, median)))
ovw <- sortOverlaps(ovw, axis=2, srt=nsrt)
plotProfiles(ovw, fname=file.path(afig.path,paste0("toAA_",SETID)),
             mai=c(.8,.5,.05,.6),
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             mtxt="Incorporated AA", ttcols=ttcols, value="median",
             rlab=LAB, llab="", ftyp=ftyp,
             ##axis2.col=aap.cols,
             vcols=acols, vbrks=abrks,
             gcols=gcols, plot.all=TRUE, ffam=ffam)

## to AA - ONLY SIGNIFICANT

ovwp <- sortOverlaps(ovw, axis=2, p.min=p.min, cut=TRUE)
nsrt <- names(sort(apply(ovwp$median, 1, median)))
ovwp <- sortOverlaps(ovwp, axis=2, srt=nsrt)

nr <- nrow(ovwp$p.value)
nc <- ncol(ovwp$p.value)
mai <- c(.8,.5,.1,.6)
nh <- nr *fh + mai[1] + mai[3]
nw <- nc *fw + mai[2] + mai[4]
fap <- rownames(ovwp$p.value)

fname <- file.path(afig.path,paste0("toAA_",SETID))
plotdev(paste0(fname,"_cut_dotplot_manual"),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
dotprofile(x=ovwp, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
box()
axis(1,  1:ncol(ovwp$p.value),
     labels=colnames(ovwp$p.value), las=2, family=ffam)
shadowtext(rep(-.1, nr), nr:1, fap, col=aap.cols[fap], xpd=TRUE,
           font=2, r=.1)
mtext("Incorporated AA", 2, 1.3, family=ffam, adj=1)
##axis(3, at=1:ncol(ovwp$num.target),
##     labels=format(ovwp$num.target[1,], big.mark=",", trim=TRUE),las=2)
axis(4, at=nrow(ovwp$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
     labels=format(ovwp$num.query[,1], big.mark=",", trim=TRUE),las=2)
dev.off()

## ROTATE
ovwr <- t(ovwp)
nr <- nrow(ovwr$p.value)
nc <- ncol(ovwr$p.value)
mai <- c(.4,.8,.6,.1)
nh <- nr *fh + mai[1] + mai[3]
nw <- nc *fw + mai[2] + mai[4]
fap <- colnames(ovwr$p.value)

fname <- file.path(afig.path,paste0("toAA_",SETID))
plotdev(paste0(fname,"_cut_dotplot_manual_rotated"),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
dotprofile(x=ovwr, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
box()
axis(2, nrow(ovwr$p.value):1,
     labels=rownames(ovwr$p.value), las=2, family=ffam)
shadowtext(1:nc, rep(0, nc), fap, col=aap.cols[fap], xpd=TRUE,
           font=2, r=.1)
mtext("Incorporated AA", 1, 1, family=ffam, adj=1)
axis(3, at=1:ncol(ovwr$num.target),
     labels=format(ovwr$num.target[1,], big.mark=",", trim=TRUE),las=2)
##axis(4, at=nrow(ovwr$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
##     labels=format(ovwr$num.query[,1], big.mark=",", trim=TRUE),las=2)
dev.off()

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
                 vcols=acols, vbrks=abrks,
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
                     vcols=acols, vbrks=abrks,
                     gcols=gcols)
    
}

## all ->G
tog <- grep(":G", rownames(ovw$p.value), value=TRUE)
ovwp <- sortOverlaps(ovw, srt=tog)

mai <- c(.8,.8,.7,.6)
if ( SETID=="tissues" ) {
    mai[1] <- 1.5
}
## calculate optimal figure height: result fields + figure margins (mai)
nr <- nrow(ovwp$p.value)
nc <- ncol(ovwp$p.value)
nh <- nr *fh + mai[1] + mai[3]
nw <- nc *fw + mai[2] + mai[4]

plotProfiles(ovwp,
             fname=file.path(afig.path,paste0("AA_",SETID,"_toG")),
             mai=mai, fw=.2, fh=.2,
             ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, llab="",  ftyp=ftyp,
             ffam=ffam,
             ##axis2.col=ft.cols,
             vcols=acols, vbrks=abrks,
             gcols=gcols, plot.all=TRUE)
fname <- file.path(afig.path,paste0("AA_",SETID,"_toG_manual"))
## combined effect size and p-value plot
plotdev(paste0(fname),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
dotprofile(x=ovwp, value="median", vbrks=abrks,
           vcols=acols, p.dot=p.dot,
           dot.sze=dot.sze, axis=1:2, xlab=NA, ylab=NA)
box()
axis(3, at=1:ncol(ovwp$num.target),
     labels=format(ovwp$num.target[1,], big.mark=",", trim=TRUE),las=2)
axis(4, at=nrow(ovwp$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
     labels=format(ovwp$num.query[,1], big.mark=",", trim=TRUE),las=2)
dev.off()




## only significant
ovwp <-ovw
ovwp <- sortOverlaps2(ovwp, axis=2, p.min=p.min, cut=TRUE)


if ( nrow(ovwp$p.value)>0 ) {

    ## sort plot
    ## by from AA and internally by to AA
    oldsrt <- rownames(ovwp$p.value)
    taa <- sub(".*:","",oldsrt)
    faa <- sub(":.*","",oldsrt)
    laa <- split(taa, faa)
    laa <- lapply(laa, function(x) x[order(match(x, aap.srt))])
    laa <- laa[order(match(names(laa),aap.srt))]
    laa <- unlist(laa)
    newsrt <- gsub("[0-9]","",paste0(names(laa),":",laa))
    ovwp <- sortOverlaps2(ovwp, axis=2, srt=newsrt)
    
    
    
    faa <- sub(":.*","",rownames(ovwp$p.value))
    ft.cols <- aa.cols[faa]
    names(ft.cols) <- rownames(ovwp$p.value)
    
    mai <- c(.8,.8,.7,.6)
    if ( SETID=="tissues" ) {
        mai[1] <- 1.5
    }
    plotProfiles(ovwp,
                 fname=file.path(afig.path,paste0("AA_",SETID,"_cut")),
                 mai=mai, fw=.2, fh=.2,
                 ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab="",  ftyp=ftyp,
                 ffam=ffam,
                 axis2.col=ft.cols,
                 vcols=acols, vbrks=abrks,
                 gcols=gcols, plot.all=TRUE)

    ## calculate optimal figure height: result fields + figure margins (mai)
    nr <- nrow(ovwp$p.value)
    nc <- ncol(ovwp$p.value)
    nh <- nr *fh + mai[1] + mai[3]
    nw <- nc *fw + mai[2] + mai[4]

    rsrt <- rownames(ovwp$p.value)
    taa <- sub(".*:","", rsrt)
    faa <- sub(":.*","", rsrt)

    fname <- file.path(afig.path,paste0("AA_",SETID,"_cut_dotplot_manual"))
    ## combined effect size and p-value plot
    plotdev(paste0(fname),
            height=nh, width=nw, res=300, type=ftyp)
    par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
    dotprofile(x=ovwp, value="median", vbrks=abrks,
               vcols=acols, p.dot=p.dot,
               dot.sze=dot.sze, axis=NA, xlab=NA, ylab=NA)
    box()
    axis(1,  1:ncol(ovwp$p.value),
         labels=colnames(ovwp$p.value), las=2, family=ffam)
    ##axis(2, length(axex):1, labels=axex, las=2, family=ffam)
    shadowtext(rep(-1.6, nr), nr:1, faa, col=aaprop.cols[AAPROP[faa]], xpd=TRUE,
               font=2, r=.1)
    text(rep(-.8,nr), nr:1, expression(""%->%""), xpd=TRUE)
    shadowtext(rep( 0, nr), nr:1, taa, col=aaprop.cols[AAPROP[taa]], xpd=TRUE,
               font=2, r=.1)
    axis(3, at=1:ncol(ovwp$num.target),
         labels=format(ovwp$num.target[1,], big.mark=",", trim=TRUE),las=2)
    axis(4, at=nrow(ovwp$num.query):1, mgp=c(1.3,.1,0), tcl=-.1,
         labels=format(ovwp$num.query[,1], big.mark=",", trim=TRUE),las=2)

    ## properties
    text(x=rep(-2.75, length(aaprop.srt)), y=c(29.5, 21.5, 13.5, 5.5), srt=90,
               labels=aaprop.srt, col=aaprop.cols[aaprop.srt], xpd=TRUE, font=2)
    
    dev.off()
 }


## TODO: plot on p-value correction
##plot(ovw$p.value, qvalue::qvalue(c(ovw$p.value))$qvalues)


## by AA->AA
source("~/work/mistrans/scripts/saap_utils.R")
for ( ds in auds ) {

    tmtd <- tmtf
    dsl <- ""
    if ( ds!="all" ) {
        dsl <- ds
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
                 mai=c(.5,.5,.6,.6), ttcols=ttcols, value="median",
                 p.min=p.min, p.txt=p.txt,
                 dot.sze=dot.sze, p.dot=p.dot,
                 rlab=LAB, llab=dsl, ftyp=ftyp,
                 vcols=acols, vbrks=abrks,
                 mtxt="Incorporated AA", 
                 mtxt1="Encoded AA", ffam=ffam,
                 axis1.las=1, fw=.2, fh=.2,
                 gcols=gcols)
}



## phypergeo of AAS types vs. Dataset/Tissues

ttmt <- tmtf[tmtf$Dataset!="Cancer",]

ovl <- clusterCluster(ttmt$fromto, ttmt$Dataset, cl2.srt=uds[uds!="Cancer"])
ovc <- sortOverlaps(ovl, p.min=p.min, cut=TRUE, axis=2)

mai <- c(.8,.8,.7,.6)
if ( SETID=="tissues" ) {
    mai[1] <- 1.5
}
## calculate optimal figure height: result fields + figure margins (mai)
nr <- nrow(ovc$p.value)
nc <- ncol(ovc$p.value)
nh <- nr *.2 + mai[1] + mai[3]
nw <- nc *.25 + mai[2] + mai[4]

plotdev(file.path(afig.path,paste0("AAStype_",SETID,"_hypergeo")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=mai, mgp=c(1.3,.3,0), tcl=-.25, family=ffam)
plotOverlaps(ovc, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab=NA, ylab=NA, text.cex = .9)
dev.off()
