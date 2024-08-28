library(segmenTools)

## attempt to get the correct path for input data;
## these should be downloaded or generated as outlined
## in README.md
proj.path <- file.path(Sys.getenv("DECODE")) 
if ( proj.path=="" ) # author's local path
    proj.path <- "/home/raim/data/mistrans"
SRC.PATH <- file.path("/home/raim/work/mistrans/decode/")


ribo.file <- file.path(proj.path, "processedData",
                       "Iwasaki19_All.RiboProElong.bed")


source("~/programs/segmenTools/R/parsers.R")
ribo <- bed2coor(ribo.file, header = c("chr", "start", "end", "score"))

## strand info by reverted start/end? would be illegal for bed, i think
any(ribo$end<ribo$start)

## fraction of multiple position values
sum(ribo$end==ribo$start)/nrow(ribo) #93% are single position values
sum(ribo$end>ribo$start)/nrow(ribo) # 6% cover several positions


## EXPAND RANGES TO FULL POSITIONS
lens <- ribo$end-ribo$start+1

## exponential decrease of covered positions, mostly 2
barplot(table(ribo$end-ribo$start+1))

## expand positions
ribo <- data.frame(chr=rep(ribo$chr, lens),
                   pos=rep(ribo$start, lens),
                   score=rep(ribo$score, lens))

head(ribo)

## number of non-consecutive runs
sum(diff(ribo$start)!=1)
