#!/bin/bash

## DATA PATHS
DECODE=/home/raim/data/mistrans/

## utils
myBWTBG=~/programs/ucsc_utils/bigWigToBedGraph

## manual download from https://gwips.ucc.ie/downloads/index.html
## mammalian > human
## Elongating Ribosomes (A Site) > Global Aggregate >
## Iwasaki19_RiboProElong_global_(child) -> DOWNLOAD

## mv ~/Downloads/Iwasaki19_All.RiboProElong.bw $DECODE/originalData

$myBWTBG $DECDE/originalData/Iwasaki19_All.RiboProElong.bw $DECDE/processedData/Iwasaki19_All.RiboProElong.bed

