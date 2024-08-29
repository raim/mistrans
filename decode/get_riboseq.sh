#!/bin/bash

## DATA PATHS
DECODE=/home/raim/data/mistrans/

## utils
myBWTBG=~/programs/ucsc_utils/bigWigToBedGraph

## manual download from https://gwips.ucc.ie/downloads/index.html
## mammalian > human

## Elongating Ribosomes (A Site) > Global Aggregate >
## Iwasaki19_RiboProElong_global_(child) -> DOWNLOAD

## Elongating Ribosomes (Footprints) > Global Aggregate >
## Iwasaki16_All.RiboCov.bw -> DOWNLOAD

## mv ~/Downloads/Iwasaki19_All.RiboProElong.bw $DECODE/originalData
## mv ~/Downloads/Iwasaki16_All.RiboCov.bw $DECODE/originalData

$myBWTBG $DECODE/originalData/Iwasaki19_All.RiboProElong.bw $DECODE/processedData/Iwasaki19_All.RiboProElong.bed
$myBWTBG $DECODE/originalData/Iwasaki19_All.RiboProElong.bw $DECODE/processedData/Iwasaki19_All.RiboProElong.bed


