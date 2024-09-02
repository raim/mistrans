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
## mv ~/Downloads/human_eIF3b_bound_40S.RiboProElong.bw $DECODE/originalData
## mv ~/Downloads/Chen20_All.RiboCov.bw $DECODE/originalData

$myBWTBG $DECODE/originalData/Iwasaki19_All.RiboProElong.bw $DECODE/processedData/Iwasaki19_All.RiboProElong.bed
$myBWTBG $DECODE/originalData/Iwasaki16_All.RiboCov.bw $DECODE/processedData/Iwasaki16_All.RiboCov.bed
$myBWTBG $DECODE/originalData/human_eIF3b_bound_40S.RiboProElong.bw $DECODE/processedData/human_eIF3b_bound_40S.RiboProElong.bed
$myBWTBG $DECODE/originalData/Chen20_All.RiboCov.bw $DECODE/processedData/Chen20_All.RiboCov.bed

## save disk space 
gzip $DECODE/processedData/Iwasaki19_All.RiboProElong.bed
gzip $DECODE/processedData/Iwasaki16_All.RiboCov.bed
gzip $DECODE/processedData/human_eIF3b_bound_40S.RiboProElong.bed
gzip $DECODE/processedData/Chen20_All.RiboCov.bed


