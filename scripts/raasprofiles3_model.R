
library(viridis)
library(segmenTools)
options(stringsAsFactors=FALSE)
library(randomForest)

## project-specific functions
source("~/work/mistrans/scripts/saap_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source("~/work/mistrans/scripts/raasprofiles3_init.R")


mdlfig.path <- file.path(fig.path,"models")
dir.create(mdlfig.path, showWarnings=FALSE)

corW <- corH <- 2.5
pmai <- c(.5,.5,.25,.25)
pmpg <- c(1.3,.3,0)

## Andrew's code

## TODO: add

rnd <- sample(1:nrow(site),floor(nrow(site)*.8))
train <- site[rnd,]
test <- site[-rnd,]

## NOTE: r=.55 when including + protein.halflife  + protein.meltingT` but
## only for about 500 proteins
model <- randomForest(RAAS.median ~  MMSeq2 + iupred3 +  flDPnn + anchor2 + DisoRDPbind + fromto + codon + protein.length + protein.intensity + protein.halflife  + protein.meltingT,
                      data = train, na.action=na.roughfix)#na.omit)
  
# Predict speed using the trained model
predictions <- predict(model, newdata = test)

plotdev(file.path(mdlfig.path,"randomForest_all"), type=ftyp,
        width=corW, height=corH, res=300)
par(mai=pmai, mgp=pmpg, tcl=-.25) 
plotCor(predictions, test$RAAS.median, xlab="prediction", ylab=xl.raas,
        title=TRUE, cor.legend=FALSE)
abline(a=0, b=1, col=5)
dev.off()


## TODO: xgboost, glm4.

## TODO: leave-on-out loop: omit each of the variables above and calculate R2
## of predicted vs. measured RAAS.

## NOTE: highest R obtained with only protein.length
