library(segmenTools)

## after Andrew Leduc's code in decode/decode_analysis/predict_raas.R

library(dplyr)
library(ggplot2)
library(xgboost)
options(stringsAsFactors=FALSE) # required for R<4

proj.path <- file.path(Sys.getenv("MISDATA")) 
if ( proj.path=="" ) {
  proj.path <- "/home/raim/data/mistrans" # author's local path
} 

fig.path <- file.path(proj.path,"figures","model")
dir.create(fig.path)

ftyp <- "png"

RAAS_data = read.delim(file.path(proj.path,'processedData',
                                 'sites_raas_unique.tsv'),
                       row.names = NULL)
RAAS_data = RAAS_data[RAAS_data$RAAS.n > 1,]

test_rows = sample(1:nrow(RAAS_data),round(nrow(RAAS_data)*.8))
train = RAAS_data[test_rows,]
test = RAAS_data[-test_rows,]

## list of features to use for random forest prediction
feat_list =  c('MMSeq2',
               ##'fromto', 'codon',
               ##'flDPnn',
               ##'DisoRDPbind',
               'iupred3',
               ##'anchor2',
               'protein.intensity',
               'protein.halflife',
               'protein.length'
               )


train_features = train %>% select(feat_list) 
test_feature = test  %>% select(feat_list) 

train_matrix <- xgb.DMatrix(data = as.matrix(train_features),
                            label = as.numeric(train$RAAS.mean))

# Set parameters for the model
params <- list(objective = "reg:squarederror", booster="gbtree",
               eval_metric="rmse")

# Train the model
model <- xgb.train(params, train_matrix, nrounds = 100)


predict_mat <- xgb.DMatrix(data = as.matrix(test_feature))

# Make predictions
predictions <- predict(model, newdata = predict_mat)

if ( interactive() ) {
    plotCor(predictions,test$RAAS.mean)

    cor(predictions,test$RAAS.mean,method = 'spearman')
}


cor_list <- c()
feat_store <- c()

for(i in feat_list){
  
  for(j in 1:100){
    
    test_rows = sample(1:nrow(RAAS_data),round(nrow(RAAS_data)*.8))
    train = RAAS_data[test_rows,]
    test = RAAS_data[-test_rows,]
    train_features = train %>% select(feat_list) 
    test_feature = test  %>% select(feat_list) 
    
    train_features <- train_features %>% select(-any_of(i))
    test_feature <- test_feature %>% select(-any_of(i))
    
    train_matrix <- xgb.DMatrix(data = as.matrix(train_features),
                                label = as.numeric(train$RAAS.mean))
    params <- list(objective = "reg:squarederror",
                   booster="gbtree", eval_metric="rmse")
    model <- xgb.train(params, train_matrix, nrounds = 100)
    predict_mat <- xgb.DMatrix(data = as.matrix(test_feature))
    predictions <- predict(model, newdata = predict_mat)
    
    cor_list <- c(cor_list,cor(predictions,test$RAAS.mean,method = 'spearman'))
    feat_store <- c(feat_store,i)
    
  }
  
}
cor_list_all <- c()
for(j in 1:100){
  
  test_rows = sample(1:nrow(RAAS_data),round(nrow(RAAS_data)*.8))
  train = RAAS_data[test_rows,]
  test = RAAS_data[-test_rows,]
  train_features = train %>% select(feat_list) 
  test_feature = test  %>% select(feat_list) 
  
  
  train_matrix <- xgb.DMatrix(data = as.matrix(train_features),
                              label = as.numeric(train$RAAS.mean))
  params <- list(objective = "reg:squarederror",
                 booster="gbtree", eval_metric="rmse")
  model <- xgb.train(params, train_matrix, nrounds = 100)
  predict_mat <- xgb.DMatrix(data = as.matrix(test_feature))
  predictions <- predict(model, newdata = predict_mat)
  
  cor_list_all <- c(cor_list_all,cor(predictions,
                                     test$RAAS.mean, method = 'spearman'))
  ##cor(predictions,test$RAAS.mean,method = 'spearman')
  
}


df_xg <- as.data.frame(cbind(cor_list, feat_store))
df_xg$cor_list <- (as.numeric(df_xg$cor_list) - median(cor_list_all)) /median(cor_list_all)*100
df_xg$cor_list <- abs(df_xg$cor_list)
df_xg$feat_store <- with(df_xg, reorder(feat_store, cor_list, median))

df_plot <- as.data.frame(predictions)
df_plot$gt <- test$RAAS.mean

p1 <- ggplot(df_plot,aes(y = predictions, x =gt )) + geom_point() + 
  xlab('Measured') + ylab('Predicted') +ggtitle('log10 RAAS')

p2 <- ggplot(df_xg,aes(y = feat_store, x = (cor_list))) + geom_boxplot() + 
    xlab('% gain of feature') + coord_cartesian(xlim = c(0,max(df_xg$cor_list))) +
    ylab('')


segmenTools::plotdev(file.path(proj.path,'figures','model',
                               'xgboost_correlation'),
                type=ftyp, height=2.5, width=2.5, res=200)
par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
segmenTools::plotCor(test$RAAS.mean, predictions, 
                     ylab="predicted", xlab="observed",
                     line.methods="ols", cor.legend=FALSE, title=TRUE)
dev.off()

write.csv(df_xg,file.path(proj.path,'processedData/xgboost_leave_one_out.csv'))
write.csv(df_plot,file.path(proj.path,'processedData/xgboost_plot.csv'))

ggsave(filename = file.path(proj.path,'figures/model/xgboost_scatter.png'),
       plot = p1, width = 3, height = 3, dpi = 300)        
ggsave(filename = file.path(proj.path,'figures/model/xgboost_features.png'),
       plot = p2, width = 4, height = 3, dpi = 300)        

          
summary(cor_list_all)
        
## PCA
docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(50)
upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(50)
ttcols <- unique(c(rev(docols), upcols))

cra <- cor(RAAS_data[,c("RAAS.median",feat_list)], use="pairwise.complete")

segmenTools::plotdev(file.path(proj.path,'figures','model',
                               'correlation'),
                type=ftyp, height=3.5, width=3.5, res=200)
par(mai=c(1.2,1.2,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
image_matrix(cra, breaks=seq(-1,1, length.out=100), col=ttcols, axis=1:2,
             xlab=NA, ylab=NA)
dev.off()
