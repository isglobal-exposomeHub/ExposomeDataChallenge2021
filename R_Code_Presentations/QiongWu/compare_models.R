setwd("/Users/qwu/Dropbox/AwardCompetitionFellow/isglobal_2021/Submit_files")
library(R.matlab)
library(pracma)
library(caret)
expos <- readMat('expos_merge.mat')
meta <- readMat('meta_merge.mat')
res <- readMat('expos_meta_res.mat')

expos_in <- expos$expos[res$t.in,]
expos_out <- expos$expos[setdiff(1:dim(expos$expos)[1],res$t.in),]
meta_in <- meta$meta[res$s.in,]
meta_out <- meta$meta[setdiff(1:dim(meta$meta)[1],res$s.in),]

## rsquared for using the full model and only exposome within the cluster
rsquared.expos.in <- rep()
for(i in 1:dim(meta_in)[1]){
  y <- meta_in[i,]
  model <- lm(y~t(expos_in))
  model.summary <- summary(model)
  rsquared.expos.in[i] <- model.summary$r.squared
}

rsquared.expos.full <- rep()
for(i in 1:dim(meta_in)[1]){
  y <- meta_in[i,]
  model <- lm(y~t(expos$expos))
  model.summary <- summary(model)
  rsquared.expos.full[i] <- model.summary$r.squared
}

c(mean(rsquared.expos.full),sd(rsquared.expos.full))
c(mean(rsquared.expos.in),sd(rsquared.expos.in))

## difference in predictive performance 
#  check for linear model
full.res <-  matrix(0,nrow = dim(meta_in)[1],ncol = 6)
sub.res <-  matrix(0,nrow = dim(meta_in)[1],ncol = 6)
names <- c("RMSE","Rsquared","MAE","RMSESD", "RsquaredSD","MAESD")
colnames(full.res) <- names
colnames(sub.res) <- names
train_control <- trainControl(method="cv", number=5)
for(i in 1:dim(meta_in)[1]){
  dat = cbind.data.frame(meta_in[i,],t(expos$expos))
  model <- train(`meta_in[i, ]`~., data=dat, trControl=train_control, method="lm")

  dat.sub = cbind.data.frame(meta_in[i,],t(expos_in))
  model.sub <- train(`meta_in[i, ]`~., data=dat.sub, trControl=train_control, method="lm")
  full.res[i,] <- as.numeric(model$results[-1])
  sub.res[i,] <- as.numeric(model.sub$results[-1])
}

colMeans(full.res)
colMeans(sub.res)
