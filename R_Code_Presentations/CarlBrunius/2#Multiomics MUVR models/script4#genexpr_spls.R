############################################################
# Exposome Challenge 
# Methylation data sPLS data reduction strategy
############################################################
#Doing sPLS analysis on the original methylation dataset
#in order to reduce the number of variables to use for MUVR modelling

library(mixOmics)
library(GEOquery)
load("genexpr.Rdata")
load("selected_exposures_partcor.RData")

x<-t(exprs(genexpr))
x<-x[order(as.integer(rownames(x))),]
genExprData<-x

for(i in 1:length(x[1,])){
  x[,i]<-scale(x[,i])
}

#Formatting the exposome dataset to only contain samples which have gene expression data
genExprExposome<-exposome[which(exposome$ID%in%rownames(x)),]
identical(rownames(genExprData), rownames(genExprExposome))
genExprExposome<-genExprExposome[,-1]
exposomeGenExpr<-genExprExposome

#Creating a number of lists into which data will be stored in the sPLS loop

spsModelList<-list()
spsLoadingsList<-list()
whichElected<-list()

#Sparsing down gene expressions using sPLS
#Loop going through every exposome variable that made it through partial correlation (n=26)
for(i in 1:length(allPartCorr)){
  y<-genExprExposome[,which(allPartCorr[i]==colnames(genExprExposome))]
  y<-scale(y)
  spsModelList[[i]]<-spls(X=x,Y=y, ncomp=4, mode="regression", keepX=c(800,600,400,200))
  spsLoadingsList[[i]]<-spsModelList[[i]]$loadings
  
  #Loop within the loop to find out which ~1000 variables that were elected by the spls
  n<-1
  whichElected[[i]]<-vector()
  
  for(l in 1:length(spsLoadingsList[[i]]$X[,1])){
    if(sum(spsLoadingsList[[i]]$X[l,])>0){
      whichElected[[i]][n]<-rownames(spsLoadingsList[[i]]$X)[l]
      n<-n+1
    }
  }
  whichElected[[i]]<-unique(whichElected[[i]])
}

save(spsModelList, file="gene_spsModelList.rda")
save(spsLoadingsList, file="gene_spsLoadingsList.rda")
save(whichElected, file="gene_whichElected.rda")


