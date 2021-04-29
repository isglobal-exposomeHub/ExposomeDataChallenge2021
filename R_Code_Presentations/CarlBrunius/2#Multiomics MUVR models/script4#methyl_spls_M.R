############################################################
# Exposome Challenge 
# Methylation data sPLS data reduction strategy
############################################################
#Doing sPLS analysis on the original gene expression dataset
#in order to reduce the number of variables to use for MUVR modelling

library(mixOmics)
library(minfi)
load("methy.Rdata")
load("selected_exposures_partcor.RData")

methy<- "getM"(methy)
x<-t(methy)
x<-x[order(as.integer(rownames(x))),]
methyData<-x

for(i in 1:length(x[1,])){
  x[,i]<-scale(x[,i])
}

#Formatting the exposome dataset to only contain samples which have gene expression data
methyExposome<-exposome[which(exposome$ID%in%rownames(x)),]
identical(rownames(methyData), rownames(methyExposome))
methyExposome<-methyExposome[,-1]
exposomemethy<-methyExposome

#Creating a number of lists into which data will be stored in the sPLS loop

spsModelList<-list()
spsLoadingsList<-list()
whichElected<-list()

#Sparsing down gene expressions using sPLS
#Loop going through every exposome variable that made it through partial correlation (n=26)
for(i in 1:length(allPartCorr)){
  y<-methyExposome[,which(allPartCorr[i]==colnames(methyExposome))]
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

save(spsModelList, file="methy_spsModelList_M.rda")
save(spsLoadingsList, file="methy_spsLoadingsList_M.rda")
save(whichElected, file="methy_whichElected_M.rda")


