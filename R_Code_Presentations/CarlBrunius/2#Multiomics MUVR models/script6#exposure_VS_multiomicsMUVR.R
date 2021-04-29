#All multiomics muvr models

library(MUVR)
library(doParallel)
load("selected_exposures_partcor.rda")
load("incommon_id.rda") 
load("multiOmicsXvars.rda")

exposome <- exposome[all_id,allPartCorr]
multiOmicsModels <- list()

###MUVRmodel for hs_bpa_madj_Log2
nCore=detectCores()-1
nRep=nCore
nOuter=8
varRatio=0.8

cl=makeCluster(nCore)
registerDoParallel(cl)

###MUVRmodel for hs_bpa_madj_Log2
x<- read.csv2(file="data_mergedhs_bpa_madj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_bpa_madj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[1]]<-MUVRModel


###MUVRmodel for hs_bupa_madj_Log2
x<- read.csv2(file="data_mergedhs_bupa_madj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_bupa_madj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[2]]<-MUVRModel


###MUVRmodel for hs_dde_cadj_Log2
x<- read.csv2(file="data_mergedhs_dde_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_dde_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[3]]<-MUVRModel


###MUVRmodel for hs_dep_cadj_Log2
x<- read.csv2(file="data_mergedhs_dep_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_dep_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[4]]<-MUVRModel


###MUVRmodel for h_Benzene_Log
x<- read.csv2(file="data_mergedh_Benzene_Log.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="h_Benzene_Log")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[5]]<-MUVRModel


###MUVRmodel for h_PM_Log
x<- read.csv2(file="data_mergedh_PM_Log.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="h_PM_Log")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[6]]<-MUVRModel


###MUVRmodel for hs_as_c_Log2
x<- read.csv2(file="data_mergedhs_as_c_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_as_c_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[7]]<-MUVRModel


###MUVRmodel for hs_cu_c_Log2
x<- read.csv2(file="data_mergedhs_cu_c_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_cu_c_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[8]]<-MUVRModel


###MUVRmodel for hs_ddt_cadj_Log2
x<- read.csv2(file="data_mergedhs_ddt_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_ddt_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[9]]<-MUVRModel


###MUVRmodel for hs_hcb_cadj_Log2
x<- read.csv2(file="data_mergedhs_hcb_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_hcb_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[10]]<-MUVRModel


###MUVRmodel for hs_mo_c_Log2
x<- read.csv2(file="data_mergedhs_mo_c_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_mo_c_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[11]]<-MUVRModel


###MUVRmodel for hs_pbde153_cadj_Log2
x<- read.csv2(file="data_mergedhs_pbde153_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pbde153_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[12]]<-MUVRModel


###MUVRmodel for hs_pcb138_cadj_Log2
x<- read.csv2(file="data_mergedhs_pcb138_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pcb138_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[13]]<-MUVRModel


###MUVRmodel for hs_pcb153_cadj_Log2
x<- read.csv2(file="data_mergedhs_pcb153_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pcb153_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[14]]<-MUVRModel


###MUVRmodel for hs_pcb170_cadj_Log2
x<- read.csv2(file="data_mergedhs_pcb170_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pcb170_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[15]]<-MUVRModel


###MUVRmodel for hs_pcb170_madj_Log2
x<- read.csv2(file="data_mergedhs_pcb170_madj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pcb170_madj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[16]]<-MUVRModel


###MUVRmodel for hs_pcb180_cadj_Log2
x<- read.csv2(file="data_mergedhs_pcb180_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pcb180_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[17]]<-MUVRModel


###MUVRmodel for hs_pfhxs_c_Log2
x<- read.csv2(file="data_mergedhs_pfhxs_c_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pfhxs_c_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[18]]<-MUVRModel


###MUVRmodel for hs_pfhxs_m_Log2
x<- read.csv2(file="data_mergedhs_pfhxs_m_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pfhxs_m_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[19]]<-MUVRModel


###MUVRmodel for hs_pfna_c_Log2
x<- read.csv2(file="data_mergedhs_pfna_c_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pfna_c_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[20]]<-MUVRModel


###MUVRmodel for hs_pfoa_c_Log2
x<- read.csv2(file="data_mergedhs_pfoa_c_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pfoa_c_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[21]]<-MUVRModel


###MUVRmodel for hs_pfoa_m_Log2
x<- read.csv2(file="data_mergedhs_pfoa_m_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pfoa_m_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[22]]<-MUVRModel


###MUVRmodel for hs_pfos_c_Log2
x<- read.csv2(file="data_mergedhs_pfos_c_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pfos_c_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[23]]<-MUVRModel


###MUVRmodel for hs_pm10_yr_hs_h_None
x<- read.csv2(file="data_mergedhs_pm10_yr_hs_h_None.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pm10_yr_hs_h_None")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[24]]<-MUVRModel


###MUVRmodel for hs_sumPCBs5_cadj_Log2
x<- read.csv2(file="data_mergedhs_sumPCBs5_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_sumPCBs5_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[25]]<-MUVRModel


###MUVRmodel for hs_sumPCBs5_cadj_Log2
x<- read.csv2(file="data_mergedhs_pcb118_cadj_Log2.csv", header=T, sep=";")
x$X <- as.character(x$X)
y<-exposome[,which(colnames(exposome)=="hs_pcb118_cadj_Log2")]
identical(x[,1],rownames(exposome))
x <- x[,-1]
MUVRModel<-MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
multiOmicsModels[[26]]<-MUVRModel

stopCluster(cl)


#Evaluating which MUVR models have Q2 values above 0.2 and listing them
partCorrAndQ2<-vector()
finalMultiOmicsmodels<-list()
finalMultiOmicsXVars<-list()
n<-1

for(i in 1:26){
  if(any(multiOmicsModels[[i]]$fitMetric$Q2>0.2)){
    partCorrAndQ2<-c(partCorrAndQ2,allPartCorr[i])
    finalMultiOmicsmodels[[n]]<-multiOmicsModels[[i]]
    finalMultiOmicsXvars[[n]]<-multiOmicsXVars[[i]]
    n<-n+1
  }
}

save.image("allMultiOmicsMUVRModels.Rdata")

