library(MUVR)
library(doParallel)
library(mixOmics)
load("exposomeEdited.Rdata")

####Setting up Exposome (x) and Phenotype (y) data frames to use with MUVR modelling
#Exposome - All variables
x<-exposome[,-1]

#Asthma - Class
y<-as.factor(phenotype$hs_asthma)

#BMI - Class
y<-as.factor(phenotype$hs_bmi_c_cat)

#BW - Regr
y<-phenotype$e3_bw

#zbmi_who - Regr
y<-phenotype$hs_zbmi_who

#hs_correct_raven - Regr
y<-phenotype$hs_correct_raven

#hs_Gen_Tot - Regr
y<-phenotype$hs_Gen_Tot

####Setting up doParallell & MUVR parameters
nCore=detectCores()-1
nRep=nCore
nOuter=8
varRatio=0.8

cl=makeCluster(nCore)
registerDoParallel(cl)

####Factor outcome models with randomized sample selection to avoid falsely low misclassification rates due to uneven group belonging of samples
####Carried out for hs_astmha

y<-as.factor(phenotype$hs_asthma)
yZero<-which(y==0)
yOne<-which(y==1)
classModel_Astma<-list()

for(i in 1:10){
  yRand<-sample(yZero,size=142)
  index<-c(yRand,yOne)

  classModel_Astma[[i]]=MUVR(X=x[index,],Y=y[index], fitness = "MISS", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
}

fitness<-data.frame()
for(i in 1:10){
  fitness[i,c(1:3)]<-classModel_Astma[[i]]$miss
}

####Carried out for hs_bmi_c
y<-as.factor(phenotype$hs_bmi_c_cat)
xTemp<-x[-which(y==1),]
y<-factor(y[-which(y==1)])
y2<-which(y==2)
y3<-which(y==3)
y4<-which(y==4)
classModel_BMICat<-list()

for(i in 1:10){
  yRand2<-sample(y2,size=131)
  yRand3<-sample(y3,size=131)
  index<-c(yRand2,yRand3, y4)

  classModel_BMICat[[i]]=MUVR(X=xTemp[index,],Y=y[index], fitness = "MISS", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
}

fitness<-data.frame()
for(i in 1:10){
  fitness[i,c(1:3)]<-classModel_BMICat[[i]]$miss
}

####All random permutations of the class models were evaluated manually and none of them had Q2 > 0.2

####Regression models don't have to be rerun several times###
regrMUVRModels<-list()

y<-as.factor(phenotype$hs_bmi_c_cat)
regrModel_BMI=MUVR(X=x,Y=y, fitness = "BER", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
regrMUVRModels[[1]]<-regrModel_BMI

y<-phenotype$e3_bw
regrModel_BW=MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
regrMUVRModels[[2]]<-regrModel_BW

y<-phenotype$hs_zbmi_who
regrModel_who=MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
regrMUVRModels[[3]]<-regrModel_BW

y<-phenotype$hs_correct_raven
regrModel_raven=MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
regrMUVRModels[[4]]<-regrModel_raven

y<-phenotype$hs_Gen_Tot
regrModel_GenTot=MUVR(X=x,Y=y, fitness = "RMSEP", method="RF", nRep=nRep, nOuter=nOuter, varRatio=varRatio)
regrMUVRModels[[5]]<-regrModel_GenTot

stopCluster(cl)

#Checking if any of the regression models had a Q2 > 0.2
#If they do the VIP-list is extracted and added to a list of vectors
n<-1
VIPlists<-list()
for(i in 1:5){
  if(any(regrMUVRModels[[i]]$fitMetric$Q2>0.2)){
    VIPlists[[n]]<-getVIP(regrMUVRModels[[i]], model="max")
    n<-n+1
  }
}

save.image("exposure_VS_phenotype_models.Rdata")
