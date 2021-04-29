############################################################
# Exposome Challenge 
# Omics Modules - Selected exposome data - Partial Correlation
############################################################

#Use MUVR & partial correlation selected exposures 
#Use exposures with omics MUVR models q2>0.2
#Put exposures and selected omics (from MUVR) in partial correlation 
#Adjust for age, sex and cohort 
#Select relevant omics variables 

#*****************************
#Partial correlation function 
#*****************************####

partcor = function(metabolites, exposure, confounders){
  p.mat<- matrix(NA, ncol(metabolites), ncol(exposure))
  corr.mat<- matrix(NA, ncol(metabolites), ncol(exposure))
  for (i in 1:ncol(metabolites)) {
    X=metabolites[, i]
    for (j in 1:ncol(exposure)) {
      tryCatch({
        Y=exposure[, j]
        fmx=lm(X~as.factor(h_cohort)+as.factor(e3_sex_None) + hs_child_age_None, data=confounders)
        fmy=lm(Y~as.factor(h_cohort)+as.factor(e3_sex_None) + hs_child_age_None, data=confounders)
        res1=resid(fmx)
        res2=resid(fmy)
        tmp <- cor.test(res1, res2, method='spearman',use='complete.obs')
        p.mat[i, j] <- tmp$p.value
        corr.mat[i, j] <- tmp$estimate
      }
      , error=function(e) {cat("ERROR :",conditionMessage(e), "\n")} )
    }}
  
  padj.mat<- matrix(NA, ncol(metabolites), ncol(exposure))
  for (i in 1:ncol(p.mat)) {
    temp=p.mat[,i]
    padj.mat[, i]=p.adjust(temp,method='fdr')
  }
  rownames(corr.mat)=colnames(metabolites)
  colnames(corr.mat)=colnames(exposure)
  rownames(padj.mat)=colnames(metabolites)
  colnames(padj.mat)=colnames(exposure)
  return(list(corr.mat, padj.mat))
}

#*****************************
#Prepare data  
#*****************************####

library(MUVR)
load("allMultiOmicsMUVRModels.Rdata")
load("exposomeEdited.Rdata")

MOBioData<-vector()
for(i in 1:length(finalMultiOmicsmodels)){
  MultiOmics=data.frame(finalMultiOmicsXVars[[i]][,intersect(colnames(finalMultiOmicsXVars[[i]]),getVIP(finalMultiOmicsmodels[[i]])$name)])
  MultiOmics$ID=as.numeric(finalMultiOmicsXVars[[i]][,1]) #rownames is ID number
  
  covariates2=merge(MultiOmics, covariates, "ID", "ID")
  exposome2=merge(MultiOmics, exposome, "ID", "ID")
  covariates2=covariates2[,-c(2:length(MultiOmics[1,]))] #remove columns of metabolites
  exposome2=exposome2[,-c(2:length(MultiOmics[1,]))]

  #sanity check identical datasets
  writeLines(paste("Identical:",identical(MultiOmics$ID, covariates2$ID, exposome2$ID)))

  #********************
  #EXPOSURES 
  #********************####

  #Loop over exposures and metabolites from MUVR model (q2>0.2 ) 
    selexp = as.data.frame(exposome2[,colnames(exposome2) %in% partCorrAndQ2])  #Selected exposures from MUVR
    metbol = as.data.frame(MultiOmics[,colnames(MultiOmics) %in% getVIP(finalMultiOmicsmodels[[i]])$name])  #Selected metabolites from MUVR  
    conf = covariates2[, c("h_cohort", "e3_sex_None", "hs_child_age_None")]  #Confounders: Age, Sex, Cohort
  
  #Do partial correlation 
  cor = as.data.frame(partcor(metbol, selexp, conf)[1])
  pval = as.data.frame(partcor(metbol, selexp, conf)[2])

  #Extract significant exposures (FDR >0.05)
  signmetbol = ifelse(pval[,1]<0.05, rownames(pval), NA)
  assign(paste("signmetbol", partCorrAndQ2, signmetbol))
  write.csv(signmetbol, file=paste("Partcor_smetbol",i,".csv",sep=""))
  MOBioData<-c(MOBioData,signmetbol)
}

MOBioData<-unique(MOBioData[!is.na(MOBioData)])

save.image("PartCorrOmicsVars.RData")
