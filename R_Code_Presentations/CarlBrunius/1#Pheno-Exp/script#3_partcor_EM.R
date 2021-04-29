############################################################
# Exposome Challenge 
# Exposure Modules - Selected exposome data - Principal Component Analysis
############################################################

#Use MUVR selected exposures (RAVEN/BMI_selexp.csv)
#Put selected exposures in partial correlation with phenotype 
#Adjust for age, sex and cohort 
#Select relevant exposures for BMI and IQ (partcor_RAVEN/BMI.csv)

load("exposomeEdited.RData")
load("exposure_VS_phenotype_models.Rdata")

#The two models which had Q2s > 0.2: who_zBMI & raven
eBMI=VIPlists[[1]]
eIQ=VIPlists[[2]]

#sanity check identical datasets
identical(phenotype$ID, covariates$ID, exposome$ID)

#*****************************
#Partial correlation function 
#*****************************####

partcor = function(exposures, phenotype, confounders){
p.mat<- matrix(NA, ncol(exposures), ncol(phenotype))
corr.mat<- matrix(NA, ncol(exposures), ncol(phenotype))
for (i in 1:ncol(exposures)) {
  X=exposures[, i]
  for (j in 1:ncol(phenotype)) {
    tryCatch({
      Y=phenotype[, j]
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

padj.mat<- matrix(NA, ncol(exposures), ncol(phenotype))
for (i in 1:ncol(p.mat)) {
  temp=p.mat[,i]
  padj.mat[, i]=p.adjust(temp,method='fdr')
}
rownames(corr.mat)=colnames(exposures)
colnames(corr.mat)=colnames(phenotype)
rownames(padj.mat)=colnames(exposures)
colnames(padj.mat)=colnames(phenotype)
return(list(corr.mat, padj.mat))
}

#********************
#BMI 
#********************####

#Get datasets
selexp = exposome[,colnames(exposome) %in% eBMI[,1]]  #Selected exposures from MUVR 
conf = covariates[, c("h_cohort", "e3_sex_None", "hs_child_age_None")]  #Confounders: Age, Sex, Cohort, others?
pt = as.data.frame(phenotype[,"hs_zbmi_who"])  #Phenotype: BMI 

#Do partial correlation 
partcor(selexp, pt, conf)
bmicor = as.data.frame(partcor(selexp, pt, conf)[1])
bmipval = as.data.frame(partcor(selexp, pt, conf)[2])

#Extract significant exposures (FDR >0.05)
signexp = ifelse(bmipval$phenotype....hs_zbmi_who..<0.05, rownames(bmipval), NA)
na.omit(signexp)
signexpBMI=signexp

write.csv(signexpBMI, file="Partcor_BMI.csv")

#********************
#IQ 
#********************#####

#Get datasets 
selexp = exposome[,colnames(exposome) %in% eIQ[,1]]  #Selected exposures from MUVR 
conf = covariates[, c("h_cohort", "e3_sex_None", "hs_child_age_None")]   #Confounders: Age, Sex, Cohort, others?
pt = as.data.frame(phenotype[,"hs_correct_raven"])   #Phenotypes: IQ (RAVEN score)

#Do partial correlation 
partcor(selexp, pt, conf)
iqcor = as.data.frame(partcor(selexp, pt, conf)[1])
iqpval = as.data.frame(partcor(selexp, pt, conf)[2])

#Extract significant exposures (FDR >0.05)
signexp = ifelse(iqpval$phenotype....hs_correct_raven..<0.05, rownames(iqpval), NA)
na.omit(signexp)
signexpIQ=signexp

write.csv(signexpIQ, file="Partcor_IQ.csv")



####Finding all unique exposures to be modelled against the multiomics dataset
allPartCorr = unique(c(signexpBMI,signexpIQ))
allPartCorr = allPartCorr[!is.na(allPartCorr)]

save.image(file="selected_exposures_partcor.Rdata")
