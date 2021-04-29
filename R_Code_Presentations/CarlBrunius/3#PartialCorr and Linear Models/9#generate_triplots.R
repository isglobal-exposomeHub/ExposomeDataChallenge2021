############################################################
# Exposome Challenge 
# Omics triplots 
############################################################

#Use selected omic variables from MUVR, partcor and lm modelling
#Do PCA on selected omic variables 
#Add correlations between PC scores and exposures 
#Add associations between PC scores and outcomes 

library(triplot)
library(ggplot2)
library(psych)
library(GPArotation)
library(MUVR)

#****************************************************************
#****************************************************************
#TRIPLOT MULTIOMICS 57(48+16) variables sign after GLM but no FDR
#****************************************************************
#****************************************************************

load("finalOmicsVars.Rdata")
load("exposomeEdited.RData")
load("metabol_serum.RData")
load("metabol_urine.RData")


####Combining all XVars into one dataframe####
#1. Combine all input-XVars into one dataframe
allXVarsCombined<-data.frame(finalMultiOmicsXVars[[1]][,1])

for(l in 1:length(finalMultiOmicsXVars)){
  allXVarsCombined<-cbind(allXVarsCombined,finalMultiOmicsXVars[[l]][,-1])
}

allXVarsCombined<-allXVarsCombined[,unique(colnames(allXVarsCombined))]
dim(allXVarsCombined)

#2. Take out the variables which made it through the filter steps
allVIPsFromModels<-vector()
for(l in 1:length(finalMultiOmicsmodels)){
  allVIPsFromModels<-c(allVIPsFromModels,getVIP(finalMultiOmicsmodels[[l]])$name)
}
allVIPsFromModels<-unique(allVIPsFromModels)

smetbol_list<-intersect(unique(unlist(PRAVENLIST)),allVIPsFromModels)
tempBMI<-intersect(unique(unlist(PBMILIST)),allVIPsFromModels)
smetbol_list<-unique(c(smetbol_list,tempBMI))
smetbol_list<-smetbol_list[order(smetbol_list)]

####Make the triplot####
#1. Preparing dataframe 
MultiOmics<-allXVarsCombined
MultiOmics$ID=as.integer(allXVarsCombined[,1])

covariates2=merge(MultiOmics, covariates, "ID", "ID")
exposome2=merge(MultiOmics, exposome, "ID", "ID")
covariates2=covariates2[,-c(2:length(MultiOmics[1,]))] #remove columns of metabolites
exposome2=exposome2[,-c(2:length(MultiOmics[1,]))]
phenotype2=merge(MultiOmics, phenotype, "ID", "ID")
phenotype2=phenotype2[,-c(2:length(MultiOmics))]

writeLines(paste("Identical:",identical(MultiOmics$ID, covariates2$ID, exposome2$ID, phenotype2$ID)))

smetbol=MultiOmics[,colnames(MultiOmics)%in%smetbol_list]

#2. Changing the metab names to real names
serumMetNames<-rownames(metabol_serum@assayData$exprs)
urineMetNames<-rownames(metabol_urine@assayData$exprs)

smetbol<-smetbol[,order(colnames(smetbol))]
#For urine
colnames(smetbol)[regexpr("\\.1",colnames(smetbol))>0 & regexpr("\\.1",colnames(smetbol))<11]<-urineMetNames[c(25,6)]

#For serum
tempNames<-vector()
for(l in 1:length(colnames(smetbol)[grep("metab",colnames(smetbol))])){
  tempNames[l]<-serumMetNames[as.integer(strsplit(colnames(smetbol)[grep("metab",colnames(smetbol))][l],"_")[[1]][2])]
}
colnames(smetbol)[grep("metab",colnames(smetbol))]<-tempNames

#3.Making the plot
pca=principal(smetbol, rotate="oblimin", nfactors = 2)
tpo=makeTPO(scores=pca$scores, loadings=pca$loadings)
triPlot(tpo, comps = c(1:2))

#CORRELATION WITH EXPOSURES LAYER 
conf = covariates2[, c("h_cohort", "e3_sex_None", "hs_child_age_None")]  #Confounders: Age, Sex, Cohort, others?
selexp = as.data.frame(exposome2[,c("hs_hcb_cadj_Log2", "hs_pcb118_cadj_Log2","hs_pcb138_cadj_Log2", "hs_pcb153_cadj_Log2", "hs_pcb170_cadj_Log2", "hs_pcb170_madj_Log2", "hs_pcb180_cadj_Log2", "hs_pfhxs_c_Log2", "hs_pfna_c_Log2", "hs_pm10_yr_hs_h_None", "hs_sumPCBs5_cadj_Log2")])  #selected 11 (MUVR, partcor and MM)

partcor(tpo$scores, selexp, conf) #from partcor script 
expcor = as.data.frame(partcor(tpo$scores, selexp, conf)[1])
tpo=addCorr(tpo, t(expcor))
rownames(tpo$corrMatrix)=c("HCB", "PCB118","PCB138", "PCB153","PCB170_c", "PCB170_m", "PCB180", "PFHxS", "PFNA", "PM10yr", "PCB5sum")

#ASSOCIATION WITH OUTCOME LAYER 
coefficients <- matrix(nrow = ncol(tpo$scores), ncol = 3)
colnames(coefficients) <- c("coefficients", "se", "p")
for (i in 1:ncol(tpo$scores)) {
  glmod<-summary(glm(phenotype2$hs_zbmi_who~tpo$scores[,i]+covariates2$h_cohort+covariates2$e3_sex_None+covariates2$hs_child_age_None,family=gaussian))
  coefficients[i, ] <- glmod$coefficients[2,c(1, 2, 4)]
}

tpo=addRisk(tpo, coefficients, name="bmi")

coefficients <- matrix(nrow = ncol(tpo$scores), ncol = 3)
colnames(coefficients) <- c("coefficients", "se", "p")
for (i in 1:ncol(tpo$scores)) {
  glmod<-summary(glm(phenotype2$hs_correct_raven~tpo$scores[,i]+covariates2$h_cohort+covariates2$e3_sex_None+covariates2$hs_child_age_None,family=gaussian))
  coefficients[i, ] <- glmod$coefficients[2,c(1, 2, 4)]
}

tpo=addRisk(tpo, coefficients, name="raven")


#full adjusted models + mediators 
coefficients <- matrix(nrow = ncol(tpo$scores), ncol = 3)
colnames(coefficients) <- c("coefficients", "se", "p")
for (i in 1:ncol(tpo$scores)) {
  glmod<-summary(glm(phenotype2$hs_zbmi_who~tpo$scores[,i]+covariates2$h_cohort+covariates2$e3_sex_None+covariates2$hs_child_age_None +
                       covariates2$h_edumc_None + covariates2$h_age_None + covariates2$h_mbmi_None + covariates2$h_parity_None + covariates2$e3_yearbir_None + 
                       covariates2$hs_wgtgain_None  + covariates2$e3_gac_None,family=gaussian))
  coefficients[i, ] <- glmod$coefficients[2,c(1, 2, 4)]
}

tpo=addRisk(tpo, coefficients, name="bmi_fulladj")

coefficients <- matrix(nrow = ncol(tpo$scores), ncol = 3)
colnames(coefficients) <- c("coefficients", "se", "p")
for (i in 1:ncol(tpo$scores)) {
  glmod<-summary(glm(phenotype2$hs_correct_raven~tpo$scores[,i]+covariates2$h_cohort+covariates2$e3_sex_None+covariates2$hs_child_age_None + 
                       covariates2$h_edumc_None + covariates2$h_age_None + covariates2$h_mbmi_None + covariates2$h_parity_None + covariates2$e3_yearbir_None + 
                       covariates2$hs_wgtgain_None  + covariates2$e3_gac_None,family=gaussian))
  coefficients[i, ] <- glmod$coefficients[2,c(1, 2, 4)]
}

tpo=addRisk(tpo, coefficients, name="raven_fulladj")


#actual triplot
checkTPO(tpo)
triPlot(tpo, comps = c(1,2), loadCut=0.2) 

#to read correlations or associations 
tpo$corrMatrix #correlations 
tpo$riskMatrix #B coefficients associations
tpo$riskP  #p-values associations
tpo$riskSE #standard errors associations

pdf(file="Triplot2comp.pdf")
triPlot(tpo, comps=c(1,2), loadCut=0.2)
dev.off()
methy@rowRanges[which(rownames(as.data.frame(methy@rowRanges))=="cg27255454"),5:6]

