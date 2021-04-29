##########################################
######## Mediation analysis with #########
##### Omics.Factor2 -- Expo.Facotr2 ######
##### Omics.Factor3 -- Expo.Facotr1 ######
##########################################
library(lavaan)
library(semPlot)


factors_OmicExpo<-readRDS("data_step_III_DNAm_0.05filtered.rds")


#setting all categorical covariates to numeric for multilevel lavaan
factors_OmicExpo$NUM_e3_sex_None<-as.numeric(factors_OmicExpo$e3_sex_None)
factors_OmicExpo$NUM_h_parity_None<-as.numeric(factors_OmicExpo$h_parity_None)
factors_OmicExpo$NUM_h_edumc_None<-as.numeric(factors_OmicExpo$h_edumc_None)


# merging level 0 and 1 of h_native_None to 1 and setting to numeric 
factors_OmicExpo$h_native_None_bin<-factors_OmicExpo$h_native_None
factors_OmicExpo$h_native_None_bin[factors_OmicExpo$h_native_None_bin == 0]<-1
factors_OmicExpo$NUM_h_native_None_bin<-as.numeric(factors_OmicExpo$h_native_None_bin)
# check
Hmisc::describe(as.factor(factors_OmicExpo$h_native_None_bin))

# performing multilevel mediation analysis

multilevelMediation <- '
level:1
hs_zbmi_who ~  e1 * Omics.Factor2 +  e2 * Omics.Factor3 + f1 * Expo.Factor1 + f2 * Expo.Factor2 + NUM_e3_sex_None + NUM_h_native_None_bin + e3_gac_None + NUM_h_parity_None + h_mbmi_None + NUM_h_edumc_None + h_age_None + hs_child_age_None
Omics.Factor2 ~ d1 * Expo.Factor2 + NUM_e3_sex_None + NUM_h_native_None_bin + e3_gac_None + NUM_h_parity_None + h_mbmi_None + NUM_h_edumc_None + h_age_None + hs_child_age_None
Omics.Factor3 ~ d2 * Expo.Factor1 + NUM_e3_sex_None + NUM_h_native_None_bin + e3_gac_None + NUM_h_parity_None + h_mbmi_None + NUM_h_edumc_None + h_age_None + hs_child_age_None

#indirect and total effects within
indirect1 := d1 * e1
indirect2 := d2 * e2
totalwith := f1 + f2 + (d1 * e1) + (d2 * e2) 

level:2
hs_zbmi_who ~~ hs_zbmi_who
hs_zbmi_who ~~ Omics.Factor2
hs_zbmi_who ~~ Omics.Factor3
Omics.Factor2 ~~ Omics.Factor2
Omics.Factor3 ~~ Omics.Factor3  
'

fit <- sem(model = multilevelMediation, data = factors_OmicExpo, cluster = "h_cohort", 
           se='standard', verbose = TRUE, missing = "ML")
summary(fit, standardized=TRUE)


tableEstimates<-parameterEstimates(fit, standardized = TRUE) 


###################################################
# Visualization of the path model
###################################################
# For visualization only take the simple model with variables of interest and manually add correct estimates/edge values 
multilevelMediationX<- '
hs_zbmi_who ~  e1 * Omics.Factor2 +  e2 * Omics.Factor3 + f1 * Expo.Factor1 + f2 * Expo.Factor2
Omics.Factor2 ~ d1 * Expo.Factor2
Omics.Factor3 ~ d2 * Expo.Factor1
'

fitX <- sem(model = multilevelMediationX, data = factors_OmicExpo)


# define graph parameters
myly<-matrix(c(2, 0,
               0,2,
                0,-2,
               -2,-2,
               -2, 2), ncol=2, byrow=TRUE)
labelsedge<-c("e1= .072***", "e2= -.068*** ","f1= .024","f2= .035","d1= .16*","d2= -.14*")
labelsedge_scaled<-c("e1= .16***", "e2= -.15*** ","f1= .062","f2= .090","d1= .19*","d2= -.17*")
labelsnode<-c("   zBMI    ","Omics.F2", "Omics.F3", " Expo.F1 "," Expo.F2 ")

# print figure unscaled estimates
semPaths(fitX,"est", intercepts = FALSE,layout=myly,exoVar = FALSE, exoCov = FALSE, residuals=FALSE,label.cex=1.7,sizeMan=8, edge.label.cex=1.3, weighted=FALSE, nodeLabels=labelsnode, edgeLabels= labelsedge)

# print figure scaled estimates
semPaths(fitX,"est", intercepts = FALSE,layout=myly,exoVar = FALSE, exoCov = FALSE, residuals=FALSE,label.cex=1.7, sizeMan=8,edge.label.cex=1.3, weighted=FALSE, nodeLabels=labelsnode,edgeLabels= labelsedge_scaled)
