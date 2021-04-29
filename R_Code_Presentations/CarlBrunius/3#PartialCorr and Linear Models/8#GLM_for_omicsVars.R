############################################################
# Exposome Challenge 
# Metabolome modules - Selected omics variables - linear modelling 
############################################################

#Use MUVR & partial correlation selected exposures 
#Use exposures with metabolite MUVR models q2>0.2
#Put exposures and metabolites in partial correlation 
#Adjust for age, sex and cohort 
#Select relevant metabolites for each exposure 
#Linear models for selected metabolites with outcomes (raven/bmi)

load("exposomeEdited.Rdata")
load("PartCorrOmicsVars.RData")


PBMILIST<-list()
PRAVENLIST<-list()

for(i in 1:length(finalMultiOmicsmodels)){
  MultiOmics=data.frame(finalMultiOmicsXVars[[i]][,intersect(colnames(finalMultiOmicsXVars[[i]]),getVIP(finalMultiOmicsmodels[[i]])$name)])
  MultiOmics$ID=as.numeric(finalMultiOmicsXVars[[i]][,1]) #rownames is ID number
  
  covariates2=merge(MultiOmics, covariates, "ID", "ID")
  exposome2=merge(MultiOmics, exposome, "ID", "ID")
  covariates2=covariates2[,-c(2:length(MultiOmics[1,]))] #remove columns of metabolites
  exposome2=exposome2[,-c(2:length(MultiOmics[1,]))]
  phenotype2=merge(MultiOmics, phenotype, "ID", "ID")
  phenotype2=phenotype2[,-c(2:length(MultiOmics))]
  
  #sanity check identical datasets, already existing matched ds?
  writeLines(paste("Identical:",identical(MultiOmics$ID, covariates2$ID, exposome2$ID, phenotype2$ID)))
  
  #********************
  #EXPOSURES 
  #********************####
  
  ###BMI
  
  smetbol_list<-MOBioData
  smetbol=MultiOmics[,colnames(MultiOmics)%in%smetbol_list ]
  
  coefficients <- matrix(nrow = ncol(smetbol), ncol = 3)
  colnames(coefficients) <- c("coefficients", "se", "p")
  for (l in 1:ncol(smetbol)) {
    glmod<-summary(glm(phenotype2$hs_zbmi_who~smetbol[,l]+covariates2$h_cohort+covariates2$e3_sex_None+covariates2$hs_child_age_None,family=gaussian))
    coefficients[l, ] <- glmod$coefficients[2,c(1, 2, 4)]
  }
  pDF<-data.frame(p=coefficients[,3])
  rownames(pDF)<-colnames(smetbol)
  signmetbolF = ifelse(pDF[,1]<0.05, rownames(pDF), NA)
  PBMILIST[[i]]<-na.omit(signmetbolF)
  
  
  ###RAVEN
  
  coefficients <- matrix(nrow = ncol(smetbol), ncol = 3)
  colnames(coefficients) <- c("coefficients", "se", "p")
  for (l in 1:ncol(smetbol)) {
    glmod<-summary(glm(phenotype2$hs_correct_raven~smetbol[,l]+covariates2$h_cohort+covariates2$e3_sex_None+covariates2$hs_child_age_None,family=gaussian))
    coefficients[l, ] <- glmod$coefficients[2,c(1, 2, 4)]
  }
  pDF<-data.frame(p=coefficients[,3])
  rownames(pDF)<-colnames(smetbol)
  signmetbolF = ifelse(pDF[,1]<0.05, rownames(pDF), NA)
  PRAVENLIST[[i]]<-na.omit(signmetbolF)
  
}

save.image("finalOmicsVars.Rdata")
