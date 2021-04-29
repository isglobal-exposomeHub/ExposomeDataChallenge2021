#Exposome Data Challenge
#Presentation Title: "Exposure Continuum Mapping for predicting health and disease in Exposome studies"
#Author: John L. Pearce
#Web: https://github.com/johnlpearce
#Date: 27 April 2021

#Bring in ECM package from GitHub
library(devtools)
devtools::install_github("johnlpearce/ECM")
library(ECM)

#Load Data in Exposome Workspace
load("exposome.RData")
#Examine data
str(covariates)
str(exposome)
str(phenotype)
str(codebook)

###################################################################################
###################################################################################
#Prepare exposome data for analysis
##Set row identifier for linkage
rownames(exposome)<-exposome$ID

##Examine codebook
names(codebook)
levels(codebook$period)
head(codebook[, c("variable_name", "labelsshort")])

##Identify prenatal data
pre_codes<-subset(codebook, period == "Pregnancy")
table(pre_codes$domain)
#subset prenatal exposures
pre_codes_expo<-subset(pre_codes, domain %in% c("Chemicals", "Outdoor exposures", "Lifestyles"))
table(pre_codes_expo$domain)
#Order data based on domain
pre_codes_expo2<-pre_codes_expo[order(pre_codes_expo$domain),]

#Build a prenatal codebook for labeling
pre_codes_tab<-data.frame(Variable_name=as.character(pre_codes_expo2$variable_name), 
                    labelsshort=as.character(pre_codes_expo2$labelsshort), 
                    domain=as.character(pre_codes_expo2$domain))
pre_codes_tab

#Subset exposome data to prenatal exposures only
pre_expo_data<-exposome[,pre_codes_tab$Variable_name]
str(pre_expo_data)
#Set colnames to short labels
colnames(pre_expo_data)<-pre_codes_tab$labelsshort
str(pre_expo_data)

##Subset prenatal exposome data to 72 continuous variables
pre_expo_data_cont<- pre_expo_data[,sapply(pre_expo_data,is.numeric)] 
str(pre_expo_data_cont)

#Specify Custom Colors for plotting based on exposure domain
library(colorspace)
custcol<-c(rev(diverge_hcl(48)), terrain_hcl(24))
custcol

#Identify variable labels
ecm_var_labs<-names(pre_expo_data_cont)

#######################################################################################
#######################################################################################
#Apply Exposure Continuum Mapping
#Step 1: Data Standardization
ecm_trn<-scale(pre_expo_data_cont)
summary(ecm_trn) #All means are now zero and units are SDs

#Step 2: Select an appropriate map size
ms<-map_size(trn_dat=ecm_trn, kmn=2, kmx=50, nstarts.=1)
ms$MAP_EVAL

#Step 3: Construct final mapping
#Identify optimal initialization values selected size
mi<-map_inits(trn_dat=ecm_trn, xdim=7, ydim=6, nstarts=5)
mi

#Construct map
expo_ecm<-map_ecm(ecm_trn, xdim=7, ydim=6, inits=mi$opt_init)

#Evaluate mapping
map_stats(expo_ecm)

#Step 4: Map Visualization
map_plot(expo_ecm, addFreq=FALSE, colormod=custcol, nodelab=TRUE, labsize=0.8)

#########################################################################
#Step 5: Exploring Health effects 
##Prep mod data
moddata<-merge(phenotype, covariates, by="ID")
##Extract ECM exposure metric
ecm_metric<-map_metric(expo_ecm)
head(ecm_metric)

#Merge moddat with exposure metric
moddata2<-merge(moddata, ecm_metric, by.x="ID", by.y="OBS")
head(moddata2)

#Apply generalized additive model for exploration of joint dose response across total exposome
library(mgcv)
#Fit model
ecm_gam<-gam(e3_bw~h_cohort + h_age_None + h_edumc_None + h_native_None+ h_parity_None + h_mbmi_None 
            + e3_sex_None + te(U, V), data=moddata2, method="REML")
#Examine Results
summary(ecm_gam)

#Visualize joint-dose response function
par(mfrow=c(1,1), pty="m", cex.lab=1.5, cex.axis=1.5)
vis.gam(ecm_gam, view=c("U", "V"), color="topo", type="link",
              ticktype="detailed", theta=45, phi=45, expand=0.75,
              zlab="Y", xlab="U-Coordinate", ylab="V-Coordinate")
title(main=" a) ECM Joint Dose-Response Function")

#Create component plane maps to assist with exploration
comptab<-data.frame(expo_ecm$codes)
dim(comptab)
library(kohonen)
par(mfrow=c(9,8), family="serif", pty="m")
for (i in 1:dim(ecm_trn)[2]){
  comp<-comptab[,i]
  plot(expo_ecm, type="property", property=comp, main=names(comptab)[i], 
       palette.name=diverging_hcl, heatkey=FALSE)
}
