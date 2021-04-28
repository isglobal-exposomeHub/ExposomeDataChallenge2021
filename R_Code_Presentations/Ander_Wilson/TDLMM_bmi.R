options(echo=TRUE)
options(stringsAsFactors = FALSE)


library(tidyverse)
library(readxl)
library(dlmtree)


setwd("../")
load("exposome_v2.RData")
load("DataPrep/covariates.rda")
load("DataPrep/time_exposure_list.rda")


# Combine responses/covariates
y <- scale(phenotype$hs_zbmi_who)
dta <- as.data.frame(cbind(y,z_base_preg,z_base_post))
# remove covariates that are already adjust for in the outcome or related to the outcome
dta <- dta[,-which(colnames(dta)%in%c("hs_child_age_None","e3_sex_Nonemale","hs_c_weight_None","hs_c_height_None"))]
# rename outcome
colnames(dta)[1] <- "y"
colnames(dta)


# Model run parameters
burn <- 100000; iter <- 500000; thin <- 10
kappa <- 0.371
trees <- 50
inttype <- "noself"

set.seed(513458)

fit <- tdlmm(y ~  ., 
             data = dta, 
             exposure.data = time_exposure_list,
             mixture.interactions = inttype, # also try 'all' and 'none'
             family = "gaussian", # change to 'gaussian' for continuous response
             n.trees = trees, n.burn = burn, n.iter = iter, n.thin = thin, 
             mix.prior = kappa, 
             shrinkage = "exposures")
sfit <- summary(fit)
sfit


save(sfit,fit, file=paste0("TDLMM/Output/TDLMM_bmi_",inttype,"_",trees,".Rdata"))



