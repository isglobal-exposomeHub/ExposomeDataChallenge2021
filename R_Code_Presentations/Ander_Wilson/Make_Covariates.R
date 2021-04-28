library(tidyverse)
library(readxl)


# load data
codebook <- read_xlsx("codebook.xlsx")
load("exposome_v2.RData")
ls()


## names of covariates and other exposure groups
covariates_preg <- codebook %>% 
  filter(domain=="Covariates" & period == "Pregnancy") %>%
  pull(variable_name) %>%
  as.character()

covariates_post <- codebook %>% 
  filter(domain=="Covariates" & period == "Postnatal") %>%
  pull(variable_name) %>%
  as.character()

## other exposure groups not being included in kernel
lifestyles_preg <- codebook %>% 
  filter(domain=="Lifestyles" & period == "Pregnancy") %>%
  pull(variable_name) %>%
  as.character()

lifestyles_post <- codebook %>% 
  filter(domain=="Lifestyles" & period == "Postnatal") %>%
  pull(variable_name) %>%
  as.character()

socioecon_post <- codebook %>%  ## there is no pregnancy version of these
  filter(domain=="Chemicals" & family=="Social and economic capital" & period == "Postnatal") %>%
  pull(variable_name) %>%
  as.character()

tobacco_preg <- codebook %>%  
  filter(domain=="Chemicals" & family=="Tobacco Smoke" & period == "Pregnancy") %>%
  pull(variable_name) %>%
  as.character()

tobacco_post <- codebook %>%  
  filter(domain=="Chemicals" & family=="Tobacco Smoke" & period == "Postnatal") %>%
  pull(variable_name) %>%
  as.character()

natspaces_preg <- codebook %>%  
  filter(domain=="Outdoor exposures" & family=="Natural Spaces" & period == "Pregnancy") %>%
  pull(variable_name) %>%
  as.character()

natspaces_post <- codebook %>%  
  filter(domain=="Outdoor exposures" & family=="Natural Spaces" & period == "Postnatal") %>%
  pull(variable_name) %>%
  as.character()

noise_preg <- codebook %>%  
  filter(domain=="Outdoor exposures" & family=="Noise" & period == "Pregnancy") %>%
  pull(variable_name) %>%
  as.character()

noise_post <- codebook %>%  
  filter(domain=="Outdoor exposures" & family=="Noise" & period == "Postnatal") %>%
  pull(variable_name) %>%
  as.character()


### get covariates and excess exposure groups

cov_base_preg <- covariates[,covariates_preg]
cov_base_post <- covariates[,covariates_post]
cov_lifestyles_preg <- exposome[,lifestyles_preg]
cov_lifestyles_post <- exposome[,lifestyles_post]
cov_socioecon_post <- exposome[,socioecon_post]
cov_tobacco_preg <- exposome[,tobacco_preg]
cov_tobacco_post <- exposome[,tobacco_post]
cov_natspaces_preg <- exposome[,natspaces_preg]
cov_natspaces_post <- exposome[,natspaces_post]
cov_noise_preg <- as.matrix(exposome[,noise_preg]) # one column
cov_noise_post <- exposome[,noise_post]

### handle specific variables
## convert year as factor (7 levels) to year as continuous var
cov_base_preg[,"e3_yearbir_None"] <- as.numeric(cov_base_preg[,"e3_yearbir_None"]) ## convert year to numeric (was a 7 level factor)
## candidate to be dropped/collapsed due to sparsity: h_fruit_preg_Ter (level 1 only has 6 observations)
cov_lifestyles_preg[,"h_fruit_preg_Ter"] <- fct_collapse(cov_lifestyles_preg[,"h_fruit_preg_Ter"], 
                                                         "(0,18.2]"=c("(0,0.6]", "(0.6,18.2]"), ## collapse first two categories
                                                         "(18.2,Inf]"=c("(18.2,Inf]"))
## Variable "e3_asmokcigd_p_None": 87% zeros, and then values from 1 to 15 resulting in an extremely skewed distribution.  could dichotomize it (uncomment the line below)
# cov_tobacco_preg[,"e3_asmokcigd_p_None"] <- as.numeric(cov_tobacco_preg[,"e3_asmokcigd_p_None"]>0) ## could dichotomize (see histogram)

## standardize continuous covariates
stdize_covs <- function(x){
  for(cc in 1:ncol(x)){
    if(is.numeric(x[,cc])){
      x[,cc] <- scale(x[,cc])
    }
  }
  return(x)
}

cov_base_preg <- stdize_covs(cov_base_preg)
cov_base_post <- stdize_covs(cov_base_post)
cov_lifestyles_preg <- stdize_covs(cov_lifestyles_preg)
cov_lifestyles_post <- stdize_covs(cov_lifestyles_post)
cov_socioecon_post <- stdize_covs(cov_socioecon_post)
cov_tobacco_preg <- stdize_covs(cov_tobacco_preg)
cov_tobacco_post <- stdize_covs(cov_tobacco_post)
cov_natspaces_preg <- stdize_covs(cov_natspaces_preg)
cov_natspaces_post <- stdize_covs(cov_natspaces_post)
cov_noise_preg <- stdize_covs(cov_noise_preg)
cov_noise_post <- stdize_covs(cov_noise_post)

## get design matrix
z_base_preg <- model.matrix(~.,data=cov_base_preg)[,-1]
z_base_post <- model.matrix(~.,data=cov_base_post)[,-1]
z_lifestyles_preg <- model.matrix(~.,data=cov_lifestyles_preg)[,-1]
z_lifestyles_post <- model.matrix(~.,data=cov_lifestyles_post)[,-1]
z_socioecon_post <- model.matrix(~.,data=cov_socioecon_post)[,-1]
z_tobacco_preg <- model.matrix(~.,data=cov_tobacco_preg)[,-1]
z_tobacco_post <- model.matrix(~.,data=cov_tobacco_post)[,-1]
z_natspaces_preg  <- model.matrix(~.,data=cov_natspaces_preg )[,-1]
z_natspaces_post  <- model.matrix(~.,data=cov_natspaces_post )[,-1]
z_noise_preg <- model.matrix(~.,data=data.frame(cov_noise_preg))[,-1]
z_noise_post <- model.matrix(~.,data=cov_noise_post)[,-1]

### combine groups for different covariate options
z_base_preg <- z_base_preg ## just "covariates" as labelled in codebook
z_base_post <- z_base_post ## just "covariates" as labelled in codebook

z_full_preg <- cbind(z_base_preg,  ## including exposure groups with mostly factors
                     z_lifestyles_preg,
                     z_tobacco_preg,
                     z_natspaces_preg,
                     z_noise_preg)
                     
z_full_post <- cbind(z_base_post,   ## including exposure groups with mostly factors
                     z_lifestyles_post,
                     z_tobacco_post,
                     z_natspaces_post,
                     z_noise_post,
                     z_socioecon_post)

save(z_base_preg,z_base_post,z_full_preg,z_full_post, file="DataPrep/covariates.rda")
