# File: A_01_prepare_exposome_data.R
# Author(s): Jaime Benavides, Lawrence Chillrud
# Date since last edit: 4/27/21

####*******************####
#### Table of Contents ####
####*******************####
####* Notes / Description
####* 0. Package imports
####* 1. Read in data
####* 2. Subset urban exposures
####* 3. Back-transform continuous (a) pre- and (b) postnatal exposures
####* 4. Save data

####*********************####
#### Notes / Description ####
####*********************####
####* This processes the given (raw) ExposomeChallenge data
####* into data we will use for our submission.

####********************####
#### 0. Package Imports ####
####********************####

print("Running script A_01...")

# 0a list packages needed
list.of.packages = c('tidyverse', 'haven','ggplot2','GGally','caret','readr',
                     'fBasics','raster', 'nortest','MASS','extrafont', 'forestplot',
                     'ggsn','survival','grid','gridExtra', 'cowplot','reshape2',
                     'naniar','rpart','rpart.plot','pROC', 'RColorBrewer', 'broom',
                     'car','olsrr','EnvStats','janitor','tableone', 'reshape', 
                     'magrittr', 'tidyselect', 'dplyr', 'tidyr', 'stringr', 'purrr', 'tibble',
                     'pracma', 'Matrix', 'foreach', 'iterators', 'parallel', 'snow', 'doSNOW',
                     'GPfit', 'heatmaply', 'kableExtra', 'mgcv', 'glmnet')

# 0b check if list of packages is installed. if not, it will install ones not yet installed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# 0c import packages
lapply(list.of.packages, require, character.only = TRUE)

####*****************####
#### 1. Read in data ####
####*****************####
load(file = paste0(raw.data.folder, "exposome_v2.RData")) # only complete cases

####***************************####
#### 2. Subset urban exposures ####
####***************************####
dom <- c("Outdoor exposures", "Indoor air")

description <- codebook %>% 
  dplyr::filter(!domain %in% c('Covariates','Phenotype')) 

# subset variables in the exposome, choose only the ones in the domains above
fam <- unique(description[description$domain %in% dom, "family"])
urban_exposome_fam <- which(description$family %in% fam)
urban_exposome <- exposome[, c("ID", rownames(description[urban_exposome_fam,]))]
urban_description <- description[urban_exposome_fam,]

####***********************************************************####
#### 3. Back-transform continuous pre- and postnatal exposures ####
####***********************************************************####

####************************####
#### 3a. Prenatal Exposures ####
####************************####
pren_urb_vars <- urban_description[urban_description$period == "Pregnancy", "variable_name"] 
pren_urb_vars <- which(urban_description$variable_name %in% pren_urb_vars)
pren_exposome <- urban_exposome[, c("ID", rownames(urban_description[pren_urb_vars,]))]
pren_description <- urban_description[pren_urb_vars, ]

# leave out categorical variables
categ_vars <- pren_description[pren_description$var_type == "factor"  |  	
                                 pren_description$transformation == "Dichotomous", "variable_name"]

rownames(pren_exposome) <- pren_exposome[, "ID"]
pren_exposome <- pren_exposome[, -1]

pren_exposome <- pren_exposome %>%
  dplyr::select(-categ_vars)
pren_description <- pren_description[-which(pren_description$variable_name %in% categ_vars) , ]

# back transform the exposure variables that have been altered to its raw shape 
# natural logarithm -> exp()
log_vars <- which(pren_description$transformation == "Natural Logarithm") 
pren_exposome[,log_vars] <- exp(pren_exposome[,log_vars]) 

# change variable names
pren_description$variable_name <- as.character(pren_description$variable_name)
pren_description[log_vars,"variable_name"] <- gsub("_Log", "", pren_description[log_vars,"variable_name"])

# square root -> ^2 
sqrt_vars <- which(pren_description$transformation == "Square root")
pren_exposome[,sqrt_vars] <- (pren_exposome[,sqrt_vars])^2

# change variable names
pren_description[sqrt_vars,"variable_name"] <- gsub("_Sqrt", "", pren_description[sqrt_vars,"variable_name"])
pren_description[,"variable_name"] <- gsub("_None", "", pren_description[,"variable_name"])

# cubic root -> ^3
cubrt_vars <- which(grepl("_pow1over3", pren_description$variable_name, fixed = TRUE))
pren_exposome[,cubrt_vars] <- (pren_exposome[,cubrt_vars])^3
pren_description[cubrt_vars,"variable_name"] <- gsub("_pow1over3", "", pren_description[cubrt_vars,"variable_name"])

# update rownames
row.names(pren_description) <- pren_description$variable_name
colnames(pren_exposome) <- row.names(pren_description)

####*************************####
#### 3b. Postnatal Exposures ####
####*************************####
post_urb_vars <- urban_description[urban_description$period == "Postnatal", "variable_name"] 
post_urb_vars <- which(urban_description$variable_name %in% post_urb_vars)
post_exposome <- urban_exposome[, c("ID", rownames(urban_description[post_urb_vars,]))]
post_description <- urban_description[post_urb_vars, ]

# leave out categorical variables
categ_vars <- post_description[post_description$var_type == "factor" |  	
                                 post_description$transformation == "Dichotomous", "variable_name"]

rownames(post_exposome) <- post_exposome[, "ID"]
post_exposome <- post_exposome[, -1]

post_exposome <- post_exposome %>%
  dplyr::select(-categ_vars)
post_description <- post_description[-which(post_description$variable_name %in% categ_vars) , ]

# choose only exposure variables with periods higher than one week   
week_day_vars <- which(grepl('week|day', post_description$labels))
post_exposome <- post_exposome[,-week_day_vars]
post_description <- post_description[-week_day_vars, ]

# back transform the exposure variables that have been altered to its raw shape 
# natural logarithm -> exp()
log_vars <- which(post_description$transformation == "Natural Logarithm")
post_exposome[,log_vars] <- exp(post_exposome[,log_vars])

# change variable names
post_description$variable_name <- as.character(post_description$variable_name)
post_description[log_vars,"variable_name"] <- gsub("_Log", "", post_description[log_vars,"variable_name"])

# square root -> ^2 
sqrt_vars <- which(post_description$transformation == "Square root")
post_exposome[,sqrt_vars] <- (post_exposome[,sqrt_vars])^2

# change variable names
post_description[sqrt_vars,"variable_name"] <- gsub("_Sqrt", "", post_description[sqrt_vars,"variable_name"])
post_description[,"variable_name"] <- gsub("_None", "", post_description[,"variable_name"])

# cubic root -> ^3
cubrt_vars <- which(grepl("_pow1over3", post_description$variable_name, fixed = TRUE))
post_exposome[,cubrt_vars] <- (post_exposome[,cubrt_vars])^3
post_description[cubrt_vars,"variable_name"] <- gsub("_pow1over3", "", post_description[cubrt_vars,"variable_name"])

# update rownames
row.names(post_description) <- post_description$variable_name
colnames(post_exposome) <- row.names(post_description)

####**************####
#### 4. Save data ####
####**************####
save(codebook, pren_exposome, pren_description, post_exposome, post_description, covariates, phenotype,
     file = paste0(generated.data.folder, "exposome_pren_postnatal.RData"))





