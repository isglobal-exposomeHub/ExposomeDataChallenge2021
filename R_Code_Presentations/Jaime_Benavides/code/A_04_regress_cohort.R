# File: A_04_regress_cohort.R
# Author(s): Jaime Benavides
# Date since last edit: 4/27/21

####*******************####
#### Table of Contents ####
####*******************####
####* Notes / Description
####* 0. Package imports
####* 1. Read in data
####* 2. Run regressions
####* 3. Save residuals

####*********************####
#### Notes / Description ####
####*********************####
####* run one linear model for each exposure variable (dependent) 
####* given cohort (independent), and obtain the residuals.

####********************####
#### 0. Package Imports ####
####********************####

print("Running script A_04...")

# 0a list packages needed:
list.of.packages = c('tidyverse', 'haven','ggplot2','GGally','caret','readr',
                     'fBasics','raster', 'nortest','MASS','extrafont', 'forestplot',
                     'ggsn','survival','grid','gridExtra', 'cowplot','reshape2',
                     'naniar','rpart','rpart.plot','pROC', 'RColorBrewer', 'broom',
                     'car','olsrr','EnvStats','janitor','tableone', 'reshape', 
                     'magrittr', 'tidyselect', 'dplyr', 'tidyr', 'stringr', 'purrr', 'tibble',
                     'pracma', 'Matrix', 'foreach', 'iterators', 'parallel', 'snow', 'doSNOW',
                     'GPfit', 'heatmaply', 'kableExtra', 'mgcv', 'glmnet')

# 0b check if list of packages is installed. if not, it will install ones not yet installed:
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# 0c import packages:
lapply(list.of.packages, require, character.only = TRUE)

####*****************####
#### 1. Read in data ####
####*****************####
load(paste0(generated.data.folder, "exposome_pren_postnatal.RData"))

# 1a check that both datasets share the same id order:
any(row.names(pren_exposome) != row.names(post_exposome))

# 1b scale data:
pren_exposome_scaled <- pren_exposome %>% scale(., center = F) %>% as.data.frame()
post_exposome_scaled <- post_exposome %>% scale(., center = F) %>% as.data.frame()

# 1c add cohort variable to exposure variables datasets:
cohorts <- covariates$h_cohort
pren_exposome_scaled$cohort <- cohorts
post_exposome_scaled$cohort <- cohorts

####********************####
#### 2. Run regressions ####
####********************####

# 2a prenatal: 
pren_exposome_scaled$id <- row.names(pren_exposome_scaled)
pren_exposome_long <- pren_exposome_scaled %>% 
  pivot_longer(
    cols = starts_with("h"),
    names_to = "exposure_name",
    values_to = "exposure_value") 

pren_exposome_res <- pren_exposome_long %>%
  group_by(exposure_name) %>% 
  do(cbind(id = .$id, lm(exposure_value ~ cohort, data = .) %>% augment)) %>%
  select(id, exposure_name, .resid) %>%
  rename(residual = .resid) %>%
  pivot_wider(names_from = exposure_name, values_from = residual) %>%
  column_to_rownames(var = "id") %>%
  as.data.frame()

# 2b postnatal:
post_exposome_scaled$id <- row.names(post_exposome_scaled)
post_exposome_long <- post_exposome_scaled %>% 
  pivot_longer(
    cols = starts_with("h"),
    names_to = "exposure_name",
    values_to = "exposure_value") 

post_exposome_res <- post_exposome_long %>%
  group_by(exposure_name) %>% 
  do(cbind(id = .$id, lm(exposure_value ~ cohort, data = .) %>% augment)) %>%
  select(id, exposure_name, .resid) %>%
  rename(residual = .resid) %>%
  pivot_wider(names_from = exposure_name, values_from = residual) %>%
  column_to_rownames(var = "id") %>%
  as.data.frame()

####*******************####
#### 3. Save residuals ####
####*******************####
save(codebook, pren_description, post_description, pren_exposome_res, post_exposome_res,
     file = paste0(generated.data.folder, "exposome_pren_postnatal_residuals.RData"))