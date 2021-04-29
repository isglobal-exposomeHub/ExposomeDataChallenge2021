rm(list=ls())
library(tidyverse)
library(readxl)

# load data
codebook <- read_xlsx("codebook.xlsx")
row.names(codebook) <- NULL
load("exposome_v2.RData")
ls()



# remove covariates and factor variables
codebook2 <- codebook %>% 
  filter(var_type=="numeric" & domain!="Covariates" & domain!="Phenotype") %>%
  filter(variable_name != "hs_sumPCBs5_madj_Log2" & variable_name != "hs_sumPCBs5_cadj_Log2") %>%   # this is a sum of other variables
  filter(variable_name != "hs_sumDEHP_madj_Log2" & variable_name != "hs_sumDEHP_cadj_Log2") %>%  # this is a sum of other variables
  select(-c("description"))
dim(codebook2)
head(codebook2)


# define indicies by unique levels of:
# domain, family, period, period_postnatal

index_identifiers <- codebook2 %>% 
  filter(period_postnatal!="Day before examination" | is.na(period_postnatal)) %>%  # remove these so that exposures do not have more than one postnatal period time
  filter(period_postnatal!="Week before examination" | is.na(period_postnatal)) %>%   # remove these so that exposures do not have more than one postnatal period time
  select(-period_postnatal) %>%
  group_by(domain, family, period) %>% 
  filter(family!="Tobacco Smoke") # included in covariates
  

index_identifiers %>% 
  group_by(family, period) %>% 
  mutate(ntimes=n()) %>%
  as.data.frame()

#----------------------
# Data for BMIM
#----------------------

dim(index_identifiers)
index_names <- index_identifiers %>% 
  select(family, period) %>%
  distinct()


index_exposure_list <- list()
for(i in 1:nrow(index_names)){
  name <- paste(str_replace_all(as.character(index_names$family[i])," ", ""),
                as.character(index_names$period[i]),
                sep="_")
  index_identifiers_temp <- index_identifiers %>% filter(family ==index_names$family[i] & period ==index_names$period[i]) %>% as.data.frame()
  
  # add index to list
  index_exposure_list[[name]] <- as.matrix(exposome[,as.character(index_identifiers_temp$variable_name)])
}

# make list of variables and indices
index_names$index_name <- names(index_exposure_list)
index_identifiers <- left_join(index_names,index_identifiers,by=c("family","domain","period"))
colnames(index_identifiers)
dim(index_identifiers)

# save data
save(index_identifiers, index_names, index_exposure_list,
     file="DataPrep/index_exposure_list.rda")

#----------------------
# Data for TDLMM
#----------------------

index_identifiers$labels <- str_replace(index_identifiers$labels, "\\(year\\)","")
index_identifiers$labels <- str_replace(index_identifiers$labels, "\\(month\\)","")

# make dataset with repeated exposures for TDLMM
time_ids <- index_identifiers %>% 
  group_by(labels) %>% 
  mutate(ntimes=n()) %>% 
  filter(ntimes==2)

time_exposure_list <- list()
for(i in unique(time_ids$labels)){
  time_ids_temp_preg  <- time_ids %>% filter(labels ==i & period=="Pregnancy") %>% as.data.frame()
  time_ids_temp_post  <- time_ids %>% filter(labels ==i & period=="Postnatal") %>% as.data.frame()
  time_exposure_list[[i]] <- as.matrix(exposome[,c(as.character(time_ids_temp_preg$variable_name),
                                                   as.character(time_ids_temp_post$variable_name))])
}

save(time_exposure_list, file="DataPrep/time_exposure_list.rda")
