##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())
source("code/tools.R")

##load data
###child exposome outdoor
load(
  "data_analysis/exposome_outdoor_data_analysis/data_preparation/expression_data"
)
load("data_analysis/exposome_outdoor_data_analysis/data_preparation/sample_info")
load("data_analysis/exposome_outdoor_data_analysis/data_preparation/variable_info")

exposome_outdoor_variable_info <-
  variable_info

exposome_outdoor_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

exposome_outdoor_expression_data =
  expression_data[, exposome_outdoor_sample_info$sample_id]

remain_idx = 
exposome_outdoor_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  }) %>% 
  `==`(0) %>% 
  which()

exposome_outdoor_variable_info = 
  exposome_outdoor_variable_info[remain_idx,]

exposome_outdoor_expression_data =
  exposome_outdoor_expression_data[remain_idx, ]

tinyTools::setwd_project()
dir.create("data_analysis/correlation_network/within_child/exposome_outdoor_vs_phenotype")
setwd("data_analysis/correlation_network/within_child/exposome_outdoor_vs_phenotype")

#####glm
exposome_outdoor_phenotype_glm = 
exposome_outdoor_sample_info[,c("Behavior", "Body.mass.index.z.score", "Intelligence.quotient")] %>% 
  purrr::map(function(x){
    temp_p = 
    exposome_outdoor_expression_data %>% 
      t() %>% 
      as.data.frame() %>% 
      purrr::map(function(y){
        #x is the phenotype
        #y is the exposome
            temp_data =
              data.frame(x = x,
                         y = y,
                         exposome_outdoor_sample_info)
            glm_reg =
              glm(
                x ~ y + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
                  Maternal.age + Child.height + Child.weight + Birthweight,
                family = gaussian,
                temp_data
              )
            
            temp =
              summary(glm_reg)$coefficients %>%
              as.data.frame()
            temp$`Pr(>|t|)`[2]
      }) %>% 
      unlist()
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

exposome_outdoor_phenotype_glm

save(exposome_outdoor_phenotype_glm, file = "exposome_outdoor_phenotype_glm")
