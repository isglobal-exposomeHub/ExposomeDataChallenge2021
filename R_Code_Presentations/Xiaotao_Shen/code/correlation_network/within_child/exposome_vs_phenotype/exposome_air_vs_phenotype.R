##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())
source("code/tools.R")

##load data
###child exposome air
load(
  "data_analysis/exposome_air_data_analysis/data_preparation/expression_data"
)
load("data_analysis/exposome_air_data_analysis/data_preparation/sample_info")
load("data_analysis/exposome_air_data_analysis/data_preparation/variable_info")

exposome_air_variable_info <-
  variable_info

exposome_air_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

exposome_air_expression_data =
  expression_data[, exposome_air_sample_info$sample_id]

exposome_air_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  })

tinyTools::setwd_project()
dir.create("data_analysis/correlation_network/within_child/exposome_air_vs_phenotype")
setwd("data_analysis/correlation_network/within_child/exposome_air_vs_phenotype")

#####glm
exposome_air_phenotype_glm = 
exposome_air_sample_info[,c("Behavior", "Body.mass.index.z.score", "Intelligence.quotient")] %>% 
  purrr::map(function(x){
    temp_p = 
    exposome_air_expression_data %>% 
      t() %>% 
      as.data.frame() %>% 
      purrr::map(function(y){
        #x is the phenotype
        #y is the exposome
            temp_data =
              data.frame(x = x,
                         y = y,
                         exposome_air_sample_info)
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

save(exposome_air_phenotype_glm, file = "exposome_air_phenotype_glm")
