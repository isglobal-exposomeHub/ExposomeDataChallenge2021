tinyTools::setwd_project()
setwd("data/")
library(tidyverse)
rm(list = ls())

load("metabol_serum.Rdata")

tinyTools::setwd_project()
setwd("data_analysis/serum_metabolome_data_analysiss/data_preparation")

expression_data = 
metabol_serum@assayData$exprs

sample_info =
  metabol_serum@phenoData@data %>%
  as.data.frame() %>%
  dplyr::rename(sample_id = ID)

variable_info =
  metabol_serum@featureData@data %>% 
  tibble::rownames_to_column(var = "variable_id")

rownames(variable_info) = NULL


sample_info = 
  sample_info %>% 
  dplyr::mutate(subject_id = sample_id) %>% 
  dplyr::mutate(sample_id = paste("child", sample_info$sample_id, sep = "_")) %>% 
  dplyr::select(sample_id, subject_id, everything())

colnames(expression_data) = sample_info$sample_id

save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")
save(sample_info, file = "sample_info")








