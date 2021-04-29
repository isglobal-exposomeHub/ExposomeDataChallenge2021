##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

##load data
####exposome
load("data_analysis/mediation_analysis/data_preparation/exposome_expression_data")
load("data_analysis/mediation_analysis/data_preparation/exposome_sample_info")
load("data_analysis/mediation_analysis/data_preparation/exposome_variable_info")

# ####internal_ome
# load("data_analysis/mediation_analysis/data_preparation/internal_ome_expression_data")
# load("data_analysis/mediation_analysis/data_preparation/internal_ome_sample_info")
# load("data_analysis/mediation_analysis/data_preparation/internal_ome_variable_info")

#####use the linear mixed model to find which exposome have association with phenotye
library(lme4)
library(lmerTest)

# ####Behavior
# behavior_glm_p = 
# exposome_expression_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(function(x) {
#     temp_data =
#       data.frame(x = x,
#                  exposome_sample_info)
#     glm_reg =
#       glm(
#         Behavior ~ x + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
#           Maternal.age + Child.height + Child.weight + Birthweight,
#         family = gaussian,
#         temp_data
#       )
#     
#     temp =
#       summary(glm_reg)$coefficients %>%
#       as.data.frame()
#     temp$`Pr(>|t|)`[2]
#   }) %>% 
#   unlist()
# 
# 
# 
# 
# ####Body.mass.index.z.score
# Body.mass.index.z.score_glm_p = 
#   exposome_expression_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(function(x) {
#     temp_data =
#       data.frame(x = x,
#                  exposome_sample_info)
#     glm_reg =
#       glm(
#         Body.mass.index.z.score ~ x + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
#           Maternal.age + Child.height + Child.weight + Birthweight,
#         family = gaussian,
#         temp_data
#       )
#     
#     temp =
#       summary(glm_reg)$coefficients %>%
#       as.data.frame()
#     temp$`Pr(>|t|)`[2]
#   }) %>% 
#   unlist()
# 
# 
# 
# ####Intelligence.quotient
# Intelligence.quotient_glm_p = 
#   exposome_expression_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(function(x) {
#     temp_data =
#       data.frame(x = x,
#                  exposome_sample_info)
#     glm_reg =
#       glm(
#         Intelligence.quotient ~ x + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
#           Maternal.age + Child.height + Child.weight + Birthweight,
#         family = gaussian,
#         temp_data
#       )
#     
#     temp =
#       summary(glm_reg)$coefficients %>%
#       as.data.frame()
#     temp$`Pr(>|t|)`[2]
#   }) %>% 
#   unlist()
# 
# exposome_phenotype_glm = 
# data.frame(behavior_glm_p, Body.mass.index.z.score_glm_p, Intelligence.quotient_glm_p) %>% 
#   tibble::rownames_to_column(var = "variable_id")
# 
# exposome_phenotype_glm$behavior_glm_p.adjust = p.adjust(exposome_phenotype_glm$behavior_glm_p, method = "BH")
# exposome_phenotype_glm$Body.mass.index.z.score_glm_p.adjust = p.adjust(exposome_phenotype_glm$Body.mass.index.z.score_glm_p, method = "BH")
# exposome_phenotype_glm$Intelligence.quotient_glm_p.adjust = p.adjust(exposome_phenotype_glm$Intelligence.quotient_glm_p, method = "BH")
# 
# sum(temp_data$behavior_glm_p.adjust < 0.05)
save(exposome_phenotype_glm, file = "exposome_phenotype_glm")

exposome_phenotype_glm$variable_id[which(exposome_phenotype_glm$behavior_glm_p.adjust < 0.05)]
exposome_phenotype_glm$variable_id[which(exposome_phenotype_glm$Body.mass.index.z.score_glm_p.adjust < 0.05)]
exposome_phenotype_glm$variable_id[which(exposome_phenotype_glm$Intelligence.quotient_glm_p.adjust < 0.05)]


# ###pie to show the how many exposome variables are associated with outcome
# data.frame(class = c("significant", "no"),
#            number = c(
#              sum(exposome_phenotype_glm$behavior_glm_p.adjust < 0.05),
#              nrow(exposome_variable_info) - sum(exposome_phenotype_glm$behavior_glm_p.adjust < 0.05)
#            )) %>%
#   ggplot(aes(x = 2, y = number, fill = class)) +
#   geom_bar(stat = "identity") +
#   coord_polar("y", start = 200)
#   theme_void() +
#   scale_fill_manual(values = c(
#     significant
#   )) +
#   xlim(.2, 2.5)




















