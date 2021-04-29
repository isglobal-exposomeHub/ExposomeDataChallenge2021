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

####internal_ome
load("data_analysis/mediation_analysis/data_preparation/internal_ome_expression_data")
load("data_analysis/mediation_analysis/data_preparation/internal_ome_sample_info")
load("data_analysis/mediation_analysis/data_preparation/internal_ome_variable_info")

#####use the linear mixed model to find which exposome have association with phenotye
library(lme4)
library(lmerTest)

dir.create("data_analysis/mediation_analysis/internal_ome_phenotype")
setwd("data_analysis/mediation_analysis/internal_ome_phenotype")

intersect_name = intersect(colnames(exposome_expression_data),
                           colnames(internal_ome_expression_data))

exposome_expression_data = exposome_expression_data[,intersect_name]
internal_ome_expression_data = internal_ome_expression_data[,intersect_name]

exposome_sample_info = 
  exposome_sample_info[match(intersect_name, exposome_sample_info$sample_id),]

internal_ome_sample_info = 
  internal_ome_sample_info[match(intersect_name, internal_ome_sample_info$sample_id),]

library(future)
library(furrr)

internal_phenotype_ome_glm_p =
  exposome_sample_info[, c("Behavior",
                           "Body.mass.index.z.score",
                           "Intelligence.quotient")] %>%
  purrr::map(function(x) {
    cat(x[1], " ")
    temp_p =
      internal_ome_expression_data %>%
      t %>%
      as.data.frame() %>%
      purrr::map(function(y) {
        #x is the phenotype
        #y is the internal ome
        temp_data =
          data.frame(x = x,
                     y = y,
                     exposome_sample_info)
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
    p.adjust(temp_p, method = "BH")
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

internal_phenotype_ome_glm_p =
  internal_phenotype_ome_glm_p %>% 
  t() %>% 
  as.data.frame()


internal_phenotype_ome_glm_p =
internal_ome_expression_data %>%
  t() %>%
  as.data.frame() %>%
  furrr::future_map(function(x) {
    temp_p = 
      exposome_sample_info[,c("Behavior", "Body.mass.index.z.score", "Intelligence.quotient")] %>% 
      purrr::map(function(y){
        temp_data =
          data.frame(x = x,
                     y = y,
                     exposome_sample_info)
        glm_reg =
          glm(
            y ~ x + Child.sex + Year.of.birth + Maternal.BMI + Gestational.age.at.birth +
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
    # p.adjust(temp_p, method = "BH")
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()



sum(temp_data$behavior_glm_p.adjust < 0.05)
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




















