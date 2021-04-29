##avoid source
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

##load data
####exposome
###child exposome chemical
load(
  "data_analysis/exposome_chemical_data_analysis/data_preparation/expression_data"
)
load("data_analysis/exposome_chemical_data_analysis/data_preparation/sample_info")
load("data_analysis/exposome_chemical_data_analysis/data_preparation/variable_info")

exposome_chemical_variable_info <-
  variable_info

exposome_chemical_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

exposome_chemical_expression_data =
  expression_data[, exposome_chemical_sample_info$sample_id]

exposome_chemical_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  })

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

exposome_outdoor_expression_data = exposome_outdoor_expression_data[remain_idx,]
exposome_outdoor_variable_info = exposome_outdoor_variable_info[remain_idx,]

dim(exposome_chemical_expression_data)
dim(exposome_air_expression_data)
dim(exposome_outdoor_expression_data)

sum(colnames(exposome_chemical_expression_data) == colnames(exposome_air_expression_data))
sum(colnames(exposome_chemical_expression_data) == colnames(exposome_outdoor_expression_data))

exposome_expression_data = rbind(exposome_chemical_expression_data, 
                                 exposome_air_expression_data,
                                 exposome_outdoor_expression_data)

exposome_variable_info = 
  rbind(
    data.frame(exposome_chemical_variable_info,
               class = "chemical"),
    data.frame(exposome_air_variable_info,
               class = "air"),
    data.frame(exposome_outdoor_variable_info,
               class = "outdoor")
  )

rownames(exposome_expression_data) == exposome_variable_info$variable_id

exposome_sample_info = exposome_air_sample_info

colnames(exposome_expression_data) == exposome_sample_info$sample_id





#####internal omics
###transcriptome
load(
  "data_analysis/transcriptome_data_analysiss/data_preparation/expression_data"
)
load("data_analysis/transcriptome_data_analysiss/data_preparation/sample_info")
load("data_analysis/transcriptome_data_analysiss/data_preparation/variable_info")

transcriptome_variable_info <-
  variable_info

transcriptome_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

transcriptome_expression_data =
  expression_data[, transcriptome_sample_info$sample_id]

transcriptome_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  })

###proteome
load(
  "data_analysis/proteome_data_analysiss/data_preparation/expression_data"
)
load("data_analysis/proteome_data_analysiss/data_preparation/sample_info")
load("data_analysis/proteome_data_analysiss/data_preparation/variable_info")

proteome_variable_info <-
  variable_info

proteome_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

proteome_expression_data =
  expression_data[, proteome_sample_info$sample_id]

proteome_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  })


###serum_metabolome
load(
  "data_analysis/serum_metabolome_data_analysiss/data_preparation/expression_data"
)
load("data_analysis/serum_metabolome_data_analysiss/data_preparation/sample_info")
load("data_analysis/serum_metabolome_data_analysiss/data_preparation/variable_info")

serum_metabolome_variable_info <-
  variable_info

serum_metabolome_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

serum_metabolome_expression_data =
  expression_data[, serum_metabolome_sample_info$sample_id]

serum_metabolome_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  })



###urine_metabolome
load(
  "data_analysis/urine_metabolome_data_analysiss/data_preparation/expression_data"
)
load("data_analysis/urine_metabolome_data_analysiss/data_preparation/sample_info")
load("data_analysis/urine_metabolome_data_analysiss/data_preparation/variable_info")

urine_metabolome_variable_info <-
  variable_info

urine_metabolome_sample_info =
  sample_info %>%
  dplyr::filter(stringr::str_detect(sample_id, pattern = "child"))

urine_metabolome_expression_data =
  expression_data[, urine_metabolome_sample_info$sample_id]

urine_metabolome_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  })

dim(transcriptome_expression_data)
dim(proteome_expression_data)
dim(serum_metabolome_expression_data)
dim(urine_metabolome_expression_data)


intersect_name = 
  Reduce(f = intersect, x = list(colnames(transcriptome_expression_data),
                                 colnames(proteome_expression_data),
                                 colnames(serum_metabolome_expression_data),
                                 colnames(urine_metabolome_expression_data)))

transcriptome_expression_data = 
  transcriptome_expression_data[,intersect_name]
transcriptome_sample_info =
  transcriptome_sample_info[match(intersect_name, transcriptome_sample_info$sample_id),]

proteome_expression_data = 
  proteome_expression_data[,intersect_name]
proteome_sample_info =
  proteome_sample_info[match(intersect_name, proteome_sample_info$sample_id),]

serum_metabolome_expression_data = 
  serum_metabolome_expression_data[,intersect_name]
serum_metabolome_sample_info =
  serum_metabolome_sample_info[match(intersect_name, serum_metabolome_sample_info$sample_id),]

urine_metabolome_expression_data = 
  urine_metabolome_expression_data[,intersect_name]
urine_metabolome_sample_info =
  urine_metabolome_sample_info[match(intersect_name, urine_metabolome_sample_info$sample_id),]

internal_ome_expression_data = rbind(transcriptome_expression_data,
                                     proteome_expression_data,
                                     serum_metabolome_expression_data,
                                     urine_metabolome_expression_data)
transcriptome_variable_info$class = "transcriptome"
proteome_variable_info$class = "proteome"
serum_metabolome_variable_info$class = "serum_metabolome"
urine_metabolome_variable_info$class = "urine_metabolome"
serum_metabolome_variable_info$CHEBI = as.character(serum_metabolome_variable_info$CHEBI)
internal_ome_variable_info = 
  dplyr::full_join(x = transcriptome_variable_info,
                   y = proteome_variable_info,
                   by = intersect(colnames(transcriptome_variable_info),
                                  colnames(proteome_variable_info))) %>% 
  dplyr::full_join(serum_metabolome_variable_info, 
                   by = intersect(colnames(serum_metabolome_variable_info),
                                  colnames(.))) %>% 
  dplyr::full_join(urine_metabolome_variable_info, 
                   by = intersect(colnames(urine_metabolome_variable_info),
                                  colnames(.)))

sum(rownames(internal_ome_expression_data) == internal_ome_variable_info$variable_id)

sum(duplicated(internal_ome_variable_info$variable_id))

head(transcriptome_variable_info)
head(proteome_variable_info)
head(serum_metabolome_variable_info)
head(urine_metabolome_variable_info)

internal_ome_variable_info %>% 
  dplyr::mutate(mol_name = 
                  case_when(
                    class == "proteome" ~ Gene_Symbol,
                    class == "serum_metabolome" ~ var,
                    class == "urine_metabolome" ~ var,
                    class == "transcriptome" ~ GeneSymbolDB2,
                    TRUE ~ NA
                  ))





internal_ome_variable_info = 
  rbind(
    data.frame(transcriptome_variable_info,
               class = "transcriptome"),
    data.frame(proteome_variable_info,
               class = "proteome"),
    data.frame(serum_metabolome_variable_info,
               class = "serum_metabolome"),
    data.frame(urine_metabolome_variable_info,
               class = "urine_metabolome")
  )

rownames(exposome_expression_data) == exposome_variable_info$variable_id

exposome_sample_info = exposome_air_sample_info

colnames(exposome_expression_data) == exposome_sample_info$sample_id



tinyTools::setwd_project()
dir.create("data_analysis/mediation_analysis")
setwd("data_analysis/mediation_analysis/")

save(exposome_expression_data, file = "exposome_expression_data")
save(exposome_sample_info, file = "exposome_sample_info")
save(exposome_variable_info, file = "exposome_variable_info")
