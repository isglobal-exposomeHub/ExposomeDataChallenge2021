tinyTools::setwd_project()
setwd("data/")
library(tidyverse)
rm(list = ls())

load("exposome.RData")

###
codebook = 
codebook %>% 
  apply(2, as.character) %>% 
  as.data.frame()
codebook$labels[codebook$labels == "DMDTP" & codebook$period == "Postnatal" & codebook$labelsshort == "DMTP"] = "DMTP"
table(codebook$domain)

#####get the sample information of all the mothers and children
dim(covariates)
dim(phenotype)

sample_info = 
  covariates %>% 
  dplyr::left_join(phenotype, by = "ID")

colnames(sample_info) = as.character(codebook$labels[match(colnames(sample_info), codebook$variable_name)])
colnames(sample_info)[1] = "subject_id"

sample_info = 
  rbind(
    data.frame(
      sample_id = paste("mother", sample_info$subject_id, sep = "_"),
      sample_info
    ),
    data.frame(
      sample_id = paste("child", sample_info$subject_id, sep = "_"),
      sample_info
    )
  ) %>% 
  dplyr::arrange(subject_id, sample_id)


#######exposome chemical
tinyTools::setwd_project()
setwd("data_analysis/exposome_chemical_data_analysis/data_preparation/")

####Covariates
codebook1 = 
  codebook %>% 
  dplyr::filter(domain == "Chemicals") %>% 
  dplyr::filter(var_type == "numeric")

codebook1$var_type
codebook1 %>% 
  dplyr::filter(var_type == "factor")

variable_name1 = as.character(codebook1$labels[codebook1$period == "Pregnancy"])
variable_name2 = as.character(codebook1$labels[codebook1$period == "Postnatal"])

variable_name = intersect(variable_name1, variable_name2)

variable1 = 
codebook1 %>% 
  dplyr::filter(period == "Pregnancy") 

variable1 = 
variable1[match(variable_name, variable1$labels),]


variable2 = 
  codebook1 %>% 
  dplyr::filter(period == "Postnatal") 

variable2 = 
  variable2[match(variable_name, variable2$labels),]


###mother data
expression_data1 = 
  exposome[,variable1$variable_name]

rownames(expression_data1) = paste("mother", exposome$ID, sep = "_")
colnames(expression_data1) = variable1$labels

###child data
expression_data2 = 
  exposome[,variable2$variable_name]

rownames(expression_data2) = paste("child", exposome$ID, sep = "_")
colnames(expression_data2) = variable2$labels


expression_data = 
rbind(expression_data1, expression_data2)

expression_data = 
expression_data %>% 
  t() %>% 
  as.data.frame()

variable_info = 
  codebook[match(rownames(expression_data), codebook$labels),] %>% 
  dplyr::mutate(variable_id = labels) %>% 
  dplyr::select(variable_id, everything()) %>% 
  dplyr::mutate(description = stringr::str_replace_all(description, " in child", "")) %>% 
  dplyr::mutate(description = stringr::str_replace_all(description, " in mother", ""))



expression_data = expression_data[,sample_info$sample_id]

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(codebook, file = "codebood")



table(codebook$domain)

#######exposome Indoor air
tinyTools::setwd_project()
setwd("data_analysis/exposome_air_data_analysis/data_preparation/")

####Covariates
codebook1 = 
  codebook %>% 
  dplyr::filter(domain == "Indoor air") %>% 
  dplyr::filter(var_type == "numeric")

variable_name1 = as.character(codebook1$labels[codebook1$period == "Pregnancy"])
variable_name2 = as.character(codebook1$labels[codebook1$period == "Postnatal"])

variable_name = variable_name2

variable1 = 
  codebook1 %>% 
  dplyr::filter(period == "Pregnancy") 

variable1 = 
  variable1[match(variable_name, variable1$labels),]


variable2 = 
  codebook1 %>% 
  dplyr::filter(period == "Postnatal") 

variable2 = 
  variable2[match(variable_name, variable2$labels),]


###mother data
expression_data1 = 
  exposome[,variable1$variable_name]

rownames(expression_data1) = paste("mother", exposome$ID, sep = "_")
colnames(expression_data1) = variable1$labels

###child data
expression_data2 = 
  exposome[,variable2$variable_name]

rownames(expression_data2) = paste("child", exposome$ID, sep = "_")
colnames(expression_data2) = variable2$labels

expression_data = 
  expression_data2

expression_data = 
  expression_data %>% 
  t() %>% 
  as.data.frame()

variable_info = 
  codebook[match(rownames(expression_data), codebook$labels),] %>% 
  dplyr::mutate(variable_id = labels) %>% 
  dplyr::select(variable_id, everything()) %>% 
  dplyr::mutate(description = stringr::str_replace_all(description, " in child", "")) %>% 
  dplyr::mutate(description = stringr::str_replace_all(description, " in mother", ""))

expression_data = expression_data[,sample_info$sample_id]

load("../../exposome_chemical_data_analysis/data_preparation/sample_info")

sample_info = sample_info[match(colnames(expression_data), sample_info$sample_id),]

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(codebook, file = "codebood")










table(codebook$domain)

#######Outdoor exposures
tinyTools::setwd_project()
setwd("data_analysis/exposome_outdoor_data_analysis/data_preparation/")

####Covariates
codebook1 = 
  codebook %>% 
  dplyr::filter(domain == "Outdoor exposures") %>% 
  dplyr::filter(var_type == "numeric")

codebook1 %>% 
  dplyr::filter(subfamily == "Humidity")

codebook1 %>% 
  dplyr::filter(labels %in% c("T", 'T(day)', "T(week)", "T(month)"))

codebook1 %>%
  dplyr::filter(period == "Postnatal" &
                  description == "Relative humidityone day before at home")

codebook1$labels[codebook1$period == "Postnatal" & codebook1$description == "Relative humidityone day before at home"] = "Hum.(day)"

variable_name1 = as.character(codebook1$labels[codebook1$period == "Pregnancy"])
variable_name2 = as.character(codebook1$labels[codebook1$period == "Postnatal"])

sum(duplicated(variable_name1))
sum(duplicated(variable_name2))

variable_name = unique(c(variable_name1, variable_name2))

variable1 = 
  codebook1 %>% 
  dplyr::filter(period == "Pregnancy") 

variable1 = 
  variable1[match(variable_name, variable1$labels),]

variable2 = 
  codebook1 %>% 
  dplyr::filter(period == "Postnatal") 

variable2 = 
  variable2[match(variable_name, variable2$labels),]

name1 = variable1$variable_name
name2 = variable2$variable_name

name1[is.na(variable1$variable_name)] = 
  name2[is.na(variable1$variable_name)]

name2[is.na(variable2$variable_name)] = 
  name1[is.na(variable2$variable_name)]


col_name1 = variable1$labels
col_name2 = variable2$labels

col_name1[is.na(col_name1)] = col_name2[is.na(col_name1)]

###mother data
expression_data1 = 
  exposome[,name1]

expression_data1[,is.na(variable1$variable_name)] = NA

rownames(expression_data1) = paste("mother", exposome$ID, sep = "_")
colnames(expression_data1) = col_name1

###child data
expression_data2 = 
  exposome[,name2]

expression_data2[,is.na(variable2$variable_name)] = NA

rownames(expression_data2) = paste("child", exposome$ID, sep = "_")
colnames(expression_data2) = col_name1

expression_data = 
  rbind(expression_data1,
        expression_data2)

expression_data = 
  expression_data %>% 
  t() %>% 
  as.data.frame()

variable_info = 
  codebook1[match(rownames(expression_data), codebook1$labels),] %>% 
  dplyr::mutate(variable_id = labels) %>% 
  dplyr::select(variable_id, everything()) %>% 
  dplyr::mutate(description = stringr::str_replace_all(description, " in child", "")) %>% 
  dplyr::mutate(description = stringr::str_replace_all(description, " in mother", ""))

expression_data = expression_data[,sample_info$sample_id]

load("../../exposome_chemical_data_analysis/data_preparation/sample_info")

sample_info = sample_info[match(colnames(expression_data), sample_info$sample_id),]

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(codebook, file = "codebood")







table(codebook$domain)

#######Lifestyles
tinyTools::setwd_project()
setwd("data_analysis/lifecycle_outdoor_data_analysis/data_preparation/")

####Lifestyles
codebook1 = 
  codebook %>% 
  dplyr::filter(domain == "Lifestyles")

variable_name1 = as.character(codebook1$labels[codebook1$period == "Pregnancy"])
variable_name2 = as.character(codebook1$labels[codebook1$period == "Postnatal"])

sum(duplicated(variable_name1))
sum(duplicated(variable_name2))

variable_name = unique(c(variable_name1, variable_name2))

variable1 = 
  codebook1 %>% 
  dplyr::filter(period == "Pregnancy") 

variable1 = 
  variable1[match(variable_name, variable1$labels),]

variable2 = 
  codebook1 %>% 
  dplyr::filter(period == "Postnatal") 

variable2 = 
  variable2[match(variable_name, variable2$labels),]

name1 = variable1$variable_name
name2 = variable2$variable_name

name1[is.na(variable1$variable_name)] = 
  name2[is.na(variable1$variable_name)]

name2[is.na(variable2$variable_name)] = 
  name1[is.na(variable2$variable_name)]

col_name1 = variable1$labels
col_name2 = variable2$labels

col_name1[is.na(col_name1)] = col_name2[is.na(col_name1)]

###mother data
expression_data1 = 
  exposome[,name1]

expression_data1[,is.na(variable1$variable_name)] = NA

rownames(expression_data1) = paste("mother", exposome$ID, sep = "_")
colnames(expression_data1) = col_name1

###child data
expression_data2 = 
  exposome[,name2]

expression_data2[,is.na(variable2$variable_name)] = NA

rownames(expression_data2) = paste("child", exposome$ID, sep = "_")
colnames(expression_data2) = col_name1

expression_data = 
  rbind(expression_data1,
        expression_data2)

expression_data = 
  expression_data %>% 
  t() %>% 
  as.data.frame()

variable_info = 
  codebook1[match(rownames(expression_data), codebook1$labels),] %>% 
  dplyr::mutate(variable_id = labels) %>% 
  dplyr::select(variable_id, everything()) %>% 
  dplyr::mutate(description = stringr::str_replace_all(description, " in child", "")) %>% 
  dplyr::mutate(description = stringr::str_replace_all(description, " in mother", ""))

expression_data = expression_data[,sample_info$sample_id]

load("../../exposome_chemical_data_analysis/data_preparation/sample_info")

sample_info = sample_info[match(colnames(expression_data), sample_info$sample_id),]

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(codebook, file = "codebood")







