# File: A_02_generate_raw_reports.R
# Author(s): Lawrence Chillrud
# Date since last edit: 4/27/21

####*******************####
#### Table of Contents ####
####*******************####
####* Notes / Description
####* 0. Package imports
####* 1. Read in data
####* 2. Define function to make raw report
####* 3. Generate reports

####*********************####
#### Notes / Description ####
####*********************####
####* From the processed ExposomeChallenge data, we generate
####* exploratory reports (in the form of .html files) of the 
####* pre- and postnatal exposome matrices, so we can compare the original
####* data to data we later recover from PCP.

####********************####
#### 0. Package Imports ####
####********************####

print("Running script A_02...")

library(pcpr)
library(PCPhelpers)
library(tidyverse)
library(ggplot2)

####*****************####
#### 1. Read in data ####
####*****************####
load(paste0(generated.data.folder, "exposome_pren_postnatal.RData"))

# 1a scale data mats by column
mat_pre <- pren_exposome %>% as.matrix() %>% scale(., center = F)
mat_post <- post_exposome %>% as.matrix() %>% scale(., center = F)

# 1b create metadata for each mat
colnames(mat_pre) <- pren_description$labels
colgroupings_pre <- pren_description$family

colnames(mat_post) <- post_description$labels
colgroupings_post <- post_description$family

rowgroupings <- covariates$h_cohort

####***************************************####
#### 2. Define function to make raw report ####
####***************************************####
make_raw_report <- function(dataset, subname, mat, scale_flag = F,
                            rowvar_name = "ID", rowvar = NULL, 
                            rowgroupings = NULL, rowgroupname = "cohort",
                            colnames = NULL, colgroupings = NULL) {
  
  filename <- paste0(paste(subname, "raw_report", sep = "_"), ".html") 
  
  rmarkdown::render(
    paste0(code.folder, "F_raw_report.Rmd"),
    params = list(
      dataset = dataset,
      subname = subname,
      mat = mat,
      scale_flag = scale_flag,
      rowvar_name = rowvar_name,
      rowvar = rowvar,
      rowgroupings = rowgroupings,
      rowgroupname = rowgroupname,
      colnames = colnames,
      colgroupings = colgroupings
    ),
    output_file = paste0(output.folder, filename)
  )
}

####*********************####
#### 3. Generate reports ####
####*********************####
make_raw_report(dataset = "exposome", subname = "prenatal", mat = mat_pre, 
                rowvar = as.numeric(rownames(mat_pre)), rowgroupings = rowgroupings, colgroupings = colgroupings_pre)

make_raw_report(dataset = "exposome", subname = "postnatal", mat = mat_post, 
                rowvar = as.numeric(rownames(mat_post)), rowgroupings = rowgroupings, colgroupings = colgroupings_post)

