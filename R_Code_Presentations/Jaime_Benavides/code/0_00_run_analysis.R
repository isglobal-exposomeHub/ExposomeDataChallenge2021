# File: 0_00_run_analysis.R
# Author(s): Lawrence Chillrud
# Date since last edit: 4/27/21

####*******************####
#### Table of Contents ####
####*******************####
####* Notes / Description
####* 
####* 0.00: Package Imports
####* 0.01: Initialize Directory Structure
####* 
####* A: Project Set Up
####* A.01: Prepare Exposome Data
####* A.02: Generate Exploratory Raw Data Reports
####* A.03: Initial PCP runs
####* A.04: Regress exposure on cohort to lessen influence of cohort
####* 
####* B: Main PCP runs
####* B.01: Conduct PCP gridsearches & render reports on residuals
####* 
####* C: Run Health Models
####* C.01: GLM Model
####* C.02: GAM Model
####* C.03: Lasso Model

####*********************####
#### Notes / Description ####
####*********************####
####* This script runs all code needed for our 
####* submission to the 2021 ExposomeDataChallenge.

####***********************####
#### 0.00: Package Imports #### 
####***********************####

####**************************************####
#### 0.01: Initialize Directory Structure #### 
####**************************************####
source(here::here("code/0_01_init_directory_structure.R"))

####*******************####
#### A: Project Set Up #### 
####*******************####

####*****************************####
#### A.01: Prepare Exposome Data #### 
####*****************************####
source(paste0(code.folder, "A_01_prepare_exposome_data.R"))

rm(list=ls()[!ls() %in% c(folder.names, "folder.names")]) # clean env.

####*********************************************####
#### A.02: Generate Exploratory Raw Data Reports #### 
####*********************************************####
source(paste0(code.folder, "A_02_generate_raw_reports.R"))

rm(list=ls()[!ls() %in% c(folder.names, "folder.names")]) # clean env.

####************************####
#### A.03: Initial PCP runs #### 
####************************####
source(paste0(code.folder, "A_03_init_pcp_runs.R"))

rm(list=ls()[!ls() %in% c(folder.names, "folder.names")]) # clean env.

####****************************************************************####
#### A.04: Regress exposure on cohort to lessen influence of cohort #### 
####****************************************************************####
source(paste0(code.folder, "A_04_regress_cohort.R"))

rm(list=ls()[!ls() %in% c(folder.names, "folder.names")]) # clean env.

####******************####
#### B: Main PCP runs #### 
####******************####

####**************************************************************####
#### B.01: Conduct PCP gridsearches & render reports on residuals #### 
####**************************************************************####
source(paste0(code.folder, "B_01_main_PCP_runs_on_res.R"))

rm(list=ls()[!ls() %in% c(folder.names, "folder.names")]) # clean env.

####**********************####
#### C: Run Health Models #### 
####**********************####

####*****************####
#### C.01: GLM Model #### 
####*****************####

print("Running script C_01...")
rmarkdown::render(paste0(code.folder, "C_01_GLM_model.Rmd"), output_file = paste0(output.folder, "GLM.html"))

####*****************####
#### C.02: GAM Model #### 
####*****************####

print("Running script C_02...")
rmarkdown::render(paste0(code.folder, "C_02_GAM_model.Rmd"), output_file = paste0(output.folder, "GAM.html"))

####*******************####
#### C.03: Lasso Model #### 
####*******************####

print("Running script C_03...")
rmarkdown::render(paste0(code.folder, "C_03_lasso_model.Rmd"), output_file = paste0(output.folder, "lasso.html"))

print("Analysis completed.")