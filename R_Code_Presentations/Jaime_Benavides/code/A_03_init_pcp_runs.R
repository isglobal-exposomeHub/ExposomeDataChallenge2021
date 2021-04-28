# File: A_03_init_pcp_runs.R
# Author(s): Lawrence Chillrud
# Date since last edit: 4/27/21

####*******************####
#### Table of Contents ####
####*******************####
####* Notes / Description
####* 0. Package imports
####* 1. Read in data
####* 2. Define function to run PCP and make corresponding report
####* 3. Run gridsearches to determine rank for init reports
####* 4. Generate reports

####*********************####
#### Notes / Description ####
####*********************####
####* This file accomplishes the following:
####* Using default lambda and mu, conduct a gridsearch over ranks 3-10
####* for each exposome data matrix (i.e. the pre- and postnatal matrices).
####* Once the ranks with the best fit for their respective data matrices 
####* are determined, generate initial PCP reports.
####* These reports can then be compared with the raw data reports.

####********************####
#### 0. Package Imports ####
####********************####

print("Running script A_03...")

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

####*************************************************************####
#### 2. Define function to run PCP and make corresponding report ####
####*************************************************************####
make_pcp_report <- function(
  dataset, subname, mat, 
  cohort, algo, param_status, rank, 
  scale_flag = F, rowvar_name = "ID", rowvar = NULL, 
  rowgroupings = NULL, rowgroupname = "cohort",
  colnames = NULL, colgroupings = NULL, run_nmf = TRUE, fa_method = "varimax") {
  
  # run PCP
  dfp <- get_pcp_defaults(mat)
  pcp <- root_pcp_noncvx_nonneg(mat, lambda = dfp$lambda, mu = dfp$mu, r = rank, verbose = T)
  
  # generate report
  filename <- paste0(paste(subname, algo, "pcp", param_status, "rank", rank, sep = "_"), ".html") 
  
  rmarkdown::render(
    paste0(code.folder, "F_pcp_report.Rmd"),
    params = list(
      algo = algo,
      dataset = dataset,
      subname = subname,
      parameters = paste(param_status, "w/rank", rank, sep = " "),
      convergence = pcp$final_iter,
      L = pcp$L,
      S = pcp$S,
      rowvar_name = rowvar_name,
      rowvar = rowvar,
      colnames = colnames,
      colgroupings = colgroupings,
      rowgroupings = rowgroupings,
      rowgroupname = rowgroupname,
      scale_flag = scale_flag,
      ranktol = 1e-04,
      sparsitytol = 1e-04,
      pcs = NULL,
      run_nmf = run_nmf,
      fa_method = fa_method),
    output_file = paste0(output.folder, filename)
  )
}

####********************************************************####
#### 3. Run gridsearches to determine rank for init reports ####
####********************************************************####

# 3a prenatal gridsearch
pre_dfp <- get_pcp_defaults(mat_pre)

pre_gs <- grid_search_cv(mat_pre, pcp_func = root_pcp_noncvx_nonnegL_na, 
                         grid_df = data.frame(r = 3:10),
                         cores = 4, perc_b = 0.2, runs = 4, 
                         lambda = pre_dfp$lambda, mu = pre_dfp$mu)

pre_best_r <- pre_gs$formatted$r[which.min(pre_gs$formatted$value)]

# 3b postnatal gridsearch
post_dfp <- get_pcp_defaults(mat_post)

post_gs <- grid_search_cv(mat_post, pcp_func = root_pcp_noncvx_nonnegL_na, 
                         grid_df = data.frame(r = 3:10),
                         cores = 4, perc_b = 0.2, runs = 4, 
                         lambda = post_dfp$lambda, mu = post_dfp$mu)

post_best_r <- post_gs$formatted$r[which.min(post_gs$formatted$value)]

##***********************####
#### 4. Generate reports ####
####*********************####

# 4a prenatal PCP report
make_pcp_report(
  dataset = "Exposome", subname = "prenatal", mat = mat_pre, 
  cohort = "ac", algo = "noncvx", 
  param_status = "default", rank = pre_best_r, 
  scale_flag = F, rowvar_name = "ID", rowvar = as.numeric(rownames(mat_pre)), 
  rowgroupings = rowgroupings, rowgroupname = "cohort",
  colnames = colnames(mat_pre), colgroupings = colgroupings_pre, fa_method = "promax"
)

# 4b postnatal PCP report
make_pcp_report(
  dataset = "Exposome", subname = "postnatal", mat = mat_post, 
  cohort = "ac", algo = "noncvx", 
  param_status = "default", rank = post_best_r, 
  scale_flag = F, rowvar_name = "ID", rowvar = as.numeric(rownames(mat_post)), 
  rowgroupings = rowgroupings, rowgroupname = "cohort",
  colnames = colnames(mat_post), colgroupings = colgroupings_post, fa_method = "promax"
)