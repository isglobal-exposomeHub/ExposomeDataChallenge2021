# File: B_01_main_PCP_runs_on_res.R
# Author(s): Lawrence Chillrud
# Date since last edit: 4/27/21

####*******************####
#### Table of Contents ####
####*******************####
####* Notes / Description
####* 0. Package imports
####* 1. Load residuals data
####* 2. Visualize residuals
####* 3. Run gridsearches to determine optimal rank
####* 4. Run PCP
####* 5. Generate PCP reports
####* 6. Pick & save patterns for health models

####*********************####
#### Notes / Description ####
####*********************####
####* Run PCP on the residuals from A_04.
####* Perform gridsearches to find the optimal rank.
####* Then output report files on the results PCP runs.
####* By running PCP on the residuals, the dominating effect
####* on the identified patterns we saw from the different cohorts will be minimized,
####* allowing us to recover consistent patterns of exposure across cohort.  

#####********************#####
##### 0. PACKAGE IMPORTS #####
#####********************#####

print("Running script B_01...")

# for PCP
library(PCPhelpers)
library(pcpr)

# for data wrangling
library(tidyverse)

# for timing:
library(tictoc)

# for FA
library(psych)

#####******************#####
##### 1. LOAD RES DATA #####
#####******************#####
load(paste0(generated.data.folder, "exposome_pren_postnatal_residuals.RData"))

# prenatal data:
pren_exposome_res <- pren_exposome_res[, pren_description$variable_name]
colnames(pren_exposome_res) <- pren_description$labels
cg_pre <- pren_description$family

# postnatal data:
post_exposome_res <- post_exposome_res[, post_description$variable_name]
colnames(post_exposome_res) <- post_description$labels
cg_post <- post_description$family

# get the cohort info:
load(paste0(generated.data.folder, "exposome_pren_postnatal.RData"))
rg <- covariates$h_cohort

mats <- list("pre" = pren_exposome_res, "post" = post_exposome_res) %>% map(as.matrix)
mats_scaled <- list("pre" = pren_exposome_res, "post" = post_exposome_res) %>% map(as.matrix) %>% map(~ scale(., center = F))
defaults <- mats %>% map(get_pcp_defaults)

#####************************#####
##### 2. VISUALIZE RESIDUALS #####
#####************************#####

mats[[1]] %>% as_tibble() %>%
  mat_to_df(., col_names = colnames(mats[[1]])) %>%
  ggplot(aes(value)) + geom_histogram() + facet_wrap(~chem) + xlim(c(-3, 3)) + ggtitle("Prenatal residuals")

mats[[2]] %>% as_tibble() %>%
  mat_to_df(., col_names = colnames(mats[[2]])) %>%
  ggplot(aes(value)) + geom_histogram() + facet_wrap(~chem) + xlim(c(-3, 3)) + ggtitle("Postnatal residuals")

mats_scaled[[1]] %>% as_tibble() %>%
  mat_to_df(., col_names = colnames(mats_scaled[[1]])) %>%
  ggplot(aes(value)) + geom_histogram() + facet_wrap(~chem) + xlim(c(-3, 3)) + ggtitle("Prenatal residuals (scaled)")

mats_scaled[[2]] %>% as_tibble() %>%
  mat_to_df(., col_names = colnames(mats_scaled[[2]])) %>%
  ggplot(aes(value)) + geom_histogram() + facet_wrap(~chem) + xlim(c(-3, 3)) + ggtitle("Postnatal residuals (scaled)")

#####****************************#####
##### 3. RUN GRIDSEARCH FOR RANK #####
#####****************************#####

print("PRE SCALED GS")
tic()
pre_gs_scaled <- grid_search_cv(mats_scaled[["pre"]], pcp_func = root_pcp_noncvx_na, 
                         grid_df = data.frame(r = 1:15),
                         cores = 4, perc_b = 0.2, runs = 4, 
                         lambda = defaults[["pre"]]$lambda, mu = defaults[["pre"]]$mu,
                         file = paste0(generated.data.folder, "pre_res_scaled_gridsearch.Rda"))
toc()

print("POST GS SCALED")
tic()
post_gs_scaled <- grid_search_cv(mats_scaled[["post"]], pcp_func = root_pcp_noncvx_na, 
                          grid_df = data.frame(r = 1:15),
                          cores = 4, perc_b = 0.2, runs = 20, 
                          lambda = defaults[["post"]]$lambda, mu = defaults[["post"]]$mu,
                          file = paste0(generated.data.folder, "post_res_scaled_gridsearch.Rda"))
toc()

#####************#####
##### 4. RUN PCP #####
#####************#####
pcp_outs <- list()
for (i in 1:length(mats_scaled)) {
  pcp_outs[[i]] <- root_pcp_noncvx(mats_scaled[[i]], lambda = defaults[[i]]$lambda, mu = defaults[[i]]$mu, r = 4, verbose = T)
}

#####*********************#####
##### 5. GENERATE REPORTS #####
#####*********************#####
make_pcp_report <- function(
  dataset, subname, mat, 
  cohort, algo, param_status, pcp, 
  scale_flag = F, rowvar_name = "ID", rowvar = NULL, 
  rowgroupings = NULL, rowgroupname = "cohort",
  colnames = NULL, colgroupings = NULL, run_nmf = TRUE, fa_method = "varimax") {
  
  filename <- paste0(paste(subname, cohort, algo, param_status, sep = "_"), ".html") 
  
  rmarkdown::render(
    paste0(code.folder, "F_pcp_report.Rmd"),
    params = list(
      algo = algo,
      dataset = dataset,
      subname = subname,
      parameters = param_status,
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

cg <- list(cg_pre, cg_post)

for (i in 1:length(pcp_outs)) {
  make_pcp_report(
    dataset = "exposome-res-scaled", subname = names(mats_scaled)[i], mat = mats_scaled[[i]], 
    cohort = "ac", algo = "noncvx", 
    param_status = "default-r4-res-scaled", pcp = pcp_outs[[i]], 
    scale_flag = F, rowvar_name = "ID", rowvar = as.numeric(rownames(mats[[i]])), 
    rowgroupings = rg, rowgroupname = "cohort",
    colnames = colnames(mats[[i]]), colgroupings = cg[[i]], run_nmf = FALSE, fa_method = "promax"
  )
}

#####**********************************#####
##### 6. PICK & SAVE PATTERNS - SCALED #####
#####**********************************#####
factors <- pcp_outs %>% map(function(p) {
  fa(p$L, n.obs = nrow(p$L), nfactors = 4, rotate = "promax", scores = "regression")
})

names(factors) <- c("pre", "post")
names(pcp_outs) <- c("pre", "post")
save(factors, file = paste0(generated.data.folder, "health_model_patterns.RData"))
save(pcp_outs, file = paste0(generated.data.folder, "health_model_pcp.RData"))

# S matrix correlations:
pcp_outs$post$S %>% GGally::ggcorr(., method = c("pairwise.complete.obs", "pearson"), 
                   limits = F, label = T, label_size = 3, label_alpha = T, 
                   size = 3, layout.exp = 1) + ggtitle("Postnatal S matrix: Pearson correlation")
