### ISGlobal Exposome Challenge 2021 - variable selection for exposures ###
### This script pre-select exposures using the Elastic-net.
### Author: Ziyue Wang


library(Biobase)
library(tidyverse)
library(glmnet)
library(ggplot2)
library(dplyr) 
library(broom)
library(gtsummary)
library(doParallel)

DATA = "~/NIEHS/Exposome_Data_Challenge/DATA/"
CODE = "~/NIEHS/Exposome_Data_Challenge/CODE/"
RESULT = "~/NIEHS/Exposome_Data_Challenge/RESULT/"


## Input Data ##
load(paste0(DATA, "exposome.rdata"))
load(paste0(DATA, "genexpr.rdata"))

expr <- exprs(genexpr)
phenoDataFrame <- pData(genexpr)
probes <- fData(genexpr)


## Pre-processing Data ##
# match samples in exposome and transcriptome data
# choose covariates: age, sex, cohort, bmi of child, birth weight and ethnicity
phenoDataFrame$ID = as.integer(phenoDataFrame$ID)
expo_cov_pheno = exposome %>%
  left_join(phenoDataFrame %>% 
              select(ID, h_ethnicity_cauc), by = "ID") %>% 
  left_join(covariates %>% 
              select(ID, e3_sex_None, hs_child_age_None, h_cohort), by = "ID") %>%
  left_join(phenotype %>% 
              select(ID, e3_bw, hs_zbmi_who, hs_asthma), by = "ID") %>%
  filter(ID %in% phenoDataFrame$ID) %>%
  mutate(hs_asthma = factor(hs_asthma)) %>%
  column_to_rownames("ID")
saveRDS(expo_cov_pheno, file = paste0(RESULT, "expo_cov_pheno.rds"))

# Sample characteristics
expo_cov_pheno %>%
  select(hs_asthma, e3_sex_None, hs_child_age_None, h_cohort, e3_bw, hs_zbmi_who, h_ethnicity_cauc) %>%
  tbl_summary(by = hs_asthma, 
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"), 
              type = list(all_dichotomous() ~ "categorical"),
              label = list(e3_sex_None ~ "Sex",
                           hs_child_age_None ~ "Age",
                           h_cohort ~ "Cohort",
                           e3_bw ~ "Birth Weight",
                           hs_zbmi_who ~ "BMI of Child",
                           h_ethnicity_cauc ~ "Ethnicity")) %>%
  add_n() %>% add_p() %>%
  modify_header(label = "**Asthma**") %>% # update the column header
  bold_labels() 


## Variable Selection (elastic net) ##
# transforms any qualitative variables into dummy variables. 
x = model.matrix(hs_asthma ~ ., data = expo_cov_pheno)[,-1]

# tuning alpha
# set a foldid, allows us to apply the same CV folds to each model
# CV for each alpha = 0, 0.1, ... , 0.9, 1.0
set.seed(123)
registerDoParallel(4)

foldid=sample(1:10, size=nrow(x), replace=TRUE)
for (i in 0:10) {
  assign(paste("fit_cv", i, sep=""), 
         cv.glmnet(x, expo_cov_pheno$hs_asthma, 
                   foldid = foldid, parallel = TRUE, 
                   alpha=i/10,
                   family="binomial"))
}
par(mfrow=c(3,4))
for(i in 0:10) {plot(get(paste("fit_cv", i, sep="")))}
plot(log(fit_cv0$lambda),fit_cv0$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=fit_cv0$name)
points(log(fit_cv2$lambda),fit_cv2$cvm,pch=19,col="grey")
points(log(fit_cv5$lambda),fit_cv5$cvm,pch=19,col="blue")
points(log(fit_cv8$lambda),fit_cv8$cvm,pch=19,col="yellow")
points(log(fit_cv10$lambda),fit_cv10$cvm,pch=19,col="green")
legend("topright",legend=c("alpha = 0 (ridge)","alpha = 0.2","alpha = 0.5","alpha = 0.8","alpha = 1 (lasso)"),
       pch=19,col=c("red","grey","blue","yellow","green"))
par(mfrow=c(1,1))


# perform 10-fold cross-validation to select lambda
set.seed(123)
registerDoParallel(4)
fit_cv <- cv.glmnet(x, expo_cov_pheno$hs_asthma, 
                    parallel = TRUE, 
                    alpha=0.2,
                    family="binomial")
plot(fit_cv, main="alpha = 0.2")
lambda_min = fit_cv$lambda.min


# fit model
fit <- glmnet(x, expo_cov_pheno$hs_asthma, family="binomial", alpha=0.2, lambda = lambda_min)

# define the selected exposures
exposure_selec_all = rownames(coef(fit))[which(coef(fit) != 0)[-1]]

exposure_selec_num = codebook %>% 
  filter(var_type == "numeric", variable_name %in% exposure_selec_all) %>%
  pull(variable_name) %>% 
  as.character()

# deal with dummy variable
exposure_selec_dummy = exposure_selec_all[!(exposure_selec_all %in% exposure_selec_num)]
exposure_selec_dummy[-grep("\\(", exposure_selec_dummy)] = c("h_pamod_t3_None",
                                                             "h_pavig_t3_None",
                                                             "hs_blueyn300_h_None",
                                                             "hs_ln_cat_h_None",
                                                             "hs_lden_cat_s_None",
                                                             "hs_lden_cat_s_None",
                                                             "hs_contactfam_3cat_num_None",
                                                             "hs_cotinine_cdich_None",
                                                             "hs_cotinine_mcat_None",
                                                             "hs_globalexp2_None",
                                                             "hs_smk_parents_None",
                                                             "hs_smk_parents_None",
                                                             "e3_sex_None")
exposure_selec_dummy[grep("\\(", exposure_selec_dummy)] = unlist(lapply(strsplit(exposure_selec_dummy[grep("\\(", exposure_selec_dummy)],"\\("), function(x) x[[1]]))

exposure_selec = c(exposure_selec_num, unique(exposure_selec_dummy))
exposure_selec = exposure_selec[exposure_selec != "e3_sex_None"]


## Save the Pre-selected Exposures ##
exposure = x %>% as.data.frame() %>% select(exposure_selec_all[exposure_selec_all != "e3_sex_Nonemale"])
coef_fit = coef(fit, s = "lambda.min")[which(coef(fit, s = "lambda.min") != 0)[-c(1,83)]]
codebook_expo_sele = codebook %>% 
  filter(variable_name %in% exposure_selec)

# covariates
covariate = expo_cov_pheno %>%
  select(e3_sex_None, hs_child_age_None, h_cohort, e3_bw, hs_zbmi_who, h_ethnicity_cauc)
# outcome - Asthma
outcome = expo_cov_pheno %>% pull(hs_asthma) 

save(exposure, coef_fit, outcome, covariate, codebook_expo_sele, 
     file = paste0(RESULT, "exposure_select.RData"))


