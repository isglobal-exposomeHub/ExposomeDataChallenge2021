### ISGlobal Exposome Challenge 2021 - mediation analysis ###
### This script is for mediation analysis for each single exposure and family-based ERS
## using HIMA package.
### Author: Ziyue Wang


library(HIMA)
library(dplyr) 
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(dplyr) 
library(tidyverse)
library(broom)
library(ggforestplot)
library(reshape2)

DATA = "~/NIEHS/Exposome_Data_Challenge/DATA/"
CODE = "~/NIEHS/Exposome_Data_Challenge/CODE/"
RESULT = "~/NIEHS/Exposome_Data_Challenge/RESULT/"


## Input Data ##
load(paste0(RESULT,"exposure_select.RData"))
load(paste0(RESULT,"DEGS_No_PCS.rdata"))
load(paste0(RESULT, "risk_score.rdata"))
load(paste0(RESULT,"gene_counts.rdata"))
expo_cov_pheno = readRDS(file = paste0(RESULT, "expo_cov_pheno.rds"))


## Define Mediators ##
mediator = gene_counts %>%
  filter(Gene_Name %in% degs$Gene_ID) %>%
  column_to_rownames("Gene_Name") %>%
  t() %>% data.frame() %>%
  rownames_to_column(var = "ID") %>%
  mutate(ID = as.integer(ID)) %>%
  arrange(ID) %>%
  column_to_rownames(var = "ID")


## Mediation Analysis (HIMA method) ##
# Option 1: single exposure
hima_out <- list()
for (i in 1:ncol(exposure)) {
  temp = hima(exposure[,i], outcome, mediator, covariate, covariate, family = "binomial")
  
  if(length(temp) > 0 ){
    hima_out[[i]] = temp
    hima_out[[i]]$gene = rownames(temp)
    hima_out[[i]]$exposure = rep(colnames(exposure)[i], nrow(temp))
  }
  else hima_out[[i]] = temp
  
  rm(temp)
}
hima_out = do.call(rbind, hima_out)
saveRDS(hima_out, file = paste0(RESULT, "hima_single_expo.rds"))

gene_med = degs %>% filter(Gene_ID %in% unique(hima_out$gene))
saveRDS(gene_med, file = paste0(RESULT, "gene_med.rds"))



# Option 2: for each exposure family-based ERS #
exposure_ERS = risk_score_data %>%
  select(Air.Pollution_rs:total_rs)

hima_out_ERS <- list()
for (i in 1:ncol(exposure_ERS)) {
  temp = hima(exposure_ERS[,i], outcome, mediator, covariate, covariate, family = "binomial")
  
  if(length(temp) > 0 ){
    hima_out_ERS[[i]] = temp
    hima_out_ERS[[i]]$gene = rownames(temp)
    hima_out_ERS[[i]]$exposre = rep(colnames(exposure_ERS)[i], nrow(temp))
  }
  else hima_out_ERS[[i]] = temp
  
  rm(temp)
}
hima_out_ERS = do.call(rbind, hima_out_ERS)
saveRDS(hima_out_ERS, file = paste0(RESULT, "hima_ERS.rds"))

gene_med_ERS = degs %>% filter(Gene_ID %in% unique(hima_out_ERS$gene))
saveRDS(gene_med_ERS, file = paste0(RESULT, "gene_med_ERS.rds"))


## Visualization ##
hima_single_expo = readRDS(file = paste0(RESULT,"hima_single_expo.rds"))
hima_ERS = readRDS(file = paste0(RESULT,"hima_ERS.rds"))

# forest plot of association effect for ERS
forestplot(risk_score_glmop %>% 
             filter(variable %in% grep("_rs",risk_score_glmop$variable, value = T)) %>%
             arrange(variable),
           name = variable,
           estimate = Estimate,
           se = Std.Error,
           pvalue = pvalue)

# mediation effect
hima_single_expo = hima_single_expo %>%
  left_join(coef_fit, by = "exposure")
ggplot(hima_single_expo[nrow(hima_single_expo):1,], 
       aes(abs(`alpha*beta`), factor(exposure, levels = unique(exposure)))) +
  geom_point(aes(colour = family, shape = factor(gene, levels = unique(gene)))) +
  scale_shape_manual(values = 0:16) + 
  labs(shape = "gene") + 
  xlab("abs(effect)") +
  ylab(NULL) +
  theme_bw()

ggplot(hima_ERS[nrow(hima_ERS):1,], 
       aes(abs(`alpha*beta`), factor(exposure, levels = unique(exposure)))) +
  geom_point(aes(shape = factor(gene))) +
  labs(shape = "gene") + 
  xlab("abs(effect)") +
  ylab(NULL) +
  theme_bw()