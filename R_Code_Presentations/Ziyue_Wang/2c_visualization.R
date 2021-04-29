### ISGlobal Exposome Challenge 2021 - visualization for exposures and ERS ###
### This script generate several plots for visualizing exposures and ERS
### Author: Ziyue Wang


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
load(paste0(RESULT, "risk_score.rdata"))
expo_cov_pheno = readRDS(file = paste0(RESULT, "expo_cov_pheno.rds"))


## Effect Size for Important Exposures ##
coef_fit = data.frame(exposure = colnames(exposure),
                      beta = coef_fit) %>%
  mutate(variable_name = exposure)

exposure_selec_dummy = coef_fit$variable_name[!coef_fit$variable_name %in% codebook_expo_sele$variable_name]
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
                                                             "hs_smk_parents_None")
exposure_selec_dummy[grep("\\(", exposure_selec_dummy)] = unlist(lapply(strsplit(exposure_selec_dummy[grep("\\(", exposure_selec_dummy)],"\\("), function(x) x[[1]]))
coef_fit$variable_name[!coef_fit$variable_name %in% codebook_expo_sele$variable_name] = exposure_selec_dummy

coef_fit = coef_fit %>%
  left_join(codebook_expo_sele, by = "variable_name")

ggplot(coef_fit[nrow(coef_fit):1,], aes(beta, factor(exposure, levels = unique(exposure))))+
  geom_point(aes(colour = family)) +
  ggtitle("Influential Exposures") +
  xlab("effect size") +
  ylab(NULL) +
  theme_bw()

ggplot(coef_fit[nrow(coef_fit):1,], aes(beta, factor(labels, levels = unique(labels))))+
  geom_point(aes(colour = family)) +
  ggtitle("Influential Exposures") +
  xlab("effect size") +
  ylab(NULL) +
  theme_bw()


## Correlation Heatmap ##
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
corr_heatmap <- function(x, order = TRUE, method = c("pearson","spearman")){
  if(!order) cormat <- cor(x, method = method)
  else cormat <- reorder_cormat(cor(x, method = method))
  melted_cormat <- melt(cormat, na.rm = TRUE)
  
  ggplot(data = melted_cormat[,], aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="correlation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank())
}

# for all exposures
x = model.matrix(hs_asthma ~ ., data = expo_cov_pheno)[,-1]
corr_heatmap(x[,ncol(x):1], order = FALSE, method = "spearman")
corr_heatmap(x[,ncol(x):1], order = TRUE, method = "spearman")

# for selected exposures
corr_heatmap(exposure[,ncol(exposure):1], order = TRUE, method = "spearman")


# for ERS
exposure_ERS = risk_score_data %>%
  select(Air.Pollution_rs:Water.DBPs_rs)
corr_heatmap(exposure_ERS, order = TRUE, method = "pearson")
