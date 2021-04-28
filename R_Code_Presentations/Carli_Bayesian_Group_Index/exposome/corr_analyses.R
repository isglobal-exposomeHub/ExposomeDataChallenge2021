load("C:/Users/Matt Carli/Desktop/IGS_Exposome_Stuff/challenge1/exposome.RDA")
library(BayesGWQS)

######################### Correlations among Covariates ######################

Z <- covariates[, -c(1:3)] # removing ID var, child sex var, and cohort var; as they are not at least ordinal
Z <- as.matrix(Z)
Z <- apply(Z, 2, as.numeric)
z_corrs <- cor(Z, method = "spearman")

# Drop year of birth, maybe child height or weight

################################### Correlations among Standardized Covars ###################################3

covars_standardized <- covariates
covars_standardized$h_mbmi_None <- (covars_standardized$h_mbmi_None - mean(covars_standardized$h_mbmi_None))/sd(covars_standardized$h_mbmi_None)
covars_standardized$hs_wgtgain_None <- (covars_standardized$hs_wgtgain_None - mean(covars_standardized$hs_wgtgain_None))/sd(covars_standardized$hs_wgtgain_None)
covars_standardized$e3_gac_None <- (covars_standardized$e3_gac_None - mean(covars_standardized$e3_gac_None))/sd(covars_standardized$e3_gac_None)
covars_standardized$h_age_None <- (covars_standardized$h_age_None - mean(covars_standardized$h_age_None))/sd(covars_standardized$h_age_None)
covars_standardized$hs_child_age_None <- (covars_standardized$hs_child_age_None - mean(covars_standardized$hs_child_age_None))/sd(covars_standardized$hs_child_age_None)
covars_standardized$hs_c_height_None <- (covars_standardized$hs_c_height_None - mean(covars_standardized$hs_c_height_None))/sd(covars_standardized$hs_c_height_None)
covars_standardized$hs_c_weight_None <- (covars_standardized$hs_c_weight_None - mean(covars_standardized$hs_c_weight_None))/sd(covars_standardized$hs_c_weight_None)

Z_stand <- covars_standardized[, -c(1:3)] # removing ID var, child sex var, and cohort var; as they are not at least ordinal
Z_stand <- as.matrix(Z_stand)
Z_stand <- apply(Z_stand, 2, as.numeric)
zstand_corrs <- cor(Z_stand, method = "spearman")

################################## Correlations of predictors with outcome ###############################

group_list <- list(
    
    c("hs_no2_yr_hs_h_Log", "hs_pm10_yr_hs_h_None", "hs_pm25_yr_hs_h_None","hs_pm25abs_yr_hs_h_Log"), # Outdoor Air
    
    c("h_Absorbance_Log", "h_NO2_Log", "h_PM_Log", "h_Benzene_Log"), # Indoor Air
    
    c("hs_as_c_Log2", "hs_cd_c_Log2", "hs_co_c_Log2", "hs_cs_c_Log2", "hs_cu_c_Log2",
      "hs_hg_c_Log2", "hs_mn_c_Log2", "hs_mo_c_Log2", "hs_pb_c_Log2"), # Metals
    
    c("hs_dde_cadj_Log2", "hs_ddt_cadj_Log2", "hs_hcb_cadj_Log2", "hs_pcb118_cadj_Log2", # Organochlorines
      "hs_pcb138_cadj_Log2", "hs_pcb153_cadj_Log2", "hs_pcb170_cadj_Log2", "hs_pcb180_cadj_Log2"),
    
    c("hs_dep_cadj_Log2", "hs_detp_cadj_Log2", "hs_dmp_cadj_Log2", "hs_dmtp_cadj_Log2"), # Organophosphate pesticides
    
    c("hs_pbde153_cadj_Log2", "hs_pbde47_cadj_Log2"), # PBDEs
    
    c("hs_pfhxs_c_Log2", "hs_pfna_c_Log2", "hs_pfoa_c_Log2", "hs_pfos_c_Log2", "hs_pfunda_c_Log2"), # PFASs
    
    c("hs_bpa_cadj_Log2", "hs_bupa_cadj_Log2", "hs_etpa_cadj_Log2", "hs_mepa_cadj_Log2",
      "hs_oxbe_cadj_Log2", "hs_prpa_cadj_Log2", "hs_trcs_cadj_Log2"), # Phenols
    
    c("hs_mbzp_cadj_Log2", "hs_mecpp_cadj_Log2", "hs_mehhp_cadj_Log2", "hs_mehp_cadj_Log2", "hs_meohp_cadj_Log2", "hs_mep_cadj_Log2",
      "hs_mibp_cadj_Log2", "hs_mnbp_cadj_Log2", "hs_ohminp_cadj_Log2", "hs_oxominp_cadj_Log2"), # Phthalates
    
    c("hs_trafload_h_pow1over3", "hs_trafnear_h_pow1over3"), # Traffic
    
    c("hs_ndvi100_h_None", "hs_ndvi100_s_None"), # Natural Spaces
    
    c("hs_hum_mt_hs_h_None", "hs_tm_mt_hs_h_None","hs_uvdvf_mt_hs_h_None"), # Meteorological
    
    c("hs_dif_hours_total_None","hs_sd_wk_None", "hs_mvpa_prd_alt_None", "hs_KIDMED_None"), # Lifestyle
    
    c("hs_accesspoints300_h_Log", "hs_builtdens300_h_Sqrt", "hs_connind300_h_Log", "hs_fdensity300_h_Log", 'hs_landuseshan300_h_None',
      'hs_popdens_h_Sqrt', 'hs_accesspoints300_s_Log', 'hs_builtdens300_s_Sqrt', 'hs_connind300_s_Log',
      'hs_fdensity300_s_Log', 'hs_landuseshan300_s_None', 'hs_popdens_s_Sqrt') # Built Environment
)

X <- make.X(exposome, 14, group_list)
Y <- phenotype$hs_asthma
mat <- as.data.frame(cbind(Y, X))
mat$hs_sd_wk_None <- -1*mat$hs_sd_wk_None
asthma_corrs <- cor(mat, method = "spearman")

################################## Correlations of predictors with outcome, by cohort ###############################

group_list <- list(
  
  c("hs_no2_yr_hs_h_Log", "hs_pm10_yr_hs_h_None", "hs_pm25_yr_hs_h_None","hs_pm25abs_yr_hs_h_Log"), # Outdoor Air
  
  c("h_Absorbance_Log", "h_NO2_Log", "h_PM_Log", "h_Benzene_Log"), # Indoor Air
  
  c("hs_as_c_Log2", "hs_cd_c_Log2", "hs_co_c_Log2", "hs_cs_c_Log2", "hs_cu_c_Log2",
    "hs_hg_c_Log2", "hs_mn_c_Log2", "hs_mo_c_Log2", "hs_pb_c_Log2"), # Metals
  
  c("hs_dde_cadj_Log2", "hs_ddt_cadj_Log2", "hs_hcb_cadj_Log2", "hs_pcb118_cadj_Log2", # Organochlorines
    "hs_pcb138_cadj_Log2", "hs_pcb153_cadj_Log2", "hs_pcb170_cadj_Log2", "hs_pcb180_cadj_Log2"),
  
  c("hs_dep_cadj_Log2", "hs_detp_cadj_Log2", "hs_dmp_cadj_Log2", "hs_dmtp_cadj_Log2"), # Organophosphate pesticides
  
  c("hs_pbde153_cadj_Log2", "hs_pbde47_cadj_Log2"), # PBDEs
  
  c("hs_pfhxs_c_Log2", "hs_pfna_c_Log2", "hs_pfoa_c_Log2", "hs_pfos_c_Log2", "hs_pfunda_c_Log2"), # PFASs
  
  c("hs_bpa_cadj_Log2", "hs_bupa_cadj_Log2", "hs_etpa_cadj_Log2", "hs_mepa_cadj_Log2",
    "hs_oxbe_cadj_Log2", "hs_prpa_cadj_Log2", "hs_trcs_cadj_Log2"), # Phenols
  
  c("hs_mbzp_cadj_Log2", "hs_mecpp_cadj_Log2", "hs_mehhp_cadj_Log2", "hs_mehp_cadj_Log2", "hs_meohp_cadj_Log2", "hs_mep_cadj_Log2",
    "hs_mibp_cadj_Log2", "hs_mnbp_cadj_Log2", "hs_ohminp_cadj_Log2", "hs_oxominp_cadj_Log2"), # Phthalates
  
  c("hs_trafload_h_pow1over3", "hs_trafnear_h_pow1over3"), # Traffic
  
  c("hs_ndvi100_h_None", "hs_ndvi100_s_None"), # Natural Spaces
  
  c("hs_hum_mt_hs_h_None", "hs_tm_mt_hs_h_None","hs_uvdvf_mt_hs_h_None"), # Meteorological
  
  c("hs_dif_hours_total_None","hs_sd_wk_None", "hs_mvpa_prd_alt_None", "hs_KIDMED_None"), # Lifestyle
  
  c("hs_accesspoints300_h_Log", "hs_builtdens300_h_Sqrt", "hs_connind300_h_Log", "hs_fdensity300_h_Log", 'hs_landuseshan300_h_None',
    'hs_popdens_h_Sqrt', 'hs_accesspoints300_s_Log', 'hs_builtdens300_s_Sqrt', 'hs_connind300_s_Log',
    'hs_fdensity300_s_Log', 'hs_landuseshan300_s_None', 'hs_popdens_s_Sqrt') # Built Environment
)

X <- make.X(exposome, 14, group_list)
Y <- phenotype$hs_asthma
Z <- covariates$h_cohort
mat <- as.data.frame(cbind(Y, X, Z))
mat$hs_sd_wk_None <- -1*mat$hs_sd_wk_None

mat1 <- mat[mat$Z==1,]
mat2 <- mat[mat$Z==2,]
mat3 <- mat[mat$Z==3,]
mat4 <- mat[mat$Z==4,]
mat5 <- mat[mat$Z==5,]
mat6 <- mat[mat$Z==6,]

mat1 <- mat1[,-c(78)]
mat2 <- mat2[,-c(78)]
mat3 <- mat3[,-c(78)]
mat4 <- mat4[,-c(78)]
mat5 <- mat5[,-c(78)]
mat6 <- mat6[,-c(78)]

asthma_corrs1 <- cor(mat1, method = "spearman")
asthma_corrs2 <- cor(mat2, method = "spearman")
asthma_corrs3 <- cor(mat3, method = "spearman")
asthma_corrs4 <- cor(mat4, method = "spearman")
asthma_corrs5 <- cor(mat5, method = "spearman")
asthma_corrs6 <- cor(mat6, method = "spearman")


