load("/home/carlimm/IGS_exposome/challenge1/exposome.RDA")
library(BayesGWQS)

group_list <- list(
  
  c("hs_no2_yr_hs_h_Log", "hs_pm25_yr_hs_h_None","hs_pm25abs_yr_hs_h_Log"), # Outdoor Air (+)
  
  c("h_NO2_Log", "h_PM_Log"), # Indoor Air (-)
  
  c("h_Absorbance_Log","h_Benzene_Log"), # Indoor Air (+)
  
  c("hs_as_c_Log2",  "hs_cs_c_Log2", "hs_hg_c_Log2"), # Metals (-)
  
  c("hs_cd_c_Log2", "hs_co_c_Log2","hs_cu_c_Log2", "hs_mn_c_Log2", "hs_mo_c_Log2", "hs_pb_c_Log2"), # Metals (+)
  
  c("hs_dde_cadj_Log2", "hs_ddt_cadj_Log2", "hs_hcb_cadj_Log2",  # Organochlorines (-)
     "hs_pcb153_cadj_Log2", "hs_pcb170_cadj_Log2", "hs_pcb180_cadj_Log2"),
  
  c("hs_pcb118_cadj_Log2","hs_pcb138_cadj_Log2"), # Organochlorines (+)
  
  c("hs_dep_cadj_Log2", "hs_detp_cadj_Log2", "hs_dmp_cadj_Log2", "hs_dmtp_cadj_Log2"), # Organophosphate pesticides (-)
  
  c("hs_pbde153_cadj_Log2", "hs_pbde47_cadj_Log2"), # PBDEs (-)
  
  c("hs_pfhxs_c_Log2", "hs_pfna_c_Log2", "hs_pfoa_c_Log2", "hs_pfos_c_Log2", "hs_pfunda_c_Log2"), # PFASs (+)
  
  c("hs_bupa_cadj_Log2", "hs_etpa_cadj_Log2", "hs_mepa_cadj_Log2","hs_prpa_cadj_Log2", "hs_trcs_cadj_Log2"), # Phenols (-)
  
  c("hs_bpa_cadj_Log2", "hs_oxbe_cadj_Log2"), # Phenols (+)
  
  c("hs_ohminp_cadj_Log2", "hs_oxominp_cadj_Log2"), # Phthalates (-)
  
  c("hs_mbzp_cadj_Log2", "hs_mecpp_cadj_Log2", "hs_mehhp_cadj_Log2", "hs_mehp_cadj_Log2", "hs_meohp_cadj_Log2", "hs_mep_cadj_Log2", "hs_mibp_cadj_Log2", "hs_mnbp_cadj_Log2"), # Phthalates (+)
  
  c("hs_trafload_h_pow1over3", "hs_trafnear_h_pow1over3"), # Traffic (mixed, no drop)
  
  c("hs_ndvi100_h_None", "hs_ndvi100_s_None"), # Natural Spaces (-)
  
  c("hs_tm_mt_hs_h_None","hs_uvdvf_mt_hs_h_None"), # Meteorological (+)
  
  c("hs_dif_hours_total_None", "hs_KIDMED_None"), # Lifestyle (-)
  
  c("hs_sd_wk_None", "hs_mvpa_prd_alt_None"), # Lifestyle (+)
  
  c("hs_accesspoints300_h_Log",  "hs_connind300_h_Log", 'hs_fdensity300_s_Log',  'hs_popdens_s_Sqrt'), # Built Environment (-)
  
  c("hs_builtdens300_h_Sqrt","hs_fdensity300_h_Log", 'hs_landuseshan300_h_None',
    'hs_popdens_h_Sqrt', 'hs_accesspoints300_s_Log', 'hs_builtdens300_s_Sqrt', 'hs_connind300_s_Log','hs_landuseshan300_s_None') # Built Environment (+)
)

exposome$hs_sd_wk_None <- -1*exposome$hs_sd_wk_None
X_temp <- make.X(exposome, 21, group_list)
x.s <- make.x.s(exposome, 21, group_list)

covars_standardized <- covariates
covars_standardized$h_mbmi_None <- (covars_standardized$h_mbmi_None - mean(covars_standardized$h_mbmi_None))/sd(covars_standardized$h_mbmi_None)
covars_standardized$hs_wgtgain_None <- (covars_standardized$hs_wgtgain_None - mean(covars_standardized$hs_wgtgain_None))/sd(covars_standardized$hs_wgtgain_None)
covars_standardized$e3_gac_None <- (covars_standardized$e3_gac_None - mean(covars_standardized$e3_gac_None))/sd(covars_standardized$e3_gac_None)
covars_standardized$h_age_None <- (covars_standardized$h_age_None - mean(covars_standardized$h_age_None))/sd(covars_standardized$h_age_None)
covars_standardized$hs_child_age_None <- (covars_standardized$hs_child_age_None - mean(covars_standardized$hs_child_age_None))/sd(covars_standardized$hs_child_age_None)
covars_standardized$hs_c_height_None <- (covars_standardized$hs_c_height_None - mean(covars_standardized$hs_c_height_None))/sd(covars_standardized$hs_c_height_None)
covars_standardized$hs_c_weight_None <- (covars_standardized$hs_c_weight_None - mean(covars_standardized$hs_c_weight_None))/sd(covars_standardized$hs_c_weight_None)
covars_standardized$e3_sex_None <- ifelse(covars_standardized$e3_sex_None=="male",1,0)

co_num <- 6 ##############

Z_temp <- covars_standardized[, c(2,3,8)] # keeping only sig covars
cohorts <- Z_temp$h_cohort
Y_temp <- phenotype$hs_asthma
df <- as.data.frame(cbind(Y_temp, X_temp, cohorts))
df <- df[df$cohorts==co_num,]

X <- df[,-c(1,ncol(df))]
X <- as.matrix(X)
X <- apply(X, 2, as.numeric)
Y <- df$Y_temp
Z <- Z_temp[Z_temp$h_cohort==co_num,]
Z <- Z[,-1]
Z <- as.matrix(Z)
Z <- apply(Z, 2, as.numeric)

work_dir <- "/home/carlimm/IGS_exposome/challenge1/childonly_benzene/corr_groups"

kidcorr_co6 <- bgwqs.fit(y = Y, x = X, z = Z, x.s = x.s, n.chains = 2, method = "quantile", n.quantiles = 4,
                     n.adapt  = 500, n.iter = 125000, n.burnin = 1000, n.thin = 1, working.dir = work_dir, DIC=T)
save(kidcorr_co6, file = "kidcorr_co6.RDA")
