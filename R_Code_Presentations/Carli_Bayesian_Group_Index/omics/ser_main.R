load("/home/carlimm/IGS_exposome/challenge1/exposome.RDA")
load("/home/carlimm/IGS_exposome/challenge3/metabol_serum.RDA")

library(dplyr)
library(Biobase)
library(BayesGWQS)

mat_off <- exprs(metabol_serum)
mat <- t(mat_off)
ID <- as.numeric(rownames(mat))
serum_mat <- data.frame(ID, mat)

featdata <- featureData(metabol_serum)
serum_family <- featdata[[1]]

# Getting Serum family groups

biogens <- paste0("metab_", which(serum_family == "biogenicamines"))

aminos <- paste0("metab_", which(serum_family == "aminoacids"))

acyls <- paste0("metab_", which(serum_family == "acylcarnitines"))

glyceros <- paste0("metab_", which(serum_family == "glycerophospholipids"))

sphingos <- paste0("metab_", which(serum_family == "sphingolipids"))

chem_list <- list(
    
    c("ID", "hs_no2_yr_hs_h_Log", "hs_pm10_yr_hs_h_None", "hs_pm25_yr_hs_h_None","hs_pm25abs_yr_hs_h_Log"), # Outdoor Air
    
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

serum_list <- list(c("ID", biogens), aminos, acyls, glyceros, sphingos)
ser_mat <- make.X(serum_mat, 5, serum_list)

exposome$hs_sd_wk_None <- -1*exposome$hs_sd_wk_None
chem_mat <- make.X(exposome, 14, chem_list)

pheno_mat <- select(phenotype, c(ID, hs_asthma))

cov_mat <- select(covariates, c(ID, h_cohort, e3_sex_None, h_age_None))

cmat1 <- merge(pheno_mat, cov_mat, by = 'ID')
cmat2 <- merge(cmat1, chem_mat, by = 'ID')
cmat <- merge(cmat2, ser_mat, by = 'ID')
cmat$h_age_None <- (cmat$h_age_None - mean(cmat$h_age_None))/sd(cmat$h_age_None)
cmat$e3_sex_None <- ifelse(cmat$e3_sex_None=="male",1,0)

corr_mat <- select(cmat, c(hs_asthma, biogens, aminos, acyls, glyceros, sphingos))

serum_corrs <- cor(corr_mat, method = "spearman")
serum_corrs <- serum_corrs[2:177,1]
biogen_neg <- names(which(serum_corrs[1:14]<0))
biogen_pos <- names(which(serum_corrs[1:14]>0))
aminos_neg <- names(which(serum_corrs[15:35]<0))
aminos_pos <- names(which(serum_corrs[15:35]>0))
acyls_neg <- names(which(serum_corrs[36:73]<0))
acyls_pos <- names(which(serum_corrs[36:73]>0))
glyceros_neg <- names(which(serum_corrs[74:162]<0))
glyceros_pos <- names(which(serum_corrs[74:162]>0))
sphingos_neg <- names(which(serum_corrs[163:176]<0))
sphingos_pos <- names(which(serum_corrs[163:176]>0))

group_list <- list(
    
    c("hs_no2_yr_hs_h_Log", "hs_pm10_yr_hs_h_None", "hs_pm25_yr_hs_h_None","hs_pm25abs_yr_hs_h_Log"), # Outdoor Air
    
    c("h_Absorbance_Log", "h_NO2_Log"), # Indoor Air (-)
    
    c("h_PM_Log", "h_Benzene_Log"), # Indoor Air (+)
    
    c("hs_as_c_Log2",  "hs_co_c_Log2", "hs_cs_c_Log2", "hs_cu_c_Log2", "hs_pb_c_Log2"), # Metals (-)
    
    c("hs_cd_c_Log2","hs_hg_c_Log2", "hs_mn_c_Log2", "hs_mo_c_Log2"),# Metals (+)
    
    c("hs_dde_cadj_Log2", "hs_ddt_cadj_Log2", "hs_hcb_cadj_Log2", # Organochlorines (-)
      "hs_pcb138_cadj_Log2", "hs_pcb153_cadj_Log2", "hs_pcb170_cadj_Log2", "hs_pcb180_cadj_Log2"),
    
    c("hs_dep_cadj_Log2", "hs_detp_cadj_Log2", "hs_dmp_cadj_Log2"), # Organophosphate pesticides (-)
    
    c("hs_pbde153_cadj_Log2", "hs_pbde47_cadj_Log2"), # PBDEs (-)
    
    c("hs_pfhxs_c_Log2", "hs_pfoa_c_Log2", "hs_pfos_c_Log2", "hs_pfunda_c_Log2"), # PFASs (+)
    
    c("hs_bpa_cadj_Log2", "hs_bupa_cadj_Log2", "hs_etpa_cadj_Log2", "hs_mepa_cadj_Log2",
      "hs_oxbe_cadj_Log2", "hs_prpa_cadj_Log2"), # Phenols (-)
    
    c("hs_mbzp_cadj_Log2", "hs_mecpp_cadj_Log2", "hs_mehhp_cadj_Log2", 
      "hs_mehp_cadj_Log2", "hs_meohp_cadj_Log2", "hs_oxominp_cadj_Log2"), # Phthalates (-)
    
    c("hs_mep_cadj_Log2","hs_mibp_cadj_Log2", "hs_mnbp_cadj_Log2", "hs_ohminp_cadj_Log2"),# Phthalates (+)
    
    c("hs_trafload_h_pow1over3", "hs_trafnear_h_pow1over3"), # Traffic (-)
    
    c("hs_ndvi100_h_None", "hs_ndvi100_s_None"), # Natural Spaces (+)
    
    c("hs_hum_mt_hs_h_None", "hs_uvdvf_mt_hs_h_None"), # Meteorological (-)
    
    c("hs_dif_hours_total_None", "hs_KIDMED_None"), # Lifestyle (-)
    
    c("hs_sd_wk_None", "hs_mvpa_prd_alt_None"), # Lifestyle (+)
    
    c("hs_accesspoints300_h_Log", "hs_builtdens300_h_Sqrt", 'hs_popdens_h_Sqrt', 'hs_builtdens300_s_Sqrt', 'hs_popdens_s_Sqrt'), # Built Environment (-)
    
    c('hs_accesspoints300_s_Log', "hs_connind300_h_Log", 'hs_connind300_s_Log', "hs_fdensity300_h_Log", 
      'hs_fdensity300_s_Log', 'hs_landuseshan300_h_None', 'hs_landuseshan300_s_None'), # Built Environment (+)
    
    biogen_neg, biogen_pos,
    
    aminos_neg, aminos_pos,
    
    acyls_neg, acyls_pos,
    
    glyceros_neg, glyceros_pos,
    
    sphingos_neg, sphingos_pos
)

X <- make.X(cmat, 29, group_list)
x.s <- make.x.s(cmat, 29, group_list)

Z <- as.matrix(select(cmat, h_cohort, e3_sex_None, h_age_None))
Z <- apply(Z, 2, as.numeric)

Y <- cmat[,2]

work_dir <- "/home/carlimm/IGS_exposome/challenge3/serum_analysis"

ser_main <- bgwqs.fit(y = Y, x = X, z = Z, x.s = x.s, n.chains = 2, method = "quantile", n.quantiles = 4,
                          n.adapt  = 500, n.iter = 100000, n.burnin = 1000, n.thin = 1, working.dir = work_dir, DIC=T)
save(ser_main, file = "ser_main.RDA")
