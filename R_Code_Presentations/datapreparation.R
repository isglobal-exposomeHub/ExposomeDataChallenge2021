#exposomeID=as.numeric(exposome$ID)
#ord=order(exposomeID) #No need to adjusted because IDs are in order.
load("D:/Grants/NIEHS/datachallenge/exposome.Rdata")
#view(codebook)
#view(exposome)
#view(covariates)
#view(phenotype)
n=dim(exposome)[1]

phenotypeID=as.numeric(phenotype$ID)
ord=order(phenotypeID) #need to be adjusted because exposome and phenotype data are not aligned
bw=phenotype$e3_bw[ord]             # birth weight
asthma=phenotype$hs_asthma[ord]     # asthma (binary)
bmizscore=phenotype$hs_zbmi_who[ord] # 6-11 yrs old bmi z-score
beh=phenotype$hs_correct_raven[ord] # behavior score
IQ=phenotype$hs_Gen_Tot[ord]        # intelligence quotient
logIQ=log(IQ)
bmicat=as.numeric(phenotype$hs_bmi_c_cat[ord])  # factor, categorized BMI


#1. Birth weight 
#   Analysis using pregnancy exposure and covariate data

# prenatal metal exposure

metalm=array(0,c(n,10))
 metalm[,1]=exposome$hs_as_m_Log2
 metalm[,2]=exposome$hs_cd_m_Log2
 metalm[,3]=exposome$hs_co_m_Log2
 metalm[,4]=exposome$hs_cs_m_Log2
 metalm[,5]=exposome$hs_cu_m_Log2
 metalm[,6]=exposome$hs_hg_m_Log2
 metalm[,7]=exposome$hs_mn_m_Log2
 metalm[,8]=exposome$hs_mo_m_Log2
 metalm[,9]=exposome$hs_pb_m_Log2
 metalm[,10]=1.0*(exposome$hs_tl_mdich_None=="Detected")

# prenatal chemical exposure

chemicalm=array(0,c(n,38))
 chemicalm[,1]=exposome$hs_dde_madj_Log2
 chemicalm[,2]=exposome$hs_ddt_madj_Log2
 chemicalm[,3]=exposome$hs_hcb_madj_Log2
 chemicalm[,4]=exposome$hs_pcb118_madj_Log2
 chemicalm[,5]=exposome$hs_pcb138_madj_Log2
 chemicalm[,6]=exposome$hs_pcb153_madj_Log2
 chemicalm[,7]=exposome$hs_pcb170_madj_Log2
 chemicalm[,8]=exposome$hs_pcb180_madj_Log2
 chemicalm[,9]=exposome$hs_sumPCBs5_madj_Log2
 chemicalm[,10]=exposome$hs_dep_madj_Log2
 chemicalm[,11]=exposome$hs_detp_madj_Log2
 chemicalm[,12]=exposome$hs_dmp_madj_Log2
 chemicalm[,13]=exposome$hs_dmtp_madj_Log2
 chemicalm[,14]=exposome$hs_pbde153_madj_Log2
 chemicalm[,15]=exposome$hs_pbde47_madj_Log2
 chemicalm[,16]=exposome$hs_pfhxs_m_Log2
 chemicalm[,17]=exposome$hs_pfna_m_Log2
 chemicalm[,18]=exposome$hs_pfoa_m_Log2
 chemicalm[,19]=exposome$hs_pfos_m_Log2
 chemicalm[,20]=exposome$hs_pfunda_m_Log2
 chemicalm[,21]=exposome$hs_bpa_madj_Log2
 chemicalm[,22]=exposome$hs_bupa_madj_Log2
 chemicalm[,23]=exposome$hs_etpa_madj_Log2
 chemicalm[,24]=exposome$hs_mepa_madj_Log2
 chemicalm[,25]=exposome$hs_oxbe_madj_Log2
 chemicalm[,26]=exposome$hs_prpa_madj_Log2
 chemicalm[,27]=exposome$hs_trcs_madj_Log2
 chemicalm[,28]=exposome$hs_mbzp_madj_Log2
 chemicalm[,29]=exposome$hs_mecpp_madj_Log2
 chemicalm[,30]=exposome$hs_mehhp_madj_Log2
 chemicalm[,31]=exposome$hs_mehp_madj_Log2
 chemicalm[,32]=exposome$hs_meohp_madj_Log2
 chemicalm[,33]=exposome$hs_mep_madj_Log2
 chemicalm[,34]=exposome$hs_mibp_madj_Log2
 chemicalm[,35]=exposome$hs_mnbp_madj_Log2
 chemicalm[,36]=exposome$hs_ohminp_madj_Log2
 chemicalm[,37]=exposome$hs_oxominp_madj_Log2
 chemicalm[,38]=exposome$hs_sumDEHP_madj_Log2

# Other prenatal exposures

airm=array(0,c(n,4))
 airm[,1]=exposome$h_abs_ratio_preg_Log 
 airm[,2]=exposome$h_no2_ratio_preg_Log 
 airm[,3]=exposome$h_pm10_ratio_preg_None 
 airm[,4]=exposome$h_pm25_ratio_preg_None 
builtm=array(0,c(n,9))
 builtm[,1]=exposome$h_accesslines300_preg_dic0
 builtm[,2]=exposome$h_accesspoints300_preg_Log
 builtm[,3]=exposome$h_builtdens300_preg_Sqrt
 builtm[,4]=exposome$h_connind300_preg_Sqrt
 builtm[,5]=exposome$h_fdensity300_preg_Log
 builtm[,6]=exposome$h_frichness300_preg_None
 builtm[,7]=exposome$h_landuseshan300_preg_None
 builtm[,8]=exposome$h_popdens_preg_Sqrt
 builtm[,9]=exposome$h_walkability_mean_preg_None
lifestylem=array(0,c(n,12))
 lifestylem[,1]=as.numeric(exposome$e3_alcpreg_yn_None)    #factor  #alcohol
 lifestylem[,2]=as.numeric(exposome$h_cereal_preg_Ter)     #factor
 lifestylem[,3]=as.numeric(exposome$h_dairy_preg_Ter)      #factor 
 lifestylem[,4]=as.numeric(exposome$h_fastfood_preg_Ter)   #factor
 lifestylem[,5]=as.numeric(exposome$h_fish_preg_Ter)       #factor 
 lifestylem[,6]=as.numeric(exposome$h_folic_t1_None)       #factor
 lifestylem[,7]=as.numeric(exposome$h_fruit_preg_Ter)      #factor
 lifestylem[,8]=as.numeric(exposome$h_legume_preg_Ter)     #factor
 lifestylem[,9]=as.numeric(exposome$h_meat_preg_Ter)       #factor
 lifestylem[,10]=as.numeric(exposome$h_pamod_t3_None)      #factor
 lifestylem[,11]=as.numeric(exposome$h_pavig_t3_None)      #factor
 lifestylem[,12]=as.numeric(exposome$h_veg_preg_Ter)       #factor
outdoorm=array(0,c(n,13))
 outdoorm[,1]=exposome$h_humidity_preg_None
 outdoorm[,2]=exposome$h_pressure_preg_None
 outdoorm[,3]=exposome$h_temperature_preg_None
 outdoorm[,4]=as.numeric(exposome$h_blueyn300_preg_None)    # factor
 outdoorm[,5]=as.numeric(exposome$h_greenyn300_preg_None)   # factor
 outdoorm[,6]=exposome$h_ndvi100_preg_None
 outdoorm[,7]=exposome$h_lden_cat_preg_None
 outdoorm[,8]=exposome$h_distinvnear1_preg_Log
 outdoorm[,9]=exposome$h_trafload_preg_pow1over3
 outdoorm[,10]=exposome$h_trafnear_preg_pow1over3
 outdoorm[,11]=exposome$h_bro_preg_Log
 outdoorm[,12]=exposome$h_clf_preg_Log
 outdoorm[,13]=exposome$h_thm_preg_Log

#2. Outcomes for 6-11 years old 
#   Analysis using pregnancy and postnatal exposure and covariate data

# post-natal metal exposure
metalc=array(0,c(n,10))
 metalc[,1]=exposome$hs_as_c_Log2
 metalc[,2]=exposome$hs_cd_c_Log2
 metalc[,3]=exposome$hs_co_c_Log2
 metalc[,4]=exposome$hs_cs_c_Log2
 metalc[,5]=exposome$hs_cu_c_Log2
 metalc[,6]=exposome$hs_hg_c_Log2
 metalc[,7]=exposome$hs_mn_c_Log2
 metalc[,8]=exposome$hs_mo_c_Log2
 metalc[,9]=exposome$hs_pb_c_Log2
 metalc[,10]=1.0*(exposome$hs_tl_cdich_None=="Detected")

#post-natal chemical exposure
chemicalc=array(0,c(n,39))
 chemicalc[,1]=exposome$hs_dde_cadj_Log2
 chemicalc[,2]=exposome$hs_ddt_cadj_Log2
 chemicalc[,3]=exposome$hs_hcb_cadj_Log2
 chemicalc[,4]=exposome$hs_pcb118_cadj_Log2
 chemicalc[,5]=exposome$hs_pcb138_cadj_Log2
 chemicalc[,6]=exposome$hs_pcb153_cadj_Log2
 chemicalc[,7]=exposome$hs_pcb170_cadj_Log2
 chemicalc[,8]=exposome$hs_pcb180_cadj_Log2
 chemicalc[,9]=exposome$hs_sumPCBs5_cadj_Log2
 chemicalc[,10]=exposome$hs_dep_cadj_Log2
 chemicalc[,11]=exposome$hs_detp_cadj_Log2
 chemicalc[,12]=exposome$hs_dmp_cadj_Log2
 chemicalc[,13]=exposome$hs_dmtp_cadj_Log2
 chemicalc[,14]=exposome$hs_pbde153_cadj_Log2
 chemicalc[,15]=exposome$hs_pbde47_cadj_Log2
 chemicalc[,16]=exposome$hs_pfhxs_c_Log2
 chemicalc[,17]=exposome$hs_pfna_c_Log2
 chemicalc[,18]=exposome$hs_pfoa_c_Log2
 chemicalc[,19]=exposome$hs_pfos_c_Log2
 chemicalc[,20]=exposome$hs_pfunda_c_Log2
 chemicalc[,21]=exposome$hs_bpa_cadj_Log2
 chemicalc[,22]=exposome$hs_bupa_cadj_Log2
 chemicalc[,23]=exposome$hs_etpa_cadj_Log2
 chemicalc[,24]=exposome$hs_mepa_cadj_Log2
 chemicalc[,25]=exposome$hs_oxbe_cadj_Log2
 chemicalc[,26]=exposome$hs_prpa_cadj_Log2
 chemicalc[,27]=exposome$hs_trcs_cadj_Log2
 chemicalc[,28]=exposome$hs_mbzp_cadj_Log2
 chemicalc[,29]=exposome$hs_mecpp_cadj_Log2
 chemicalc[,30]=exposome$hs_mehhp_cadj_Log2
 chemicalc[,31]=exposome$hs_mehp_cadj_Log2
 chemicalc[,32]=exposome$hs_meohp_cadj_Log2
 chemicalc[,33]=exposome$hs_mep_cadj_Log2
 chemicalc[,34]=exposome$hs_mibp_cadj_Log2
 chemicalc[,35]=exposome$hs_mnbp_cadj_Log2
 chemicalc[,36]=exposome$hs_ohminp_cadj_Log2
 chemicalc[,37]=exposome$hs_oxominp_cadj_Log2
 chemicalc[,38]=exposome$hs_sumDEHP_cadj_Log2
 chemicalc[,39]=1.0*(exposome$hs_dmdtp_cdich_None=="detected")

# post-natal air exposures

airc=array(0,c(n,12))
 airc[,1]=exposome$hs_no2_dy_hs_h_Log
 airc[,2]=exposome$hs_no2_wk_hs_h_Log
 airc[,3]=exposome$hs_no2_yr_hs_h_Log
 airc[,4]=exposome$hs_pm10_dy_hs_h_None
 airc[,5]=exposome$hs_pm10_wk_hs_h_None
 airc[,6]=exposome$hs_pm10_yr_hs_h_None
 airc[,7]=exposome$hs_pm25_dy_hs_h_None
 airc[,8]=exposome$hs_pm25_wk_hs_h_None
 airc[,9]=exposome$hs_pm25_yr_hs_h_None
 airc[,10]=exposome$hs_pm25abs_dy_hs_h_Log
 airc[,11]=exposome$hs_pm25abs_wk_hs_h_Log
 airc[,12]=exposome$hs_pm25abs_yr_hs_h_Log

#post-natal built exposures  
builtc=array(0,c(n,15))
 builtc[,1]=exposome$hs_accesslines300_h_dic0
 builtc[,2]=exposome$hs_accesspoints300_h_Log
 builtc[,3]=exposome$hs_builtdens300_h_Sqrt
 builtc[,4]=exposome$hs_connind300_h_Log
 builtc[,5]=exposome$hs_fdensity300_h_Log
 builtc[,6]=exposome$hs_landuseshan300_h_None
 builtc[,7]=exposome$hs_popdens_h_Sqrt
 builtc[,8]=exposome$hs_walkability_mean_h_None
 builtc[,9]=exposome$hs_accesslines300_s_dic0
 builtc[,10]=exposome$hs_accesspoints300_s_Log
 builtc[,11]=exposome$hs_builtdens300_s_Sqrt
 builtc[,12]=exposome$hs_connind300_s_Log
 builtc[,13]=exposome$hs_fdensity300_s_Log
 builtc[,14]=exposome$hs_landuseshan300_s_None
 builtc[,15]=exposome$hs_popdens_s_Sqrt

# post-natal indoor air exposures
indoorc=array(0,c(n,5))
 indoorc[,1]=exposome$h_Absorbance_Log
 indoorc[,2]=exposome$h_Benzene_Log
 indoorc[,3]=exposome$h_NO2_Log
 indoorc[,4]=exposome$h_PM_Log
 indoorc[,5]=exposome$h_TEX_Log

# post-natal life style exposures
lifestylec=array(0,c(n,27))
 lifestylec[,1]=as.numeric(exposome$h_bfdur_Ter)              #factor  #alcohol
 lifestylec[,2]=as.numeric(exposome$hs_bakery_prod_Ter)       #factor
 lifestylec[,3]=as.numeric(exposome$hs_beverages_Ter)         #factor 
 lifestylec[,4]=as.numeric(exposome$hs_break_cer_Ter)         #factor
 lifestylec[,5]=as.numeric(exposome$hs_caff_drink_Ter)        #factor 
 lifestylec[,6]=as.numeric(exposome$hs_dairy_Ter)             #factor
 lifestylec[,7]=as.numeric(exposome$hs_fastfood_Ter)          #factor
 lifestylec[,8]=as.numeric(exposome$hs_KIDMED_None)           #factor
 lifestylec[,9]=as.numeric(exposome$hs_mvpa_prd_alt_None)     #factor
 lifestylec[,10]=as.numeric(exposome$hs_org_food_Ter)         #factor
 lifestylec[,11]=as.numeric(exposome$hs_pet_cat_r2_None)      #factor
 lifestylec[,12]=as.numeric(exposome$hs_pet_dog_r2_None)      #factor
 lifestylec[,13]=as.numeric(exposome$hs_pet_None)             #factor 
 lifestylec[,14]=as.numeric(exposome$hs_proc_meat_Ter)        #factor
 lifestylec[,15]=as.numeric(exposome$hs_readymade_Ter)        #factor 
 lifestylec[,16]=as.numeric(exposome$hs_sd_wk_None)           #factor
 lifestylec[,17]=as.numeric(exposome$hs_total_bread_Ter)      #factor
 lifestylec[,18]=as.numeric(exposome$hs_total_cereal_Ter)     #factor
 lifestylec[,19]=as.numeric(exposome$hs_total_fish_Ter)       #factor
 lifestylec[,20]=as.numeric(exposome$hs_total_fish_Ter)       #factor
 lifestylec[,21]=as.numeric(exposome$hs_total_lipids_Ter)     #factor  
 lifestylec[,22]=as.numeric(exposome$hs_total_meat_Ter)       #factor
 lifestylec[,23]=as.numeric(exposome$hs_total_potatoes_Ter)   #factor 
 lifestylec[,24]=as.numeric(exposome$hs_total_sweets_Ter)     #factor
 lifestylec[,25]=as.numeric(exposome$hs_total_veg_Ter)        #factor 
 lifestylec[,26]=as.numeric(exposome$hs_total_yog_Ter)        #factor
 lifestylec[,27]=as.numeric(exposome$hs_dif_hours_total_None) #factor

# post natal outdoor exposures
outdoorc=array(0,c(n,19))
 outdoorc[,1]=exposome$hs_hum_mt_hs_h_None
 outdoorc[,2]=exposome$hs_tm_mt_hs_h_None
 outdoorc[,3]=exposome$hs_uvdvf_mt_hs_h_None
 outdoorc[,4]=as.numeric(exposome$hs_hum_dy_hs_h_None)    # factor
 outdoorc[,5]=as.numeric(exposome$hs_hum_wk_hs_h_None)   # factor
 outdoorc[,6]=exposome$hs_tm_dy_hs_h_None
 outdoorc[,7]=exposome$hs_tm_wk_hs_h_None
 outdoorc[,8]=exposome$hs_uvdvf_dy_hs_h_None
 outdoorc[,9]=exposome$hs_uvdvf_wk_hs_h_None
 outdoorc[,10]=exposome$hs_blueyn300_s_None
 outdoorc[,11]=exposome$hs_greenyn300_s_None
 outdoorc[,12]=exposome$hs_blueyn300_h_None
 outdoorc[,13]=exposome$hs_greenyn300_h_None
 outdoorc[,14]=exposome$hs_ndvi100_h_None
 outdoorc[,15]=exposome$hs_ln_cat_h_None
 outdoorc[,16]=exposome$hs_ln_cat_h_None
 outdoorc[,17]=exposome$hs_lden_cat_s_None
 outdoorc[,18]=exposome$hs_trafload_h_pow1over3
 outdoorc[,19]=exposome$hs_trafnear_h_pow1over3

#special chemicals (post-natal)
schemc=array(0,c(n,7))
 schemc[,1]=exposome$FAS_cat_None
 schemc[,2]=exposome$hs_contactfam_3cat_num_None
 schemc[,3]=exposome$hs_hm_pers_None
 schemc[,4]=exposome$hs_participation_3cat_None
 schemc[,5]=exposome$hs_cotinine_cdich_None
 schemc[,6]=exposome$hs_globalexp2_None
 schemc[,7]=exposome$hs_smk_parents_None

## Covariate data for both pregnancy and postnatal

covariate=array(0,c(n,10))
 covariate[,1]=covariates$h_mbmi_None                 # maternal BMI
 covariate[,2]=covariates$hs_wgtgain_None             # weight gains during pregnancy
 covariate[,3]=covariates$e3_gac_None                 # Gestational age 
 covariate[,4]=as.numeric(covariates$e3_sex_None)     # factor (child sex)
 covariate[,5]=as.numeric(covariates$e3_yearbir_None) # factor (year of birth) 
 covariate[,6]=covariates$h_age_None                  # maternal age
 covariate[,7]=as.numeric(covariates$h_cohort)        # factor (cohort)
 covariate[,8]=as.numeric(covariates$h_edumc_None)    # factor (maternal education level)
 covariate[,9]=as.numeric(covariates$h_native_None)   # factor (native)
 covariate[,10]=as.numeric(covariates$h_parity_None)  # factor (parity)

#special covariates for children (post-natal)
 scovariate=array(0,c(n,3))
 scovariate[,1]=covariates$hs_c_height_None
 scovariate[,2]=covariates$hs_c_weight_None
 scovariate[,3]=covariates$hs_child_age_None


