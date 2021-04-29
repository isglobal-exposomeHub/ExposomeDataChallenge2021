#######################################################################################################################
# ISGlobal Exposome Project
# Title: Air pollution and childhood cognition: a g-computation approach 
# to assess mediation by a mixture of metals

# Created by: Charlie Roscoe, Hari Iyer, Huichu Li
# Contact: Hari_Iyer@dfci.harvard.edu
# April 18, 2021

# Objective: Create analysis dataset for quantile g-computation analysis to estimate the extent to which association
# between air pollution and cognition may be mediated by metals mixture

## 1) Prepare covariates, PM25 exposure, mixtures
## 2) Run table 1, correlation, boxplots
## 3) Fit outcome models, mediator models
## 4) Run quantile g-computation
## 5) Get controlled direct effects for association between air pollution and raven under interventions to change metals
## mixture

## Notes on notation: We sometimes refer to A, M, Y where A=exposure (air pollution),
## M=mediator (metals mixture), Y=outcome (Raven score)
#######################################################################################################################

# the following function checks whether the required packages have been installed
# if the package has been installed, then it will be loaded into the environment
# if not, then R will install and load the package

packages<-c("tidyverse",
            "magrittr",
            "qgcomp",
            "sandwich",
            "Amelia",
            "table1",
            "corrplot",
            "gtools",
            "haven",
            "gridExtra",
            "car",
            "boot",
            "foreign",
            "lmtest",
            "sandwich"
)

pkg.check<-lapply(packages,
                  FUN = function(x){
                    if(!require(x,character.only = TRUE)){
                      install.packages(x,dependencies = T)
                      library(x,character.only = T)
                    }
                  }      
)


# data
load("C:/Users/hphiy/Documents/Charlie Exposome Project/exposome.RData") ## change to yours

View(codebook)
colnames(exposome)
head(exposome)
colnames(phenotype)
head(phenotype)


#selecting variables for this analysis: covarites from children 
 df<-subset(exposome, 
            select=c("ID",
                     "hs_pm25_yr_hs_h_None", # exposure: PM2.5
                     "hs_as_c_Log2","hs_cd_c_Log2","hs_co_c_Log2","hs_cs_c_Log2","hs_cu_c_Log2","hs_hg_c_Log2","hs_mn_c_Log2","hs_mo_c_Log2","hs_pb_c_Log2",
                     "hs_tl_cdich_None","hs_smk_parents_None",
                     "hs_participation_3cat_None","FAS_cat_None",# SES; social participation and family affluence score
                     "hs_ndvi100_h_None", #greenness
                     "e3_asmokcigd_p_None", #maternal active smoking
                     "hs_tm_mt_hs_h_None", # temperature
                     "hs_total_fish_Ter", # fish intake, total (in tertiles)
                     "h_fish_preg_Ter", # fish intake preg mom
                     "hs_mvpa_prd_alt_None", # physical activity
                     "hs_popdens_h_Sqrt", # population density
                     "hs_pm10_yr_hs_h_None","hs_no2_yr_hs_h_Log", #PM`0 and NO2 in the previous year of examination`
                     "h_pm25_ratio_preg_None"))  
 
 
#THALLIUM IN METALS IS A FACTOR detected/undetected (NOT NUMERIC) 
#table(df$hs_tl_cdich_None)
#Detected Undetected 
#102       1199 

df_covs<-subset(covariates, select=c("ID", "h_age_None", "h_edumc_None", "hs_child_age_None", "e3_sex_None", "h_cohort", "h_mbmi_None","hs_c_height_None","hs_c_weight_None")) 
df_phen<-subset(phenotype, select=c("ID", "hs_correct_raven"))

#inner join 
df<-merge(df, df_covs)
df<-merge(df, df_phen)

summary(df)

## create quintiles of air pollution for data description
df$pm25.f<-factor(quantcut(df$hs_pm25_yr_hs_h_None,5),levels=c("[4.83,10.1]","(10.1,12.1]","(12.1,13.7]","(13.7,15.7]","(15.7,21.9]"),labels=c("Q1","Q2","Q3","Q4","Q5"))

df$e3_sex_None<-as.factor(df$e3_sex_None)
df$hs_tl_cdich_None<-as.factor(df$hs_tl_cdich_None)
df$h_cohort<-as.factor(df$h_cohort)
df$h_edumc_None<-factor(df$h_edumc_None,levels=c(1:3),labels=c("Primary school","Secondary school","University degree or above"))

# check missingness
isTRUE(is.na(df)) 
# returns FALSE, meaning no missingness...


label(df$hs_pm25_yr_hs_h_None)       <- "Annual average PM2.5"
label(df$hs_as_c_Log2)       <- "Arsenic"
label(df$hs_cd_c_Log2)       <- "Cadmium"
label(df$hs_co_c_Log2)       <- "Cobalt"
label(df$hs_cs_c_Log2)       <- "Caesium"
label(df$hs_cu_c_Log2)       <- "	Copper"
label(df$hs_hg_c_Log2)       <- "Mercury"
label(df$hs_mn_c_Log2)       <- "Manganese"
label(df$hs_mo_c_Log2)       <- "Molybdenum"
label(df$hs_pb_c_Log2)       <- "Lead"
label(df$hs_tl_cdich_None)       <- "Thallium"
label(df$hs_correct_raven)       <- "RAVEN score"
label(df$hs_child_age_None) <- "Age at assessment"
label(df$pm25.f) <- "Quintile of PM2.5"
label(df$e3_sex_None) <- "Sex"
label(df$h_cohort) <- "Cohort no."
label(df$h_mbmi_None) <- "child BMI"
label(df$hs_ndvi100_h_None) <- "NDVI at Home"
label(df$e3_asmokcigd_p_None) <- "Parent Smoking"
label(df$FAS_cat_None) <- "Family Affluence"
label(df$hs_tm_mt_hs_h_None)<-"Monthly temperature"
label(df$hs_total_fish_Ter)<-"Total fish consumption"
label(df$hs_mvpa_prd_alt_None)<-"Moderate to Vigorous physical activity (min/day)"
label(df$hs_popdens_h_Sqrt)<-"Population density"
label(df$h_age_None)<-"Maternal age"
label(df$h_edumc_None)<-"Maternal education"


# scale PM2.5 by 10 ug/m3 
df$pm25_10ug <- df$hs_pm25_yr_hs_h_None/10
  
######################################################################
# TABLE 1
######################################################################

table1(~hs_pm25_yr_hs_h_None+ 
       hs_child_age_None + e3_sex_None +h_cohort + h_mbmi_None +
         hs_ndvi100_h_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter +
         hs_mvpa_prd_alt_None + h_age_None + h_edumc_None + hs_correct_raven | pm25.f, data=df)

hist(df$hs_correct_raven)
df$ln_raven <- log(df$hs_correct_raven)

##################################################################################
# CORRELATIONS
##################################################################################

df_numeric<-df[,c("hs_pm25_yr_hs_h_None", "hs_as_c_Log2","hs_cd_c_Log2","hs_co_c_Log2","hs_cs_c_Log2","hs_cu_c_Log2","hs_hg_c_Log2","hs_mn_c_Log2","hs_mo_c_Log2","hs_pb_c_Log2")]
colnames(df_numeric)<-c("PM2.5","As","Cd","Co","Cs","Cu","Hg","Mn","Mo","Pb")
M <- cor(df_numeric,method="spearman")
corrplot(M, type="upper", method = "circle",tl.col="black",tl.srt=45)

#child age is -ve assoc with AP, strange, should be mindful of this in adj. ests


####################################################################################
## Step 1: Models for X-Y associaiton

## Note: adjusting for site obliterates PM2.5 > raven association. For this analysis,
## we will not adjust for site. This probably is a tricky issue because variability 
## in air pollution exposure is likely strongly correlated with site


## Checked normality - it's not great, but range is wide enough and histogram 
## isn't too bad. Does not meet Poisson (mean ne var, distribution is more normal) 

##> var(df$hs_correct_raven)
##[1] 41.53608
##> mean(df$hs_correct_raven)
##[1] 26.28824
##> median(df$hs_correct_raven)
##[1] 27
####################################################################################

## Univariable model
mod<-glm(formula = hs_correct_raven ~ pm25_10ug, data=df)
summary(mod)
confint(mod)

## Adjust for age and sex
mod2<-glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None, data=df)
summary(mod2)
confint(mod2)


## Additonally adjust for maternal BMI, smoking, SES, age, education, and residential NDVI
mod3<-glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
           hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + h_age_None + h_edumc_None, data=df)
summary(mod3)
confint(mod3)

## Additonally adjust for temperature, fish consumption, physical activity, and population density
mod4<-glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
           hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
            hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None, data=df)
summary(mod4)


anova(mod, mod4, test="Chisq") ## model with all covariates is better

## WITH COHORT
mod5<-glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None + h_cohort +
            hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
            hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None, data=df)
summary(mod5) ## not statistically significant at 0.05
confint(mod5)

######################################
# FINAL MODEL: MODEL 4
######################################

mod4<-glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
            hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
            hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None, data=df)
summary(mod4)
confint(mod4)


## Coefficients:
##   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                             7.006243   1.992633   3.516 0.000453 ***
##  pm25_10ug                              -2.048009   0.473817  -4.322 1.66e-05 ***
##  hs_child_age_None                       2.341523   0.108587  21.564  < 2e-16 ***
##  e3_sex_Nonemale                        -0.165106   0.275302  -0.600 0.548794    
##  h_mbmi_None                             0.059303   0.025174   2.356 0.018635 *  
##  hs_ndvi100_h_None                       0.435172   1.033773   0.421 0.673859    
##  e3_asmokcigd_p_None                    -0.040497   0.071652  -0.565 0.572041    
## FAS_cat_NoneMiddle                      1.286438   0.435345   2.955 0.003184 ** 
## FAS_cat_NoneHigh                        1.076131   0.452676   2.377 0.017587 *  
##  hs_tm_mt_hs_h_None                     -0.115370   0.019077  -6.048 1.93e-09 ***
##  hs_total_fish_Ter(1.5,3]                0.010238   0.316747   0.032 0.974220    
## hs_total_fish_Ter(3,Inf]                0.481662   0.332144   1.450 0.147259    
## hs_mvpa_prd_alt_None                   -0.011435   0.006496  -1.760 0.078598 .  
## hs_popdens_h_Sqrt                       0.013158   0.003248   4.051 5.40e-05 ***
##  h_age_None                             -0.009961   0.026788  -0.372 0.710076    
## h_edumc_NoneSecondary school            1.182750   0.404248   2.926 0.003496 ** 
##  h_edumc_NoneUniversity degree or above  2.419217   0.412070   5.871 5.51e-09 ***


#########################################################################################
##
## Step 2: M-Y association
##
#########################################################################################

## save name of mixture vriables in the variable "Xnm"

Xnm <- c('hs_as_c_Log2','hs_cd_c_Log2','hs_co_c_Log2',
         'hs_cs_c_Log2','hs_cu_c_Log2','hs_hg_c_Log2',
         'hs_mn_c_Log2','hs_mo_c_Log2','hs_pb_c_Log2')

covars = c('hs_pm25_yr_hs_h_None', 'hs_child_age_None','e3_sex_None','h_mbmi_None','hs_ndvi100_h_None',
           'e3_asmokcigd_p_None','FAS_cat_None','hs_tm_mt_hs_h_None','hs_total_fish_Ter','h_fish_preg_Ter',
           'hs_mvpa_prd_alt_None','hs_popdens_h_Sqrt','h_age_None','h_edumc_None')

nonessentialXnm <- c('hs_as_c_Log2','hs_cd_c_Log2','hs_co_c_Log2','hs_cs_c_Log2','hs_hg_c_Log2',"hs_mo_c_Log2","hs_pb_c_Log2")
essentialXnm <- c('hs_cu_c_Log2','hs_mn_c_Log2')

## just mixture and outcome: the default is by quartile (q=4)
qc.fit <- qgcomp.noboot(hs_correct_raven ~.,dat=df[,c(Xnm, 'hs_correct_raven')], family=gaussian())

qc.fit

## The overall mixture effect from quantile g-computation (psi1) is interpreted as the effect on the outcome of increasing every exposure 
## by one quantile, possibly conditional on covariates. 
## Given the overall exposure effect, the weights are considered fixed and so do not have confidence intervals or p-values.

## Scaled effect size (positive direction, sum of positive coefficients = 3.36)
## hs_cs_c_Log2 hs_hg_c_Log2 hs_as_c_Log2 hs_cd_c_Log2 hs_mn_c_Log2 
## 0.7062       0.1691       0.0568       0.0355       0.0325 

##  Scaled effect size (negative direction, sum of negative coefficients = -2.6)
## hs_cu_c_Log2 hs_co_c_Log2 hs_mo_c_Log2 hs_pb_c_Log2 
## 0.307        0.287        0.210        0.196 

## Mixture slope parameters (Delta method CI):
  
##  Estimate Std. Error  Lower CI Upper CI t value Pr(>|t|)
## (Intercept) 25.16675    0.52694 24.133962  26.1995 47.7599   <2e-16
## psi1         0.76255    0.33839  0.099305   1.4258  2.2534   0.0244
##
##

## model with covariates and air pollution

qcfit.2 <- qgcomp.noboot(hs_correct_raven ~ hs_pm25_yr_hs_h_None + hs_child_age_None + e3_sex_None + h_mbmi_None + hs_ndvi100_h_None +
              e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
              hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + 
              hs_as_c_Log2 + hs_cd_c_Log2 + hs_co_c_Log2 +
              hs_cs_c_Log2 + hs_cu_c_Log2 + hs_hg_c_Log2 +
              hs_mn_c_Log2 + hs_mo_c_Log2 + hs_pb_c_Log2,
              expnms=Xnm,
              df, family=gaussian(), q=4)

qcfit.2
plot(qcfit.2)


## switched sign

## Scaled effect size (positive direction, sum of positive coefficients = 0.537)
## hs_cs_c_Log2 hs_as_c_Log2 
## 0.701        0.299 

## Scaled effect size (negative direction, sum of negative coefficients = -1.25)
## hs_cu_c_Log2 hs_mo_c_Log2 hs_co_c_Log2 hs_pb_c_Log2 hs_hg_c_Log2 hs_cd_c_Log2 hs_mn_c_Log2 
## 0.2624       0.1982       0.1666       0.1504       0.1094       0.0785       0.0346 

## Mixture slope parameters (Delta method CI):
  
##  Estimate Std. Error Lower CI Upper CI t value  Pr(>|t|)
## (Intercept)  8.50210    2.07706   4.4311 12.57305  4.0933 4.516e-05
## psi1        -0.71008    0.28761  -1.2738 -0.14638 -2.4689   0.01368

## Interpretation - for a one quantile increase jointly across metals, there is a -0.74986 unit change (95% CI: -1.314, -0.18576) in
## raven score - so increasing metals mixture leads to reduction in raven score (as expected)


## split based on mixture

# 40/60% training/validation split
set.seed(1231124)
trainidx <- sample(1:nrow(df), round(nrow(df)*0.4))
valididx <- setdiff(1:nrow(df),trainidx)
traindata <- df[trainidx,]
validdata <- df[valididx,]
dim(traindata) # 181 observations = 40% of total


splitres <- qgcomp.partials(fun="qgcomp.noboot", f=hs_correct_raven~., q=4, 
                            traindata=traindata[,c(Xnm, covars, "hs_correct_raven")],validdata=validdata[,c(Xnm, covars, "hs_correct_raven")], expnms=Xnm)
splitres

## Variables with positive effect sizes in training data: hs_hg_c_Log2, hs_cs_c_Log2
## Variables with negative effect sizes in training data: hs_pb_c_Log2, hs_mo_c_Log2, hs_cu_c_Log2, hs_cd_c_Log2, hs_co_c_Log2, hs_as_c_Log2, hs_mn_c_Log2
## Partial effect sizes estimated in validation data
## Positive direction Mixture slope parameters (Delta method CI):
  
##  Estimate Std. Error Lower CI Upper CI t value Pr(>|t|)
## (Intercept) 19.78598    7.29631  5.48548 34.08647  2.7118  0.00684
##psi1         0.10317    0.25824 -0.40297  0.60931  0.3995  0.68961

## Negative direction Mixture slope parameters (Delta method CI):
  
##  Estimate Std. Error Lower CI  Upper CI t value Pr(>|t|)
## (Intercept)  6.92690    2.76832   1.5011 12.352697  2.5022  0.01255
## psi1        -0.76895    0.34503  -1.4452 -0.092707 -2.2287  0.02612

plot(splitres$pos.fit)
plot(splitres$neg.fit)


## NON-ESSENTIAL ONLY

qc.fit.nonessential <- qgcomp.noboot(hs_correct_raven~.,dat=df[,c(Xnm, covars, 'hs_correct_raven')], expnms=nonessentialXnm, family=gaussian())
qc.fit.nonessential
## Scaled effect size (positive direction, sum of positive coefficients = 0.545)
## hs_cs_c_Log2 hs_as_c_Log2 
## 0.695        0.305 

## Scaled effect size (negative direction, sum of negative coefficients = -0.854)
## hs_mo_c_Log2 hs_co_c_Log2 hs_pb_c_Log2 hs_hg_c_Log2 hs_cd_c_Log2 
## 0.285        0.227        0.209        0.166        0.112 

## Mixture slope parameters (Delta method CI):
  
##  Estimate Std. Error Lower CI Upper CI t value  Pr(>|t|)
## (Intercept) 22.90552    5.77050 11.59554  34.2155  3.9694 7.601e-05
## psi1        -0.30889    0.27918 -0.85608   0.2383 -1.1064    0.2688

## ESSENTIAL

qc.fit.essential <- qgcomp.noboot(hs_correct_raven~.,dat=df[,c(Xnm, covars, 'hs_correct_raven')], expnms=essentialXnm, family=gaussian())
qc.fit.essential
## Scaled effect size (positive direction, sum of positive coefficients = 0)
## None

## Scaled effect size (negative direction, sum of negative coefficients = -0.378)
## hs_cu_c_Log2 hs_mn_c_Log2 
## 0.854        0.146 

## Mixture slope parameters (Delta method CI):
  
##  Estimate Std. Error Lower CI  Upper CI t value  Pr(>|t|)
## (Intercept)  7.37020    2.11895  3.21714 11.523265  3.4782 0.0005214
## psi1        -0.37777    0.15422 -0.68004 -0.075496 -2.4495 0.0144374

## This entire mixture effect is driven primarily by adverse effects of copper

plot(qc.fit.essential)
plot(qc.fit.nonessential)

############################################################################################
# exposure-outcome association with and without adjusting for mediators
# Mediators were transformed to quartile 
# We computed both conditional X-Y association (holding all covariates constant) and marginal estimates using g-computation
############################################################################################

## conditional estimates: using the quartile exposure matrix from qgcomp
mixq <- as.data.frame(qcfit.2$qx)

df <- df %>%
  mutate(
    surrID = row_number()
  )


mixq <- mixq %>%
  mutate(
    surrID = row_number()
  )

df_q <- df %>%
  full_join(mixq, by="surrID")


# pm2.5 and raven score no mixture
mod_cde_nm<-glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
              hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
              hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None, data=df)
summary(mod_cde_nm)

## pm25_10ug                              -2.062137   0.472303  -4.366 1.37e-05 ***

## with mixture
mod_cde_m<-glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
                  hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
                  hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + hs_as_c_Log2_q +
                 hs_cd_c_Log2_q + hs_co_c_Log2_q + hs_cs_c_Log2_q +
                 hs_cu_c_Log2_q + hs_hg_c_Log2_q + hs_mn_c_Log2_q + hs_mo_c_Log2_q + hs_pb_c_Log2_q, data=df_q)
summary(mod_cde_m) 

## pm25_10ug                              -1.592215   0.508219  -3.133  0.00177 **  

## there is evidence of mediation by metals mixture because the point estimate has attenuated.

## G-computation for marginal estimates of the X-Y association

#################################################################################
# Exposure-mediator interactions - are 9 interactions significant?
# Likelihood ratio test comparing models with or without interaction
#################################################################################

## M-Y #1 (no interactions between mixture Qs and outcome)
mod_my1 <- glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
                 hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
                 hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + 
                 hs_as_c_Log2_q + hs_cd_c_Log2_q + hs_co_c_Log2_q +
                 hs_cs_c_Log2_q + hs_cu_c_Log2_q + hs_hg_c_Log2_q +
                 hs_mn_c_Log2_q + hs_mo_c_Log2_q + hs_pb_c_Log2_q, 
               data=onesample)
summary(mod_my1)

## M-Y #2 (interactions between mixture Qs and outcome)
mod_my2 <- glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
                 hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
                 hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + 
                 hs_as_c_Log2_q + hs_cd_c_Log2_q + hs_co_c_Log2_q +
                 hs_cs_c_Log2_q + hs_cu_c_Log2_q + hs_hg_c_Log2_q +
                 hs_mn_c_Log2_q + hs_mo_c_Log2_q + hs_pb_c_Log2_q +
                 hs_as_c_Log2_q*pm25_10ug + hs_cd_c_Log2_q*pm25_10ug + hs_co_c_Log2_q*pm25_10ug +
                 hs_cs_c_Log2_q*pm25_10ug + hs_cu_c_Log2_q*pm25_10ug + hs_hg_c_Log2_q*pm25_10ug +
                 hs_mn_c_Log2_q*pm25_10ug + hs_mo_c_Log2_q*pm25_10ug + hs_pb_c_Log2_q*pm25_10ug, 
               data=onesample)
summary(mod_my2)

lrtest(mod_my1, mod_my2) ## 9 df test
##   #Df  LogLik Df  Chisq Pr(>Chisq)  
## 1  29 -3744.1                       
## 2  38 -3733.8  9 20.626    0.01442 *

## Suggests we should include interactions


########################################################################
# Estimate CDE with bootstrap
#
## So here, we are doing the following:
## Effect of shifting air pollution from Q1 to Q3 of PM2.5
## Note: these are scaled to be consistent with our effect corresponding to 10ug/m3 change in PM2.5
# 
########################################################################

summary(df$pm25_10ug)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.4829  1.0410  1.3110  1.2916  1.5122  2.1917 

########################################################################
## Step 1. Estimate effect of changing PM2.5 from 4q to 1q on raven score. No intervention on the intermediators 
## function to calculate difference in means 
########################################################################

standardization_obs <- function(data, indices) {
  # create a dataset with 3 copies of each subject
  d <- data[indices,] # 1st copy: equal to original one`
  d$interv <- -1
  d0 <- d # 2nd copy: treatment set to 0, outcome to missing
  d0$interv <- 0
  d0$pm25_10ug <- 1.0410 ## PM2.5 exposure for 25p in study population (scaled by 10)
  d0$hs_correct_raven <- NA
  d1 <- d # 3rd copy: treatment set to 1, outcome to missing
  d1$interv <- 1
  d1$pm25_10ug <- 1.5122 ## PM2.5 exposure for 75p in study population (scaled by 10)
  d1$hs_correct_raven <- NA
  d.onesample <- rbind(d, d0, d1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (interv= -1)
  # parameter estimates are used to predict mean outcome for observations with set 
  # treatment (interv=0 and interv=1)
  fit <- glm(formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
               hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
               hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None, data=d.onesample)
  
  d.onesample$predicted_meanY <- predict(fit, d.onesample)
  
  # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
  return(c(mean(d.onesample$predicted_meanY[d.onesample$interv==-1]),
           mean(d.onesample$predicted_meanY[d.onesample$interv==0]),
           mean(d.onesample$predicted_meanY[d.onesample$interv==1]),
           mean(d.onesample$predicted_meanY[d.onesample$interv==1])-
             mean(d.onesample$predicted_meanY[d.onesample$interv==0])))
}

# bootstrap to get 95% confidence intervals
set.seed(123)
results <- boot(data=df_q, statistic=standardization_obs, R=1000)

getCI <- function(indat){
  # generating confidence intervals
  se <- c(sd(indat$t[,1]), sd(indat$t[,2]), 
          sd(indat$t[,3]), sd(indat$t[,4]))
  Estimate <- indat$t0
  ll <- Estimate - qnorm(0.975)*se
  ul <- Estimate + qnorm(0.975)*se
  
  V1<-c("Observed", "PM25 = 25p", "PM25 = 75p", 
        "Diff PM25 75p-25p")
  
  f.out<- data.frame(cbind.data.frame(V1, Estimate, se, ll, ul))
  return(f.out)
}

no_intervention<-getCI(results)

no_intervention

## V1           Estimate                se                ll                ul
## 1          Observed   26.2882398155265 0.178953646367018   25.937497113745   26.638982517308
## 2        PM25 = 25p   26.8049843613093 0.214932645722387  26.3837241165915  27.2262446060271
## 3        PM25 = 75p   25.8333053133702 0.193902504264449  25.4532633884997  26.2133472382406
## 4 Diff PM25 75p-25p -0.971679047939137 0.226194443362292 -1.41501201043231 -0.52834608544596


######################################################################################################
## Step 2. Estimate CDEs for change in PM2.5, fixing mixtures to each quartile
## with and without PM2.5*mixtures interaction
######################################################################################################

# creating a function using g-computation to estimate CDE.
# Adding the 'formula' statement so the user can define their own model
# Can incoporate glm family (?)
# Can incorporate models w/ and w/o X-M interaction (user defined)
# Can create quantile data (user defined)

gcomp_cde<-function(formula=NULL, # glm formula
                    data, #input dataset
                    indices, # do not define anything...
                    family=gaussian, #family statement for GLM, guassian by default
                    predict.type="link", # type pf prediction to be passed to 'predict()' function. Alternatives are 'response' and 'terms'
                    quantile=F, #indicator of whether the exposure mixture has been converted into quantiless, if=FALSE then we need to make quantiles in this function
                    q=4,# quantiles
                    m_names=NULL, # names of mixture (as mediators)
                    x_name=NULL, #variable name of the exposure
                    y_name=NULL, #name of the outcome variable
                    m_fixed_quantile = 1, #quantile level of the mediator to be fixed at (min=1)
                    x_intervention=c(0,1), #a vector of contrast for exposure
                    seed=123,# set random seed, default 123
                    R=1000){
  set.seed(seed)
  #make quantile data
  if(!isTRUE(quantile)){
    data[,m_names]<-data.frame(apply(data[,m_names],2,function(x) as.numeric(quantcut(x,q,na.rm=T))))
  }
  f<-formula(formula)
  # starting the model
  fun1<-function(data,indices){
    # create a dataset with 3 copies of each subject
    d <- data[indices,] # 1st copy: equal to original one`
    d$interv <- -1
    
    d0 <- d # 2nd copy: treatment set to 0, outcome to missing
    d0$interv <- 0
    d0[,x_name]<- x_intervention[1]
    d0[,m_names]<-m_fixed_quantile
    d0[,y_name] <- NA
    
    d1 <- d # 3rd copy: treatment set to 1, outcome to missing
    d1$interv <- 1
    d1[,x_name]<- x_intervention[2]
    d1[,m_names]<-m_fixed_quantile
    d1[,y_name] <- NA
    
    d.onesample <- rbind(d, d0, d1) # combining datasets
    
    fit <- glm(formula = f,data=d.onesample,family=family)
    
    d.onesample$predicted_meanY <- predict(fit, d.onesample, type=predict.type)
    
    # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
    return(c(mean(d.onesample$predicted_meanY[d.onesample$interv==-1]),
             mean(d.onesample$predicted_meanY[d.onesample$interv==0]),
             mean(d.onesample$predicted_meanY[d.onesample$interv==1]),
             mean(d.onesample$predicted_meanY[d.onesample$interv==1])-
               mean(d.onesample$predicted_meanY[d.onesample$interv==0])))
  }
  boot_out<-boot(data=data,statistic = fun1,R=R)
  
  getCI <- function(indat){
    # generating confidence intervals
    se <- c(sd(indat$t[,1]), sd(indat$t[,2]), 
            sd(indat$t[,3]), sd(indat$t[,4]))
    Estimate <- indat$t0
    ll <- Estimate - qnorm(0.975)*se
    ul <- Estimate + qnorm(0.975)*se
    
    V1<- c("Observed", "PM25 = 25p", "PM25 = 75p", "Diff PM25 75p-25p")
    f.out<- data.frame(cbind.data.frame(V1, Estimate, se, ll, ul))
    
    return(f.out)
  }
  cde_estimates<-getCI(boot_out)
  list.out<-list(cde_estimates,boot_out)
  names(list.out)<-c("cde_estimates","boot_out")
  return(list.out)
} 


#######################################################################################
# CDE - No PM2.5*mixtures interaction
#######################################################################################

## usage: computing difference of raven score if we reduce PM2.5 from the 75th to the 25th percentile, while keeping blood metal levels at Q1-Q4, respectively
## CDE without interaction: should be the same for all quantiles as we assumed no interactions

cde_1<-gcomp_cde(data=df,
                 formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
                 hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
                 hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + 
                 hs_as_c_Log2 + hs_cd_c_Log2 + hs_co_c_Log2 + 
                 hs_cs_c_Log2 + hs_cu_c_Log2 + hs_hg_c_Log2 + 
                 hs_mn_c_Log2 + hs_mo_c_Log2 + hs_pb_c_Log2,
                 indices,
                 m_names=Xnm,
                 x_name="pm25_10ug",
                 y_name="hs_correct_raven",
                 m_fixed_quantile = 1, #quantile level of the mediator to be fixed at 
                 x_intervention=quantile(df$pm25_10ug,probs=c(0.25,0.75)), #a vector of contrast for exposure
                 R=1000)
cde_1$cde_estimates


########################################################################################################
## CDE with PM2.5*mixtures interaction
########################################################################################################

## CDE with interaction
cde_int_1<-gcomp_cde(data=df,
                     formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
                       hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
                       hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + 
                       hs_as_c_Log2 + hs_cd_c_Log2 + hs_co_c_Log2 + 
                       hs_cs_c_Log2 + hs_cu_c_Log2 + hs_hg_c_Log2 + 
                       hs_mn_c_Log2 + hs_mo_c_Log2 + hs_pb_c_Log2 +
                       hs_as_c_Log2*pm25_10ug + hs_cd_c_Log2*pm25_10ug + hs_co_c_Log2*pm25_10ug +
                       hs_cs_c_Log2*pm25_10ug + hs_cu_c_Log2*pm25_10ug + hs_hg_c_Log2*pm25_10ug +
                       hs_mn_c_Log2*pm25_10ug + hs_mo_c_Log2*pm25_10ug + hs_pb_c_Log2*pm25_10ug,
                     indices,
                     m_names=Xnm,
                     x_name="pm25_10ug",
                     y_name="hs_correct_raven",
                     m_fixed_quantile = 1, #quantile level of the mediator to be fixed at 
                     x_intervention=quantile(df$pm25_10ug,probs=c(0.25,0.75)), #a vector of contrast for exposure
                     R=1000)
cde_int_1$cde_estimates


cde_int_2<-gcomp_cde(data=df,
                     formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
                       hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
                       hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + 
                       hs_as_c_Log2 + hs_cd_c_Log2 + hs_co_c_Log2 + 
                       hs_cs_c_Log2 + hs_cu_c_Log2 + hs_hg_c_Log2 + 
                       hs_mn_c_Log2 + hs_mo_c_Log2 + hs_pb_c_Log2 +
                       hs_as_c_Log2*pm25_10ug + hs_cd_c_Log2*pm25_10ug + hs_co_c_Log2*pm25_10ug +
                       hs_cs_c_Log2*pm25_10ug + hs_cu_c_Log2*pm25_10ug + hs_hg_c_Log2*pm25_10ug +
                       hs_mn_c_Log2*pm25_10ug + hs_mo_c_Log2*pm25_10ug + hs_pb_c_Log2*pm25_10ug,
                     indices,
                     m_names=Xnm,
                     x_name="pm25_10ug",
                     y_name="hs_correct_raven",
                     m_fixed_quantile = 2, #quantile level of the mediator to be fixed at 
                     x_intervention=quantile(df$pm25_10ug,probs=c(0.25,0.75)), #a vector of contrast for exposure
                     R=1000)
cde_int_2$cde_estimates


cde_int_3<-gcomp_cde(data=df,
                     formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
                       hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
                       hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + 
                       hs_as_c_Log2 + hs_cd_c_Log2 + hs_co_c_Log2 + 
                       hs_cs_c_Log2 + hs_cu_c_Log2 + hs_hg_c_Log2 + 
                       hs_mn_c_Log2 + hs_mo_c_Log2 + hs_pb_c_Log2 +
                       hs_as_c_Log2*pm25_10ug + hs_cd_c_Log2*pm25_10ug + hs_co_c_Log2*pm25_10ug +
                       hs_cs_c_Log2*pm25_10ug + hs_cu_c_Log2*pm25_10ug + hs_hg_c_Log2*pm25_10ug +
                       hs_mn_c_Log2*pm25_10ug + hs_mo_c_Log2*pm25_10ug + hs_pb_c_Log2*pm25_10ug,
                     indices,
                     m_names=Xnm,
                     x_name="pm25_10ug",
                     y_name="hs_correct_raven",
                     m_fixed_quantile = 3, #quantile level of the mediator to be fixed at 
                     x_intervention=quantile(df$pm25_10ug,probs=c(0.25,0.75)), #a vector of contrast for exposure
                     R=1000)
cde_int_3$cde_estimates


cde_int_4<-gcomp_cde(data=df,
                     formula = hs_correct_raven ~ pm25_10ug + hs_child_age_None + e3_sex_None + h_mbmi_None +
                       hs_ndvi100_h_None + e3_asmokcigd_p_None + FAS_cat_None + hs_tm_mt_hs_h_None + hs_total_fish_Ter + h_fish_preg_Ter +
                       hs_mvpa_prd_alt_None + hs_popdens_h_Sqrt + h_age_None + h_edumc_None + 
                       hs_as_c_Log2 + hs_cd_c_Log2 + hs_co_c_Log2 + 
                       hs_cs_c_Log2 + hs_cu_c_Log2 + hs_hg_c_Log2 + 
                       hs_mn_c_Log2 + hs_mo_c_Log2 + hs_pb_c_Log2 +
                       hs_as_c_Log2*pm25_10ug + hs_cd_c_Log2*pm25_10ug + hs_co_c_Log2*pm25_10ug +
                       hs_cs_c_Log2*pm25_10ug + hs_cu_c_Log2*pm25_10ug + hs_hg_c_Log2*pm25_10ug +
                       hs_mn_c_Log2*pm25_10ug + hs_mo_c_Log2*pm25_10ug + hs_pb_c_Log2*pm25_10ug,
                     indices,
                     m_names=Xnm,
                     x_name="pm25_10ug",
                     y_name="hs_correct_raven",
                     m_fixed_quantile = 4, #quantile level of the mediator to be fixed at 
                     x_intervention=quantile(df$pm25_10ug,probs=c(0.25,0.75)), #a vector of contrast for exposure
                     R=1000)
cde_int_4$cde_estimates


######################################################################################################################
# Plotting CDEs
######################################################################################################################

## some housekeeping

prepPlot <- function(indat, lab){
  indat$cde_type <- lab
  indat<-indat[c(2,3),-c(3)]
  return(indat)
}


labvec <- c("NoInterv","NoIntx","IntxQ1","IntxQ2","IntxQ3","IntxQ4")


plotlist <- list(no_intervention,cde_1[[1]],cde_int_1[[1]],cde_int_2[[1]],cde_int_3[[1]],cde_int_4[[1]])
for(i in seq_along(plotlist)){
  plotlist[[i]]<-prepPlot(plotlist[[i]],labvec[i])
}

ls_plot<-dplyr::bind_rows(plotlist)

ls_plot$int.f <- factor(ls_plot$cde_type, levels=c("NoInterv","NoIntx","IntxQ1","IntxQ2","IntxQ3","IntxQ4"),
                        labels=c("NoInterv", "NoMetalsIntx","IntxQ1","IntxQ2","IntxQ3","IntxQ4"))

## Plot

dodge <- position_dodge(width=0.5)
p1 <- ggplot(data=ls_plot, aes(x=int.f, y=Estimate, ymin=ll, ymax=ul, colour=V1)) + 
  geom_pointrange(position=dodge) +
  xlab("Hypothetical Interventions on Metals") + ylab("Counterfactual Raven Score under Interventions") +
  labs(color="Intervention on PM2.5") +
  theme_bw() +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size = 15), legend.position = "bottom") 

#######################################################################################################
# Table with full CDE results
#######################################################################################################

## create list of all outputs

plotlist_f <- list(no_intervention,cde_1[[1]],cde_int_1[[1]],cde_int_2[[1]],cde_int_3[[1]],cde_int_4[[1]])
ls_tab <- dplyr::bind_rows(plotlist_f)


## it looks like there were interaction between PM2.5 and metals 
## Greater impaovement for kids with the higher blood metal levels if we intervene PM2.5

## calcuating PE (proportion elimiated)
# PE=(difference w/o intervention-difference w/ intervention)/(difference w/o intervention)

#difference without intervention was calculated in line 764 using the standardization_obs funciton

#PE point estimate
pe_est<-((results$t0[4]-cde_1$boot_out$t0[4])/results$t0[4])*100

# 95% confidence interval: do we need this or not? what about those negative values
# pe_nointeract<-data.frame(cbind(results$t[,4],cde_1$boot_out$t[,4]))
# names(pe_nointeract)<-c("no_intervention","with_intervention")
# pe_nointeract$pe<-(pe_nointeract[,1]-pe_nointeract[,2])/pe_nointeract[,1]
# 
# quantile(pe_nointeract$pe,probs=c(0.025,0.975))

# PE by metal quantile and with PM-metal interaction: not sure if this is right or meaningful...
pe_est_q1<-((results$t0[4]-cde_int_1$boot_out$t0[4])/results$t0[4])*100
pe_est_q2<-((results$t0[4]-cde_int_2$boot_out$t0[4])/results$t0[4])*100
pe_est_q3<-((results$t0[4]-cde_int_3$boot_out$t0[4])/results$t0[4])*100
pe_est_q4<-((results$t0[4]-cde_int_4$boot_out$t0[4])/results$t0[4])*100

######################################################################
# Supplemental slides
# Air Pollutants by Quintile of PM2.5
######################################################################

p_as <- df %>%
  ggplot(aes(x=pm25.f, y=hs_as_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none") + xlab("Quintiles of PM2.5") + ylab("Arsenic")


p_cd <- df %>%
  ggplot(aes(x=pm25.f, y=hs_cd_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none")+ xlab("Quintiles of PM2.5") + ylab("Cadmium")


p_co <- df %>%
  ggplot(aes(x=pm25.f, y=hs_co_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none") + xlab("Quintiles of PM2.5") + ylab("Cobalt")


p_cs <- df %>%
  ggplot(aes(x=pm25.f, y=hs_cs_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none") + xlab("Quintiles of PM2.5") + ylab("Cesium")


p_cu <- df %>%
  ggplot(aes(x=pm25.f, y=hs_cu_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none") + xlab("Quintiles of PM2.5") + ylab("Copper")


p_hg <- df %>%
  ggplot(aes(x=pm25.f, y=hs_hg_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none") + xlab("Quintiles of PM2.5") + ylab("Mercury")


p_mn <- df %>%
  ggplot(aes(x=pm25.f, y=hs_mn_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none") + xlab("Quintiles of PM2.5") + ylab("Manganese")


p_mo <- df %>%
  ggplot(aes(x=pm25.f, y=hs_mo_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none") + xlab("Quintiles of PM2.5") + ylab("Molybdenum")


p_pb <- df %>%
  ggplot(aes(x=pm25.f, y=hs_pb_c_Log2, fill=pm25.f)) +
  geom_boxplot() + theme(legend.position="none") + xlab("Quintiles of PM2.5") + ylab("Lead")


grid.arrange(p_as, p_cd, p_co, p_cs, p_cu, p_hg, p_mn, p_mo, p_pb, ncol=3, nrow=3)

## plot air pollution levels by cohort

p_ap_site <- df %>%
  ggplot(aes(x=h_cohort, y=hs_pm25_yr_hs_h_None, fill=h_cohort)) +
  geom_boxplot() + theme(legend.position="none")







