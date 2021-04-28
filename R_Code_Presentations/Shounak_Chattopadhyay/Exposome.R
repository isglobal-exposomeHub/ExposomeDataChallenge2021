setwd("C:/Users/shoun/Google Drive/Research_Duke/Shape Constrained Inference/Exposome Data Challenge/Data")
load("C:/Users/shoun/Google Drive/Research_Duke/Shape Constrained Inference/Exposome Data Challenge/Data/Exposome Data/exposome.RData")

# res_var = covariates$e3_gac_None*7 + rnorm(n = length(covariates$e3_gac_None),
#                                          mean = 0, sd = 0.5)

#exposure_matrix = NULL

n = 1301

res_var = phenotype$e3_bw

# exposure_matrix = cbind(exposure_matrix, 
#                         exposome$hs_pfoa_m_Log2,  exposome$hs_as_m_Log2,
#                         exposome$hs_cu_m_Log2, 
#                         exposome$hs_pb_m_Log2,
#                         exposome$hs_cd_m_Log2, exposome$hs_co_m_Log2, exposome$hs_hg_m_Log2,
#                         exposome$hs_pfhxs_m_Log2, exposome$hs_pfna_m_Log2,
#                         exposome$hs_pfos_m_Log2, exposome$hs_pfunda_m_Log2)


exposure_matrix = cbind((exposome$hs_sumPCBs5_madj_Log2 + rnorm(n, 0, 0.1)),
                        (exposome$hs_dde_madj_Log2 + rnorm(n, 0, 0.1)),
                        exposome$hs_pfoa_m_Log2)

covariate_matrix = cbind(covariates$h_mbmi_None, covariates$hs_wgtgain_None,
                         covariates$e3_gac_None, covariates$h_age_None,
                         (as.numeric(exposome$e3_alcpreg_yn_None) - rep(1,1301)),
                         ifelse(exposome$e3_asmokcigd_p_None > 0, 1, 0))

whole_data = as.data.frame(cbind(res_var, covariate_matrix, exposure_matrix))

whole_data = whole_data[which(covariates$e3_gac_None >= 37),]

colnames(whole_data) = c("BirthWt","PrePregnancyWt",
                         "WtGain","GestPeriod","AgeofMother","Alcohol","Smoking",
                         "SumofPCBs","DDE","PFOA")

n = dim(whole_data)[1]
p = dim(whole_data[,8:10])[2]
p_cov = dim(whole_data[,2:7])[2]

X_stand = matrix(0, nrow = n, ncol = p)
#X = matrix(0, nrow = n, ncol = p)
X = whole_data[,8:10]
y = rep(0, n)
Z = matrix(0, nrow = n, ncol = p_cov)

for(i in 1:n)
{
  
  y[i] = (whole_data$BirthWt[i] - mean(whole_data$BirthWt))/
              sd(whole_data$BirthWt)
  
  # for(j in 1:p)
  # {
  #   
  #   X_stand[i,j] = (exposure_matrix[i,j] - mean(exposure_matrix[,j]))/sd(exposure_matrix[,j])
  #   
  # }
  
}

for(i in 1:n)
{
  
  for(j in 1:p)
  {
    
    # X[i,j] = (exposure_matrix[i,j] - min(exposure_matrix[,j]))/
    #           (max(exposure_matrix[,j]) - min(exposure_matrix[,j]))
    
    # X_stand[i,j] = exp(exposure_matrix[i,j])/(1 + exp(exposure_matrix[i,j]))
    
    X_stand[i,j] = X[i,j]
    
  }
  
}

for(i in 1:n)
{
  for(j in 1:p)
  {
    
    X[i,j] = (X_stand[i,j] - min(X_stand[,j]))/(max(X_stand[,j]) - min(X_stand[,j]))
    
  }
}

covariate_matrix = whole_data[,2:7]

for(i in 1:n)
{
  
  for(j in 1:6)
  {
    
    if((j == 5) || (j == 6))
    {
      Z[i,j] = covariate_matrix[i,j]
    }else{
      Z[i,j] = (covariate_matrix[i,j] - mean(covariate_matrix[,j]))/sd(covariate_matrix[,j])
      
    }
    
  }
  
  
}

M = 6
N = 4
eps_HMC = rep(0.2, choose(p,2))
#eps_HMC = 0.1
MC = 30000
coef_vec = rep(-1, choose(p,2))
nu_MALA = 0.1
zero_ind = rep(1, choose(p,2))























