rm(list=ls())

# install the meshed package -- requires compilation
install.packages("R_Code_Presentations/Michele_Peruzzi/meshed_0.21.04.29.tar.gz", repos=NULL, type="source")

library(meshed)
library(tidyverse)
library(magrittr)
library(scico)

load(url("https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/exposome_NA.RData?raw=true"))

set.seed(20210429)

# select some exposures
exposomeNA_num <- exposomeNA[,exposomeNA %>% as.list() %>% lapply(is.numeric) %>% unlist()]
exposomeNA_num %<>% 
  dplyr::select(ID, contains("pcb"),
                hs_KIDMED_None, 
                hs_mvpa_prd_alt_None,
                hs_dif_hours_total_None,
                hs_ndvi100_h_None, 
                hs_ndvi100_s_None)

colnames(exposomeNA_num)[-1] <- paste0("X_", colnames(exposomeNA_num)[-1])


# covariates with missing data
Ytildes <- covariatesNA %>% 
  dplyr::select(ID, 
                hs_wgtgain_None, 
                h_mbmi_None,  
                hs_c_height_None, 
                hs_c_weight_None)
colnames(Ytildes)[-1] <- paste0("Yt_", colnames(Ytildes)[-1])

# covariates & confounders
covariatesNA_sub <- covariatesNA %>%
  dplyr::filter(complete.cases(h_edumc_None, h_native_None, h_age_None, h_parity_None)) %>%
  dplyr::select(-h_mbmi_None, 
                -hs_wgtgain_None, 
                -hs_c_height_None, 
                -hs_c_weight_None) %>%
  mutate(e3_sex_None = 1*(e3_sex_None == "female"))
covariatesNA_sub <- model.matrix(~. -1, data=covariatesNA_sub) %>% as.data.frame()
colnames(covariatesNA_sub)[-1] <- paste0("Z_", colnames(covariatesNA_sub)[-1])

# outcomes
phenotypeNA_sub <- phenotypeNA %>% 
  dplyr::select(ID, 
                e3_bw,
                hs_zbmi_who,
                hs_bmi_c_cat) %>% 
  mutate(overweight = 1*(hs_bmi_c_cat %in% c(3,4))) %>% 
  dplyr::select(-hs_bmi_c_cat)
colnames(phenotypeNA_sub)[-1] <- paste0("Y_", colnames(phenotypeNA_sub)[-1])

dataf <- covariatesNA_sub %>% 
  left_join(phenotypeNA_sub) %>% 
  left_join(Ytildes) %>%
  left_join(exposomeNA_num)


Z <- dataf %>% 
  dplyr::select(contains("Z_", ignore.case=F)) %>%
  as.matrix() 
# gelman : http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf
Z[,c("Z_e3_gac_None", "Z_h_age_None", "Z_hs_child_age_None")] %<>%
  apply(2, function(x) (x-mean(x))/(2*sd(x))) 

X <- dataf %>% 
  dplyr::select(contains("X_", ignore.case=F)) %>% 
  as.matrix() %>% apply(2, function(x) scale(x) %>% as.numeric())

Y <- dataf %>% 
  dplyr::select(contains("Yt_", ignore.case=F),
                contains("Y_", ignore.case=F)) %>% 
  as.matrix()

# types
data_types <- c("poisson", "gaussian", "gaussian", 
                "gaussian", "gaussian", 
                "gaussian", "binomial")

# center gaussian outcomes
Y[,data_types == "gaussian"] <- Y[,data_types=="gaussian"] %>% apply(2, function(x) (x-mean(x, na.rm=T))/sd(x,na.rm=T))


# # # # # # # # # # # # # # # # # # # # # #
#           run piMeshedGP                #
# # # # # # # # # # # # # # # # # # # # # #

meshed_out <- pimeshed(Y, X, Z, k=4, proj_dim=2,
                       block_size = 30,
                       n_samples = 1000,
                       n_burnin = 5000, #1000,
                       n_thin = 1,
                       n_threads = 10,
                       print_every = 20,
                       family = data_types,
                       debug = list(sample_beta=T, sample_tausq=T, 
                                    sample_theta=T, sample_w=T, sample_lambda=T,
                                    verbose=F, debug=F))

tcp_code <- '
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace std;

//[[Rcpp::export]]
arma::cube tcp(const arma::cube& x){
  arma::cube result = arma::zeros(x.n_rows, x.n_rows, x.n_slices);
  for(int i=0; i<x.n_slices; i++){
    result.slice(i) = x.slice(i) * arma::trans(x.slice(i));
  }
  return result;
}
'
Rcpp::sourceCpp(code=tcp_code)

LLt <- meshed_out$lambda_mcmc %>% tcp() %>% apply(1:2, mean)

LLt_norm <- diag(1/sqrt(diag(LLt))) %*% LLt %*% diag(1/sqrt(diag(LLt)))
LLt_norm %<>% as.data.frame()
colnames(LLt_norm) <- rownames(LLt_norm) <- 
  colnames(Y) %>% 
  sapply(function(x) gsub("Y_", "", x)) %>%
  sapply(function(x) gsub("Yt_", "", x))

LLt_norm %<>% tibble::rownames_to_column("Outcome_row") %>%
  gather(Outcome_col, Value, -Outcome_row)

# correlations across outcomes
ggplot(LLt_norm, aes(Outcome_row, Outcome_col, fill=Value)) +
  geom_tile() +
  scale_fill_scico() +
  labs(x=NULL, y=NULL, fill="Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1))


# confounders on outcomes
Bmean <- meshed_out$beta_mcmc %>% 
  apply(1:2, mean) %>% 
  apply(2, function(x) x/max(abs(x)))
Bmean %<>% as.data.frame()
colnames(Bmean) <- 
  colnames(Y) %>% 
  sapply(function(x) gsub("Y_", "", x)) %>%
  sapply(function(x) gsub("Yt_", "", x))
rownames(Bmean) <- 
  colnames(Z) %>% 
  sapply(function(x) gsub("Z_", "", x))
Bmean %<>% tibble::rownames_to_column("Covariate") %>%
  gather(Outcome, Value, -Covariate)

# correlations across outcomes
ggplot(Bmean, aes(Outcome, Covariate, fill=Value)) +
  geom_tile() +
  scale_fill_scico(palette="broc") +
  labs(x=NULL, y=NULL, fill="Reg. coeff.\n(norm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1))


# exposures params 
# first row is used for parameter expansion and is included in lambda already
thetamean <- meshed_out$theta_mcmc[-1,,] %>% 
  apply(1:2, mean) %>% 
  apply(2, function(x) x/max(abs(x)))
thetamean %<>% as.data.frame()
colnames(thetamean) <- paste0("comp_", 1:4)
rownames(thetamean) <- 
  colnames(X) %>% 
  sapply(function(x) gsub("X_", "", x))
thetamean %<>% tibble::rownames_to_column("Exposure") %>%
  gather(Outcome, Value, -Exposure)

# correlations across outcomes
ggplot(thetamean, aes(Outcome, Exposure, fill=Value)) +
  geom_tile() +
  scale_fill_scico(palette="roma") +
  labs(x=NULL, y=NULL, fill="Exposure\nparam") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1))


# bkmr
library(bkmr)

# bkmr uses Z as the object for exposures

bkmr_time <- 
  system.time({
  fitkm <- kmbayes(y = Y[,6], Z = X, X = Z, iter = 2000, verbose = T, varsel = T) })


unimeshed_time <- 
  system.time({
  unimeshed_out <- pimeshed(Y[,6,drop=F], X, Z, k=1, proj_dim=2,
                       block_size = 30,
                       n_samples = 1000,
                       n_burnin = 1000,
                       n_thin = 1,
                       n_threads = 10,
                       print_every = 20,
                       family = "gaussian", #data_types,
                       debug = list(sample_beta=T, sample_tausq=T, 
                                    sample_theta=T, sample_w=T, sample_lambda=T,
                                    verbose=F, debug=F)) })



