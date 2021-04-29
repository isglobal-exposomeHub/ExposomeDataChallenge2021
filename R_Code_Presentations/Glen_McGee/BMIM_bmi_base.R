options(echo=TRUE)
options(stringsAsFactors = FALSE)

library(tidyverse)
# devtools::install_github("glenmcgee/bsmim2")
library(bsmim2)


setwd("../")
load("exposome_v2.RData")
load("DataPrep/covariates.rda")
load("DataPrep/index_exposure_list.rda")

#################################
## 2-way interaction plots

## function to predict two-way interaction curves
pred_twoway <- function(obj,qtls=c(0.1,0.9),qtl_lims=c(0.01,0.99),pts=20,include_median=TRUE){
  
  ## get predictions at each level (skip 0.5 since it is implicitly computed)
  res <- list()
  for(mm in 1:ncol(obj$rho)){
    res_mm <- list()
    for(qq in 1:length(qtls)){
      res_mm[[qq]] <- predict_hnew_indexwise2(obj,crossM=mm,qtl=qtls[[qq]],qtl_lims=qtl_lims,points=pts)
    }
    names(res_mm) <- qtls  ## label
    res[[mm]] <- res_mm
  }
  
  ## combine predictions into dataframe for plotting
  df_var1 <- df_var2 <- df_grid <- df_quantile <- df_est <- c()
  for(xx in 1:ncol(obj$rho)){
    for(yy in 1:ncol(obj$rho)){
      if(xx==yy){
        next
      }
      for(qq in 1:length(qtls)){
        df_var1 <- c(df_var1,rep(paste0("Index ",xx),pts))
        df_var2 <- c(df_var2,rep(paste0("Index ",yy),pts))
        df_grid <- c(df_grid,res[[yy]][[qq]]$grid[(xx-1)*pts+1:pts])
        df_est <- c(df_est,res[[yy]][[qq]]$mean[(xx-1)*pts+1:pts])
        df_quantile <- c(df_quantile,rep(qtls[qq],pts))
      }
      if(include_median==TRUE){ ## implicitly predicted above
        df_var1 <- c(df_var1,rep(paste0("Index ",xx),pts))
        df_var2 <- c(df_var2,rep(paste0("Index ",yy),pts))
        df_grid <- c(df_grid,res[[xx]][[qq]]$grid[(xx-1)*pts+1:pts]) ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_est <- c(df_est,res[[xx]][[qq]]$mean[(xx-1)*pts+1:pts])   ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_quantile <- c(df_quantile,rep(0.5,pts))
      }
      
    }
  }
  
  pred_df <- data.frame(var1=df_var1,var2=df_var2,grid=df_grid,quantile=df_quantile,est=df_est)
  
  return(pred_df)
}



z <- cbind(z_base_preg,z_base_post)
# remove covariates that are already adjust for in the outcome or related to the outcome
z <- z[,-which(colnames(z)%in%c("hs_child_age_None","e3_sex_Nonemale","hs_c_weight_None","hs_c_height_None"))]

nit <- 20000
thin = 20

  fit <- bsmim2(y=scale(phenotype$hs_zbmi_who), 
                x=index_exposure_list, 
                z =z , 
                spike_slab=FALSE, ## variable selection
                gauss_prior=TRUE, ## gaussian slab on thetastar
                prior_theta_slab_sd=0.2, ## sd of gaussian slab
                stepsize_theta=0.2, ## stepsize for random walk 
                num_theta_steps=10,
                draw_h=FALSE,
                niter=nit, 
                nthin=thin,
                ) 
  save(fit, file=paste0("BMIM/Output/bmim_bmi_base.Rdata"))
  index_names <- colnames(index_exposure_list)
  pred_assoc <- predict_hnew_assoc2(fit)
  pred_overall <- predict_hnew_assoc2(fit,overall = TRUE)
  pred_ind <- predict_hnew_indexwise2(fit)
  sum_theta <- summarize_thetas(fit)
  pred_inter <- pred_twoway(fit)
  save(fit, index_names, pred_assoc, pred_overall, pred_ind, sum_theta, pred_inter,
       file=paste0("BMIM/Output/bmim_bmi_base.Rdata"))
