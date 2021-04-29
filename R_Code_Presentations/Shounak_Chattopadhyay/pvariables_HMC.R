library(bayesplot)
library(MASS)
library(splines)
library(gam)
library(splines2)
library(tmg)
library(truncnorm)
library(mvtnorm)

#### Note: April 19, 2021 - code has no intercept to ensure identifiability. #### 

##M+4 = #splines for main effects, N+3 = #splines for interaction tensor products

penmatt<-function(M)
{
  
  matt = matrix(0, nrow = M, ncol = M)
  
  diag(matt)[1:(M-1)] = 2
  diag(matt)[M] = 1
  
  for(i in 1:M)
  {
    for(j in 1:M)
    {
      if(i == j+1)
      {
        matt[i,j] = -1
      }
      else if(j == i+1)
      {
        matt[i,j] = -1
      }
    }
  }
  
  return(matt)
  
}

SIDsampler<-function(y, X, Z, M, N, MC, coef_vec, eps_HMC){
  
  n = dim(X)[1]
  p = dim(X)[2]
  p_cov = dim(Z)[2]
  
  #### Main Effects Splines ####
  
  ME_nspl = M+4
  
  ME_list = array(0, dim = c(n, ME_nspl, p))
  ind = 1
  
  for(ind in 1:p)
  {
    
    ### Integral constraint ###
    
    me_knots = seq(1/M, 1-(1/M), length.out = M)
    me_spl = bSpline(X[,ind], knots = me_knots, intercept = TRUE)

    Xmat = sweep(me_spl, 2, colMeans(me_spl))
    
    ME_list[,,ind] = Xmat
    
  }
  
  ME_mat = NULL
  for(ind in 1:p)
  {
    
    ME_mat = cbind(ME_mat, ME_list[,,ind])
    
  }
  
  #### Interaction Effects Splines ####
  
  IE_nspl = N+3
  
  IE_list = array(0, dim = c(n, (IE_nspl)*(IE_nspl), choose(p, 2)))
  precond_list = array(0, dim = c((IE_nspl)^2, (IE_nspl)^2, choose(p,2)))
  k = 1
  
  for(k in 1:choose(p,2))
  {
    
    #Obtain inverse index (u,v)
    
    quad_dis = (2*p - 1)^2 - 8*k
    u = ceiling(0.5*((2*p-1) - quad_dis^0.5))
    v = p + k - (u*(p - 0.5*(u+1)))
    
    # Now store matrix
    
    int_knots_1 = seq(1/(N+1), 1-(1/(N+1)), length.out = N)
    int_knots_2 = seq(1/(N+1), 1-(1/(N+1)), length.out = N)

    int_spl1 = bSpline(X[,u], knots = int_knots_1, intercept = FALSE)
    int_spl2 = bSpline(X[,v], knots = int_knots_2, intercept = FALSE)
    
    # int_spl1 = bSpline(X[,u], knots = int_knots_1, intercept = TRUE)
    # int_spl2 = bSpline(X[,v], knots = int_knots_2, intercept = TRUE)
    # 
    # int_spl1[,1] = int_spl1[,1] - rep(1, n)
    # int_spl2[,1] = int_spl2[,1] - rep(1, n)
    
    for(i in 1:n)
    {
      IE_list[i,,k] = kronecker(int_spl1[i,], int_spl2[i,])
    }
    
    #precond_list[,,k] = (0.1^2)*diag(IE_nspl^2)
    #precond_list[,,k] = (t(IE_list[,,k]) %*% IE_list[,,k])/n
    precond_list[,,k] = (nu_MALA^2)*diag(IE_nspl^2)
    
  }
  
  ##### Define the penalty matrices #####
  
  SigmaME_inv = penmatt(ME_nspl)
  SigmaME = solve(SigmaME_inv)

  SigmaInt_inv = kronecker(penmatt(IE_nspl), penmatt(IE_nspl))
  SigmaInt = solve(SigmaInt_inv)
  
  #### HMC Potential Energy and Gradient Functions ####
  
  pot_fn_HMC <- function(Gam, sigsq, lamb, tausq, Siginv, M, R, coef)
  {
    
    mainpart = M %*% Gam
    
    part1 = sum((R - (coef*mainpart^2))^2)/(2*sigsq)
    #part1 = (sum((mainpart^2)^2) - 2*sum((mainpart^2)*R))/(2*sigsq)
    
    part2_qf = as.numeric(t(Gam) %*% Siginv %*% Gam)
    part2 = (lamb*part2_qf)/(2*sigsq*tausq)
    
    pot_output = part1 + part2
    return(pot_output)
    
  }
  
  grad_fn_HMC <- function(Gam, sigsq, lamb, tausq, Siginv, M, R, coef)
  {
    
    mainpart = M %*% Gam
    
    part1 = lamb*(Siginv %*% Gam)/(sigsq*tausq)
    
    part2_coeff = 2/sigsq
    part2_vec = (mainpart^2 - coef*R) * mainpart
    
    grad_output = part1 + (part2_coeff * (t(M) %*% part2_vec))
    
    return(as.numeric(grad_output))
    
  }
  
  # pot_fn_precond <- function(Gam, sigsq, lamb, tausq, Siginv, M, R, 
  #                                precond)
  # {
  #   
  #   return(pot_fn_HMC(precond %*% Gam, sigsq, lamb, tausq, Siginv, M, R))
  #   
  # }
  # 
  # grad_fn_precond <- function(Gam, sigsq, lamb, tausq, Siginv, M, R,
  #                             precond)
  # {
  #   
  #   return(precond %*% grad_fn_HMC(precond %*% Gam, sigsq, lamb, tausq, 
  #                                  Siginv, M, R))
  #   
  # }
  
  #### Define parameter storages ####
  
  lambda_stor = matrix(1, nrow = MC, ncol = p)
  delta_stor = matrix(1, nrow = MC, ncol = choose(p,2))
  nu_stor = matrix(1, nrow = MC, ncol = choose(p,2))
  tausq_stor = rep(1, MC)
  sigmasq_stor = rep(0.1, MC)
  alpha_stor = rep(0, MC)
  a_stor = matrix(1, nrow = MC, ncol = p)
  cov_effect = matrix(0, nrow = MC, ncol = p_cov)
  
  ls_mat = solve(t(Z) %*% Z)
  
  ME_coeff_stor = matrix(0, nrow = MC, ncol = p*(ME_nspl))
  IE_coeff_stor = array(0, dim = c(MC, IE_nspl^2, choose(p,2)))
  
  rho_stor = matrix(0, nrow = IE_nspl^2, ncol = choose(p,2))
  
  R0 = y - (Z %*% ((ls_mat %*% t(Z)) %*% y))
  
  for(k in 1:choose(p,2))
  {
    
    # IE_coeff_stor[1,,k] = mvrnorm(n = 1, mu = rep(0, IE_nspl^2),
    #         Sigma = sigmasq_stor[1]*SigmaInt*tausq_stor[1]/delta_stor[1,k])
    
    #IE_coeff_stor[1,,k] = rep(1, IE_nspl^2)
    
    IE_coeff_stor[1,,k] = solve(t(IE_list[,,k]) %*% IE_list[,,k] + 0.05*diag(IE_nspl^2)) %*% 
                          (t(IE_list[,,k]) %*% (abs(R0))^(1/2))/sqrt(choose(p,2))
    
  }
  
  accept_HMC = matrix(0, nrow = MC, ncol = choose(p,2))
  accept_HMC[1,] = rep(1, choose(p,2))
  
  int_select = matrix(0, nrow = MC, ncol = choose(p,2))
  
  #### Start Sampling ####
  
  int_effects = matrix(0, nrow = n, ncol = choose(p,2))
  for(k in 1:choose(p,2))
  {
    
    int_effects[,k] = coef_vec[k]*(IE_list[,,k] %*% IE_coeff_stor[1,,k])^2
    
  }
  
  m = 2
  
  t1 = proc.time()
  
  for(m in 2:MC)
  {
    
    #### y, X, M, N, eps_HMC, MC
    
    sum_inteff = rowSums(int_effects)
    
    #### Covariate effects ####
    
    R_cov = y - (ME_mat %*% ME_coeff_stor[m-1,]) - sum_inteff
    
    ls_est = ls_mat %*% (t(Z) %*% R_cov)
    
    cov_effect[m,] = mvrnorm(n = 1, mu = ls_est, 
                             Sigma = sigmasq_stor[m-1]*ls_mat)
    
    #### Sample alpha (don't need to do this anymore in new model) ####
    
    # R_alpha = y - (ME_mat %*% ME_coeff_stor[m-1,]) - (sum_inteff)
    # 
    # # alpha_stor[m] = rnorm(n = 1, mean = mean(R_alpha),
    # #                       sd = sqrt(sigmasq_stor[m-1]/n))
    # 
    # alpha_stor[m] = 0
    
    #### Sample beta (main effects coeff) ####
    
    R_beta = y - (alpha_stor[m]*rep(1,n)) - (sum_inteff) - (Z %*% cov_effect[m,])
    
    Delta_mat_inv = diag(1/lambda_stor[m-1,])
    Psi_mat = kronecker(Delta_mat_inv, SigmaME)
    
    proxy_mat = solve((t(ME_mat) %*% ME_mat) + solve(Psi_mat))
    cov_beta = sigmasq_stor[m-1]*proxy_mat
    mean_beta = proxy_mat %*% (t(ME_mat) %*% R_beta)
    
    ME_coeff_stor[m,] = mvrnorm(n = 1, mu = mean_beta, 
                                Sigma = cov_beta)
    
    Svec = rep(0, p)
    for(ME_ind in 1:p)
    {
      
      which_main = (1 + (ME_ind - 1)*(ME_nspl)):((ME_nspl)*(ME_ind))
      
      Svec[ME_ind] = as.numeric(t(ME_coeff_stor[m,which_main]) %*% 
                                  SigmaME_inv %*%
                                  (ME_coeff_stor[m,which_main]))
      
      #### Sample lambda_j ####
      
      # shape_par = ((ME_nspl) + 1)/2
      # rate_par = 0.5 + (0.5*Svec[ME_ind]/sigmasq_stor[m-1])
      # 
      # lambda_stor[m,ME_ind] = rgamma(n = 1, shape = shape_par,
      #                                rate = rate_par)
      
      shape_par = ((ME_nspl) + 1)/2
      rate_par = a_stor[m-1,ME_ind] + (0.5*Svec[ME_ind]/sigmasq_stor[m-1])

      lambda_stor[m, ME_ind] = rgamma(n = 1, shape = shape_par, rate = rate_par)
      
      #### Sample a_j ####
      
      a_stor[m, ME_ind] = rgamma(n = 1, shape = 1, rate = lambda_stor[m, ME_ind] + 1)
      a_stor[m, ME_ind] = 0.5
      
    }
    
    #### Sample gammas (interaction coefficients) and deltas ####
    
    Qvec = rep(0, choose(p,2))
    
    k=1
    for(k in 1:choose(p,2))
    {
      
      if(coef_vec[k] != 0)
      
      {
        #### Sample gamma_k = gamma_{uv} ####
        
        R_gamma = y - (alpha_stor[m]*rep(1,n)) - (ME_mat %*% ME_coeff_stor[m,]) -
          (rowSums(int_effects[,-k])) - (Z %*% cov_effect[m,])
        
        ##Use MALA with preconditioning. Use eps_HMC = d^(-1/2) - Yuansi et al paper.
        
        tau_MALA = (eps_HMC[k]^2)/2
        
        Mmat = precond_list[,,k]
        
        current_gamma = IE_coeff_stor[m-1,,k]
        new_rho = mvrnorm(n = 1, mu = rep(0, IE_nspl^2), Sigma = 2*tau_MALA*Mmat)
        
        current_rho = new_rho
        
        ## Orthogonal rho idea (highly experimental) ##
        
        # current_gamma = IE_coeff_stor[m-1,,k]
        # 
        # if(accept_HMC[m-1,k] == 1)
        # {
        #   new_rho = mvrnorm(n = 1, mu = rep(0, IE_nspl^2),
        #                     Sigma = 2*tau_MALA*Mmat)
        # }else{
        #   rho_tilde = mvrnorm(n = 1, mu = rep(0, IE_nspl^2),
        #                       Sigma = 2*tau_MALA*Mmat)
        # 
        #   new_rho = rho_tilde - ((sum(rho_tilde*rho_stor[,k])/
        #                             sum(rho_stor[,k]^2))*rho_stor[,k])
        # }
        # 
        # current_rho = new_rho
        
        new_gamma = current_gamma - (tau_MALA*(Mmat %*% grad_fn_HMC(current_gamma, 
                                                                    sigmasq_stor[m-1], delta_stor[m-1,k], 
                                                                    tausq_stor[m-1], SigmaInt_inv, IE_list[,,k],
                                                                    R_gamma, coef_vec[k]))) + new_rho
        
        current_U = pot_fn_HMC(current_gamma, sigmasq_stor[m-1], 
                               delta_stor[m-1,k], tausq_stor[m-1], SigmaInt_inv, 
                               IE_list[,,k], R_gamma, coef_vec[k])
        
        current_K = dmvnorm(x = current_gamma, mean = new_gamma -
                              tau_MALA*(Mmat %*% grad_fn_HMC(new_gamma, sigmasq_stor[m-1],
                                                             delta_stor[m-1,k],tausq_stor[m-1],
                                                             SigmaInt_inv, IE_list[,,k], R_gamma, coef_vec[k])),
                            sigma = 2*tau_MALA*Mmat, log = TRUE)
        
        # current_K = sum((new_gamma - current_gamma + 
        #             (tau_MALA*grad_fn_precond(current_gamma, 
        #               sigmasq_stor[m-1], delta_stor[m-1,k], 
        #               tausq_stor[m-1], SigmaInt_inv, IE_list[,,k],
        #               R_gamma, Mmat)))^2)/(4*tau_MALA)
        
        new_U = pot_fn_HMC(new_gamma, sigmasq_stor[m-1], delta_stor[m-1,k],
                           tausq_stor[m-1], SigmaInt_inv, IE_list[,,k],
                           R_gamma, coef_vec[k])
        
        v1 = as.numeric(new_gamma)
        v2 = current_gamma - tau_MALA*(Mmat %*% grad_fn_HMC(current_gamma, 
                                                            sigmasq_stor[m-1], delta_stor[m-1,k],
                                                            tausq_stor[m-1], SigmaInt_inv, 
                                                            IE_list[,,k], R_gamma, coef_vec[k]))
        
        new_K = dmvnorm(x = v1, mean = v2, sigma = 2*tau_MALA*Mmat, log = TRUE)
        
        U1 = runif(1)
        
        energy_diff = exp(current_U - new_U + current_K - new_K)
        
        if(is.finite(energy_diff) == T)
        {
          if(U1 <= energy_diff)
          {
            IE_coeff_stor[m,,k] = new_gamma
            accept_HMC[m,k] = 1
          }else
          {
            IE_coeff_stor[m,,k] = IE_coeff_stor[m-1,,k]
            accept_HMC[m,k] = 0
          }
          #print(paste("No Issues with MALA #", k, sep=""))
        }else
        {
          IE_coeff_stor[m,,k] = IE_coeff_stor[m-1,,k]
          accept_HMC[m,k] = 0
          #print(paste("HMC Diverged for MALA #", k, sep=""))
        }
        
        #rho_stor[,k] = new_rho
        
        ## Update interaction effects matrix and var-selection coefficient##
        
        int_effects[,k] = coef_vec[k]*(IE_list[,,k] %*% IE_coeff_stor[m,,k])^2
        
        int_select[m,k] = sum(int_effects[,k])/n
        
        #print(int_effects[,k])
        
        Qvec[k] = as.numeric(t(IE_coeff_stor[m,,k]) %*% SigmaInt_inv %*%
                               IE_coeff_stor[m,,k])
        
        # #### Update delta_{uv} = delta_k ####
        # 
        # shape_del = 0.5*(1 + (IE_nspl)*(IE_nspl))
        # rate_del = 0.5 + (0.5*Qvec[k]/(sigmasq_stor[m-1]*tausq_stor[m-1]))
        # 
        # delta_stor[m,k] = rgamma(n = 1, shape = shape_del, 
        #                          rate = rate_del)
        
        #### Update delta_{uv} = delta_k ####
        
        shape_del = 0.5 + 0.5*((IE_nspl)*(IE_nspl))
        rate_del = nu_stor[m-1,k] + (0.5*Qvec[k]/(sigmasq_stor[m-1]*tausq_stor[m-1]))
        
        delta_stor[m,k] = rgamma(n = 1, shape = shape_del, 
                                 rate = rate_del)
        
        #### Update nu_{uv} = nu_k ####
        
        nu_stor[m,k] = rgamma(n = 1, shape = 1, rate = delta_stor[m,k] + 1)
      }
      
    }
    
    #### Sample tausq ####
    
    tau_sum = sum(delta_stor[m,]*Qvec)
    tau_shape = 0.5 + (0.5 * choose(p,2) * (IE_nspl) * (IE_nspl))
    tau_rate = 0.5 + (0.5*tau_sum/sigmasq_stor[m-1])

    tausq_stor[m] = 1/rgamma(n = 1, shape = tau_shape,
                             rate = tau_rate)
    
    #### Sample sigmasq ####
    
    R_sigsq = y - (alpha_stor[m]*rep(1,n)) - (ME_mat %*% ME_coeff_stor[m,]) -
                  (rowSums(int_effects)) - (Z %*% cov_effect[m,])
    s1_me = sum(lambda_stor[m,]*Svec)
    
    sigsq_shape = 0.5*(n + p*(ME_nspl) + (IE_nspl)*(IE_nspl)*choose(p,2) + 1)
    sigsq_rate = 0.5*(sum(R_sigsq^2) + s1_me + (tau_sum/tausq_stor[m]) + 1)
    
    sigmasq_stor[m] = 1/rgamma(n = 1, shape = sigsq_shape, 
                               rate = sigsq_rate)
    
    
      print("MCMC Iterate")
      print(m)
      if(m > 500)
      {
        print("Accept Ratio")
        print(colMeans(accept_HMC[(m-500):m,]))
      }
      print("Overall Accept Ratio")
      print(colMeans(accept_HMC[1:m,]))
      print("h(1,1)'s")
      print(IE_coeff_stor[m,IE_nspl^2,]^2)
      print("Error variance")
      print(sigmasq_stor[m])
    #}
    
  }
  
  t2 = proc.time()
  
  print("Time Elapsed")
  print((t2-t1)[3])
  
}

####### Plotting etc #######

#good_index = (MC - 500):MC
#good_index = seq(6005, 10000, by = 5)

good_index = 20000:30000

# me1 = ME_coeff_stor[good_index,1:14] %*% t(ME_list[,,1])
# plot(X[,1][order(X[,1])], colMeans(me1)[order(X[,1])], type = "l",
#      main = "ME1", xlab = "x1", ylab = "f_1(x_1)")
# lines(X[,1][order(X[,1])], apply(me1, 2, quantile, 0.025)[order(X[,1])],
#       col = "blue", lty = 3)
# lines(X[,1][order(X[,1])], apply(me1, 2, quantile, 0.975)[order(X[,1])],
#       col = "blue", lty = 3)
# lines(X[,1][order(X[,1])], true_main_mat[,1][order(X[,1])], col="red")
# 
# me2 = ME_coeff_stor[good_index,11:20] %*% t(ME_list[,,2])
# plot(X[,2][order(X[,2])], colMeans(me2)[order(X[,2])], type = "l",
#      main = "ME2", xlab = "x2", ylab = "f_2(x_2)")
# lines(X[,2][order(X[,2])], apply(me2, 2, quantile, 0.025)[order(X[,2])],
#       col = "blue", lty = 3)
# lines(X[,2][order(X[,2])], apply(me2, 2, quantile, 0.975)[order(X[,2])],
#       col = "blue", lty = 3)
# lines(X[,2][order(X[,2])], true_main_mat[,2][order(X[,2])], col="red")


#### Plot main effects ####

par(mfrow = c(1,1))

for(j in 1:p)
{
  
  main_ind = (1 + (j - 1)*(ME_nspl)):((ME_nspl)*j)
  
  memod = ME_coeff_stor[good_index,main_ind] %*% t(ME_list[,,j])
  plot(X[,j][order(X[,j])], colMeans(memod)[order(X[,j])], type = "l",
       main = paste("ME",j,sep=""), xlab = paste("x",j,sep=""), ylab = "ME", ylim = c(-1,1.5))
  lines(X[,j][order(X[,j])], apply(memod, 2, quantile, 0.025)[order(X[,j])],
        col = "blue", lty = 3)
  lines(X[,j][order(X[,j])], apply(memod, 2, quantile, 0.975)[order(X[,j])],
        col = "blue", lty = 3)
  #lines(X[,j][order(X[,j])], true_main_mat[,j][order(X[,j])], col="red")
  
  #legend("topright", fill = c("Black", "Red"), legend = c("Estimate", "True"))
  
}

#### Plot interaction effects ####

par(mfrow=c(2,3))

for(k in 1:choose(p,2))
{
  
  quad_dis = (2*p - 1)^2 - 8*k
  u = ceiling(0.5*((2*p-1) - quad_dis^0.5))
  v = p + k - (u*(p - 0.5*(u+1)))
  
  ## Given u, v, k = 
  
  x1_plot = seq(0.01, 0.99, length.out = 30)
  x2_plot = seq(0.01, 0.99, length.out = 30)
  
  true_int_plot = matrix(0, nrow = 30, ncol = 30)
  #new_true_int_plot = matrix(0, nrow = 30, ncol = 30)
  
  if(zero_ind[k] == 1)
  {
    
    for(i in 1:30)
    {
      for(j in 1:30)
      {
        
        true_int_plot[i,j] = -(x1_plot[i]*x2_plot[j])^2
        #new_true_int_plot[i,j] = 4*(x1_plot[i]*x2_plot[j])*(x1_plot[i]-x2_plot[j])^2
        
      }
    }
    
  }
  
  
  spline1 = bSpline(x = x1_plot, knots = seq(1/(N+1), 1-(1/(N+1)), length.out = N),
                    intercept = FALSE)
  spline2 = bSpline(x = x2_plot, knots = seq(1/(N+1), 1-(1/(N+1)), length.out = N),
                    intercept = FALSE)
  
  # beta_mean_vec = colMeans(model_pars$beta_vec_int)
  # beta_mean_mat = matrix(beta_mean_vec, nrow = dim(spline1)[2], ncol = dim(spline2)[2], byrow = T)
  
  est_int_plot = matrix(0, nrow = 30, ncol = 30)
  low_int_plot = matrix(0, nrow = 30, ncol = 30)
  high_int_plot = matrix(0, nrow = 30, ncol = 30)
  
  #int_const = diag(betamat[good_index,] %*% A %*% t(betamat[good_index,]))
  
  for(i in 1:30)
  {
    for(j in 1:30)
    {
      v_ij = kronecker(spline1[i,], spline2[j,])
      
      est_int_plot[i,j] = mean((IE_coeff_stor[good_index,,k] %*% v_ij)^2)
      #est_int_plot[i,j] = median((IE_coeff_stor[good_index,,k] %*% v_ij)^2)
      
      low_int_plot[i,j] = quantile((IE_coeff_stor[good_index,,k] %*% v_ij)^2, 0.025)
      high_int_plot[i,j] = quantile((IE_coeff_stor[good_index,,k] %*% v_ij)^2, 0.975)
    }
  }
  
  if(zero_ind[k] == 0)
  {
    
    persp(x = x1_plot, y = x2_plot, z = true_int_plot, theta = -50, 
          main= paste("True Surface", c(u,v), sep=""), xlim = c(0,1), ylim = c(0,1),
          zlim = c(0,1),
          ticktype = "detailed", xlab = "", ylab = "", zlab = "")
    persp(x = x1_plot, y = x2_plot, z = est_int_plot, theta = -50, 
          main = paste("Estimated Surface", c(u,v), sep=""), xlim = c(0,1), ylim = c(0,1),
          zlim = c(0,1),
          ticktype = "detailed", xlab = "", ylab = "", zlab = "")
    
  }else
  {
    
    # persp(x = x1_plot, y = x2_plot, z = true_int_plot, theta = -50, 
    #       main= paste("True Surface", c(u,v), sep=""), xlim = c(0,1), ylim = c(0,1),
    #       ticktype = "detailed", xlab = "", ylab = "", zlab = "")
    # persp(x = x1_plot, y = x2_plot, z = new_true_int_plot, theta = -50, 
    #       main= paste("New True Surface", c(u,v), sep=""), xlim = c(0,1), ylim = c(0,1),
    #       ticktype = "detailed", xlab = "", ylab = "", zlab = "")
    persp(x = x1_plot, y = x2_plot, z = est_int_plot, theta = -50, 
          main = paste("Estimated Surface", c(u,v), sep=""), xlim = c(0,1), ylim = c(0,1),
          ticktype = "detailed", xlab = "x1", ylab = "x2", zlab = "")
    #, zlim = c(0,1)
    persp(x = x1_plot, y = x2_plot, z = low_int_plot, theta = -50,
          main = paste("2.5% Surface", c(u,v), sep=""), xlim = c(0,1), ylim = c(0,1),
          ticktype = "detailed", xlab = "", ylab = "", zlab = "")
    persp(x = x1_plot, y = x2_plot, z = high_int_plot, theta = -50,
          main = paste("97.5% Surface", c(u,v), sep=""), xlim = c(0,1), ylim = c(0,1),
          ticktype = "detailed", xlab = "", ylab = "", zlab = "")
    print("Max Int")
    print(max(est_int_plot))
    print("Mean Int")
    print(mean(est_int_plot))
    
  }
  
}

par(mfrow = c(1,1))

#dev.off()

mean(alpha_stor[good_index])
quantile(alpha_stor[good_index], c(0.025,0.975))
mean(sigmasq_stor[good_index])  
quantile(sigmasq_stor[good_index], c(0.025,0.975))
colMeans(cov_effect[good_index,])
apply(cov_effect[good_index,], 2, quantile, 0.025)
apply(cov_effect[good_index,], 2, quantile, 0.975)

### Plot acceptance rate ###

accept_rate = rep(0, MC+1-1000)

for(i in 1000:MC)
{
  
  accept_rate[i-1000+1] = mean(colMeans(accept_HMC[(i-500):i,]))
  
}
  
plot(1000:MC, accept_rate, type = "l", main = "Accept Rate")
  
### Which interactions are zero? ###

colMeans(int_select[good_index,])
apply(int_select[good_index,],2,quantile,0.025)
apply(int_select[good_index,],2,quantile,0.975)





