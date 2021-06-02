
################################################################################################
#Non linear mediation analysis with high-dimensional mediators whose causal structure is unkown#
#Loh WW, Moerkerke B, Loeys T, Vansteelandt S#
#https://github.com/wwloh/interventional-hdmed#
################################################################################################

rm(list=ls())
library(data.table)
library(boot)

################################
####load and organize dataset###
################################

A<-read.table("/Users/eciliberto/Documents/progetti/athlete/interventional_hdmed_master/challenge/exposure.csv",sep = ",", header = TRUE)
A<-A[,-1]
M<-read.table("/Users/eciliberto/Documents/progetti/athlete/interventional_hdmed_master/challenge/mediators.csv",sep = ",", header = TRUE)
M<-M[,-1] 
#M<-M[,-c(1,6,15:26,36,40:41,82:83)] #only continuous mediators
Y<-read.table("/Users/eciliberto/Documents/progetti/athlete/interventional_hdmed_master/challenge/outcome.csv",sep = ",", header = TRUE)
Y<-Y[,-1]
Y[Y<2500]<-1 #binary outcome
Y[Y>=2500]<-0
L<-read.table("/Users/eciliberto/Documents/progetti/athlete/interventional_hdmed_master/challenge/covariates.csv",sep = ",", header = TRUE)
L<-L[,-1]

id<-c(1:length(A)) # id to preserve initial ordering of data
colnames(M) <- paste0("M",1:length(M))

###############################################################
#generate variables for functions
###############################################################

t <- ncol(M) # number of mediators
L <- as.matrix(L)
Lnames<-colnames(L)
A <- as.numeric(A)
M <- as.matrix(M)
Mnames<-colnames(M)
Y <- as.numeric(Y)
id<-as.numeric(id)
#fitY.family<-"binomial"
#fitM.family<-"gaussian"

#only continuous mediators and binary outcome
data<-data.frame(cbind(id,A,L,M,Y))
names(data)<-c("id","A",Lnames, Mnames, "Y")
completenames<-c("A",Lnames, Mnames)
fitY.form <- as.formula(paste0("Y~",paste(completenames, collapse="+"))) 
fitA.form <- as.formula(paste0("A~",paste(Lnames,collapse="+"))) 

###############################################################
# helper function to create duplicated data for one individual#
###############################################################

Dupdata <- function(t) {
  out <- diag(1,nrow=t,ncol=t) # matrix with diagonal of 1: for indirect effects via each mediator
  out <- rbind(0,out,1)  # first row all 0s, last row all 1s
  out <- cbind(0,out) # first column all zeroes for a0
  out <- rbind(0,out[nrow(out),],1,out) # for direct, joint, and mutual effects
  colnames(out) <- paste0("a",0:t)
  out <- cbind("Marg"=c(rep(0,3),rep(1,nrow(out)-3)),"a.i"=1:nrow(out),out)
  return(out)
}

#################################################################################
# function to estimate direct and indiretc effects without confounders selection#
#################################################################################

estimator_noselection <- function(data,mc_draws,fitY.family,pt_est=TRUE) {
  res <- list()

  # fit outcome model =========================================================
  if (fitY.family=="gaussian") {
    fitY <- glm(fitY.form, family = gaussian("identity"), data = data)  
  } else {
    fitY <- glm(fitY.form, family = binomial("logit"), data = data)  
  }
  
  
  # sampling probabilities for each individual's mediator values ==============
  fitA <- glm(fitA.form, family = binomial("logit"), data = data)
  pA1hat <- predict(fitA,type="response") # predicted prob. of treatment
  pA0hat <- 1-pA1hat # predicted prob. of control
  wt_A1 <- 1/pA1hat # inverse prob. of treatment weight
  wt_A0 <- 1/pA0hat # inverse prob. of control weight
  ind_A1 <- which(data$A==1) # indices for individuals in treatment group
  ind_A0 <- which(data$A==0) # indices for individuals in control group
  wt_A1.raw <- wt_A1[ind_A1]
  wt_A1 <- wt_A1.raw/sum(wt_A1.raw)
  wt_A0.raw <- wt_A0[ind_A0]
  wt_A0 <- wt_A0.raw/sum(wt_A0.raw)
  
  # mediator column names
    ## observed mediator values in each treatment group
  Mobs <- list(data[data$A==0, Mnames],data[data$A==1, Mnames])
  n_A0 <- sum(data$A==0)
  n_A1 <- sum(data$A==1)
  
     
  # helper function to predict Y for different outcome models
  PredictY <- function(onedat) {
    ## first column is always treatment
    colnames(onedat)[1] <- "A"
    Ya <- predict.glm(fitY, type="response", newdata=onedat)
    return(Ya)
  }
  
  # duplicated data for each individual =======================================
  dat <- data.table(data[,c("id",Lnames)]) #### keep only id, L to reduce memory
  setkey(dat)
  alevels <- dat[,as.data.table(Dupdata(t)),by=id] ###for each subject generate the expansion
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels)
  
   ## sample *joint* mediator values
    data.joint<-joint_means<-list()
    for (i in 1:max(id)){
  	mydt <- dat[id==i]
    mydt_mc <- mydt[rep(1:nrow(mydt),each=mc_draws)] #ripete numero righe di id=1 mc_draws volte
    mydt_mc <- cbind(mydt_mc,"mc"=rep(1:mc_draws,times=nrow(mydt)))
    setkey(mydt_mc)
    
    Mtilde <- data.frame(matrix(NA,nrow=sum(mydt_mc$Marg==0),ncol=t))
    colnames(Mtilde) <- Mnames
    
    ## sample *joint* mediator values, predict potential outcomes using sampled mediator value and average on MC drawss########
    #### a1==0
    a1_0 <- mydt_mc$Marg==0 & mydt_mc$a1==0
    Mtilde[a1_0,] <- Mobs[[1]][sample(n_A0, size=sum(a1_0), replace=TRUE,prob=wt_A0),] #sample from subjects with A=0 
    rm(a1_0)
    #### a1==1
    a1_1 <- mydt_mc$Marg==0 & mydt_mc$a1==1
    Mtilde[a1_1,] <- Mobs[[2]][sample(n_A1, size=sum(a1_1), replace=TRUE,prob=wt_A1),]
    rm(a1_1)
    data.joint[[i]]<-cbind(mydt_mc[Marg==0,],Mtilde)
    names(data.joint[[i]])<-c(colnames(dat),"mc",Mnames)
    rm(Mtilde)
    rm(mydt_mc)
    data.joint[[i]][, "Y.a" := PredictY(onedat=data.joint[[i]][, c("a0",Lnames,Mnames), with=FALSE])]
    setkey(data.joint[[i]])
    ## average over MC draws
    joint_means[[i]] <- data.joint[[i]][,lapply(.SD, mean), by=c("Marg","a.i",paste0("a",0:t)),.SDcols="Y.a"]
    setkey(joint_means[[i]])
    }
    joint_means<-do.call(rbind,joint_means)
    ## average over all individuals
    mu_hat_joint <- joint_means[, lapply(.SD,mean), by=c("Marg","a.i",paste0("a",0:t)),.SDcols="Y.a"]
    setnames(mu_hat_joint,old="Y.a",new="Y")
    setkey(mu_hat_joint)
       
    ## sample *marginal* mediator values, predict potential outcomes using sampled mediator value and average on MC drawss########
    data.marginal<-marginal_means<-list()
    for (i in 1:max(id)){
  	mydt <- dat[id==i]
    mydt_mc <- mydt[rep(1:nrow(mydt),each=mc_draws)] #ripete numero righe di id=1 mc_draws volte
    mydt_mc <- cbind(mydt_mc,"mc"=rep(1:mc_draws,times=nrow(mydt)))
    setkey(mydt_mc)
    
    Mtilde <- data.frame(matrix(NA,nrow=sum(mydt_mc$Marg==1),ncol=t))
    colnames(Mtilde) <- Mnames
    
    M.as <- lapply(1:t, function(s) {
      a_s <- unlist(mydt_mc[Marg==1,paste0("a",s),with=FALSE]) 
      Ms.as <- rep(NA,length(a_s))
      Ms.as[a_s==0] <- Mobs[[1]][sample(n_A0, size=sum(a_s==0), replace=TRUE,prob=wt_A0),s]
      Ms.as[a_s==1] <- Mobs[[2]][sample(n_A1, size=sum(a_s==1), replace=TRUE,prob=wt_A1),s]
      return(Ms.as)
    })
    Mtilde <- do.call(cbind,M.as)
    rm(M.as)
    data.marginal[[i]]<-cbind(mydt_mc[Marg==1,],Mtilde)
    names(data.marginal[[i]])<-c(colnames(dat),"mc",Mnames)
    rm(Mtilde)
    rm(mydt_mc)
    data.marginal[[i]][, "Y.a" := PredictY(onedat=data.marginal[[i]][, c("a0",Lnames,Mnames), with=FALSE])]
    setkey(data.marginal[[i]])
    ## average over MC draws
    marginal_means[[i]] <- data.marginal[[i]][,lapply(.SD, mean), by=c("Marg","a.i",paste0("a",0:t)),.SDcols="Y.a"]
    setkey(marginal_means[[i]])
   
    }
    marginal_means<-do.call(rbind,marginal_means)
    ## average over all individuals
    mu_hat_marginal <- marginal_means[, lapply(.SD,mean), by=c("Marg","a.i",paste0("a",0:t)),.SDcols="Y.a"]
    setnames(mu_hat_marginal,old="Y.a",new="Y")
    setkey(mu_hat_marginal)
    
  
  # Delta transformation matrix for effects
  mu_hat<-rbind(mu_hat_joint,mu_hat_marginal)
  eff_mat <- matrix(0,nrow=t+4,ncol=nrow(mu_hat))
  eff_mat[1,2:3] <- c(-1,1) # direct effect
  eff_mat[2,1:2] <- c(-1,1) # joint indirect effect
  for (s in 1:(t+1)) {
    eff_mat[s+2,c(4,s+4)] <- c(-1,1) # indirect effect via Ms
  }
  # indirect effect via mutual
  eff_mat[nrow(eff_mat),] <- eff_mat[2,] - eff_mat[nrow(eff_mat)-1,]
  
  if (fitY$family$family=="binomial") {
    mu_hat[, Y := log(Y/(1-Y))] # logit function
  }
  est <- (eff_mat %*% mu_hat[,Y])[,1]
  names(est) <- c(paste0("g",0:1),paste0("t",1:t),"t_margsum","mu")
  
  if (pt_est==TRUE) {
    var_est <- NULL
  }
  return(est)
}

 
#########################################################################################
#############################bootstrap###################################################
#########################################################################################
 
set.seed(13133)

#bootstrap of the estimation procedure
effect <- function(data, index) { 
  boot.data <- as.data.frame(data[index, ])
  
  L <- as.matrix(cbind(boot.data$e3_sex_None,boot.data$h_age_None,boot.data$h_cohort,boot.data$h_native_None, boot.data$h_parity_None))
  A <- as.numeric(boot.data$A)
  M <- as.matrix(boot.data[,8:(ncol(boot.data)-1)])
  Y <- as.numeric(boot.data$Y)
  boot.data<-boot.data[,-1] #ridefinire id variable
  boot.data <- cbind("id"=1:nrow(boot.data), boot.data)
 
 est<-estimator_noselection(boot.data,mc_draws=50,fitY.family="binomial",pt_est=TRUE)
 return(est)
}

replicates=50
options(nwarnings = 200) # per incrementare il numero max di items
timestamp()
bootSE <- boot(data=data, statistic=effect, R=replicates)
timestamp()
