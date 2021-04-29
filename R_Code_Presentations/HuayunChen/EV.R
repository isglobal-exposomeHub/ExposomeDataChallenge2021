############################################################################
#  Function EV compute explained variation of the outcome
#    by the exposure adjusted for the covariate.
#      1. outcome: a vector contains data on the outcome
#      2. exposure: a matrix contains data on the multiple exposures
#      3. covariate: a mtrix contains data on the covariates to be adjusted
#    Function interaction create pairwise interaction variables, which
#       is needed only when interactions are considered.
#    Function EV calls function R2estequ for the computation.
#############################################################################
EV=function(outcome,exposure,covariate){
  # Estimate variation of outcome explained by exposure 
  # after adjusted for covariates
  # outcome[n]
  # exposure[n,p]
  # covariate[n,q]
  n=length(outcome)
  p=dim(exposure)[2]
  q=dim(covariate)[2]

#1. Adjust for covariates
  if(q>0){
    for(j in 1:p){
      fit=lm(exposure[,j]~covariate)
      exposure[,j]=fit$residual
      }
  }
#2. Decorrelation by linear transformation
 corexp=cor(exposure)
 Esvd=svd(t(exposure))
 tran=Esvd$u%*%diag(1/Esvd$d)%*%t(Esvd$u)/sqrt(n)
 x=exposure%*%tran

# Analyze the main effect of exposure

   for(j in 1:p){
     mu=mean(x[,j])
     sdx=sd(x[,j])
     x[,j]=(x[,j]-mu)/sdx
     }
   sdy=sd(outcome)
   y=(outcome-mean(outcome))/sdy


  fit1=R2EstEqu(y,x,lam=0.1,niter=5)
  #fit2=EigenPrismFull(y,x)  
  #fit3=R2GCTA(y,x,lam=0.2,bt='T', btn=2000)
  #return(list(fit1,fit2,fit3))
  return(fit1)
}

interaction=function(x){
   # create interaction terms for x
  n=dim(x)[1]
  p=dim(x)[2]

  xx=array(0,c(n,(p-1)*(p-2)/2))
  ii=0
  for(i in 1:(p-2)){for(j in (i+1):(p-1)){
    ii=ii+1
    xx[,ii]=x[,i]*x[,j]
   }}
  list(xx)
}

R2EstEqu=function(y,x,lam,niter=1){
   # Fast estimate r^2 without additional weighting
   # y==outcome
   # x==covariates
   # lam==parameter adjusting the format of weighting matrix

   n=dim(x)[1]
   p=dim(x)[2]
   
#1. Standardization

   for(j in 1:p){
     mu=mean(x[,j])
     sdx=sd(x[,j])
     x[,j]=(x[,j]-mu)/sdx
     }
   sdy=sd(y)
   y=(y-mean(y))/sdy

#2. Sigular value decomposition
   Xsvd=svd(x,nv=0) #nv=0 means not computing v matrix 
     # singular value decomposition
     # $u%*%diag($d)%*%t($v)=X, t($u)%*%$u=I, t($v)%*%$v=I 
   Mev=Xsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.  

 for(ii in 1:niter){ #iteration to update lambda                   
   Dev=(Mev-1)/(1+lam*Mev)^2  #vector of eigenvalues of weight matrix 

#3. Compute the estimators 
   uy=t(Xsvd$u)%*%y
   u1=t(Xsvd$u)%*%rep(1,n)
   if(n>=p){
     com=sum(u1^2*(Dev+1))/n-1
     num=sum(uy^2*(Dev+1))-sum(y^2)-sum(Dev)+n-p
     den=sum(Dev*(Mev-1))+n-p
   }else{
     com=sum(u1^2*Dev)/n
     num=sum(uy^2*Dev)-sum(Dev)
     den=sum(Dev*(Mev-1))  
    }

   r2=(num+com)/(den+com)
   if(r2>1){r2=1}
   if(r2<0){r2=0}

  lam=r2
  }
   # estimate variance

   Wev=Mev*Dev
   if(n>=p){
     W=Xsvd$u%*%diag(Dev+1)%*%t(Xsvd$u)-diag(rep(1,n))
     vest=2*p*var(Wev)*r2^2
     vest=vest+4*r2*(1-r2)*sum(Dev^2*Mev)
     vest=vest+2*(1-r2)^2*(sum(Dev^2)+n-p)
   }else{
     W=Xsvd$u%*%diag(Dev)%*%t(Xsvd$u)
     vest=2*n*var(Wev)*r2^2 #because the length of Wev is n now.
     vest=vest+4*r2*(1-r2)*sum(Dev^2*Mev)
     vest=vest+2*(1-r2)^2*sum(Dev^2)   
   }
   M=Xsvd$u%*%diag(Mev)%*%t(Xsvd$u)
   add=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)#-r2^2
   add=add-r2^2*(sum(Mev)/n)^2
   vest1=vest+sum(diag(W^2))*(add*(add>0)-2*(1-r2)^2)
   #vest1=vest1-r2^2*(2*sum(Mev^2)-sum(diag(M^2)))*mean(diag(W^2))

   vest1=vest1/den^2
   vest=vest/den^2

   ci=r2+sqrt(vest)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
   ci1=r2+sqrt(vest1)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
   ci[2*c(1:3)-1]=ci[2*c(1:3)-1]*(ci[2*c(1:3)-1]>0)
   ci[2*c(1:3)]=ci[2*c(1:3)]*(ci[2*c(1:3)]<1)+1.0*(ci[2*c(1:3)]>=1)
   ci1[2*c(1:3)-1]=ci1[2*c(1:3)-1]*(ci1[2*c(1:3)-1]>0)
   ci1[2*c(1:3)]=ci1[2*c(1:3)]*(ci1[2*c(1:3)]<1)+1.0*(ci1[2*c(1:3)]>=1)

#5. output result 

   list(c(r2,vest,vest1),ci,ci1)
   #[[1]]==> estimate, estimated variance (under normal), estimated variance.
   #[[2]]==>confidence interval under normal
   #[[3]]==>confidence interval without normal 
}
############################################################################################
# Output of EV has three parts which are
#    [[1]]=(EV estimate, variance estimate under normal, variance estimate without normal)
#    [[2]]=(90% CI, 95% CI, 99% CI) under normality assumption.
#    [[3]]=(90% CI, 95% CI, 99% CI) without normality assumption.
############################################################################################ 

