#####################################################
#     1. Analysis of birth weight
#
#####################################################
noadj=matrix(rep(1,length(bw)),ncol=1)
  ##for use in the place of covariate when no covariate is adjusted.

#1.1 Linear model with variable selection
 
tcomb=cbind(metalm,chemicalm,airm,builtm,lifestylem,outdoorm,covariate)
ecomb=cbind(metalm,chemicalm,airm,builtm,lifestylem,outdoorm)

summary(lm(bw~tcomb))
 acomb=cbind(metalm[,4],chemicalm[,7],chemicalm[,10],chemicalm[,13],
   chemicalm[,24],chemicalm[,32],builtm[,6],lifestylem[,4],outdoorm[,11],
   covariate[,1:4],covariate[,8],covariate[,10])
summary(lm(bw~acomb))
 bcomb=cbind(chemicalm[,13],outdoorm[,11],
   covariate[,1:4],covariate[,8],covariate[,10])
summary(lm(bw~bcomb))

summary(lm(bw~ecomb))
 xx=ecomb
 ccomb=cbind(xx[,4],xx[,16],xx[,17],xx[,23],xx[,34],xx[,44],xx[,51],xx[,54],
             xx[,58],xx[,61],xx[,65],xx[,67],xx[,84],xx[,86])
summary(lm(bw~ccomb))
 dcomb=cbind(xx[,51],xx[,84],xx[,86])
summary(lm(bw~dcomb))

#1.2. EV analysis without interaction
 #1.2.1. no covariate adjustment

EV(bw,covariate,noadj)
EV(bw,ecomb,noadj)
EV(bw,dcomb,noadj)
EV(bw,tcomb,noadj)
EV(bw,bcomb,noadj)

 #1.2.2. covariate adjusted analyses

EV(bw,ecomb,covariate)
EV(bw,metalm,covariate)
EV(bw,chemicalm,covariate)
EV(bw,lifestylem,covariate)

#1.3. EV analysis of interactions
  
zz=interaction(metalm)[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)  
zz=interaction(chemicalm)[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)
zz=interaction(lifestylem)[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)
zz=interaction(builtm)[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)
zz=interaction(cbind(airm,outdoorm))[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)

#############################################################
#         2. Analysis of intelligence quotient
#
#############################################################
noadj=matrix(rep(1,length(IQ)),ncol=1)
      
#2.1 Linear model with variable selection
             
iq=log(1.0+IQ)
tchemicalc=cbind(chemicalc[,1:38],schemc)
tcovariate=cbind(covariate,scovariate)
lifestylec=cbind(lifestylec[,1:19],lifestylec[,21:27])
outdoorc=cbind(outdoorc[,1:15],outdoorc[,17:19])
 
 tcomb=cbind(metalc,tchemicalc,airc,builtc,lifestylec,outdoorc,indoorc,tcovariate,bw)
 ecomb=cbind(metalc,tchemicalc,airc,builtc,lifestylec,outdoorc,indoorc)

summary(lm(iq~tcomb))
 xx=tcomb 
 fcomb=cbind(xx[,3],xx[,9],xx[,18],xx[,25],xx[,27],xx[,29],xx[,32],xx[,38],xx[,48],
         xx[,58],xx[,60],xx[,64],xx[,86],xx[,97],xx[,101],xx[,108],xx[,112],
         xx[,119],xx[,123],xx[,128],xx[,135],xx[,140],xx[,142],xx[,143])
summary(lm(iq~fcomb))
 gcomb=cbind(xx[,9],xx[,18],xx[,29],xx[,38],
         xx[,58],xx[,97],xx[,108],xx[,112],
         xx[,128],xx[,140],xx[,143])
summary(lm(iq~gcomb))
 hcomb=cbind(xx[,9],xx[,18],xx[,29],xx[,38],
         xx[,58],xx[,112],xx[,128],xx[,140])
summary(lm(iq~hcomb))

#2.2. EV analysis without interaction
 #2.2.1. no covariate adjustment

EV(iq,tcovariate,noadj)
EV(iq,ecomb,noadj)
EV(iq,tcomb,noadj)
EV(iq,xxx,noadj) 

#2.2.2. covariate adjusted analyses

EV(iq,ecomb,tcovariate)
EV(iq,metalc,tcovariate)
EV(iq,tchemicalc,tcovariate)
EV(iq,lifestylec,tcovariate)
 
#2.3. EV analysis of interactions
 
zz=interaction(metalc)[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)  
zz=interaction(tchemicalc)[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)  
zz=interaction(lifestylec)[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)  
zz=interaction(builtc)[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)  
zz=interaction(cbind(airc,outdoorc,indoorc))[[1]]
EV(bw,zz,ecomb)
EV(bw,zz,tcomb)  



 ecomb=cbind(metalc,tchemicalc,airc,builtc,lifestylec,outdoorc,indoorc)
summary(lm(iq~ecomb))
 xx=ecomb 
 fcomb=cbind(xx[,3],xx[,9],xx[,18],xx[,25],xx[,27],xx[,29],xx[,32],xx[,33],
             xx[,38],xx[,48],xx[,58],xx[,60],xx[,61],xx[,64],xx[,72],xx[,86],
             xx[,97],xx[,101],xx[,108],xx[,112])
summary(lm(iq~fcomb))
 gcomb=cbind(xx[,3],xx[,9],xx[,18],xx[,29],xx[,32],xx[,33],
             xx[,38],xx[,58],xx[,61],
             xx[,97],xx[,108],xx[,112])
summary(lm(iq~gcomb))

EV(iq,ecomb,noadj)
EV(iq,gcomb,noadj)
 
#################################
#     3. Other outcomes
#
#################################
EV(bmizscore,tcovariate,noadj)
EV(bmizscore,tcomb,noadj)
EV(bmizscore,ecomb,noadj)
EV(bmizscore,ecomb,tcovariate)

EV(beh,tcovariate,noadj)
EV(beh,tcomb,noadj)
EV(beh,ecomb,noadj)
EV(beh,ecomb,tcovariate)

