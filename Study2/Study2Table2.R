#Common data-generating function

#popsize set at 1500000, can increase in large ssize is needed
#ssize set at 500, may have trouble when not enough patient available

#Defaults (scenario 1) true marg RR = 
#OR_C<-3 #effect of vaccine on covid
#OR_WI<-1 #no effect of vaccine on W from other infection
#OR_WC<-5 #effect of vaccine on covid symptoms
#OR_H<-1.5 #effect of vaccine on hospitalization among those with symptoms
datagen<-function(seed=sample(1:1000000,size=1),ssize=500,popsize=1500000,OR_C=1.5,OR_WI=1,OR_WC=3.5,OR_H=1.5,em=0,cfV0=F,cfV1=F,return_full=F){
  set.seed(seed)
  
  #generate data
  C<-runif(n=popsize, -3, 3)
  U1<-rbinom(n=popsize,size=1,prob=0.5) #affects both
  U2<-rbinom(n=popsize,size=1,prob=0.5) #affects covid
  
  if(cfV0==T) V=rep(0,popsize); if(cfV1==T) V=rep(1,popsize);
  if(cfV0==F&cfV1==F){
    V<-rbinom(prob=plogis(1.5 + 0.3*C - abs(C) - sin(pi*C) ),size=1,n=popsize) #prevalence is around 0.61
  }
  #Infection (with something) has some common risk factors U1 and C
  Infec<-rbinom(prob=plogis(0.5*C-5+0.5*U1),size=1,n=popsize) #current prevalence around 0.007
  
  #Infected with COVID
  Infec_COVID<- rbinom(prob=plogis( -log(OR_C)*V -4 + 2*C - 0.15*exp(C)+log(2)*U2*(1.5-V)-2*U1), size=1,n=popsize) #0.009
  #symptoms based on infection
  #can come from either or both infections, if present
  W=rep(0,popsize)
  W[Infec==1]<-rbinom(prob=plogis(2+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
  W[Infec_COVID==1]<-rbinom(prob=plogis(-5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
  #mean(W[Infec==1|Infec_COVID==1]) #25%
  #mean(W[Infec_COVID==1]) #39%
  #mean(W[Infec==1]) #12%
  #mean(W[Infec_COVID==1&V==1]) #22%
  
  #hospitalization, only possible if symptoms present
  H=rep(0,popsize)
  H[W==1]<-rbinom(prob=plogis(1+0.5*C[W==1]+log(OR_H)*V[W==1]-0.5*U1[W==1]),size=1,n=sum(W==1))
  #mean(H[W==1]) #83% with severe symptoms go to hospital
  
  #selection on outcome for testing (does not condition on infectious status, just being in the hospital)
  R<-sample(which(H==1),ssize) #sample randomly from those in hospital to control the study size
  
  if(return_full==F){
    dat<-as.data.frame(cbind(Y=Infec_COVID,V=V,C=C)[R,])
  } else{dat<-as.data.frame(cbind(Infec_COVID=Infec_COVID,Infec=Infec,H=H,W=W,V=V,C=C))}
  return(dat)
}

ssize=5000
dat<-datagen(ssize= ssize, em=0)
sum(dat$Y)/ssize

######################################################################
###Run Simulation Study
# IPW estimator
library(sandwich)
# IPW correct
mod_IPW_c <- function(dat){
  TNDmod_g_col<-glm(V ~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  ipw <- ifelse(dat$V == 1, 1/g1_cont, 1/(1 - g1_cont))
  modY.ipw <- glm(Y ~ V, family=binomial(link = "logit"), weights = ipw, data=dat)
  est.ipw <- exp(modY.ipw$coefficients[2])
  se.ipw <- sqrt(vcovHC(modY.ipw)[2,2])
  CI_l <- est.ipw *exp(- 1.96 * se.ipw )
  CI_u <- est.ipw *exp( 1.96 * se.ipw )
  return(list(est = est.ipw, se = se.ipw, CI =  c(CI_l, CI_u)))
}
# IPW wrong
mod_IPW_w <- function(dat){
  TNDmod_g_col<-glm(V ~ 1 ,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  ipw <- ifelse(dat$V == 1, 1/g1_cont, 1/(1 - g1_cont))
  modY.ipw <- glm(Y ~ V, family=binomial(link = "logit"), weights = ipw, data=dat)
  est.ipw <- exp(modY.ipw$coefficients[2])
  se.ipw <- sqrt(vcovHC(modY.ipw)[2,2])
  CI_l <- est.ipw *exp(- 1.96 * se.ipw )
  CI_u <- est.ipw *exp( 1.96 * se.ipw )
  return(list(est = est.ipw, se = se.ipw, CI =  c(CI_l, CI_u)))
}

# outcome regression correct
mod_OR_c1 <- function(dat){
  TNDmod<-glm(Y~ V + C + exp(C),family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  
  TNDmod1<-glm(Y ~ C + exp(C),family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  preY <- pmax(pmin(preY, 1), 0)
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  est <- mean(Q1)/mean(Q0)
  return(list(est = est))
}

mod_OR_c <- function(dat){
  est <- mod_OR_c1(dat)$est
  nbs <- 50
  bsest<-rep(NA,nbs)
  for(i in 1:nbs){
    resamps<-sample(1:nrow(dat),size=nrow(dat),replace=T)
    datk<-dat[resamps,]
    bsest[i] <- mod_OR_c1(datk)$est
  }
  bs_var <- var(bsest)
  CI = quantile(bsest,c(0.025,0.975))
  return( list(est = est, CI = as.vector(CI)) )
}


#  outcome regression wrong
mod_OR_w1 <- function(dat){
  TNDmod<-glm(Y ~ V,family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  TNDmod1<-glm(Y ~ 1,family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  preY <- pmax(pmin(preY, 1), 0)
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  (est <- mean(Q1)/mean(Q0))
  return(list(est = est))
}

mod_OR_w <- function(dat){
  est <- mod_OR_w1(dat)$est
  nbs <- 50
  bsest<-rep(NA,nbs)
  for(i in 1:nbs){
    resamps<-sample(1:nrow(dat),size=nrow(dat),replace=T)
    datk<-dat[resamps,]
    bsest[i] <- mod_OR_w1(datk)$est
  }
  bs_var <- var(bsest)
  CI = quantile(bsest,c(0.025,0.975))
  return( list(est = est, CI = as.vector(CI)) )
}

######################################################################
# EIF estimator 1 (Equation 10):
# S1) both are correct
modEIF1a<-function(dat){
  TNDmod_g_col<-glm(V~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)));
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  g0_cont <- 1 - g1_cont
  TNDmod1<-glm(Y ~ C + exp(C),family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  TNDmod<-glm(Y~ V + C + exp(C),family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  preY <- pmax(pmin(preY, 1), 0)
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - preY))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - preY))
  psi.1 <- mean(dat$Y*dat$V/g1_cont - Q1*A1)
  psi.0 <- mean(dat$Y*(1-dat$V)/g0_cont - Q0*A0)
  mod_eif <- psi.1/psi.0
  eifln <-  ((dat$Y*dat$V/g1_cont - Q1*A1)/psi.1) - ((dat$Y*(1-dat$V)/g0_cont - Q0*A0)/ psi.0)
  varln <-  var(eifln)/nrow(dat)
  
  CI_l1 <- exp(log(mod_eif) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif) + 1.96 * sqrt(varln) )
  
  eifpsi <- (dat$Y*dat$V/g1_cont - Q1*A1)/psi.0 - (psi.1/psi.0)*(dat$Y*(1-dat$V)/g0_cont - Q0*A0)/psi.0
  var <- var(eifpsi)/nrow(dat)
  CI_l2 <- mod_eif - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif + 1.96 * sqrt(var) 
  return(list( est = mod_eif, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}

modEIF1a(dat)
# S2) Out is correct, but PS is wrong
modEIF1b<-function(dat){
  TNDmod_g_col<-glm(V~  1,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  g0_cont <- 1 - g1_cont
  TNDmod1<-glm(Y ~ C + exp(C),family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  
  TNDmod<-glm(Y~ V + C + exp(C),family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  preY <- pmax(pmin(preY, 1), 0)
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - preY))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - preY))
  psi.1 <- mean(dat$Y*dat$V/g1_cont - Q1*A1)
  psi.0 <- mean(dat$Y*(1-dat$V)/g0_cont - Q0*A0)
  mod_eif <- psi.1/psi.0
  eifln <-  ((dat$Y*dat$V/g1_cont - Q1*A1)/psi.1) - ((dat$Y*(1-dat$V)/g0_cont - Q0*A0)/ psi.0)
  varln <-  var(eifln)/nrow(dat)
  
  CI_l1 <- exp(log(mod_eif) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif) + 1.96 * sqrt(varln) )
  
  eifpsi <- (dat$Y*dat$V/g1_cont - Q1*A1)/psi.0 - (psi.1/psi.0)*(dat$Y*(1-dat$V)/g0_cont - Q0*A0)/psi.0
  var <- var(eifpsi)/nrow(dat)
  CI_l2 <- mod_eif - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif + 1.96 * sqrt(var) 
  return(list( est = mod_eif, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}
# S3) PS is correct, but Out is wrong
modEIF1c<-function(dat){
  TNDmod_g_col<-glm(V~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  g0_cont <- 1 - g1_cont
  TNDmod1<-glm(Y ~ 1,family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  
  TNDmod<-glm(Y~ V,family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  preY <- pmax(pmin(preY, 1), 0)
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - preY))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - preY))
  psi.1 <- mean(dat$Y*dat$V/g1_cont - Q1*A1)
  psi.0 <- mean(dat$Y*(1-dat$V)/g0_cont - Q0*A0)
  mod_eif <- psi.1/psi.0
  eifln <-  ((dat$Y*dat$V/g1_cont - Q1*A1)/psi.1) - ((dat$Y*(1-dat$V)/g0_cont - Q0*A0)/ psi.0)
  varln <-  var(eifln)/nrow(dat)
  
  CI_l1 <- exp(log(mod_eif) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif) + 1.96 * sqrt(varln) )
  
  eifpsi <- (dat$Y*dat$V/g1_cont - Q1*A1)/psi.0 - (psi.1/psi.0)*(dat$Y*(1-dat$V)/g0_cont - Q0*A0)/psi.0
  var <- var(eifpsi)/nrow(dat)
  CI_l2 <- mod_eif - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif + 1.96 * sqrt(var) 
  return(list( est = mod_eif, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}
# S4) Both are wrong
modEIF1d<-function(dat){
  TNDmod_g_col<-glm(V~  1,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  g0_cont <- 1 - g1_cont
  TNDmod1<-glm(Y ~ 1,family=binomial(),data=dat)
  preY=predict(TNDmod1,type="response")
  
  TNDmod<-glm(Y~ V,family=binomial(),data=dat)
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  preY <- pmax(pmin(preY, 1), 0)
  w1 = (1 - preY)/(1 - mu1)
  w0 = (1 - preY)/(1 - mu0)
  Q1 = w1 * mu1
  Q0 = w0 * mu0
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - preY))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - preY))
  psi.1 <- mean(dat$Y*dat$V/g1_cont - Q1*A1)
  psi.0 <- mean(dat$Y*(1-dat$V)/g0_cont - Q0*A0)
  mod_eif <- psi.1/psi.0
  eifln <-  ((dat$Y*dat$V/g1_cont - Q1*A1)/psi.1) - ((dat$Y*(1-dat$V)/g0_cont - Q0*A0)/ psi.0)
  varln <-  var(eifln)/nrow(dat)
  
  CI_l1 <- exp(log(mod_eif) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif) + 1.96 * sqrt(varln) )
  
  eifpsi <- (dat$Y*dat$V/g1_cont - Q1*A1)/psi.0 - (psi.1/psi.0)*(dat$Y*(1-dat$V)/g0_cont - Q0*A0)/psi.0
  var <- var(eifpsi)/nrow(dat)
  CI_l2 <- mod_eif - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif + 1.96 * sqrt(var) 
  return(list( est = mod_eif, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}


######################################################################
### EIF 2 (Equation 11)
# EIF estimator 2:S1) both are correct
modEIF2a<-function(dat){
  TNDmod_g_col<-glm(V~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V + C + exp(C),family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  psi.1 <- mean(dat$Y*dat$V/g1_cont - mu1*A1)
  psi.0 <- mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  mod_eif2 <- psi.1/psi.0
  eifln <-  ((dat$Y*dat$V/g1_cont - mu1*A1)/psi.1) - ((dat$Y*(1-dat$V)/g0_cont - mu0*A0)/ psi.0)
  varln <-  var(eifln)/nrow(dat)
  
  CI_l1 <- exp(log(mod_eif2) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif2) + 1.96 * sqrt(varln) )
  
  eifpsi <- (dat$Y*dat$V/g1_cont - mu1*A1)/psi.0 - (psi.1/psi.0)*(dat$Y*(1-dat$V)/g0_cont - mu0*A0)/psi.0
  var <- var(eifpsi)/nrow(dat)
  CI_l2 <- mod_eif2 - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif2 + 1.96 * sqrt(var) 
  return(list( est = mod_eif2, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}

modEIF2a(dat)
modEIF2b<-function(dat){
  TNDmod_g_col<-glm(V~ 1,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V + C + exp(C),family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  psi.1 <- mean(dat$Y*dat$V/g1_cont - mu1*A1)
  psi.0 <- mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  mod_eif2 <- psi.1/psi.0
  eifln <-  ((dat$Y*dat$V/g1_cont - mu1*A1)/psi.1) - ((dat$Y*(1-dat$V)/g0_cont - mu0*A0)/ psi.0)
  varln <-  var(eifln)/nrow(dat)
  
  CI_l1 <- exp(log(mod_eif2) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif2) + 1.96 * sqrt(varln) )
  
  eifpsi <- (dat$Y*dat$V/g1_cont - mu1*A1)/psi.0 - (psi.1/psi.0)*(dat$Y*(1-dat$V)/g0_cont - mu0*A0)/psi.0
  var <- var(eifpsi)/nrow(dat)
  CI_l2 <- mod_eif2 - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif2 + 1.96 * sqrt(var) 
  return(list( est = mod_eif2, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}

modEIF2c<-function(dat){
  TNDmod_g_col<-glm(V~  C + abs(C) +  sin(pi*C),family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V,family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  psi.1 <- mean(dat$Y*dat$V/g1_cont - mu1*A1)
  psi.0 <- mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  mod_eif2 <- psi.1/psi.0
  eifln <-  ((dat$Y*dat$V/g1_cont - mu1*A1)/psi.1) - ((dat$Y*(1-dat$V)/g0_cont - mu0*A0)/ psi.0)
  varln <-  var(eifln)/nrow(dat)
  
  CI_l1 <- exp(log(mod_eif2) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif2) + 1.96 * sqrt(varln) )
  
  eifpsi <- (dat$Y*dat$V/g1_cont - mu1*A1)/psi.0 - (psi.1/psi.0)*(dat$Y*(1-dat$V)/g0_cont - mu0*A0)/psi.0
  var <- var(eifpsi)/nrow(dat)
  CI_l2 <- mod_eif2 - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif2 + 1.96 * sqrt(var) 
  return(list( est = mod_eif2, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}

modEIF2d<-function(dat){
  TNDmod_g_col<-glm(V~ 1,family=binomial(),data=dat,subset=(dat$Y==0))
  g1_cont<-predict(TNDmod_g_col,type="response",newdata=as.data.frame(cbind(C=dat$C,V=dat$V,Y=dat$Y)))
  g1_cont <- pmax(pmin(g1_cont, 1), 0)
  g0_cont <- 1 - g1_cont
  
  TNDmod<-(glm(Y~ V,family=binomial(),data=dat)) 
  mu1=predict(TNDmod,newdata=as.data.frame(cbind(V=1,C=dat$C)),type="response")
  mu0=predict(TNDmod,newdata=as.data.frame(cbind(V=0,C=dat$C)),type="response")
  mu1 <- pmax(pmin(mu1, 1), 0)
  mu0 <- pmax(pmin(mu0, 1), 0)
  A1 <- (1 - dat$Y)*(dat$V - g1_cont)/(g1_cont* (1 - mu1))
  A0 <- (1 - dat$Y)*((1-dat$V) - g0_cont)/(g0_cont* (1 - mu0))
  
  psi.1 <- mean(dat$Y*dat$V/g1_cont - mu1*A1)
  psi.0 <- mean(dat$Y*(1-dat$V)/g0_cont - mu0*A0)
  mod_eif2 <- psi.1/psi.0
  eifln <-  ((dat$Y*dat$V/g1_cont - mu1*A1)/psi.1) - ((dat$Y*(1-dat$V)/g0_cont - mu0*A0)/ psi.0)
  varln <-  var(eifln)/nrow(dat)
  
  CI_l1 <- exp(log(mod_eif2) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif2) + 1.96 * sqrt(varln) )
  
  eifpsi <- (dat$Y*dat$V/g1_cont - mu1*A1)/psi.0 - (psi.1/psi.0)*(dat$Y*(1-dat$V)/g0_cont - mu0*A0)/psi.0
  var <- var(eifpsi)/nrow(dat)
  CI_l2 <- mod_eif2 - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif2 + 1.96 * sqrt(var) 
  return(list( est = mod_eif2, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}

CI1=CI2=CI3=CI3.2=CI4=CI4.2=CI5=CI5.2=CI6=CI6.2=CI7=CI7.2=CI8=CI8.2=CI9=CI9.2=CI10=CI10.2=CI11 = CI12 = c(0,0)

for (i in 1:1000){
  dat<-datagen(ssize=250, em=0)
  #######################################################
  # Marginal RR 
  # IPW ps correct
  est1<-mod_IPW_c(dat)$est
  CI1 <- mod_IPW_c(dat)$CI
  # IPW ps wrong
  est2<-mod_IPW_w(dat)$est
  CI2 <- mod_IPW_w(dat)$CI
  #######################################################
  # EIF estimator 1
  #Case 1: All correct
  res3 <- modEIF1a(dat)
  est3 <- res3$est
  CI3 <-  res3$CI1
  CI3.2 <-  res3$CI2
  #######################################################
  #Case 2: Outcome model is right, but PS is not
  res4 <- modEIF1b(dat)
  est4 <- res4$est
  CI4 <-  res4$CI1
  CI4.2 <-  res4$CI2
  #######################################################
  # Case 3: (PS is correct, but the outcome model is not correct)
  res5 <- modEIF1c(dat)
  est5 <- res5$est
  CI5 <-  res5$CI1
  CI5.2 <-  res5$CI2
  #######################################################
  # Both are not correct
  res6 <- modEIF1d(dat)
  est6 <- res6$est
  CI6 <-  res6$CI1
  CI6.2 <-  res6$CI2
  #######################################################
  # EIF estimator 2
  #Case 1: All correct
  res7 <- modEIF2a(dat)
  est7 <- res7$est
  CI7 <-  res7$CI1
  CI7.2 <-  res7$CI2
  #######################################################
  #Case 2: Outcome model is right, but PS is not
  res8 <- modEIF2b(dat)
  est8 <- res8$est
  CI8 <-  res8$CI1
  CI8.2 <-  res8$CI2
  #######################################################
  # Case 3: (PS is correct, but the outcome model is not correct)
  res9 <- modEIF2c(dat)
  est9 <- res9$est
  CI9 <-  res9$CI1
  CI9.2 <-  res9$CI2
  #######################################################
  # Both are not correct
  res10 <- modEIF2d(dat)
  est10 <- res10$est
  CI10 <-  res10$CI1
  CI10.2 <-  res10$CI2
  #######################################################
  # OR estimators
  est11 <- mod_OR_c(dat)$est
  CI11 <- mod_OR_c(dat)$CI
  est12 <- mod_OR_w(dat)$est
  CI12 <- mod_OR_w(dat)$CI
  write(c(i,est1,CI1, est2,CI2, est3,CI3,CI3.2,est4,CI4,CI4.2, est5,CI5,CI5.2, 
          est6,CI6,CI6.2, est7,CI7,CI7.2,est8,CI8, CI8.2,est9,CI9, CI9.2,est10,CI10,CI10.2, est11, CI11, est12, CI12),file="Study2results250_0808CI.txt",ncolumns=70,append=T)
}




res1<-read.table("Study2results250_0808CI.txt",header=F)
head(res1)
res1 <- na.omit(res1)
dim(res1)
#truth mRR
psi = 0.128
colnames(res1)[c(2, 5, seq(8, 44,5), 48,51)] <- c("IPWCorrect", "IPWrong", "BothCorrect1", "OutCorrect1", "PSCorrect1", "BothWrong1", "BothCorrect2", "OutCorrect2", "PSCorrect2", "BothWrong2",  "ORCorrect", "ORWrong")
(apply(res1[, c(2, 5, seq(8, 44,5), 48,51)], 2, mean) - psi)
(apply(res1[, c(2, 5, seq(8, 44,5), 48,51)], 2, median) - psi)
apply(res1[,c(2, 5, seq(8, 44,5), 48,51)], 2, sd)


#CIs
mean(psi<=res1$V4 & psi>=res1$V3) # correct ipw
mean(psi<=res1$V7 & psi>=res1$V6) # incorrect ipw
mean(psi<=res1$V10 & psi>=res1$V9) # bothcorrect ief 1
mean(psi<=res1$V30 & psi>=res1$V29) # bothcorrect tnddr
mean(psi<=res1$V35 & psi>=res1$V34) # out cor tnddr
mean(psi<=res1$V40 & psi>=res1$V39) # ps cor tnddr
mean(psi<=res1$V45 & psi>=res1$V44) # bothcorrect itnddr


mean(psi<=res1$V50 & psi>=res1$V49) # out correct
mean(psi<=res1$V53 & psi>=res1$V52) # out wrong

