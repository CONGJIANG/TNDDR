#ssize set at 500, may have trouble when not enough patient available
datagen<-function(seed=sample(1:1000000,size=1),ssize=500,popsize=1500000,OR_C=1.5,OR_WI=1,OR_WC=3.5,OR_H=1.5,em=1,cfV0=F,cfV1=F,return_full=F){
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
  Infec_COVID<- rbinom(prob=plogis( -log(OR_C)*V -4 + 2*C - 0.15*exp(C)+ em*V*C + log(2)*U2*(1.5-V)-2*U1), size=1,n=popsize) #0.009
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


TNDdat = datagen(ssize=5000, em=1)
ssize=5000
sum(TNDdat$Y)/ssize  # 0.1308

library(earth)
library(dplyr)

element.earth <- function(TNDdf_train){
  # Data splitting
  s <- sample(1:length(TNDdf_train$Y), length(TNDdf_train$Y)/2)
  TNDdf_train1 <- TNDdf_train[s,]
  TNDdf_train2 <-  TNDdf_train[-s,]
  
  TNDdf_train_ctr1 <- subset(TNDdf_train1, Y==0)
  # training
  mod_g1_ctr <- earth(
    V ~ .,  
    data = subset(TNDdf_train_ctr1, select = -Y),
    glm=list(family=binomial)
  )
  
  TNDdf_train_ctr2 <- subset(TNDdf_train2, Y==0)
  # training
  mod_g2_ctr <- earth(
    V ~ .,  
    data = subset(TNDdf_train_ctr2, select = -Y),
    glm=list(family=binomial)
  )
  g1_cont <- TNDdf_train$V
  g1_cont[-s]<-predict(mod_g1_ctr,type="response",newdata=as.data.frame(cbind(C=TNDdf_train2$C,V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y)))
  g1_cont[s]<-predict(mod_g2_ctr,type="response",newdata=as.data.frame(cbind(C=TNDdf_train1$C,V=rep(1, nrow(TNDdf_train1)),Y=TNDdf_train1$Y)))
  
  
  # mu
  Out_mu1 <- earth(
    Y ~ .,  
    data = TNDdf_train1,
    glm=list(family=binomial)
  )
  Out_mu2 <- earth(
    Y ~ .,  
    data = TNDdf_train2,
    glm=list(family=binomial)
  )
  
  mu1 <- TNDdf_train$Y; mu0 <- TNDdf_train$Y;
  mu1[-s] <- predict(Out_mu1,newdata=as.data.frame(cbind(V=1,C = select(TNDdf_train2, !c(V,Y)) )),type="response")
  mu1[s] <- predict(Out_mu2,newdata=as.data.frame(cbind(V=1,C = select(TNDdf_train1, !c(V,Y)) )),type="response")               
  
  mu0[-s] <- predict(Out_mu1,newdata=as.data.frame(cbind(V=0,C = select(TNDdf_train2, !c(V,Y)) )),type="response")
  mu0[s] <- predict(Out_mu2,newdata=as.data.frame(cbind(V=0,C = select(TNDdf_train1, !c(V,Y)) )),type="response") 
  
  # m0
  Out_m1 <- earth(
    Y ~ .,  
    data = subset(TNDdf_train1, select = -V),
    glm=list(family=binomial)
  )
  Out_m2 <- earth(
    Y ~ .,  
    data = subset(TNDdf_train2, select = -V),
    glm=list(family=binomial)
  )
  m0 <- TNDdf_train$Y
  m0[-s] <- 1- predict(Out_m1,newdata = select(TNDdf_train2, !c(V,Y)), type="response")
  m0[s] <- 1- predict(Out_m2,newdata = select(TNDdf_train1, !c(V,Y)), type="response")
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 =g1_cont, g0 = 1- g1_cont, w1 = m0/(1 - mu1), w0 = m0/(1 - mu0)))
}

element <- function(TNDdf_train){
  # Data splitting
  s <- sample(1:length(TNDdf_train$Y), length(TNDdf_train$Y)/2)
  TNDdf_train1 <- TNDdf_train[s,]
  TNDdf_train2 <-  TNDdf_train[-s,]
  
  TNDdf_train_ctr1 <- subset(TNDdf_train1, Y==0)
  # training
  mod_g1_ctr <- fit_hal(
    X = select(TNDdf_train_ctr1, !c(V,Y)),
    Y = TNDdf_train_ctr1$V,
    family = "binomial"
  )
  
  TNDdf_train_ctr2 <- subset(TNDdf_train2, Y==0)
  # training
  mod_g2_ctr <- fit_hal(
    X = select(TNDdf_train_ctr2, !c(V,Y)),
    Y = TNDdf_train_ctr2$V,
    family = "binomial"
  )
  g1_cont <- TNDdf_train$V
  g1_cont[-s]<-predict(mod_g1_ctr,type="response",new_data=as.data.frame(cbind(C=TNDdf_train2$C,V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y)))
  g1_cont[s]<-predict(mod_g2_ctr,type="response",new_data=as.data.frame(cbind(C=TNDdf_train1$C,V=rep(1, nrow(TNDdf_train1)),Y=TNDdf_train1$Y)))
  
  #PS full model
  mod_g1_full <- fit_hal(
    X = select(TNDdf_train1, !c(V,Y)),
    Y = TNDdf_train1$V,
    family = "binomial"
  )
  
  mod_g2_full <- fit_hal(
    X = select(TNDdf_train2, !c(V,Y)),
    Y = TNDdf_train2$V,
    family = "binomial"
  )
  g1_full <- TNDdf_train$V
  g1_full[-s]<-predict(mod_g1_full,type="response",new_data=as.data.frame(cbind(C=TNDdf_train2$C,V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y)))
  g1_full[s]<-predict(mod_g2_full,type="response",new_data=as.data.frame(cbind(C=TNDdf_train1$C,V=rep(1, nrow(TNDdf_train1)),Y=TNDdf_train1$Y)))
  
  
  # mu
  Out_mu1 <- fit_hal(X = select(TNDdf_train1, !Y),
                     Y = TNDdf_train1$Y,
                     family = "binomial"
  )
  Out_mu2 <- fit_hal(X = select(TNDdf_train2, !Y),
                     Y = TNDdf_train2$Y,
                     family = "binomial"
  )
  mu1 <- TNDdf_train$Y; mu0 <- TNDdf_train$Y;
  mu1[-s] <- predict(Out_mu1,new_data=as.data.frame(cbind(V=1,C = select(TNDdf_train2, !c(V,Y)) )),type="response")
  mu1[s] <- predict(Out_mu2,new_data=as.data.frame(cbind(V=1,C = select(TNDdf_train1, !c(V,Y)) )),type="response")               
  
  mu0[-s] <- predict(Out_mu1,new_data=as.data.frame(cbind(V=0,C = select(TNDdf_train2, !c(V,Y)) )),type="response")
  mu0[s] <- predict(Out_mu2,new_data=as.data.frame(cbind(V=0,C = select(TNDdf_train1, !c(V,Y)) )),type="response") 
  
  # m0
  Out_m1 <- fit_hal(X = select(TNDdf_train1, !c(V,Y)),
                    Y = TNDdf_train1$Y,
                    family = "binomial"
  )
  Out_m2 <- fit_hal(X = select(TNDdf_train2, !c(V,Y)),
                    Y = TNDdf_train2$Y,
                    family = "binomial"
  )
  m0 <- TNDdf_train$Y
  m0[-s] <- 1- predict(Out_m1,new_data = select(TNDdf_train2, !c(V,Y)), type="response")
  m0[s] <- 1- predict(Out_m2,new_data = select(TNDdf_train1, !c(V,Y)), type="response")
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 =g1_cont, g0 = 1- g1_cont, w1 = m0/(1 - mu1), w0 = m0/(1 - mu0), w1.ps = g1_full/g1_cont, w0.ps = (1 - g1_full)/(1 - g1_cont) ))
}


#####################################

library(sandwich)
res <- element(TNDdat)
mod_IPW <- function(TNDdat, res){
  TNDdat$ipw <- ifelse(TNDdat$V == 1, 1/res$g1, 1/res$g0)
  modY.ipw <- glm(Y ~ V, family=binomial(link = "log"), weights = ipw, data=TNDdat)
  est.ipw <- exp(modY.ipw$coefficients[2])
  se.ipw <- sqrt(vcovHC(modY.ipw)[2,2])
  
  CI_l <- est.ipw *exp(- 1.96 * se.ipw )
  CI_u <- est.ipw *exp( 1.96 * se.ipw )
  return(list(est = est.ipw, se = se.ipw, CI =  c(CI_l, CI_u)))
}

mod_IPW(TNDdat, res)

mod_EIF1OUT <- function(TNDdat, res){
  # proposed eif estimator 1 with Out ratio weights
  A.1 <- res$w1*res$mu1*((1 - TNDdat$Y)*(TNDdat$V - res$g1))/(res$g1* res$m0)
  A.0 <- res$w0*res$mu0*((1 - TNDdat$Y)*((1-TNDdat$V) - res$g0))/(res$g0* res$m0)
  psi.1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - A.1)
  psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - A.0)
  mod_eif <- psi.1/psi.0
  eifln <-  ((TNDdat$Y*TNDdat$V/res$g1 - A.1)/psi.1) - ((TNDdat$Y*(1-TNDdat$V)/res$g0 - A.0)/ psi.0)
  varln <-  var(eifln)/nrow(TNDdat)
  
  CI_l1 <- exp(log(mod_eif) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif) + 1.96 * sqrt(varln) )
  
  eifpsi <- (TNDdat$Y*TNDdat$V/res$g1 - A.1)/psi.0 - (psi.1/psi.0)*(TNDdat$Y*(1-TNDdat$V)/res$g0 - A.0)/psi.0
  var <- var(eifpsi)/nrow(TNDdat)
  CI_l2 <- mod_eif - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif + 1.96 * sqrt(var) 
  return(list( est = mod_eif, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}

mod_EIF1OUT(TNDdat, res)

mod_EIF1PS <- function(TNDdat, res){
  # proposed eif estimator 1 with PS ratio weights
  A.1 <- res$w1.ps*res$mu1*((1 - TNDdat$Y)*(TNDdat$V - res$g1))/(res$g1* res$m0)
  A.0 <- res$w0.ps*res$mu0*((1 - TNDdat$Y)*((1-TNDdat$V) - res$g0))/(res$g0* res$m0)
  psi.1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - A.1)
  psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - A.0)
  mod_eif <- psi.1/psi.0
  eifln <-  ((TNDdat$Y*TNDdat$V/res$g1 - A.1)/psi.1) - ((TNDdat$Y*(1-TNDdat$V)/res$g0 - A.0)/ psi.0)
  varln <-  var(eifln)/nrow(TNDdat)
  
  CI_l1 <- exp(log(mod_eif) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif) + 1.96 * sqrt(varln) )
  
  eifpsi <- (TNDdat$Y*TNDdat$V/res$g1 - A.1)/psi.0 - (psi.1/psi.0)*(TNDdat$Y*(1-TNDdat$V)/res$g0 - A.0)/psi.0
  var <- var(eifpsi)/nrow(TNDdat)
  CI_l2 <- mod_eif - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif + 1.96 * sqrt(var) 
  return(list( est = mod_eif, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}

mod_EIF1PS(TNDdat, res)

mod_EIF2 <- function(TNDdat, res){
  # proposed eif estimator 2: SAME as the previous method, i.e.,eif
  A.1 <- ((1 - TNDdat$Y)*(TNDdat$V - res$g1))/(res$g1* (1 -res$mu1))
  A.0 <- ((1 - TNDdat$Y)*((1-TNDdat$V) - res$g0))/(res$g0* (1 - res$mu0))
  psi.1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)
  psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)
  mod_eif2 <- psi.1/psi.0
  eifln <-  ((TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)/psi.1) - ((TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)/ psi.0)
  varln <-  var(eifln)/nrow(TNDdat)
  
  CI_l1 <- exp(log(mod_eif2) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif2) + 1.96 * sqrt(varln) )
  
  eifpsi <- (TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)/psi.0 - (psi.1/psi.0)*(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)/psi.0
  var <- var(eifpsi)/nrow(TNDdat)
  CI_l2 <- mod_eif2 - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif2 + 1.96 * sqrt(var) 
  return(list( est = mod_eif2, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}

mod_EIF2(TNDdat, res)

mod_OutReg<- function(TNDdat, res){
  mod_OR1 <- mean(res$mu1 * res$w1)/mean(res$mu0 * res$w0)
  mod_OR2 <- mean(res$mu1 * res$w1.ps)/mean(res$mu0 * res$w0.ps)
  return(list(est1 = mod_OR1, est2 = mod_OR2 ))
}

mod_OutReg(TNDdat, res)
library(hal9001)
for (i in 1:1000){
  TNDdat = datagen(ssize=5000, em=0.5)
  res <- element(TNDdat)
  #######################################################
  #IPW with controls to fit g
  est.1<-mod_IPW(TNDdat, res)
  est1<- est.1$est
  se1 <- est.1$se
  CI1<- est.1$CI
  #######################################################
  # Outcome regression with propensity ratio
  est2<-mod_OutReg(TNDdat, res)$est1
  est3<-mod_OutReg(TNDdat, res)$est2
  #######################################################
  ## Method 3: The proposed estimator that is based on the EIF
  # eif with Out-ratio
  est.4 <- mod_EIF1OUT(TNDdat, res)
  est4<- est.4$est
  OUTvarln <- est.4$varln
  OUTvar <- est.4$var
  CI4<- est.4$CI1
  # eif with PS-ratio
  est.5 <- mod_EIF1PS(TNDdat, res)
  est5<- est.5$est
  PSvarln <- est.5$varln
  PSvar <- est.5$var
  CI5<- est.5$CI1
  # eif with PS-ratio
  est.6 <- mod_EIF2(TNDdat, res)
  est6<- est.6$est
  varln <- est.6$varln
  var <- est.6$var
  CI6<- est.6$CI1
  CI7<- est.6$CI2
  write(c(i,est1,CI1,est2,est3, est4, CI4,est5, CI5, est6, CI6, CI7, se1, OUTvarln, OUTvar, PSvarln, PSvar, varln, var),file="FULLhal5000_em05.txt",ncolumns=50,append=T)
}

res1<-read.table("FULLhal5000_em05.txt",header=F)
res1 <- na.omit(res1)
head(res1)
dim(res1)
colnames(res1)[c(1, 2, 5, 6, 7, 10, 13, 18, 23,24)] <- c("i", "IPW", "out-Out", "out-PS", "EIF-out", "EIF-ps","EIF", "Se-ipw","Varln-EIF","Var-EIF")
summary(res1, na.rm = FALSE)
mRR = 0.166
mRR = 0.2837254
apply(res1, 2, mean) - mRR
apply(res1, 2, median) - mRR
apply(res1, 2, mean)
apply(res1, 2, median)
apply(res1, 2, sd)


#CIs
psi = mRR
mean(psi<=res1$V4 & psi>=res1$V3) #
mean(psi<=res1$V9 & psi>=res1$V8) #
mean(psi<=res1$V12 & psi>=res1$V11)
mean(psi<=res1$V15 & psi>=res1$V14)
mean(psi<=res1$V17 & psi>=res1$V16)

res1<-read.table("/Users/congjiang/Downloads/Fullhal5000_em2cross.txt",header=F)
res1 <- na.omit(res1)
head(res1)
dim(res1)
colnames(res1)[c(1, 2, 5, 6, 7, 10, 13, 18, 23,24)] <- c("i", "IPW", "out-Out", "out-PS", "EIF-out", "EIF-ps","EIF", "Var-ipw","Varln-EIF","Var-EIF")
summary(res1, na.rm = FALSE)
mRR = 0.263
apply(res1, 2, mean) - mRR
apply(res1, 2, median) - mRR
apply(res1, 2, mean)
apply(res1, 2, sd)


#CIs
psi = mRR
mean(psi<=res1$V4 & psi>=res1$V3) #
mean(psi<=res1$V9 & psi>=res1$V8) #
mean(psi<=res1$V12 & psi>=res1$V11)
mean(psi<=res1$V15 & psi>=res1$V14)
mean(psi<=res1$V17 & psi>=res1$V16)

res1<-read.table("/Users/congjiang/Downloads/Fullhal5000_em3cross.txt",header=F)
res1 <- na.omit(res1)
head(res1)
dim(res1)
colnames(res1)[c(1, 2, 5, 6, 7, 10, 13, 18, 23,24)] <- c("i", "IPW", "out-Out", "out-PS", "EIF-out", "EIF-ps","EIF", "Var-ipw","Varln-EIF","Var-EIF")
summary(res1, na.rm = FALSE)
mRR = 0.415
apply(res1, 2, mean) - mRR
apply(res1, 2, median) - mRR
apply(res1, 2, mean)
apply(res1, 2, sd)


#CIs
psi = mRR
mean(psi<=res1$V4 & psi>=res1$V3) #
mean(psi<=res1$V9 & psi>=res1$V8) #
mean(psi<=res1$V12 & psi>=res1$V11)
mean(psi<=res1$V15 & psi>=res1$V14)
mean(psi<=res1$V17 & psi>=res1$V16)
