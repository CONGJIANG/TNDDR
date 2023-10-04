set.seed(202308)
datagen<-function(seed=sample(1:1000000,size=1),ssize=500,popsize=1500000,OR_C=1.5,OR_WI=1,OR_WC=3.5,OR_H=1.5,em=0.45,cfV0=F,cfV1=F,return_full=F){
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

N <- 5000
nsim <- 100
TndDAT <- datagen(ssize = N, em = 0.45)
sum(TndDAT$Y)/N
# Define the sample size
n <- seq(300, N, by = 200)
# Define an empty matrix to store the estimated means
estIPW <- matrix(NA, nsim, length(n))
estOR <- matrix(NA, nsim, length(n))
estEIF <- matrix(NA, nsim, length(n))
for (j in 1:nsim) {
  TndDAT <- datagen(ssize = N, em = 0.45)
  for (i in 1:length(n)) {
    TNDdat <- TndDAT[1:n[i], ]
    res <- element(TNDdat)
    est.1<-mod_IPW(TNDdat, res)
    estIPW[j,i] <- est.1$est
    estOR[j,i] <-mod_OutReg(TNDdat, res)$est1
    est.2 <- mod_EIF2(TNDdat, res)
    estEIF[j,i] <- est.2$est
  }
}

mRR.em1 = 0.2646 # em = 0.45 
mRR.em1 = 0.300 # em = 0.55 

# Compute the biases
bias.ipw <- (apply(na.omit(estIPW), 2, mean) - mRR.em1) 
bias.or <- (apply(na.omit(estOR), 2, mean) - mRR.em1) 
bias.eif <- (apply(na.omit(estEIF), 2, mean) - mRR.em1)

res <- as.data.frame(cbind(n = n, bias.ipw = bias.ipw,bias.or = bias.or, bias.eif = bias.eif))
res <- na.omit(res[2:25,])
# Plot the biases as a function of sample size
par(mar = c(5, 5, 4, 2) + 0.1)  # Adjust the margin settings
plot(res$n, res$bias.or, type = "b", pch = 4, col = "blue",
     xlab = "Sample size (n)", ylab = expression(paste(hat(psi)[mRR]  - psi[mRR])), ylim = c(-0.2, 0.2))
lines(res$n, res$bias.eif, type = "b",pch = 19, col = "red")
lines(res$n, res$bias.ipw, type = "b", pch = 2, col = "black")
lines(res$n, 1/res$n^(1/2), pch = 0,  col = "brown")
lines(res$n, -1/sqrt(res$n), pch = 0,  col = "brown")
abline(h = 0, lty = 2)
legend("bottomright", legend = c("IPW","EIF", "OutReg", expression(1/sqrt(n))), pch = c(2, 19, 4, NA), col = c("black","red", "blue", "brown"), lty = 1)
title(expression(paste("Root-n Consistency of Proposed Estimators (", beta[em], "= 0.45)")))

write.csv(res,file="TNDStudy1halrootnV3em450831.csv")
res <- read.csv("/Users/congjiang/Downloads/TNDStudyRootN450209.txt")
dim(res)



