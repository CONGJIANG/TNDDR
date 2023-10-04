## @ Sara : 
## Il y a deux endroits ou j'ai indique que tu devais faire quelque chose de specifique,
## L'une est pour indiquer le chemin d'acces aux donnees, l'autre pour les variables
## categorielles (mais le 2e est juste pour valider)
## Autrement, tu devrais pouvoir rouler tout le code d'un coup

## Preparation
#Before running this VE analyse code,  the TND dataset should be well prepared: 
#Outcome is binary (SARS-Cov-2 infection or not), Covid19 Vaccination status is binary, e.g., unvaccinated (V=0) vs fully vaccinated + 14 days (V=1)
#Covariates that we will use to adjust is selected, such as ages group (e.g., < 18, 18~60, or >=60), gender, race, 
#and the calendar month of a test-positive subject's first positive COVID test and a test-negative subject's last COVID test, etc. 

# need the following packages:
require("dplyr") # to use select function
require("sandwich")
require("haven")
require("sas7bdat")
library(earth)


# Load a SAS dataset
## @Sara : Ci-dessous, tu dois indiquer le chemin d'acces aux donnees
TNDdat <- read.sas7bdat("path/to/TND_sas_dataset.sas7bdat")
## TNDdat <- read.sas7bdat("C:\\Users\\Denis Talbot\\Dropbox\\Travail\\Supervision\\CongJiang\\exemple_r.sas7bdat")
## TNDdat = rbind(TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat, TNDdat)

# head(TNDdat)
# convert all categorical variables to factors
## @Sara : Ci-dessous, il faut ecrire la liste des noms des variables categorielles
##         j'ai mis sex et multimorb meme si elles sont binaires, mais ce n'est pas ncessaire
factors <- c("sex", "epiweek_test", "multimorb")
# factors <- c("sex", "quintmat", "quitsoc")
TNDdat[factors] <- lapply(TNDdat[factors], as.factor)
baselinevars <- names(dplyr::select(TNDdat, !c(V,Y)))
form = as.formula(paste0("TNDdat$Y ~ -1 + ", paste(baselinevars, collapse = " + "), sep = ""));

# Rename variables (accordingly)
TNDdat2 = data.frame(TNDdat$Y, TNDdat$V, model.matrix(form, data = TNDdat));
names(TNDdat2) = c("Y", "V", paste0("X", 1:(ncol(TNDdat2)-2)));
TNDdat = TNDdat2;

################### Creation des fonctions necessaires pour l'estimation


######## fit nuisance models: HAL method w/ data splitting
# TNDdf_train = TNDdat;
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
  g1_cont[-s]<-predict(mod_g1_ctr,type="response",newdata=as.data.frame(cbind(select(TNDdf_train2, !c(V,Y)),V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y)))
  g1_cont[s]<-predict(mod_g2_ctr,type="response",newdata=as.data.frame(cbind(select(TNDdf_train1, !c(V,Y)),V=rep(1, nrow(TNDdf_train1)),Y=TNDdf_train1$Y)))
  
  
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
  mu1[-s] <- predict(Out_mu1,newdata=as.data.frame(cbind(V=1, select(TNDdf_train2, !c(V,Y)) )),type="response")
  mu1[s] <- predict(Out_mu2,newdata=as.data.frame(cbind(V=1, select(TNDdf_train1, !c(V,Y)) )),type="response")               
  
  mu0[-s] <- predict(Out_mu1,newdata=as.data.frame(cbind(V=0, select(TNDdf_train2, !c(V,Y)) )),type="response")
  mu0[s] <- predict(Out_mu2,newdata=as.data.frame(cbind(V=0, select(TNDdf_train1, !c(V,Y)) )),type="response") 
  
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
######## Methods for VE
### IPW estimator
mod_IPW <- function(TNDdat, res){
  TNDdat$ipw <- ifelse(TNDdat$V == 1, 1/res$g1, 1/res$g0)
  modY.ipw <- glm(Y ~ V, family=binomial(link = "log"), weights = ipw, data=TNDdat)
  est.ipw <- exp(modY.ipw$coefficients[2])
  se.ipw <- sqrt(vcovHC(modY.ipw)[2,2])
  
  CI_l <- est.ipw *exp(- 1.96 * se.ipw )
  CI_u <- est.ipw *exp( 1.96 * se.ipw )
  return(list(est = est.ipw, se = se.ipw, CI =  c(CI_l, CI_u)))
}
### EIF based estimator with OUTCOME ratio debiasing weights
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
### EIF based estimator with PS ratio debiasing weights
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

#### proposed TNDDR, SAME as eif1OUT
mod_EIF2 <- function(TNDdat, res){
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

mod_OutReg<- function(TNDdat, res){
  mod_OR1 <- mean(res$mu1 * res$w1)/mean(res$mu0 * res$w0)
  mod_OR2 <- mean(res$mu1 * res$w1.ps)/mean(res$mu0 * res$w0.ps)
  
  nbs <- 50
  bsest<-rep(NA,nbs)
  for(i in 1:nbs){
    resamps<-sample(1:nrow(TNDdat),size=nrow(TNDdat),replace=T)
    datk<-TNDdat[resamps,]
    res <- element.earth(datk)
    bsest[i] <- mean(res$mu1 * res$w1)/mean(res$mu0 * res$w0)
  }
  bs_var <- var(bsest)
  CI = quantile(bsest,c(0.025,0.975))
  return(list(est1 = mod_OR1,CI.OR1 = CI, est2 = mod_OR2 ))
}

################### FIN : Creation des fonctions necessaires pour l'estimation



############## Estimation avec les differentes methodes (on obtient des RR)
a = Sys.time();
res <- element.earth(TNDdat)
b = Sys.time();
b - a;

summary(res$mu0);
summary(res$mu1);
summary(res$m0);
summary(res$g1);
summary(res$g0);

mod_IPW(TNDdat, res)
mod_EIF1OUT(TNDdat, res)
mod_EIF1PS(TNDdat, res)
mod_EIF2(TNDdat, res)
mod_OutReg(TNDdat, res)
