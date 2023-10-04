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

# Load a SAS dataset
## @Sara : Ci-dessous, tu dois indiquer le chemin d'acces aux donnees
TNDdat <- read.sas7bdat("path/to/TND_sas_dataset.sas7bdat")
## TNDdat <- read.sas7bdat("C:\\Users\\Denis Talbot\\Dropbox\\Travail\\Supervision\\CongJiang\\exemple_r.sas7bdat")

# convert all categorical variables to factors
## @Sara : Ci-dessous, il faut ecrire la liste des noms des variables categorielles
##         j'ai mis sex et multimorb meme si elles sont binaires, mais ce n'est pas ncessaire
factors <- c("sex", "epiweek_test", "multimorb")
TNDdat[factors] <- lapply(TNDdat[factors], as.factor)
# Rename variables (accordingly)

# isolate the names of baseline covariates
baselinevars <- names(dplyr::select(TNDdat, !c(V,Y)))
baselinevars


################### Creation des fonctions necessaires pour l'estimation


######## fit nuisance models:
element <- function(TNDdat){
  baselinevars <- names(dplyr::select(TNDdat, !c(V,Y)))
  ######## fit models
  # adjust the exposure variable (primary interest) + covariates
  # fit mu model
  mu.formula <- as.formula(paste("Y~ V +", 
                                 paste(baselinevars, 
                                       collapse = "+")))
  
  mu.fit <- glm(mu.formula,family="binomial", data=TNDdat)
  mu1 <- predict(mu.fit,newdata=as.data.frame(cbind(V=1,select(TNDdat, !c(V,Y)))),type="response")
  mu0 <- predict(mu.fit,newdata=as.data.frame(cbind(V=0,select(TNDdat, !c(V,Y)))),type="response")
  # fit m model
  m.formula <- as.formula(paste("Y~", 
                                paste(baselinevars, 
                                      collapse = "+")))
  m.fit <- glm(mu.formula,family="binomial", data=TNDdat)
  m0 <- 1 - predict(m.fit,type="response")
  # debiasing weights
  #w1 <-m0/(1 - mu1)
  #w0 <-m0/(1 - mu0)
  # fit propensity models using control data
  ps.formula <- as.formula(paste("V ~",
                                 paste(baselinevars,
                                       collapse = "+")))
  ps.fit <- glm(ps.formula,family="binomial", data=TNDdat, subset=(TNDdat$Y==0))
  g1_cont<-predict(ps.fit,type="response",newdata= TNDdat);
  g0_cont <- 1 - g1_cont
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 =g1_cont, g0 = g0_cont, w1 = m0/(1 - mu1), w0 = m0/(1 - mu0), mu.fit = mu.fit))
}

######## Methods for VE
# 1. IPW estimate and SE
library(sandwich)
estIPW <- function(TNDdat, res){
  TNDdat$ipw <- ifelse(TNDdat$V == 1, 1/res$g1, 1/res$g0)
  # modY.ipw <- glm(Y ~ V, family=binomial(link = "log"), weights = ipw, data=TNDdat)
  # if does not work try:
  modY.ipw <- glm(Y ~ V, family=binomial(link = "log"), weights = ipw, data=TNDdat)
  est.ipw <- exp(modY.ipw$coefficients[2])
  se.ipw <- sqrt(vcovHC(modY.ipw)[2,2])
  CI_l <- est.ipw * exp(- 1.96*se.ipw)
  CI_u <- est.ipw * exp( 1.96*se.ipw) 
  return(list( est.ipw = est.ipw, CI.ipw = c(CI_l, CI_u) ))
}

# 2. proposed eif estimator
estEIF1 <- function(TNDdat, res){
  A1 <- (1 - TNDdat$Y)*(TNDdat$V - res$g1)/(res$g1* res$m0)
  A0 <- (1 - TNDdat$Y)*((1-TNDdat$V) - res$g0)/(res$g0* res$m0)
  psi1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$w1* res$mu1*A1)
  psi0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$w0* res$mu0*A0)
  mod_eif <- psi1/psi0
  
  eifln <-  (TNDdat$Y*TNDdat$V/res$g1 - res$w1* res$mu1*A1)/psi1 - (TNDdat$Y*(1-TNDdat$V)/res$g0 - res$w0* res$mu0*A0)/ psi0
  varln <-  var(eifln)/nrow(TNDdat)
  CI_l1 <- mod_eif *exp(- 1.96 * sqrt(varln) )
  CI_u1 <- mod_eif *exp( 1.96 * sqrt(varln) )
  
  eifpsi <- (TNDdat$Y*TNDdat$V/res$g1 - res$w1* res$mu1*A1)/psi0 - (TNDdat$Y*(1-TNDdat$V)/res$g0 - res$w0* res$mu0*A0)/psi0
  var <- var(eifpsi)/nrow(TNDdat)
  CI_l2 <- mod_eif - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif + 1.96 * sqrt(var) 
  return(list( est = mod_eif, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}


# 3 proposed eif estimator 2: SAME as the previous method, i.e.,eif
estEIF2 <- function(TNDdat, res){
  A.1 <- (1 - TNDdat$Y)*(TNDdat$V - res$g1)/(res$g1* (1 -res$mu1))
  A.0 <- (1 - TNDdat$Y)*((1-TNDdat$V) - res$g0)/(res$g0* (1 - res$mu0))
  psi.1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)
  psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)
  mod_eif2 <- psi.1/psi.0
  eifln <-  (TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)/psi.1 - (TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)/ psi.0
  varln <-  var(eifln)/nrow(TNDdat)
  CI_l1 <- mod_eif2 *exp(- 1.96 * sqrt(varln) )
  CI_u1 <- mod_eif2 *exp( 1.96 * sqrt(varln) )
  eifpsi <- (TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)/psi.0 - (psi.1/psi.0)*(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)/psi.0
  var <- var(eifpsi)/nrow(TNDdat)
  CI_l2 <- mod_eif2 - 1.96 * sqrt(var) 
  CI_u2 <- mod_eif2 + 1.96 * sqrt(var) 
  return(list( est = mod_eif2, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2)))
}



# 4. Outcome regression estimator
estOR <- function(TNDdat, res){
  mod_OR <- mean(res$mu1 * res$w1)/mean(res$mu0 * res$w0)
  nbs <- 50
  bsest<-rep(NA,nbs)
  for(i in 1:nbs){
    resamps<-sample(1:nrow(TNDdat),size=nrow(TNDdat),replace=T)
    datk<-TNDdat[resamps,]
    res <- element(datk)
    bsest[i] <- mean(res$mu1 * res$w1)/mean(res$mu0 * res$w0)
  }
  bs_var <- var(bsest)
  CI = quantile(bsest,c(0.025,0.975))
  return( list(est.OR = mod_OR, CI.OR = CI) )
}

################### FIN : Creation des fonctions necessaires pour l'estimation



############## Estimation avec les differentes methodes (on obtient des RR)
a = Sys.time();
res <- element(TNDdat)
b = Sys.time();
b - a;
summary(res$mu0);
summary(res$mu1);
summary(res$m0);
summary(res$g1);
summary(res$g0);
summary(res$mu.fit);
estIPW(TNDdat, res);
estEIF1(TNDdat, res)
estEIF2(TNDdat, res)
estOR(TNDdat, res)



