# Install and load required packages
library(sandwich)
library(randomForest)
library(nnet)
library(dplyr)
library(ranger)
library(hal9001)
library(earth)

# Define function using Random Forest to predict probabilities
element.randomForest <- function(TNDdf_train) {
  # Ensure the variables V and Y are treated as factors for classification
  TNDdf_train$V <- as.factor(TNDdf_train$V)
  TNDdf_train$Y <- as.factor(TNDdf_train$Y)
  
  # Data splitting: split data into two equal parts randomly
  s <- sample(1:nrow(TNDdf_train), nrow(TNDdf_train) / 2)
  TNDdf_train1 <- TNDdf_train[s, ]
  TNDdf_train2 <- TNDdf_train[-s, ]
  
  # Subset for controls (Y == 0) in first split
  TNDdf_train_ctr1 <- subset(TNDdf_train1, Y == 0)
  
  mod_g1_ctr <- ranger(
    V ~ .,
    data = subset(TNDdf_train_ctr1, select = -Y),
    num.trees = 200,  # Set number of trees to 200
    mtry = 1,         # Set mtry to 1 for 2-3 covariates
    min.node.size = 60,  # Regularization
    sample.fraction = 0.33,
    splitrule = "extratrees",  # More randomness in splits
    probability = TRUE
  )
  
  # Subset for controls in second split
  TNDdf_train_ctr2 <- subset(TNDdf_train2, Y == 0)
  # Train Random Forest on second control subset using ranger
  mod_g2_ctr <- ranger(
    V ~ .,
    data = subset(TNDdf_train_ctr2, select = -Y),
    num.trees = 200,  # Set number of trees to 200
    mtry = 1,         # Set mtry to 1 for 2-3 covariates
    min.node.size = 60,  # Regularization
    sample.fraction = 0.33,
    splitrule = "extratrees",  # More randomness in splits
    probability = TRUE
  )
  
  # Initialize g1_cont to store probabilities instead of factor levels
  g1_cont <- rep(NA, nrow(TNDdf_train))  # Initialize as numeric vector
  
  # Ensure newdata has the same structure as training data
  newdata_train2 <- TNDdf_train2 %>%
    mutate(V = factor(rep(1, nrow(TNDdf_train2)), levels = levels(TNDdf_train$V)))
  
  g1_cont[-s] <- predict(mod_g1_ctr, data = newdata_train2)$predictions[, 2]
  
  newdata_train1 <- TNDdf_train1 %>%
    mutate(V = factor(rep(1, nrow(TNDdf_train1)), levels = levels(TNDdf_train$V)))
  
  g1_cont[s] <- predict(mod_g2_ctr, data = newdata_train1)$predictions[, 2]
  
  # Train Random Forest models for predicting Y (mu1 and mu0) using ranger
  Out_mu1 <- ranger(
    Y ~ .,
    data = TNDdf_train1,
    num.trees = 50,  
    mtry = 1,      
    probability = TRUE
  )
  
  Out_mu2 <- ranger(
    Y ~ .,
    data = TNDdf_train2,
    num.trees = 50,  
    mtry = 1,        
    probability = TRUE
  )
  
  # Initialize vectors to store mu1 and mu0 predictions
  mu1 <- rep(NA, nrow(TNDdf_train))
  mu0 <- rep(NA, nrow(TNDdf_train))
  
  # Predict mu1: the probability of Y = 1 when V = 1
  newdata_mu1_train2 <- TNDdf_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdf_train2)), levels = levels(TNDdf_train$V)))
  
  mu1[-s] <- predict(Out_mu1, data = newdata_mu1_train2)$predictions[, 2]
  
  newdata_mu1_train1 <- TNDdf_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdf_train1)), levels = levels(TNDdf_train$V)))
  
  mu1[s] <- predict(Out_mu2, data = newdata_mu1_train1)$predictions[, 2]
  
  # Predict mu0: the probability of Y = 1 when V = 0
  newdata_mu0_train2 <- TNDdf_train2 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdf_train2)), levels = levels(TNDdf_train$V)))
  
  mu0[-s] <- predict(Out_mu1, data = newdata_mu0_train2)$predictions[, 2]
  
  newdata_mu0_train1 <- TNDdf_train1 %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdf_train1)), levels = levels(TNDdf_train$V)))
  
  mu0[s] <- predict(Out_mu2, data = newdata_mu0_train1)$predictions[, 2]
  
  # Train Random Forest models for m0 (1 - Y, or probability of Y = 0) using ranger
  Out_m1 <- ranger(
    Y ~ .,
    data = subset(TNDdf_train1, select = -V),
    num.trees = 500,  # Set number of trees to 50
    mtry = 1,         # Set mtry to 1 for 2-3 covariates
    probability = TRUE
  )
  
  Out_m2 <- ranger(
    Y ~ .,
    data = subset(TNDdf_train2, select = -V),
    num.trees = 500,  # Set number of trees to 50
    mtry = 1,         # Set mtry to 1 for 2-3 covariates
    probability = TRUE
  )
  
  # Initialize m0 to store predictions
  m0 <- rep(NA, nrow(TNDdf_train))
  
  # Predict m0: probability of Y = 0
  m0[-s] <- 1 - predict(Out_m1, data = select(TNDdf_train2, -c(V, Y)))$predictions[, 2]
  m0[s] <- 1 - predict(Out_m2, data = select(TNDdf_train1, -c(V, Y)))$predictions[, 2]
  
  # Return all necessary elements
  mu1 <- pmin(pmax(mu1, 0.0000001), 0.9999999)
  mu0 <- pmin(pmax(mu0, 0.0000001), 0.9999999)
  m0 <- pmin(pmax(m0, 0.0000001), 0.9999999)
  g1 <- pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1, g0 = g0, w1 = m0 / (1 - mu1), w0 = m0 / (1 - mu0)))
}


# Define function using Neural Networks to predict probabilities
element.nnet <- function(TNDdf_train){
  
  # Ensure the variables V and Y are treated as factors for classification
  TNDdf_train$V <- as.factor(TNDdf_train$V)
  TNDdf_train$Y <- as.factor(TNDdf_train$Y)
  
  # Data splitting: split data into two equal parts randomly
  s <- sample(1:nrow(TNDdf_train), nrow(TNDdf_train) / 2)
  TNDdf_train1_scaled <- TNDdf_train[s, ]
  TNDdf_train2_scaled <- TNDdf_train[-s, ]
  
  # Subset for controls (Y == 0) in first split
  TNDdf_train_ctr1 <- subset(TNDdf_train1_scaled, Y == 0)
  # Train Neural Network on first control subset
  mod_g1_ctr <- nnet(
    V ~ .,
    data = subset(TNDdf_train_ctr1, select = -Y),
    size = 5, maxit = 50
  )
  
  # Subset for controls in second split
  TNDdf_train_ctr2 <- subset(TNDdf_train2_scaled, Y == 0)
  # Train Neural Network on second control subset
  mod_g2_ctr <- nnet(
    V ~ .,
    data = subset(TNDdf_train_ctr2, select = -Y),
    size = 5, maxit = 50
  )
  
  # Initialize g1_cont to store probabilities
  g1_cont <- rep(NA, nrow(TNDdf_train))
  
  # Predict probabilities for controls
  newdata_train2 <- TNDdf_train2_scaled %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdf_train2_scaled)), levels = levels(TNDdf_train$V)))
  
  g1_cont[-s] <- predict(mod_g1_ctr, newdata = newdata_train2, type = "raw")
  
  newdata_train1 <- TNDdf_train1_scaled %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdf_train1_scaled)), levels = levels(TNDdf_train$V)))
  
  g1_cont[s] <- predict(mod_g2_ctr, newdata = newdata_train1, type = "raw")
  
  # Train Neural Network models for predicting Y (mu1 and mu0)
  Out_mu1 <- nnet(
    Y ~ .,
    data = TNDdf_train1_scaled,
    size = 5, maxit = 50
  )
  
  Out_mu2 <- nnet(
    Y ~ .,
    data = TNDdf_train2_scaled,
    size = 5, maxit = 50
  )
  
  # Initialize vectors to store mu1 and mu0 predictions
  mu1 <- rep(NA, nrow(TNDdf_train))
  mu0 <- rep(NA, nrow(TNDdf_train))
  
  # Predict mu1: the probability of Y = 1 when V = 1
  newdata_mu1_train2 <- TNDdf_train2_scaled %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdf_train2_scaled)), levels = levels(TNDdf_train$V)))
  
  mu1[-s] <- predict(Out_mu1, newdata = newdata_mu1_train2, type = "raw")
  
  newdata_mu1_train1 <- TNDdf_train1_scaled %>%
    select(-Y) %>%
    mutate(V = factor(rep(1, nrow(TNDdf_train1_scaled)), levels = levels(TNDdf_train$V)))
  
  mu1[s] <- predict(Out_mu2, newdata = newdata_mu1_train1, type = "raw")
  
  # Predict mu0: the probability of Y = 1 when V = 0
  newdata_mu0_train2 <- TNDdf_train2_scaled %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdf_train2_scaled)), levels = levels(TNDdf_train$V)))
  
  mu0[-s] <- predict(Out_mu1, newdata = newdata_mu0_train2, type = "raw")
  
  newdata_mu0_train1 <- TNDdf_train1_scaled %>%
    select(-Y) %>%
    mutate(V = factor(rep(0, nrow(TNDdf_train1_scaled)), levels = levels(TNDdf_train$V)))
  
  mu0[s] <- predict(Out_mu2, newdata = newdata_mu0_train1, type = "raw")
  
  # Train Neural Network models for m0 (1 - Y, or probability of Y = 0)
  Out_m1 <- nnet(
    Y ~ .,
    data = subset(TNDdf_train1_scaled, select = -V),
    size = 5, maxit = 50
  )
  
  Out_m2 <- nnet(
    Y ~ .,
    data = subset(TNDdf_train2_scaled, select = -V),
    size = 5, maxit = 50
  )
  
  # Initialize m0 to store predictions
  m0 <- rep(NA, nrow(TNDdf_train))
  
  # Predict m0: probability of Y = 0
  m0[-s] <- 1 - predict(Out_m1, newdata = select(TNDdf_train2_scaled, -c(V, Y)), type = "raw")
  m0[s] <- 1 - predict(Out_m2, newdata = select(TNDdf_train1_scaled, -c(V, Y)), type = "raw")
  
  # Return all necessary elements
  mu1 <- pmin(pmax(mu1, 0.0000001), 0.9999999)
  mu0 <- pmin(pmax(mu0, 0.0000001), 0.9999999)
  m0 <- pmin(pmax(m0, 0.0000001), 0.9999999)
  g1 <- pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1,g0 = g0, w1 = m0 / (1 - mu1),w0 = m0 / (1 - mu0)))
}




######## fit nuisance models: HAL method w/ data splitting
element.HAL <- function(TNDdf_train){
  # Data splitting
  s <- sample(1:length(TNDdf_train$Y), length(TNDdf_train$Y)/2)
  TNDdf_train1 <- TNDdf_train[s,]
  TNDdf_train2 <- TNDdf_train[-s,]
  
  TNDdf_train_ctr1 <- subset(TNDdf_train1, Y==0)
  #training
  mod_g1_ctr <- fit_hal(
    X = select(TNDdf_train_ctr1, !c(V,Y)),
    Y = TNDdf_train_ctr1$V,
    family = "binomial"
  )
  
  TNDdf_train_ctr2 <- subset(TNDdf_train2, Y==0)
  #training
  mod_g2_ctr <- fit_hal(
    X = select(TNDdf_train_ctr2, !c(V,Y)),
    Y = TNDdf_train_ctr2$V,
    family = "binomial"
  )
  g1_cont <- TNDdf_train$V
  g1_cont[-s]<-predict(mod_g1_ctr,type="response",new_data=as.data.frame(cbind(C=TNDdf_train2$C,V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y )))
  g1_cont[s]<-predict(mod_g2_ctr,type="response",new_data=as.data.frame(cbind(C=TNDdf_train1$C,V=rep(1, nrow(TNDdf_train1)),Y=TNDdf_train1$Y) ))
  
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
  g1_full[-s]<-predict(mod_g1_full,type="response",new_data=as.data.frame(cbind(C=TNDdf_train2$C,V=rep(1, nrow(TNDdf_train2)),Y=TNDdf_train2$Y )))
  g1_full[s]<-predict(mod_g2_full,type="response",new_data=as.data.frame(cbind(C=TNDdf_train1$C,V=rep(1, nrow(TNDdf_train1)),Y=TNDdf_train1$Y) ))
  
  
  #mu
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
  mu1 <- pmin(pmax(mu1, 0.0000001), 0.9999999)
  mu0 <- pmin(pmax(mu0, 0.0000001), 0.9999999)
  m0 <- pmin(pmax(m0, 0.0000001), 0.9999999)
  g1 <- pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1,g0 = g0, w1 = m0 / (1 - mu1),w0 = m0 / (1 - mu0)))
}



element.earth <- function(TNDdf_train){
  # Data splitting
  s <- sample(1:length(TNDdf_train$Y), length(TNDdf_train$Y)/2)
  TNDdf_train1 <- TNDdf_train[s,]
  TNDdf_train2 <- TNDdf_train[-s,]
  
  TNDdf_train_ctr1 <- subset(TNDdf_train1, Y==0)
  #training
  mod_g1_ctr <- earth(
    V ~ .,
    data = subset(TNDdf_train_ctr1, select = -Y),
    glm=list(family=binomial)
  )
  
  TNDdf_train_ctr2 <- subset(TNDdf_train2, Y==0)
  #training
  mod_g2_ctr <- earth(
    V ~ .,
    data = subset(TNDdf_train_ctr2, select = -Y),
    glm=list(family=binomial)
  )
  g1_cont <- TNDdf_train$V
  g1_cont[-s]<-predict(mod_g1_ctr,type="response",newdata=as.data.frame(cbind(select(TNDdf_train2, !c(V,Y)),V=rep(1, nrow(TNDdf_train2) ),Y=TNDdf_train2$Y)))
  g1_cont[s]<-predict(mod_g2_ctr,type="response",newdata=as.data.frame(cbind(select(TNDdf_train1, !c(V,Y)),V=rep(1, nrow(TNDdf_train1)) ,Y=TNDdf_train1$Y)))
  
  
  #mu
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
  
  mu1 <- pmin(pmax(mu1, 0.0000001), 0.9999999)
  mu0 <- pmin(pmax(mu0, 0.0000001), 0.9999999)
  m0 <- pmin(pmax(m0, 0.0000001), 0.9999999)
  g1 <- pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.0000001), 0.9999999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1,g0 = g0, w1 = m0 / (1 - mu1),w0 = m0 / (1 - mu0)))
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

mod_IPW1 <- function(TNDdat, res, bootstrap_CI = TRUE, method = "RandomForest"){
  IPW.est <- mean(TNDdat$Y*TNDdat$V/res$g1)/mean(TNDdat$Y*(1-TNDdat$V)/res$g0)
  if (bootstrap_CI) {
    nbs <- 50
    bsest <- rep(NA, nbs)
    
    for(i in 1:nbs){
      resamps <- sample(1:nrow(TNDdat), size = nrow(TNDdat), replace = TRUE)
      datk <- TNDdat[resamps,]
      
      # Fit model based on the chosen method
      if (method == "HAL") {
        res <- element.HAL(datk)
      } else if (method == "Earth") {
        res <- element.earth(datk)
      } else if (method == "RandomForest") {
        res <- element.randomForest(datk)
      } else if (method == "NN") {
        res <- element.nnet(datk)
      } else {
        stop("Unsupported method. Please choose from 'HAL', 'Earth', 'RandomForest', or 'NN'.")
      }
      
      bsest[i] <- mean(datk$Y*datk$V/res$g1)/mean(datk$Y*(1-datk$V)/res$g0)
    }
    
    bs_var <- var(bsest)
    CI <- quantile(bsest, c(0.025, 0.975), na.rm = TRUE)
  } else {
    CI <- NA
  }
  return(list(est = IPW.est, CI = CI))
}

#### proposed TNDDR
mod_TNDDR <- function(TNDdat, res){
  A.1 <- ((1 - TNDdat$Y)*(TNDdat$V - res$g1))/(res$g1* (1 -res$mu1))
  A.0 <- ((1 - TNDdat$Y)*((1-TNDdat$V) - res$g0))/(res$g0* (1 - res$mu0))
  psi.1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)
  psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)
  mod_eif2 <- pmin(pmax(psi.1/psi.0, 0.001), 0.999)
  
  eifln <-  ((TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1 - psi.1)/psi.1) - ((TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0 - psi.0)/ psi.0)
  varln <-  var(eifln)/nrow(TNDdat)
  
  CI_l1 <- exp(log(mod_eif2) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif2) + 1.96 * sqrt(varln) )
  
  eifpsi <- (TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1 - psi.1)/psi.0 - mod_eif2*(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0 - psi.0)/psi.0
  var <- var(eifpsi)/nrow(TNDdat)
  CI_l2 <- mod_eif2 - 1.96 * sqrt(var)
  CI_u2 <- mod_eif2 + 1.96 * sqrt(var)
  
  varn2 <- mean((mod_eif2* TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1 - TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)^2)
  denJ <- psi.0^2
  var2 <- varn2/(denJ * nrow(TNDdat))
  CI_l3 <- mod_eif2 - 1.96 * sqrt(var2)
  CI_u3 <- mod_eif2 + 1.96 * sqrt(var2)
  return(list( est = mod_eif2, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2), CI3 =  c(CI_l3, CI_u3)))
}

#TNDdat <- datagen(ssize = 10000, em = -0.25)
#res <- fit_methods[["Earth"]](TNDdat)
#mod_EIF2(TNDdat, res)

mod_OutReg <- function(TNDdat, res, bootstrap_CI = TRUE, method = "RandomForest"){
  mod_OR1 <- mean(res$mu1 * res$w1) / mean(res$mu0 * res$w0)
  if (bootstrap_CI) {
    nbs <- 50
    bsest <- rep(NA, nbs)
    
    for(i in 1:nbs){
      resamps <- sample(1:nrow(TNDdat), size = nrow(TNDdat), replace = TRUE)
      datk <- TNDdat[resamps,]
      
      # Fit model based on the chosen method
      if (method == "HAL") {
        res <- element.HAL(datk)
      } else if (method == "Earth") {
        res <- element.earth(datk)
      } else if (method == "RandomForest") {
        res <- element.randomForest(datk)
      } else if (method == "NN") {
        res <- element.nnet(datk)
      } else {
        stop("Unsupported method. Please choose from 'HAL', 'Earth', 'RandomForest', or 'NN'.")
      }
      
      bsest[i] <- mean(res$mu1 * res$w1) / mean(res$mu0 * res$w0)
    }
    
    bs_var <- var(bsest)
    CI <- quantile(bsest, c(0.025, 0.975), na.rm = TRUE)
  } else {
    CI <- NA
  }
  return(list(est = mod_OR1, CI.OR = CI))
}

