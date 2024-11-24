# Load necessary libraries
library(dplyr)
library(randomForest)
library(nnet)
library(earth)
require("dplyr") # to use select function
require("sandwich")
require("haven")
require("sas7bdat")
require("hal9001")

source("TrueValue.R")

element_glm <- function(TNDdf_train, ps_model, out_model1, out_model2) {
  # Data splitting
  s <- sample(1:nrow(TNDdf_train), nrow(TNDdf_train) / 2)
  TNDdf_train1 <- TNDdf_train[s,]
  TNDdf_train2 <- TNDdf_train[-s,]
  
  TNDdf_train_ctr1 <- subset(TNDdf_train1, Y==0)
  TNDdf_train_ctr2 <- subset(TNDdf_train2, Y==0)
  # Training glm models for the treatment effect
  mod_g1_ctr <- glm(ps_model, data = subset(TNDdf_train_ctr1, select = -Y), family = binomial())
  mod_g2_ctr <- glm(ps_model, data = subset(TNDdf_train_ctr2, select = -Y), family = binomial())
  
  g1_cont <- TNDdf_train$V
  g1_cont[-s] <- predict(mod_g1_ctr, newdata = as.data.frame(cbind(select(TNDdf_train2, -c(V, Y)), V = rep(1, nrow(TNDdf_train2)), Y = TNDdf_train2$Y)), type = "response")
  g1_cont[s] <- predict(mod_g2_ctr, newdata = as.data.frame(cbind(select(TNDdf_train1, -c(V, Y)), V = rep(1, nrow(TNDdf_train1)), Y = TNDdf_train1$Y)), type = "response")
  
  # Training glm models for the outcome
  Out_mu1 <- glm(out_model1, data = TNDdf_train1, family = binomial())
  Out_mu2 <- glm(out_model1, data = TNDdf_train2, family = binomial())
  
  mu1 <- TNDdf_train$Y
  mu0 <- TNDdf_train$Y
  
  mu1[-s] <- predict(Out_mu1, newdata = as.data.frame(cbind(V = 1, select(TNDdf_train2, -c(V, Y)))), type = "response")
  mu1[s] <- predict(Out_mu2, newdata = as.data.frame(cbind(V = 1, select(TNDdf_train1, -c(V, Y)))), type = "response")
  
  mu0[-s] <- predict(Out_mu1, newdata = as.data.frame(cbind(V = 0, select(TNDdf_train2, -c(V, Y)))), type = "response")
  mu0[s] <- predict(Out_mu2, newdata = as.data.frame(cbind(V = 0, select(TNDdf_train1, -c(V, Y)))), type = "response")
  
  # Training glm models for m0
  Out_m1 <- glm(out_model2, data = subset(TNDdf_train1, select = -V), family = binomial)
  Out_m2 <- glm(out_model2, data = subset(TNDdf_train2, select = -V), family = binomial)
  
  m0 <- TNDdf_train$Y
  m0[-s] <- 1 - predict(Out_m1, newdata = select(TNDdf_train2, -c(V, Y)), type = "response")
  m0[s] <- 1 - predict(Out_m2, newdata = select(TNDdf_train1, -c(V, Y)), type = "response")
  
  summary(Out_m1)
  mu1 <- pmin(pmax(mu1, 0.001), 0.999)
  mu0 <- pmin(pmax(mu0, 0.001), 0.999)
  m0 <- pmin(pmax(m0, 0.001), 0.999)
  g1 <- pmin(pmax(g1_cont, 0.001), 0.999)
  g0 <- 1 - pmin(pmax(g1_cont, 0.001), 0.999)
  
  return(list(mu1 = mu1, mu0 = mu0, m0 = m0, g1 = g1,g0 = g0, w1 = m0 / (1 - mu1),w0 = m0 / (1 - mu0)))
}




######## Methods for VE
### IPW estimator
mod_IPW1GLM <- function(TNDdat, res){
  TNDdat$ipw <- ifelse(TNDdat$V == 1, 1/res$g1, 1/res$g0)
  modY.ipw <- glm(Y ~ V, family = poisson(link = "log"), weights = ipw, data = TNDdat)
  
  est.ipw <- exp(modY.ipw$coefficients[2])
  se.ipw <- sqrt(vcovHC(modY.ipw, type = "HC0")[2,2])
  
  CI_l <- est.ipw *exp(- 1.96 * se.ipw )
  CI_u <- est.ipw *exp( 1.96 * se.ipw )
  return(list(est = est.ipw, se = se.ipw, CI =  c(CI_l, CI_u)))
}

mod_IPW <- function(TNDdat, res, bootstrap_CI = TRUE, ps_model, out_model1, out_model2){
  IPW.est <- mean(TNDdat$Y*TNDdat$V/res$g1)/mean(TNDdat$Y*(1-TNDdat$V)/res$g0)
  if (bootstrap_CI) {
    nbs <- 200
    bsest <- rep(NA, nbs)
    
    for(i in 1:nbs){
      resamps <- sample(1:nrow(TNDdat), size = nrow(TNDdat), replace = TRUE)
      datk <- TNDdat[resamps,]
      res <- element_glm(datk, ps_model, out_model1, out_model2)
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
mod_EIF2 <- function(TNDdat, res, bootstrap_CI = TRUE){
  A.1 <- ((1 - TNDdat$Y)*(TNDdat$V - res$g1))/(res$g1* (1 -res$mu1))
  A.0 <- ((1 - TNDdat$Y)*((1-TNDdat$V) - res$g0))/(res$g0* (1 - res$mu0))
  psi.1 <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)
  psi.0 <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)
  mod_eif2 <- pmin(pmax(psi.1/psi.0, 0.0001), 0.9999)
  eifln <-  ((TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1 - psi.1)/psi.1) - ((TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0-psi.0)/ psi.0)
  varln <-  var(eifln)/nrow(TNDdat)
  
  CI_l1 <- exp(log(mod_eif2) - 1.96 * sqrt(varln) )
  CI_u1 <- exp(log(mod_eif2) + 1.96 * sqrt(varln) )
  
  eifpsi <- ((TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1 - psi.1)/psi.0) - ((psi.1/psi.0)*(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0 - psi.0)/psi.0)
  var <- var(eifpsi)/nrow(TNDdat)
  CI_l2 <- mod_eif2 - 1.96 * sqrt(var)
  CI_u2 <- mod_eif2 + 1.96 * sqrt(var)
  
  
  varn2 <- mean((mod_eif2* TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1 - TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)^2)
  denJ <- psi.0^2
  var2 <- varn2/(denJ * nrow(TNDdat))
  CI_l3 <- mod_eif2 - 1.96 * sqrt(var2)
  CI_u3 <- mod_eif2 + 1.96 * sqrt(var2)
  
  if (bootstrap_CI) {
    nbs <- 500
    bsest <- rep(NA, nbs)
    for(i in 1:nbs){
      resamps <- sample(1:nrow(TNDdat), size = nrow(TNDdat), replace = TRUE)
      datk <- TNDdat[resamps,]
      res <- element_glm(datk, ps_model, out_model1, out_model2)
      A.1[i] <- ((1 - TNDdat$Y)*(TNDdat$V - res$g1))/(res$g1* (1 -res$mu1))
      A.0[i] <- ((1 - TNDdat$Y)*((1-TNDdat$V) - res$g0))/(res$g0* (1 - res$mu0))
      psi.1[i] <- mean(TNDdat$Y*TNDdat$V/res$g1 - res$mu1*A.1)
      psi.0[i] <- mean(TNDdat$Y*(1-TNDdat$V)/res$g0 - res$mu0*A.0)
      bsest[i] <- pmin(pmax(psi.1[i]/psi.0[i], 0.001), 0.999)
    }
    bs_var <- var(bsest)
    CI.b <- quantile(bsest, c(0.025, 0.975), na.rm = TRUE)
  } else {
    CI.b <- NA
  }
  
  return(list( est = mod_eif2, varln = varln, var = var, CI1 =  c(CI_l1, CI_u1), CI2 =  c(CI_l2, CI_u2), CI3 =  c(CI_l3, CI_u3), CI.boots = CI.b))
}

mod_OutRegGLM <- function(TNDdat, res, bootstrap_CI = TRUE, ps_model, out_model1, out_model2){
  mod_OR1 <- mean(res$mu1 * res$w1) / mean(res$mu0 * res$w0)
  if (bootstrap_CI) {
    nbs <- 200
    bsest <- rep(NA, nbs)
    
    for(i in 1:nbs){
      resamps <- sample(1:nrow(TNDdat), size = nrow(TNDdat), replace = TRUE)
      datk <- TNDdat[resamps,]
      res <- element_glm(datk, ps_model, out_model1, out_model2)
      bsest[i] <- mean(res$mu1 * res$w1) / mean(res$mu0 * res$w0)
    }
    bs_var <- var(bsest)
    CI <- quantile(bsest, c(0.025, 0.975), na.rm = TRUE)
  } else {
    CI <- NA
  }
  return(list(est = mod_OR1, CI.OR = CI))
}


output_dir <- ""


TNDresGLM <- function(B = 1000, ssize = 1000, em = -0.15, bootstrap_CI = FALSE, ps_model, out_model1, out_model2) {
  # Construct the output file name based on the method
  output_file <- paste0(output_dir, "Nov5GLMBotht", ssize, "em_", em, ".txt")
  # Main loop to run simulations and apply chosen methods
  for (i in 1:B) {
    # Generate data
    TNDdat <- datagen(ssize = ssize, em = em)
    res <- element_glm(TNDdat, ps_model = ps_model, out_model1 = out_model1, out_model2 = out_model2)
    
    est1 <- mod_IPW(TNDdat, res, bootstrap_CI = bootstrap_CI, ps_model, out_model1, out_model2)
    est1_val <- est1$est
    CI1 <- est1$CI
    
    est2 <- mod_OutRegGLM(TNDdat, res, bootstrap_CI = bootstrap_CI, ps_model, out_model1, out_model2)
    est2_val <- est2$est
    CI2 <- est2$CI.OR
    
    est3 <- mod_EIF2(TNDdat, res, bootstrap_CI = F)
    est3_val <- est3$est
    TNDCI1 <- est3$CI1
    TNDCI2 <- est3$CI2
    TNDCI3 <- est3$CI3
    
    # Save results
    write(
      c(i, est1_val, CI1, est2_val, CI2, est3_val,TNDCI1, TNDCI2, TNDCI3),
      file = output_file,
      ncolumns = 50,
      append = TRUE
    )
  }
}

# C + log(C) + sin(pi * C)
# V + C + V * C + + 0.5 * exp(C)* (1 + 0.15 * cos(C))
ps_model <- V ~ C + log(C) + sin(pi * C)
out_model1 <- Y ~ C + V + I(V * C) + exp(C) + exp(C) * cos(C)
out_model2 <- Y ~ C + exp(C) + exp(C) * cos(C)
TNDresGLM(B = 500, ssize = 8000, em = -0.25, bootstrap_CI = TRUE, ps_model = ps_model, out_model1 = out_model1, out_model2 = out_model2)









####################################################################################
####### Logistic regression
#######
####################################################################################


mu.fit <- glm(out_model1,family="binomial", data=TNDdat)
exp(coef(mu.fit)["V"])

mod_logistic <- function(TNDdat, bootstrap_CI = TRUE, out_model){
  mu.fit <- glm(out_model,family="binomial", data=TNDdat)
  res.logistic <- exp(coef(mu.fit)["V"])
  if (bootstrap_CI) {
    nbs <- 200
    bsest <- rep(NA, nbs)
    
    for(i in 1:nbs){
      resamps <- sample(1:nrow(TNDdat), size = nrow(TNDdat), replace = TRUE)
      datk <- TNDdat[resamps,]
      mu.fit <- glm(out_model,family="binomial", data=datk)
      bsest[i] <- exp(coef(mu.fit)["V"])
    }
    bs_var <- var(bsest)
    CI <- quantile(bsest, c(0.025, 0.975), na.rm = TRUE)
  } else {
    CI <- NA
  }
  return(list(est = res.logistic, CI.logistic = CI))
}

mod_logistic(TNDdat, out_model = out_model1)

TNDresGLM <- function(B = 1000, ssize = 1000, em = -0.15, bootstrap_CI = FALSE, ps_model, out_model1) {
  # Construct the output file name based on the method
  output_file <- paste0(output_dir, "logistic", ssize, "em_", em, ".txt")
  # Main loop to run simulations and apply chosen methods
  for (i in 1:B) {
    # Generate data
    TNDdat <- datagen(ssize = ssize, em = em)
    
    est0 <- mod_logistic(TNDdat, out_model = out_model1)
    est0_val <- est0$est
    CI0 <- est0$CI.logistic
    # Save results
    write(
      c(i, est0_val, CI0),
      file = output_file,
      ncolumns = 50,
      append = TRUE
    )
  }
}

out_model1 <- Y ~ C + V
#TNDresGLM(B = 500, ssize = 1000, em = 0.25, bootstrap_CI = TRUE, ps_model = ps_model, out_model1 = out_model1, out_model2 = out_model2)


library(geex)

psglm <- glm(V ~ C, data = subset(data,data$Y==0), family = binomial)

# Define the estimating function based on the system of equations
my_estfun <- function(data) {
  Y <- data$Y
  V <- data$V
  C <- as.model.matrix(subset(data, select = -c(Y, V))) # assuming C is a matrix of covariates
  
  # Define the estimating function that depends on theta
  function(theta) {
    alpha <- theta[1:length(C[1, ])]  # extracting alpha (logistic regression coefficients)
    psi_v <- theta[length(C[1, ]) + 1]  # psi_v parameter
    psi_v0 <- theta[length(C[1, ]) + 2]  # psi_v0 parameter
    psi_mRR <- theta[length(C[1, ]) + 3]  # psi_mRR parameter
    
    expit_C_alpha <- 1 / (1 + exp(-C %*% alpha))  # expit(C_i^T * alpha)
    
    # The system of equations
    c(
      (Y == 0) * ((V - expit_C_alpha) %*% C),  # First equation
      (q0 * Y * V) / expit_C_alpha - psi_v,  # Second equation
      (q0 * Y * (1 - V)) / (1 - expit_C_alpha) - psi_v0,  # Third equation
      (psi_v / psi_v0) - psi_mRR  # Fourth equation
    )
  }
}

# Example data
set.seed(123)
n <- 100
data <- data.frame(
  Y = rbinom(n, 1, 0.5),      # Binary outcome Y
  V = rbinom(n, 1, 0.5),      # Binary or continuous variable V
  C = matrix(rnorm(n * 3), n)  # Matrix of covariates C (3 covariates in this example)
)

# Define the q0 constant (assumed constant here, adjust as needed)
q0 <- 1

# Use m_estimate to estimate the parameters
results <- m_estimate(
  estFUN = my_estfun,
  data = data,
  root_control = setup_root_control(start = c(rep(0, 3), 1, 1, 1))  # starting values for theta (alpha, psi_v, psi_v0, psi_mRR)
)


alpha <- rep(0,3)
# Check the results
summary(results)
