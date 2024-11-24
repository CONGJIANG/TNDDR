library(dplyr)

# Check significance of V and output message
sig_msg <- function(model, cond) {
  p_value <- summary(model)$coefficients["V", 4]
  if (p_value < 0.05) {
    return(paste(cond, ": V is significant (p =", round(p_value, 4), ")"))
  } else {
    return(paste(cond, ": V is not significant (p =", round(p_value, 4), ")"))
  }
}


coinf_ctrl_exch <- function(ssize = 5000, popsize = 15000000, OR_C = 3,
                            OR_W1 = 1, OR_W2 = 2.5, em = 0.15,
                            cfV0 = FALSE, cfV1 = FALSE, return_full = FALSE, co_inf_para = 0.00001) {
  # Generate confounders and covariates
  C <- runif(n = popsize, 0.1, 3)  # Continuous confounder
  U1 <- rbinom(n = popsize, size = 1, prob = 0.5)  # Unmeasured binary covariate
  U2 <- rbinom(n = popsize, size = 1, prob = 0.5)  # Unmeasured binary covariate
  
  # Generate Vaccination Status
  if (cfV0 == TRUE) {
    V <- rep(0, popsize)
  } else if (cfV1 == TRUE) {
    V <- rep(1, popsize)
  } else {
    V <- rbinom(n = popsize, size = 1, prob = plogis(0.5 * (0.5 + 1.5 * C - log(C) - 2.5 * sin(pi * C))))
  }
  
  # Convert to log-odds for logistic regression
  bsl_I1 <- qlogis(co_inf_para)
  bsl_I2 <- qlogis(co_inf_para)
  
  # Generate independent infections
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(0.35 * C + bsl_I1 + 6.5 * U1))  # Infection with other viruses
  lambda_covid <- bsl_I2 + 0.15 * C + 0.5 * exp(C)* (1 + 0.15 * cos(C)) -log(OR_C)*V + em * V * C + log(1.2) * U2 * (1.5 - V) - 2 * U1
  I2 <- rbinom(n = popsize, size = 1, prob = plogis(lambda_covid))  # Infection with SARS-CoV-2
  
  # Calculate the percentage of co-infections
  co_inf <- sum(I1 == 1 & I2 == 1)
  per_co_inf <- co_inf / popsize * 100
  
  # Generate symptoms W1 and W2
  W1 <- rep(0, popsize)
  W1[I1 == 1] <- rbinom(
    n = sum(I1 == 1),
    size = 1,
    prob = plogis(-0.5 + 0.5 * C[I1 == 1] - log(OR_W1) * V[I1 == 1] - 0.5 * U1[I1 == 1])
  )
  
  W2 <- rep(0, popsize)
  W2[I2 == 1] <- rbinom(
    n = sum(I2 == 1),
    size = 1,
    prob = plogis(-3.75 + 2 * C[I2 == 1] - log(OR_W2) * V[I2 == 1] - 1 * U1[I2 == 1] + 0.5 * U2[I2 == 1] * (1 - V[I2 == 1]))
  )
  
  # Combine symptoms into a unified status
  W <- pmax(W1, W2)
  
  #hospitalization, only possible if symptoms present
  H=rep(0,popsize)
  H[W==1]<-rbinom(prob=plogis(-1.5+0.5*C[W==1] -1.5*U1[W==1]),size=1,n=sum(W==1))
  #mean(H[W==1]) # with severe symptoms go to hospital
  # Create the data frame
  dat <- data.frame(C, U1, U2, V, I1, I2, W1, W2, W, H)
  # Controls conditions
  dat$Con1 <- dat$H == 1 & dat$I2 == 0
  dat$Con2 <- dat$I1* dat$H == 1
  
  # Fit logistic regression models for conditions
  model_Con1 <- glm(Con1 ~ V + C, data = dat, family = binomial)
  model_Con2 <- glm(Con2 ~ V + C, data = dat, family = binomial)
  
  msg_Con1 <- sig_msg(model_Con1, "Con1")
  msg_Con2 <- sig_msg(model_Con2, "Con2")
  return(list(
    per_co_inf = per_co_inf,
    significance_messages = list(msg_Con1, msg_Con2)
  ))
}

coinf_ctrl_exch(em = 0.15)
coinf_ctrl_exch(em = 0.25)



datagen <- function(ssize = 5000, popsize = 15000000, OR_C = 3,
                    OR_W1 = 1, OR_W2 = 2.5, em = 0.15,
                    cfV0 = FALSE, cfV1 = FALSE, return_full = FALSE, co_inf_para = 0.00001) {
  # Generate confounders and covariates
  C <- runif(n = popsize, 0.1, 3)  # Continuous confounder
  U1 <- rbinom(n = popsize, size = 1, prob = 0.5)  # Unmeasured binary covariate
  U2 <- rbinom(n = popsize, size = 1, prob = 0.5)  # Unmeasured binary covariate
  
  # Generate Vaccination Status
  if (cfV0 == TRUE) {
    V <- rep(0, popsize)
  } else if (cfV1 == TRUE) {
    V <- rep(1, popsize)
  } else {
    V <- rbinom(n = popsize, size = 1, prob = plogis(0.5 * (0.5 + 1.5 * C - log(C) - 2.5 * sin(pi * C))))
  }
  
  # Convert to log-odds for logistic regression
  bsl_I1 <- qlogis(co_inf_para)
  bsl_I2 <- qlogis(co_inf_para)
  
  # Generate independent infections
  I1 <- rbinom(n = popsize, size = 1, prob = plogis(0.35 * C + bsl_I1 + 6.5 * U1))  # Infection with other viruses
  lambda_covid <- bsl_I2 + 0.15 * C + 0.5 * exp(C)* (1 + 0.15 * cos(C)) -log(OR_C)*V + em * V * C + log(1.2) * U2 * (1.5 - V) - 2 * U1
  I2 <- rbinom(n = popsize, size = 1, prob = plogis(lambda_covid))  # Infection with SARS-CoV-2
  
  # Calculate the percentage of co-infections
  #co_inf <- sum(I1 == 1 & I2 == 1)
  #per_co_inf <- co_inf / popsize * 100
  
  # Generate symptoms W1 and W2
  W1 <- rep(0, popsize)
  W1[I1 == 1] <- rbinom(
    n = sum(I1 == 1),
    size = 1,
    prob = plogis(-0.5 + 0.5 * C[I1 == 1] - log(OR_W1) * V[I1 == 1] - 0.5 * U1[I1 == 1])
  )
  
  W2 <- rep(0, popsize)
  W2[I2 == 1] <- rbinom(
    n = sum(I2 == 1),
    size = 1,
    prob = plogis(-3.75 + 2 * C[I2 == 1] - log(OR_W2) * V[I2 == 1] - 1 * U1[I2 == 1] + 0.5 * U2[I2 == 1] * (1 - V[I2 == 1]))
  )
  
  # Combine symptoms into a unified status
  W <- pmax(W1, W2)
  
  #hospitalization, only possible if symptoms present
  H=rep(0,popsize)
  H[W==1]<-rbinom(prob=plogis(-1.5+0.5*C[W==1] -1.5*U1[W==1]),size=1,n=sum(W==1))
  #mean(H[W==1]) # with severe symptoms go to hospital
  # Selection on outcome for testing (hospitalization)
  R <- sample(which(H == 1), ssize, replace = TRUE)  # Sample randomly from those hospitalized
  
  if (return_full == FALSE) {
    dat <- as.data.frame(cbind(Y = I2, V = V, C = C)[R, ])
  } else {
    dat <- as.data.frame(cbind(Infec_COVID = I2, Infec = I1, H = H, W = W, V = V, C = C))
  }
  
  return(dat)
}



# Example of using the integrated function
test_data <- datagen(ssize = 1000, em = 0.25)
head(test_data)
sum(test_data$Y)/nrow(test_data)
sum(test_data$V)/nrow(test_data)

test.full <- datagen(ssize = 1000, em = 0.25, OR_C = 3,
                     OR_W1 = 1, OR_W2 = 2.5, return_full = T)
test.tnd <- test.full[test.full$H == 1,]
co_inf <- sum(test.tnd$Infec_COVID == 1 & test.tnd$Infec == 1)
(per_co_inf <- co_inf / nrow(test.tnd) * 100)

sum(test.tnd$Infec_COVID)/nrow(test.tnd)
sum(test.tnd$V)/nrow(test.tnd)

