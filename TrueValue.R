#Common data-generating function
#popsize set at 15000000, can increase in large ssize is needed
#ssize set at 5000, may have trouble when not enough patient available

# Note: The true values vary for each setting with different parameters, 
# such as the values of em (coefficient of modification).

datagen <- function(seed=sample(1:1000000,size=1), ssize = 5000, popsize = 15000000, OR_C = 3,
                      OR_W1 = 1, OR_W2 = 2.5, em = 0.25,
                      cfV0 = FALSE, cfV1 = FALSE, return_full = FALSE, co_inf_para = 0.00001) {
    set.seed(seed)
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

dat<-datagen(ssize=5000)


#true marginal RR
orvect<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect[j]<-mean(datfull1$H*datfull1$Infec_COVID)/mean(datfull0$H*datfull0$Infec_COVID)
}
(psi = mean(orvect)) #0.04209
hist(orvect)

#true marginal RR for Infec_COVID
orvect1<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect1[j]<-mean(datfull1$Infec_COVID)/mean(datfull0$Infec_COVID)
}
mean(orvect1) 
hist(orvect1)

#true marginal RR for Infec_COVID * W (severe COVID disease)
orvect2<-rep(NA,50)
for (j in 1:50){
  datfull1<-datagen(cfV1=T,return_full=T)
  datfull0<-datagen(cfV0=T,return_full=T)
  orvect2[j]<-mean(datfull1$W*datfull1$Infec_COVID)/mean(datfull0$W*datfull0$Infec_COVID)
}
mean(orvect2) 
hist(orvect2)






####################################################################################
####### DGP parameters' setting
#######
####################################################################################



# Cong's updated data generating code
datagen_c <- function(ssize = 5000, popsize = 15000000, OR_C = 3,
                      OR_W1 = 1, OR_W2 = 2.5, em = 0.5,
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


prevalence_c_fun = function(effmod, nn, ssize=5000, OR_C = 3,
                            OR_W1 = 1, OR_W2 = 2.5, co_inf_para = 1e-3) {
  stats = matrix(NA, nrow = nn, ncol = 18)
  for (j in 1:nn) {
    da = datagen_c(return_full=T, ssize=ssize, em= effmod, OR_C = OR_C,
                   OR_W2 = OR_W2, co_inf_para = co_inf_para)
    tnd = sample(which(da$H==1), ssize)
    stats[j, ] = c(mean(da$V), mean(da$Infec), 
                   mean(da$Infec_COVID), mean(da$Infec==1 & da$Infec_COVID==1), 
                   mean(da$W==1), 
                   mean(da$W[da$Infec_COVID==1]), 
                   mean(da$W[da$Infec==1|da$Infec_COVID==1]),
                   mean(da$W[da$Infec_COVID==1& da$V==1]), mean(da$W[da$Infec_COVID==1& da$V==0]), 
                   mean(da$H), mean(da$H[da$Infec_COVID==1]), mean(da$H[da$W==1]), 
                   mean(da$H[da$Infec_COVID==1& da$V==1]), mean(da$H[da$Infec_COVID==1& da$V==0]), 
                   
                   mean(da[tnd,]$Infec_COVID), mean(da[tnd,]$Infec),
                   sum(rowSums(da[tnd, c("Infec","Infec_COVID")])==2)/ssize, mean(da[tnd,]$V))
  }
  return(stats= stats)
}

prevalences_em1 = prevalence_c_fun(effmod=0.25, OR_C = 3,OR_W2 = 2.5, nn=30, co_inf_para = 0.00001)
prevalences_em2 = prevalence_c_fun(effmod=0.25, OR_C = 3,OR_W2 = 2.5, nn=30, co_inf_para = 0.001)
prevalences_em2 = prevalence_c_fun(effmod=0.25, OR_C = 3,OR_W2 = 2.5, nn=30, co_inf_para = 0.1)


prev_c = round(data.frame(rbind(colMeans(prevalences_em1), colMeans(prevalences_em2), colMeans(prevalences_em3))), digits =4)
colnames(prev_c)= c("p_V","p_I1","p_I2", "p_[I1=I2=1]", "p_W", "p_W[I2=1]", "p_W[I1=1|I2=1]",  
                    "p_W[V=1,I2=1]","p_W[V=0,I2=1]",
                    "p_H","p_H[I2=1]","p_H[W=1]","p_H[V=1,I2=1]","p_H[V=0,I2=1]",
                    "p_TND_I2", "p_TND_I1", "p_TND_[I1=I2=1]", "p_TND_V")
rownames(prev_c)= c("co_inf_para = 0.00001","co_inf_para = 0.001", "co_inf_para = 0.1")
t(prev_c)



########### compute true values of cRR and mRR
tru_c_fun <- function(effmod, nn, ssize=5000, co_inf_para = 1e-5) {
  tru_RR <- matrix(NA, nrow = nn, ncol = 2)
  for (j in 1:nn){
    datfull <- datagen_c(return_full=T, em= effmod, ssize=ssize, co_inf_para = co_inf_para)
    Qmod <- glm(I(Infec_COVID*W*H) ~ V+C, data=datfull, family = binomial)
    tru_RR[j, 1] <- exp(Qmod$coef[2])
    
    datfull1 <- datagen_c(cfV1=T, return_full=T, em= effmod, ssize=ssize, co_inf_para = co_inf_para)
    datfull0 <- datagen_c(cfV0=T, return_full=T, em= effmod, ssize=ssize, co_inf_para = co_inf_para)
    tru_RR[j, 2] <- mean(datfull1$H*datfull1$Infec_COVID)/mean(datfull0$H*datfull0$Infec_COVID)
  }
  return(colMeans(tru_RR, na.rm = T))
}
tru_c_fun(effmod = 1, nn=2)
tru_c_fun(effmod = 0.25, nn=2)


