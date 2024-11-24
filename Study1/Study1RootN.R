source("ElementFunction.R")
source("TrueValue.R")

# Define a function map for model fitting methods
fit_methods <- list(
  HAL = element.HAL,
  Earth = element.earth,
  RandomForest = element.randomForest,
  NN = element.nnet
)

RootN <- function(ssize = 5000, nsim = 300,  method = "NN", bootstrap_CI = FALSE){
  n <- seq(300, ssize, by = 200)
  estIPW <- matrix(NA, nsim, length(n))
  estOR <- matrix(NA, nsim, length(n))
  estEIF <- matrix(NA, nsim, length(n))
  methods <- c("HAL", "Earth", "RandomForest", "NN")
  # Check if the method is supported
  if (!method %in% names(fit_methods)) {
    stop("Unsupported method. Please choose from 'HAL', 'Earth', 'RandomForest', or 'NN'.")
  }
  
  for (j in 1:nsim) {
    TndDAT <- datagen(ssize = ssize, em = 0.25)
    for (i in 1:length(n)) {
      TNDdat <- TndDAT[1:n[i], ]
      res <- fit_methods[[method]](TNDdat)
      
      est1 <- mod_IPW1(TNDdat, res, bootstrap_CI = bootstrap_CI, method = method)
      estIPW[j,i] <- est1$est
      
      
      est2 <- mod_OutReg(TNDdat, res, bootstrap_CI = bootstrap_CI, method = method)
      estOR[j,i] <- est2$est
      
      
      est3 <- mod_TNDDR(TNDdat, res)
      estEIF[j,i] <- est3$est
    }
  }
  return(list(IPW = estIPW, OR = estOR, TNDDR = estEIF))
}

root.res <- RootN(ssize = 10000, nsim = 300,  method = "HAL")


mRR.em1 = 0.508 # em = -0.25 
estIPW <- root.res$IPW;  estOR <- root.res$OR;  estEIF <- root.res$TNDDR;  
# Compute the biases
bias.ipw <- (apply(na.omit(estIPW), 2, median) - mRR.em1) 
bias.or <- (apply(na.omit(estOR), 2, median) - mRR.em1) 
bias.eif <- (apply(na.omit(estEIF), 2,  median) - mRR.em1)




n = seq(300, 10000, by = 200)
res <- as.data.frame(cbind(n = n, bias.ipw = bias.ipw,bias.or = bias.or, bias.eif = bias.eif))
# Plot the biases as a function of sample size
par(mar = c(5, 5, 4, 2) + 0.4)  # Adjust the margin settings
plot(res$n, res$bias.or, type = "b", pch = 4, col = "blue",
     xlab = "Sample size (n)", ylab = expression(paste(hat(psi)[mRR]  - psi[mRR])), ylim = c(-0.2, 0.2))
lines(res$n, res$bias.eif, type = "b",pch = 19, col = "red")
lines(res$n, res$bias.ipw, type = "b", pch = 2, col = "black")
lines(res$n, 1/res$n^(1/2), pch = 0,  col = "brown")
lines(res$n, -1/sqrt(res$n), pch = 0,  col = "brown")
abline(h = 0, lty = 2)
legend("bottomright", legend = c("IPW","TNDDR", "OutReg", expression(1/sqrt(n))), pch = c(2, 19, 4, NA), col = c("black","red", "blue", "brown"), lty = 1)
title(expression(paste("Root-n Consistency of TNDDR (NN) (em = 0.25)")))

write.csv(res, file = "Nov7TNDStudy1HAL1rootn(em0.25).csv")
res <- read.csv("Nov7TNDStudy1Earthrootn(em0.25).csv")

