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

source("ElementFunction.R")
source("TrueValue.R")
output_dir <- ""
# Define a function map for model fitting methods
fit_methods <- list(
  HAL = element.HAL,
  Earth = element.earth,
  RandomForest = element.randomForest,
  NN = element.nnet
)


TNDres <- function(B = 1000, ssize = 1000, co_inf_para = 0.00001, em = 0.25, method = "NN", bootstrap_CI = FALSE){
  #method = "RandomForest"
  # Define a list of methods to apply
  methods <- c("HAL", "Earth", "RandomForest", "NN")
  # Check if the method is supported
  if (!method %in% names(fit_methods)) {
    stop("Unsupported method. Please choose from 'HAL', 'Earth', 'RandomForest', or 'NN'.")
  }
  
  # Construct the output file name based on the method
  output_file <- paste0(output_dir, "Nov5_", method, ssize, "em_", em, ".txt")
  # Main loop to run simulations and apply chosen methods
  for (i in 1:B) {
    # Generate data
    TNDdat <- datagen(ssize = ssize, em = em, co_inf_para = co_inf_para)
    # Fit model using the corresponding function from the function map
    res <- fit_methods[[method]](TNDdat)
    # Apply estimation functions and collect results
    #est0 <- mod_IPW(TNDdat, res)
    #est0_val <- est0$est
    #CI0 <- est0$CI
    
    est1 <- mod_IPW1(TNDdat, res, bootstrap_CI = bootstrap_CI, method = method)
    est1_val <- est1$est
    CI1 <- est1$CI
    
    est2 <- mod_OutReg(TNDdat, res, bootstrap_CI = bootstrap_CI, method = method)
    est2_val <- est2$est
    CI2 <- est2$CI.OR
    
    est3 <- mod_TNDDR(TNDdat, res)
    est3_val <- est3$est
    TNDCI1 <- est3$CI1
    TNDCI2 <- est3$CI2
    TNDCI3 <- est3$CI3
    
    # Save results
    write(
      c(i, est1_val, CI1,est2_val, CI2, est3_val, TNDCI1, TNDCI2, TNDCI3),
      file = output_file,
      ncolumns = 50,
      append = TRUE
    )
  }
}





TNDres(B = 500, ssize = 8000, em = 0.25, co_inf_para = 0.00001,  method = "HAL", bootstrap_CI = F)

