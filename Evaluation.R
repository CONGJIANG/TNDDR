compute_metrics <- function(data, theta_true, method_cols, ci_bounds_pairs_list) {
  # Initialize a nested list to store results for each method and CI pair
  results_list <- list()
  
  # Loop over each specified method column
  for (method_name in method_cols) {
    # Extract estimates for the current method
    estimates <- data[[method_name]]
    
    # Calculate the bias
    bias <- estimates - theta_true
    
    # Calculate metrics
    mean_bias <- mean(bias, na.rm = TRUE)
    median_bias <- median(bias, na.rm = TRUE)
    rmse <- sqrt(mean(bias^2, na.rm = TRUE))
    mc_standard_error <- sd(estimates, na.rm = TRUE) / sqrt(length(estimates))
    
    # Initialize a list to store CI coverage results for different CI pairs
    ci_coverage_results <- list()
    
    # Loop through each set of CI bounds for the current method
    for (i in seq_along(ci_bounds_pairs_list[[method_name]])) {
      ci_bounds <- ci_bounds_pairs_list[[method_name]][[i]]
      lower_bounds <- data[[ci_bounds[1]]]
      upper_bounds <- data[[ci_bounds[2]]]
      
      # Calculate CI coverage
      ci_covered <- sum(theta_true >= lower_bounds & theta_true <= upper_bounds, na.rm = TRUE)
      ci_coverage <- ci_covered / length(estimates)
      
      # Give each CI coverage a unique identifier
      ci_identifier <- paste("CI", i, sep = "")
      ci_coverage_results[[ci_identifier]] <- ci_coverage
    }
    
    # Store results in a list
    results_list[[method_name]] <- list(
      mean_bias = mean_bias,
      median_bias = median_bias,
      rmse = rmse,
      mc_standard_error = mc_standard_error,
      ci_coverages = ci_coverage_results
    )
  }
  
  return(results_list)
}



output_file <- "/Users/congjiang/Downloads/rstudio-export/HAL5000em_0.25.txt"
res <- read.table(output_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
colnames(res)[c(1, 2, 4, 6)] <- c("i", "IPW", "Out", "EIF")
summary(res, na.rm = FALSE)

# Usage example

theta_true <- 0.507  # Replace with the actual known true parameter value
method_columns <- c("IPW", "Out", "EIF")  # Columns for method estimates

# For EIF, specify all three different CI bounds
ci_bounds_pairs_list <- list(
  "EIF"  = list(c("V7", "V8"),
                c("V9", "V10"))       # CI bounds for Out   
)

results <- compute_metrics(res, theta_true, method_columns, ci_bounds_pairs_list)

# Print results for each method
print(results)


output_file <- "/Users/congjiang/Downloads/rstudio-export/NN5000em_0.25.txt"
res <- read.table(output_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
colnames(res)[c(1, 2, 5, 8)] <- c("i", "IPW", "Out", "EIF")
summary(res, na.rm = FALSE)

# Usage example

theta_true <- 0.507  # Replace with the actual known true parameter value
method_columns <- c("IPW", "Out", "EIF")  # Columns for method estimates

# For EIF, specify all three different CI bounds
ci_bounds_pairs_list <- list(
  "IPW" = list(c("V3", "V4")),          # CI bounds for IPW1
  "Out" = list(c("V6", "V7")),          # CI bounds for IPW2
  "EIF"  = list(c("V9", "V10"),
                c("V11", "V12"))       # CI bounds for Out   
)

results <- compute_metrics(res, theta_true, method_columns, ci_bounds_pairs_list)

# Print results for each method
print(results)



# Usage example
output_file <- "/Users/congjiang/Downloads/rstudio-export/GLMBotht5000em_0.25.txt"
res <- read.table(output_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
colnames(res)[c(1, 2, 5, 8)] <- c("i", "IPW", "Out", "EIF")
res <- res[!(res$EIF < 0.00100000 | res$EIF > 0.999), ]
summary(res, na.rm = FALSE)

theta_true <- 0.124  # Replace with the actual known true parameter value
method_columns <- c("IPW", "Out", "EIF")  # Columns for method estimates

# For EIF, specify all three different CI bounds
ci_bounds_pairs_list <- list(
  "IPW" = list(c("V3", "V4")),          # CI bounds for IPW1
  "Out" = list(c("V6", "V7")),          # CI bounds for IPW2
  "EIF"  = list(c("V9", "V10"),         # CI bounds for Out
                c("V11", "V12"),         # CI1 bounds for EIF
                c("V13", "V14"))         # CI3 bounds for EIF
)

results <- compute_metrics(res, theta_true, method_columns, ci_bounds_pairs_list)

# Print results for each method
print(results)

