### Simulation of p-hacking ###
# Date: 11.05.2024

## save R info 
# Rsession <- sessionInfo()
# saveRDS(Rsession, file='Rsession.rds')


## libraries
library('MASS')
library('dplyr')
library('ggplot2')
library('tidyverse')
library('parallel')
library('microbenchmark') #check speed of functions 
library('data.table') # more flexible data structure
library('mvtnorm')


## function for DGP of p-hacking 
sim_phacking <- function(m = 5, rho = 0, sd = 1, mhu = 0) {
  corr_matrix <- matrix(rho, nrow = m, ncol = m)
  diag(corr_matrix) <- 1
  
  std_devs <- rep(sd, m)
  custom_covariance <- corr_matrix * outer(std_devs, std_devs)
  diag(custom_covariance) <- 1 
  
  # set.seed(42)
  data_mvrnorm <- mvrnorm(m, rep(mhu, m), custom_covariance)
  set.seed(99)
  data_mvrnorm_zr <- mvrnorm(m, rep(mhu, m), custom_covariance)
  data_mvrnorm <- as.data.frame(data_mvrnorm)
  
  data_all <- data.frame()
  for (i in 1:nrow(data_mvrnorm)) {
    position_max <- which(abs(data_mvrnorm[i,]) == max(abs(data_mvrnorm[i,])))[1]
    zo_value <- data_mvrnorm[i, position_max]
    zr_value <- data_mvrnorm_zr[i, position_max]
    
    po <- (2 * (1 - pnorm(abs(zo_value))))
    pr <- (2 * (1 - pnorm(abs(zr_value))))
    
    sig_10p_zo <- ifelse(po <= 0.10, 1, 0)
    sig_5p_zo <- ifelse(po <= 0.05, 1, 0)
    sig_1p_zo <- ifelse(po <= 0.01, 1, 0)
    sig_05p_zo <- ifelse(po <= 0.005, 1, 0)
    
    sig_10p_zr <- ifelse(pr <= 0.10, 1, 0)
    sig_5p_zr <- ifelse(pr <= 0.05, 1, 0)
    sig_1p_zr <- ifelse(pr <= 0.01, 1, 0)
    sig_05p_zr <- ifelse(pr <= 0.005, 1, 0)
    
    results <- data.frame(
      df = paste("sim_m", m, "_r", rho, "_ew", mhu, sep = ""),
      position_zo = position_max,
      zo = zo_value,
      zr = zr_value,
      po = po,
      pr = pr,
      sig_10p_zo = sig_10p_zo,
      sig_5p_zo = sig_5p_zo,
      sig_1p_zo = sig_1p_zo,
      sig_05p_zo = sig_05p_zo,
      sig_10p_zr = sig_10p_zr,
      sig_5p_zr = sig_5p_zr,
      sig_1p_zr = sig_1p_zr,
      sig_05p_zr = sig_05p_zr,
      p_hacking_10p = ifelse(sig_10p_zo == 1 & sig_10p_zr == 0, 1, 0),
      p_hacking_5p = ifelse(sig_5p_zo == 1 & sig_5p_zr == 0, 1, 0),
      p_hacking_1p = ifelse(sig_1p_zo == 1 & sig_1p_zr == 0, 1, 0),
      p_hacking_05p = ifelse(sig_05p_zo == 1 & sig_05p_zr == 0, 1, 0),
      param_m = m, 
      param_rho = rho, 
      param_mhu = mhu
    )
    
    data_all <- rbind(data_all, results)
  }
  
  return(list(data_all = data_all, data_mvrnorm = data_mvrnorm))
}


# Function to generate datasets with different settings in parallel
sim_datasets_phacking <- function(m = 5, rho = 0, mhu = 0, niter = 10) {
  # Generate parameter combinations
  param_combi <- expand.grid(niter = 1:niter, m = m, rho = rho, mhu = mhu)
  
  # Set seed for reproducibility 
  set.seed(42)
  param_combi$seeds <- sample(1:1000000, nrow(param_combi), replace = FALSE)
  
  # Set number of cores to use 
  n_cores <- detectCores() / 2
  
  # Function to generate data for each combination
  generate_data <- function(i) {
    set.seed(param_combi$seeds[i])  # Random seed for each iteration
    data <- sim_phacking(m = param_combi$m[i], rho = param_combi$rho[i], mhu = param_combi$mhu[i])
    data$data_all <- data$data_all %>% mutate(seed = param_combi$seeds[i], niter = param_combi$niter[i])
    return(data)
  }
  
  # Use mclapply for parallel execution
  data_all_list <- mclapply(1:nrow(param_combi), generate_data, mc.cores = n_cores)
  
  # Combine all data frames
  data_all_combined <- bind_rows(lapply(data_all_list, function(x) x$data_all))
  data_mvrnorm_combined <- bind_rows(lapply(data_all_list, function(x) x$data_mvrnorm))
  
  # Return the combined data frame
  return(list(data_all = data_all_combined, data_mvrnorm = data_mvrnorm_combined))
}

# Example usage
# result <- sim_datasets_phacking(m = c(2, 5), rho = c(0, 0.5), niter = 10)
# head(result$data_all)
# head(result$data_mvrnorm)
# tail(result$data_all)
# tail(result$data_mvrnorm)



## Function for p-value adjustment 
adjusted_pm <- function(zo, m, rho, mean = 0) {
  #Set z_values to absolute values
  lower_zo <- rep(-abs(zo), m)
  upper_zo <- rep(abs(zo), m)
  
  # Create correlation matrix
  corr_matrix <- matrix(rho, nrow = m, ncol = m)
  diag(corr_matrix) <- 1 
  
  padj <- pmvnorm(lower = lower_zo, upper = upper_zo, mean = rep(mean, m), sigma = corr_matrix, keepAttr=F)
  
  pm <- 1-padj
  return(pm)
}

##check function 
# adjusted_pm(zo = sim_phacking_test$zo[1:m] , m = 5, rho = 0.5, mean=1)


## find correlation for range of pm 
find_max_correlation <- function(corr_matrix) {
  # Convert the correlation matrix to a data frame for easier manipulation
  corr_df <- as.data.frame(as.table(corr_matrix))
  
  # Filter out the perfect correlations (diagonal elements)
  corr_df <- corr_df[corr_df$Freq != 1, ]
  
  # Find the maximum correlation
  max_corr <- max(corr_df$Freq)
  
  return(max_corr)
}

