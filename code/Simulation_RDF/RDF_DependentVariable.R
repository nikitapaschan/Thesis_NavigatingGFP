######################################################
# RDF Selective reporting of dependent variables
#####################################################
setwd('/Users/nikitapaschan/Thesis_NavigatingGFP/')


# Function to generate correlated data
generate_correlated_matrix <- function(n = 100, k = 3, r = 0, mu = 0) {
  mean_vector <- rep(mu, k)
  correlation_matrix <- matrix(r, nrow = k, ncol = k)
  diag(correlation_matrix) <- 1
  
  correlated_matrix <- mvrnorm(n = n, mu = mean_vector, Sigma = correlation_matrix)
  return(correlated_matrix)
}

# Function to perform statistical tests and collect results
perform_tests <- function(n = 100, k = 5, r = 0.5, mu = 0, seed = 123, iteration = 1) {
  set.seed(seed)
  
  # Generate data
  data_matrix <- generate_correlated_matrix(n = n, k = k, r = r, mu = mu)
  
  # Create grouping variable
  g1 <- sample(c(rep(0, n/2), rep(1, n/2)), n, replace = FALSE)
  data_testing <- as.data.frame(cbind(data_matrix, g1))
  
  # Exclude the last column from the data frame for tests
  data_test <- data_testing[, -ncol(data_testing)]
  
  # Create an empty data frame to store the results
  result_testDV <- data.frame()
  
  # Iterate through each column and perform tests
  for (col in names(data_test)) {
    
    # T-test
    ttest <- t.test(data_testing[g1 == 1, col], data_testing[g1 == 0, col])
    
    
    # Capture results in data set
    result_testDV <- rbind(result_testDV,
                           data.frame(
                             iteration = iteration,
                             n = n,
                             k = k,
                             r = r,
                             mu = mu,
                             seed = seed,
                             Variable = col,
                             t_test = ttest$statistic,
                             pvalue_ttest = ttest$p.value,
                             sig_ttest_5p = ifelse(ttest$p.value < 0.05, 1, 0)
                             
                           ))
  }
  
  return(result_testDV)
}

# Example usage
set.seed(142)
k_values <- c(3, 5, 10) #add one setting with 1000 to see extreme 
r_values <- c(0, 0.3, 0.8)
niter <- 10000

param_combi_testing <- expand.grid(niter = 1:niter, k = k_values, r = r_values)
param_combi_testing$seed <- sample(1:1000000, nrow(param_combi_testing), replace = FALSE)

# Create a list to store results
all_results <- list()

# Perform tests for each combination of parameters
for (i in 1:nrow(param_combi_testing)) {
  params <- param_combi_testing[i, ]
  results <- perform_tests(n = 100, k = params$k, r = params$r, mu = 0, seed = params$seed, iteration = params$niter)
  all_results[[i]] <- results
}

# Combine all results into a single data frame
result_testDV <- bind_rows(all_results)

# Print the final results
head(result_testDV)

```

```{r}
# Summarize the data to get the maximum absolute t_test for each parameter combination
phacked_testDV <- result_testDV %>%
  group_by(iteration, n, k, r, seed) %>%
  summarize(
    t_test = max(abs(t_test), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  inner_join(result_testDV, by = c("iteration", "n", "k", "r", "seed", "t_test")) %>%
  select(c(iteration, n, k, r, seed, Variable, t_test, pvalue_ttest, sig_ttest_5p))

head(phacked_testDV)
```

```{r}
# Calculate correlations

result_testDV <- as.data.table(result_testDV)
cor_testDV <- data.table()

for (k_value in unique(result_testDV$k)) { #unique(result_testDV$k)
  for (r_value in unique(result_testDV$r)) {
    variables <- unique(result_testDV[k == k_value]$Variable)
    for (i in 1:(length(variables) - 1)) {
      for (j in (i + 1):length(variables)) {
        var1 <- variables[i]
        var2 <- variables[j]
        cor_value <- cor(result_testDV[k == k_value & r == r_value & Variable == var1]$t_test,
                         result_testDV[k == k_value & r == r_value & Variable == var2]$t_test)
        cor_testDV <- rbind(cor_testDV, data.table(
          k = k_value,
          r = r_value,
          var1 = var1,
          var2 = var2,
          # Variables = paste(var1, var2, sep = "-"),
          correlation = cor_value
        ))
      }
    }
  }
}

# head the correlation results
view(cor_testDV)
# write_csv(cor_testDV, file = '20240709_cor_testDV.csv')

phacked_testDV$pm_minr <- NA
phacked_testDV$pm_maxr <- NA


# Iterate through the unique combinations of k and r
for (k_value in unique(phacked_testDV$k)) {
  for (r_value in unique(phacked_testDV$r)) {
    
    # Filter the correlation data based on the current k and r values
    filtered_data <- cor_testDV[k == k_value & r == r_value]
    
    # Calculate min and max correlations
    min_rho_testDV <- min(filtered_data$correlation)
    max_rho_testDV <- max(abs(filtered_data$correlation))  # Taking absolute max to find the highest magnitude correlation
    
    # Apply the adjusted_pm function to each row of phacked_testDV for the current k and r
    phacked_testDV$pm_minr[phacked_testDV$k == k_value & phacked_testDV$r == r_value] <- 
      sapply(phacked_testDV$t_test[phacked_testDV$k == k_value & phacked_testDV$r == r_value], 
             function(zo) adjusted_pm(zo, k_value, min_rho_testDV))
    
    phacked_testDV$pm_maxr[phacked_testDV$k == k_value & phacked_testDV$r == r_value] <- 
      sapply(phacked_testDV$t_test[phacked_testDV$k == k_value & phacked_testDV$r == r_value], 
             function(zo) adjusted_pm(zo, k_value, max_rho_testDV))
  }
}
view(phacked_testDV)


phacked_testDV <- phacked_testDV %>%
  mutate(sig_flg_pm_maxr = ifelse(pm_maxr <= 0.05, 1, 0)
         , sig_flg_change_maxr = ifelse(sig_ttest_5p != sig_flg_pm_maxr, 1, 0)
         , diff_po_pm_maxr = pm_maxr - pvalue_ttest
         , sig_flg_pm_minr = ifelse(pm_minr <= 0.05, 1, 0)
         , sig_flg_change_minr = ifelse(sig_ttest_5p != sig_flg_pm_minr, 1, 0)
         , diff_po_pm_minr = pm_minr - pvalue_ttest)

table(phacked_testDV$sig_flg_change_maxr)
table(phacked_testDV$sig_flg_change_minr)

phacked_testDV <- as.data.table(phacked_testDV)
head(phacked_testDV[phacked_testDV$sig_flg_change_maxr ==1], 10)
table(phacked_testDV$sig_ttest_5p, phacked_testDV$sig_flg_pm_maxr) 

table(phacked_testDV$sig_ttest_5p, phacked_testDV$sig_flg_pm_minr) 

mean(phacked_testDV$diff_po_pm_minr)
mean(phacked_testDV$diff_po_pm_maxr)

plot(phacked_testDV$pm_maxr, phacked_testDV$pvalue_ttest)
plot(density(result_testDV$t_test))
plot(density(phacked_testDV$pm_maxr))


head(phacked_testDV)
table(phacked_testDV$test_name)/sum(table(phacked_testDV$test_name))
table(phacked_testDV$k, phacked_testDV$sig_ttest_5p)
