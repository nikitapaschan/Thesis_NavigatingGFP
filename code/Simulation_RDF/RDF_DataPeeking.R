######################################################
# RDF Data peeking
#####################################################
setwd('/Users/nikitapaschan/Thesis_NavigatingGFP/')

# Define the sample size sequence
n <- c(20, 30, 40, 50, 100)

# Initialize data frames to store results
data_df <- data.table()
result_peeking <- data.table()

set.seed(123)
seeds <- sample(1:100000, 10000, replace = FALSE)

# Set seed for reproducibility
for (s in seeds) {
  
  continue_collecting <- 1  # Initialize the flag for data collection
  
  # Loop over the sample sizes
  for (i in n) {
    # if (continue_collecting == 1) {
    # Generate data
    set.seed(s)
    x <- rnorm(n = i, mean = 0.000001, sd = 1) 
    set.seed(s+1)
    y <- rnorm(n = i, mean = 0, sd = 1)
    
    # Combine data into a single data frame
    data_df <- rbind(data_df, data.frame(group = c(rep("x", i), rep("y", i)), value = c(x, y)))
    
    # Perform the two-sample t-test
    t_test <- t.test(value ~ group, data = data_df)
    
    # Determine if data collection should continue
    continue_collecting <- ifelse(t_test$p.value >= 0.05, 1, 0)
    
    # Store the results
    result_peeking <- rbind(result_peeking, data.frame(seed = s
                                                       , n = i
                                                       , t_test = t_test$statistic
                                                       , p_value = t_test$p.value
                                                       , sig_ttest_5p = ifelse(t_test$p.value < 0.05, 1, 0)
                                                       , continue_collecting = continue_collecting
    ))
    
    # If the p-value is significant, break the inner loop
    # if (continue_collecting == 0) {
    #   break
    # }
  }
}
# }

# Print the final results
head(result_peeking)
table(result_peeking$continue_collecting)
table(result_peeking$n, result_peeking$continue_collecting)
table(result_peeking$n)
# Correlation 
cor_peeking <- data.frame()

# Get unique n values
n_values <- unique(result_peeking$n)

# Calculate correlations for unique combinations
#NP: funktioniert nicht, da dimensionen 
for (i in 1:(length(n_values) - 1)) {
  for (j in (i + 1):length(n_values)) {
    n_value1 <- n_values[i]
    n_value2 <- n_values[j]
    #
    cor_value <- cor(result_peeking[n == n_value1]$t_test,
                     result_peeking[n == n_value2]$t_test, use = 'complete.obs' )
    #
    cor_peeking <- rbind(cor_peeking, data.table(
      n_value1 = n_value1,
      n_value2 = n_value2,
      correlation = cor_value
    ))
  }
}
#
view(cor_peeking)  


#save prelimary results
# write_csv(cor_peek, file = '20240712_cor_peek.csv')

boxplot(cor_peek$correlation ~ cor_peek$n_value1)

#example adjusted p-value
adjusted_pm(zo = 5.248976, m = 5, rho = 0.2306337)