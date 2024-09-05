######################################################
# RDF Selective reporting of independent variables
#####################################################
setwd('/Users/nikitapaschan/Thesis_NavigatingGFP/')

generate_correlated_data <- function(n, m, correlation, seed) {
  set.seed(seed)
  sigma <- matrix(correlation, nrow = m, ncol = m)
  diag(sigma) <- 1
  data <- mvrnorm(n, mu = rep(0, m), Sigma = sigma)
  data_df <- as.data.frame(data)
  return(data_df)
}

trial_data_regIV <- generate_correlated_data(100, 5, 0.3, 123)
trial_data_regIV$y <- rnorm(100)

summary(lm(y ~ V1, data = trial_data_regIV))

# Set parameters
n_values <- c(100)
m_values <- c(3, 5, 10)
niter <- 10000
correlation_values <- c(0, 0.3, 0.8)

# Generate seeds for reproducibility
set.seed(123)
seeds <- sample(1:1000000, niter * length(n_values) * length(m_values) * length(correlation_values), replace = FALSE)
seeds[10000]
tail(seeds, 5)
# Initialize 
results_regIV <- data.frame()

# Index for accessing the correct seed
seed_index <- 1

# Total number of iterations for progress bar
total_iterations <- length(n_values) * length(m_values) * length(correlation_values) * niter
progress <- 0

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

# Nested for loops to iterate over all parameter combinations
for (n in n_values) {
  for (m in m_values) {
    for (corr in correlation_values) {
      for (iteration in 1:niter) {
        
        # Set the seed for reproducibility
        set.seed(seeds[seed_index])
        seed_index <- seed_index + 1
        
        # Generate response vector y
        regIV_y <- rnorm(n)
        
        # Generate predictor matrix x with specified n, m, and correlation
        regIV_x <- generate_correlated_data(n, m, corr, seeds[seed_index])
        
        # Loop through each predictor (column of x)
        for (i in 1:m) {
          
          # Fit the linear model
          summary_model <- summary(lm(y ~ x, data = data.frame(y = regIV_y, x = regIV_x[, i])))
          
          # Extract beta value and p-value
          beta_value <- coef(summary_model)[2, "Estimate"]
          t_value <- coef(summary_model)[2, "t value"]
          p_value <- coef(summary_model)[2, "Pr(>|t|)"]
          sig_flag <- ifelse(p_value < 0.05, 1, 0)
          
          # Store the results
          results_regIV <- rbind(results_regIV, data.frame(
            iteration = iteration,
            n = n,
            m = m,
            correlation = corr,
            column_used = paste0("x", i),
            beta = beta_value,
            tvalue = t_value,
            p_value = p_value,
            sig_flag_5p = sig_flag
          ))
        }
        
        # Update progress bar after each iteration
        progress <- progress + 1
        setTxtProgressBar(pb, progress)
      }
    }
  }
}

# Close the progress bar
close(pb)

tail(results_regIV)


# Modify the naming convention
results_regIV <- results_regIV %>%
  mutate(name = paste(n, m, correlation, sep = "_"))

tail(results_regIV, 20)

# MaxZ
phacked_regIV <- results_regIV %>%
  group_by(iteration, n, m, correlation) %>%
  filter(abs(tvalue) == max(abs(tvalue))) %>%
  ungroup()

head(phacked_regIV)


# Function to calculate the correlation matrix for each parameter combination
calculate_correlations <- function(df) {
  beta_matrix <- df %>%
    select(iteration, column_used, tvalue) %>%
    pivot_wider(names_from = column_used, values_from = tvalue, names_prefix = "tvalue_")
  
  cor_matrix <- cor(beta_matrix[, -1], use = "complete.obs")
  return(cor_matrix)
}

table(phacked_regIV$column_used)

# Split results_regIV by parameter combinations and calculate correlation matrices
results_split <- split(results_regIV, list(results_regIV$n, results_regIV$m, results_regIV$correlation))

cor_matrices <- lapply(results_split, calculate_correlations)
cor_matrices[[5]]

# Create the correlation plot
cor_matrix_regIV <- ggcorrplot(
  cor_matrices$'100.5.0.3', 
  hc.order = TRUE, 
  type = "lower",
  lab = TRUE,
  lab_size = 8,  # Increase the size of the text inside the cells
  colors = brewer.pal(n = 3, name = "PRGn")
) + 
  ggtitle("Correlations of Regression Coefficients between Independent Variables") +  
  theme(
    axis.text.x = element_text(size = 19, angle = 45, hjust = 1),  # Increase font size and rotate x-axis labels
    axis.text.y = element_text(size = 19),  # Increase font size for y-axis labels
    plot.title = element_text(size = 20, hjust = 0.5, face = 'bold'),  # Increase font size and center the title
    legend.text = element_text(size = 15),  # Increase font size of the legend text
    legend.title = element_text(size = 17),  # Increase font size of the legend title
    plot.margin = margin(10, 10, 10, 10)  # Add some margin around the plot
  ) 


cor_matrix_regIV

ggsave(cor_matrix_regIV, 
       file = '/Users/nikitapaschan/Master-Thesis/plots/cor_matrix/cor_matrix_regIV.pdf', 
       width = 10, height = 8, units = "in")


# Function to calculate mean correlation excluding the diagonal
mean_correlation <- function(cor_matrix) {
  # Extract the upper triangle without the diagonal
  upper_triangle <- cor_matrix[upper.tri(cor_matrix)]
  
  # Calculate and return the mean of the upper triangle correlations
  return(mean(upper_triangle))
}

# Initialize
mean_correlations_list <- list()

# Generate mean correlation for each parameter combination
for (name in names(cor_matrices)) {
  # Extract the components of the name
  components <- unlist(strsplit(name, "\\."))
  
  # mean correlation for the current matrix
  mean_corr <- mean_correlation(cor_matrices[[name]])
  
  # Append the results
  mean_correlations_list[[name]] <- c(
    mean = mean_corr,
    n = as.numeric(components[1]),
    m = as.numeric(components[2])
    ,correlation = as.numeric(components[4])
  )
}

mean_correlations_df <- do.call(rbind, mean_correlations_list)
mean_correlations_df <- as.data.frame(mean_correlations_df)

# Print the final data frame
mean_correlations_df
# Adjust correlation values 
mean_correlations_df$correlation[is.na(mean_correlations_df$correlation)] <- 0.0
mean_correlations_df$correlation[mean_correlations_df$correlation == 3] <- 0.3
mean_correlations_df$correlation[mean_correlations_df$correlation == 8] <- 0.8

boxplot(mean_correlations_df$mean ~ mean_correlations_df$correlation, xlab = "Correlation Value", ylab = "Mean Correlation", main = "Mean Correlation by Correlation Value")


# Adjusted p-value  
phacked_regIV$pm <- pbmapply(function(zo, m, rho) adjusted_pm(zo, m, rho), 
                                phacked_regIV$tvalue, phacked_regIV$m, phacked_regIV$correlation)

phacked_regIV <- phacked_regIV %>%
  mutate(sig_flg_pm = ifelse(pm <= 0.05, 1, 0)
         , sig_flg_change = ifelse(sig_flag_5p != sig_flg_pm, 1, 0)
         , diff_po_pm_maxr = pm - p_value)

table(phacked_regIV$sig_flg_change)
table(phacked_regIV$sig_flag_5p)/nrow(phacked_regIV)
head(phacked_regIV)
