######################################################
# RDF Opportunistic variable transformation
#####################################################
setwd('/Users/nikitapaschan/Thesis_NavigatingGFP/')

# Function to apply transformations and fit models
qrp_transformation <- function(n, seed) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Generate normal data
  data_trans <- data.frame(
    x = rnorm(n, mean = 0, sd = 1),
    y = rnorm(n, mean = 0, sd = 1)
  )
  
  # Handle cases where transformations might fail
  epsilon <- .Machine$double.eps
  
  # Apply transformations with adjustments
  data_trans <- data_trans %>%
    mutate(
      x_log = log(x + 5), #to avoid nans 
      y_log = log(y + 5),  
      x_reciprocal = 1 / (x + epsilon),# to avoid division by zero
      y_reciprocal = 1 / (y + epsilon),
      x_sqrt = sqrt(abs(x)),
      y_sqrt = sqrt(abs(y))
    )
  
  load('20240808_np_qrp.RData')
  
  # Function to extract coefficients and p-values
  extract_coefs <- function(model) {
    coefs <- summary(model)$coefficients
    data.frame(
      beta = coefs[2, "Estimate"],
      t_value = coefs[2, "t value"],
      p_value = coefs[2, "Pr(>|t|)"]
    )
  }
  
  # Create linear models
  lm_notrans <- lm(y ~ x, data = data_trans)
  lm_logx <- lm(y ~ x_log, data = data_trans)
  lm_logy <- lm(y_log ~ x, data = data_trans)
  lm_logxy <- lm(y_log ~ x_log, data = data_trans)
  
  lm_recix <- lm(y ~ x_reciprocal, data = data_trans)
  lm_reciy <- lm(y_reciprocal ~ x, data = data_trans)
  lm_recixy <- lm(y_reciprocal ~ x_reciprocal, data = data_trans)
  
  lm_sqrtx <- lm(y ~ x_sqrt, data = data_trans)
  lm_sqrty <- lm(y_sqrt ~ x, data = data_trans)
  lm_sqrtxy <- lm(y_sqrt ~ x_sqrt, data = data_trans)
  
  # Extract results for each model
  results <- list(
    extract_coefs(lm_notrans),
    extract_coefs(lm_logx),
    extract_coefs(lm_logy),
    extract_coefs(lm_logxy),
    extract_coefs(lm_recix),
    extract_coefs(lm_reciy),
    extract_coefs(lm_recixy),
    extract_coefs(lm_sqrtx),
    extract_coefs(lm_sqrty),
    extract_coefs(lm_sqrtxy)
  )
  
  # Perform Shapiro-Wilk test on residuals of lm_notrans
  shapiro_test <- shapiro.test(residuals(lm_notrans))
  shapiro_flag <- ifelse(shapiro_test$p.value < 0.05, 1, 0)
  
  # Back-transform betas where y was transformed
  beta_logy_backtrans <- exp(results[[3]]$beta)
  beta_logxy_backtrans <- exp(results[[4]]$beta)
  beta_reciy_backtrans <- 1 / results[[6]]$beta
  beta_recixy_backtrans <- 1 / results[[7]]$beta
  beta_sqrty_backtrans <- results[[9]]$beta^2
  beta_sqrtxy_backtrans <- results[[10]]$beta^2
  
  # Combine results into a single data frame row
  result_row <- data.frame(
    iteration = i,  # Iteration number
    n = n,
    beta_notrans = round(results[[1]]$beta, 5),
    tvalue_notrans = round(results[[1]]$t_value, 5),
    pvalue_notrans = round(results[[1]]$p_value, 5),
    sig_5p_notrans = ifelse(results[[1]]$p_value <= 0.05, 1, 0),
    shapiro_flag = shapiro_flag,
    
    # Log transformation
    beta_logx = round(results[[2]]$beta, 5),
    tvalue_logx = round(results[[2]]$t_value, 5),
    pvalue_logx = round(results[[2]]$p_value, 5),
    sig_5p_logx = ifelse(results[[2]]$p_value <= 0.05, 1, 0),
    beta_logy = round(results[[3]]$beta, 5),
    beta_logy_backtrans = round(beta_logy_backtrans, 5),  
    tvalue_logy = round(results[[3]]$t_value, 5),
    pvalue_logy = round(results[[3]]$p_value, 5),
    sig_5p_logy = ifelse(results[[3]]$p_value <= 0.05, 1, 0),
    beta_logxy = round(results[[4]]$beta, 5), 
    beta_logxy_backtrans = round(beta_logxy_backtrans, 5),  
    tvalue_logxy = round(results[[4]]$t_value, 5),
    pvalue_logxy = round(results[[4]]$p_value, 5),
    sig_5p_logxy = ifelse(results[[4]]$p_value <= 0.05, 1, 0),
    
    # Reciprocal transformation
    beta_recix = round(results[[5]]$beta, 5), 
    tvalue_recix = round(results[[5]]$t_value, 5),
    pvalue_recix = round(results[[5]]$p_value, 5), 
    sig_5p_recix = ifelse(results[[5]]$p_value <= 0.05, 1, 0),
    beta_reciy = round(results[[6]]$beta, 5),
    beta_reciy_backtrans = round(beta_reciy_backtrans, 5),  
    tvalue_reciy = round(results[[6]]$t_value, 5),
    pvalue_reciy = round(results[[6]]$p_value, 5),
    sig_5p_reciy = ifelse(results[[6]]$p_value <= 0.05, 1, 0),
    beta_recixy = round(results[[7]]$beta, 5),
    beta_recixy_backtrans = round(beta_recixy_backtrans, 5),  
    tvalue_recixy = round(results[[7]]$t_value, 5),
    pvalue_recixy = round(results[[7]]$p_value, 5),
    sig_5p_recixy = ifelse(results[[7]]$p_value <= 0.05, 1, 0),
    
    # Square root transformation
    beta_sqrtx = round(results[[8]]$beta, 5),
    tvalue_sqrtx = round(results[[8]]$t_value, 5),
    pvalue_sqrtx = round(results[[8]]$p_value, 5),
    sig_5p_sqrtx = ifelse(results[[8]]$p_value <= 0.05, 1, 0),
    beta_sqrty = round(results[[9]]$beta, 5),
    beta_sqrty_backtrans = round(beta_sqrty_backtrans, 5),  
    tvalue_sqrty = round(results[[9]]$t_value, 5),
    pvalue_sqrty = round(results[[9]]$p_value, 5),
    sig_5p_sqrty = ifelse(results[[9]]$p_value <= 0.05, 1, 0),
    beta_sqrtxy = round(results[[10]]$beta, 5),
    beta_sqrtxy_backtrans = round(beta_sqrtxy_backtrans, 5),  
    tvalue_sqrtxy = round(results[[10]]$t_value, 5),
    pvalue_sqrtxy = round(results[[10]]$p_value, 5),
    sig_5p_sqrtxy = ifelse(results[[10]]$p_value <= 0.05, 1, 0)
  )
  
  return(result_row)
}


# Number of iterations and sample size
niter <- 10000
n <- 100
set.seed(42)
seeds_trans <- sample(1:1000000, niter, replace = FALSE)  # Different seed for each iteration

results_trans <- do.call(rbind, pblapply(seeds_trans, function(seed) qrp_transformation(n, seed), cl = detectCores() - 1))


# head the combined results
head(results_trans)
results_trans <- as.data.table(results_trans)
colnames(results_trans)

# Find values of selective reporting 
results_trans_betas <- results_trans[, .(tvalue_notrans, tvalue_logx, tvalue_logy, tvalue_logxy, 
                                         tvalue_recix, tvalue_reciy, tvalue_recixy, 
                                         tvalue_sqrtx, tvalue_sqrty, tvalue_sqrtxy)]

trans_names <- colnames(results_trans_betas)

# Initialize
phacked_trans <- data.table()

# Loop through each row of the data
for (i in 1:nrow(results_trans_betas)) {
  # Find the index of the column with the maximum absolute t-value
  position_max_beta <- which.max(abs(as.numeric(results_trans_betas[i, ])))
  
  # Extract test name and zo value
  trans_name_beta <- trans_names[[position_max_beta]]
  zo_value_beta <- results_trans_betas[i, ..position_max_beta][[1]]
  
  # Append to the final data table
  phacked_trans <- rbind(phacked_trans, data.table(
    position_beta = position_max_beta,
    trans_name_beta = trans_name_beta,
    zo_test = zo_value_beta
  ), use.names = TRUE, fill = TRUE)
}

table(phacked_trans$trans_name_beta) # tvalue_reciy is the most common transformation

# Correlation matrix
results_trans_betas
correlation_matrix_trans <- cor(results_trans_betas)
correlation_matrix_trans
mean(correlation_matrix_trans)

# Create the correlation plot
cor_matrix_trans <- ggcorrplot(
  correlation_matrix_trans, 
  hc.order = TRUE, 
  type = "lower",
  lab = TRUE,
  lab_size = 3.5,  
  colors = brewer.pal(n = 3, name = "PRGn")
) + 
  ggtitle("Correlations of Regression Coefficients by Variable Transformation") +  
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 16),  
    plot.title = element_text(size = 18, hjust = 0.5, face = 'bold'), 
    legend.text = element_text(size = 13),  
    legend.title = element_text(size = 15), 
    plot.margin = margin(13, 13, 13, 13)  
  ) 

# Create text box plot
textbox_cor_matrix_trans <- ggdraw() +
  draw_text("Short names:\ntvalue = T-test of coefficient \nx/ y/ xy = Transformed variable \nnotrans = No transformation \nsqrt = Square root \nreci = Reciprocal\nlog = Logarithm", 
            x = 0.9, y = -0.9, hjust = 0.8, vjust = -1, size = 14)


# Combine the correlation plot and the text box plot
cor_matrix_trans <- plot_grid(cor_matrix_trans, textbox_cor_matrix_trans, ncol = 1, rel_heights = c(2, 0.5))

cor_matrix_trans

# ggsave(cor_matrix_trans, 
       # file = 'plots/cor_matrix/cor_matrix_trans.pdf', 
       # width = 10, height = 8, units = "in")



# ggsave(cor_matrix_trans, 
#        file = 'plots/cor_matrix/cor_matrix_trans.pdf', 
#        width = 10, height = 8, units = "in")

head(phacked_trans)

# adjusted p-value: 
m_trans <- 10
min_rho_trans <- min(correlation_matrix_trans)
max_rho_trans <-find_max_correlation(correlation_matrix_trans) 

# Apply adjusted_pm 
phacked_trans$pm_minr <- sapply(phacked_trans$zo_test, function(zo) adjusted_pm(zo, m_trans, min_rho_trans))
phacked_trans$pm_maxr <- sapply(phacked_trans$zo_test, function(zo) adjusted_pm(zo, m_trans, max_rho_trans))


phacked_trans <- phacked_trans %>%
  mutate(sig_flg_pm_maxr = ifelse(pm_maxr <= 0.05, 1, 0)
         , sig_flg_change_maxr = ifelse(sig_flg_trans != sig_flg_pm_maxr, 1, 0)
         , diff_po_pm_maxr = pm_maxr - po_trans
         , sig_flg_pm_minr = ifelse(pm_minr <= 0.05, 1, 0)
         , sig_flg_change_minr = ifelse(sig_flg_trans != sig_flg_pm_minr, 1, 0)
         , diff_po_pm_minr = pm_minr - po_trans)

phacked_trans <- as.data.table(phacked_trans)
head(phacked_trans)
table(phacked_trans$test_name)/sum(table(phacked_trans$test_name))

# save(results_trans, file = '20240902_np_rdf_transformation.RData')
# load('20240902_np_rdf_transformation.RData')

