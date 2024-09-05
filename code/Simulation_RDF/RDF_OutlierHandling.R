######################################################
# RDF Outlier handling
#####################################################
setwd('/Users/nikitapaschan/Thesis_NavigatingGFP/')

set.seed(123)
n <- 100
niter <- 10000

# Generating the data
data_outlier <- data.table(matrix(rcauchy(n * niter, location = 0, scale = 1), nrow = n, ncol = niter))

# Adding grouping and response variable
data_outlier[, `:=`(g1 = sample(c(rep(0, n/2), rep(1, n/2)), n, replace = FALSE),
                    y = rnorm(n, mean = 0, sd = 1))]

# Initialize
result_outlier <- data.table()

handle_outliers_and_tests <- function(col) {
  # Identify outliers using boxplot method
  out_boxplot <- boxplot.stats(data_outlier[[col]])$out
  outlier_flag_col_box <- paste0('outl_box_', col)
  data_outlier[, (outlier_flag_col_box) := fifelse(get(col) %in% out_boxplot, 1, 0)]
  
  # Winsorize the data
  data_outlier_winsorized <- copy(data_outlier)
  data_outlier_winsorized[[col]] <- Winsorize(data_outlier[[col]], probs = c(0.05, 0.95))
  
  # MAD outlier detection
  mad_value <- mad(data_outlier[[col]])
  median_value <- median(data_outlier[[col]])
  out_mad <- which(abs(data_outlier[[col]] - median_value) > 3 * mad_value)
  outlier_flag_col_mad <- paste0('outl_mad_', col)
  data_outlier[, (outlier_flag_col_mad) := fifelse(seq_len(nrow(data_outlier)) %in% out_mad, 1, 0)]
  
  # Rosner test for outliers
  rosner_test <- rosnerTest(data_outlier[[col]], k = 5)
  out_rosner <- rosner_test$all.stats$Obs.Num[rosner_test$all.stats$Outlier]
  outlier_flag_col_rosner <- paste0('outl_rosner_', col)
  data_outlier[, (outlier_flag_col_rosner) := fifelse(seq_len(nrow(data_outlier)) %in% out_rosner, 1, 0)]
  
  # T-tests
  formula_t <- as.formula(paste(col, "~ g1"))
  ttest <- t.test(formula_t, data = data_outlier)
  ttest_box <- t.test(formula_t, data = data_outlier[get(outlier_flag_col_box) == 0])
  ttest_winsor <- t.test(formula_t, data = data_outlier_winsorized)
  ttest_mad <- t.test(formula_t, data = data_outlier[get(outlier_flag_col_mad) == 0])
  ttest_rosner <- t.test(formula_t, data = data_outlier[get(outlier_flag_col_rosner) == 0])
  
  # Linear models
  formula_lm <- as.formula(paste("y ~ ", col))
  lm_all <- lm(formula_lm, data = data_outlier)
  lm_box <- lm(formula_lm, data = data_outlier[get(outlier_flag_col_box) == 0])
  lm_winsor <- lm(formula_lm, data = data_outlier_winsorized)
  lm_mad <- lm(formula_lm, data = data_outlier[get(outlier_flag_col_mad) == 0])
  lm_rosner <- lm(formula_lm, data = data_outlier[get(outlier_flag_col_rosner) == 0])
  
  # Cook's distance for outliers
  cooks_d <- cooks.distance(lm_all)
  out_cooks <- which(cooks_d > (4 / n))
  outlier_flag_col_cooks <- paste0('outl_cooks_', col)
  data_outlier[, (outlier_flag_col_cooks) := fifelse(seq_len(nrow(data_outlier)) %in% out_cooks, 1, 0)]
  lm_cooks <- lm(formula_lm, data = data_outlier[get(outlier_flag_col_cooks) == 0])
  
  # Collect results
  return(data.table(Variable = col,
                    # No outlier handling
                    t_test_all = ttest$statistic, 
                    pvalue_ttest_all = ttest$p.value, 
                    sig_ttest_all_5p = as.integer(ttest$p.value <= 0.05),
                    # Boxplot outlier handling
                    t_test_box = ttest_box$statistic, 
                    pvalue_ttest_box = ttest_box$p.value, 
                    sig_ttest_box_5p = as.integer(ttest_box$p.value <= 0.05),
                    cnt_outlier_box = sum(data_outlier[[outlier_flag_col_box]] == 1),
                    # Winsorizing
                    t_test_winsor = ttest_winsor$statistic, 
                    pvalue_ttest_winsor = ttest_winsor$p.value, 
                    sig_ttest_winsor_5p = as.integer(ttest_winsor$p.value <= 0.05),
                    # MAD
                    t_test_mad = ttest_mad$statistic, 
                    pvalue_ttest_mad = ttest_mad$p.value, 
                    sig_ttest_mad_5p = as.integer(ttest_mad$p.value <= 0.05),
                    cnt_outlier_mad = sum(data_outlier[[outlier_flag_col_mad]] == 1),
                    # Rosner test
                    t_test_rosner = ttest_rosner$statistic, 
                    pvalue_ttest_rosner = ttest_rosner$p.value, 
                    sig_ttest_rosner_5p = as.integer(ttest_rosner$p.value <= 0.05),
                    cnt_outlier_rosner = sum(data_outlier[[outlier_flag_col_rosner]] == 1), 
                    ### Outlier Handling on Regression ###
                    # No outlier handling
                    lm_all = summary(lm_all)$coefficients[2, 1], 
                    pvalue_lm_all = summary(lm_all)$coefficients[2, 4], 
                    sig_lm_all_5p = as.integer(summary(lm_all)$coefficients[2, 4] <= 0.05),
                    # Boxplot outlier handling
                    lm_box = summary(lm_box)$coefficients[2, 1],
                    pvalue_lm_box = summary(lm_box)$coefficients[2, 4],
                    sig_lm_box_5p = as.integer(summary(lm_box)$coefficients[2, 4] <= 0.05),
                    # Winsorizing
                    lm_winsor = summary(lm_winsor)$coefficients[2, 1],
                    pvalue_lm_winsor = summary(lm_winsor)$coefficients[2, 4],
                    sig_lm_winsor_5p = as.integer(summary(lm_winsor)$coefficients[2, 4] <= 0.05),
                    # MAD
                    lm_mad = summary(lm_mad)$coefficients[2, 1],
                    pvalue_lm_mad = summary(lm_mad)$coefficients[2, 4],
                    sig_lm_mad_5p = as.integer(summary(lm_mad)$coefficients[2, 4] <= 0.05),
                    # Rosner test
                    lm_rosner = summary(lm_rosner)$coefficients[2, 1],
                    pvalue_lm_rosner = summary(lm_rosner)$coefficients[2, 4],
                    sig_lm_rosner_5p = as.integer(summary(lm_rosner)$coefficients[2, 4] <= 0.05),
                    # Cook's distance
                    lm_cooks = summary(lm_cooks)$coefficients[2, 1],
                    pvalue_lm_cooks = summary(lm_cooks)$coefficients[2, 4],
                    sig_lm_cooks_5p = as.integer(summary(lm_cooks)$coefficients[2, 4] <= 0.05),
                    cnt_outlier_cooks = sum(data_outlier[[outlier_flag_col_cooks]] == 1)))
}

# Process in Chunks to Manage Memory
chunk_size <- 500 # Adjust as needed
n_chunks <- ceiling(niter / chunk_size)

for (i in 1:n_chunks) {
  cat("Processing chunk", i, "of", n_chunks, "\n")
  start_col <- (i - 1) * chunk_size + 1
  end_col <- min(i * chunk_size, niter)
  
  chunk_results <- rbindlist(lapply(names(data_outlier)[start_col:end_col], handle_outliers_and_tests))
  result_outlier <- rbind(result_outlier, chunk_results)
  
  # Run garbage collection to free up memory
  gc()
}


# head or Analyze the Results
head(result_outlier)
#save as rdata 
# save(result_outlier, file = "datasets/result_outlier.RData")

####################
#### maxZ / phacked 
######################
# Convert result_outlier to data.table for efficient operations
result_outlier <- as.data.table(result_outlier)

# Select the columns with t-test statistics
result_outlier_teststat <- result_outlier %>%
  select(t_test_all, t_test_box, t_test_winsor
         , t_test_mad, t_test_rosner)

# Get column names
test_names <- colnames(result_outlier_teststat)
result_outlier_teststat <- as.data.frame(result_outlier_teststat)

# Initialize empty data table
phacked_outlier_test <- data.table()

# Loop through each row of the data
for (i in 1:nrow(result_outlier_teststat)) {
  # Find the index of the column with the maximum absolute value (vectorized)
  position_max <- which.max(abs(result_outlier_teststat[i, ]))
  
  # Extract test name and zo value using vectorization
  test_name <- test_names[position_max]
  zo_value <- result_outlier_teststat[i, position_max]
  
  # Append the result to the final data table (efficiently)
  phacked_outlier_test <- rbind(phacked_outlier_test, data.table(position = position_max,
                                                       test_name = test_name
                                                       ,zo_test = zo_value
  ))
}

head(phacked_outlier_test)


table(phacked_outlier_test$test_name)

# Initialize 
po_outlier_test <- rep(NA, nrow(result_outlier))
sig_flg_outlier_test <- rep(NA, nrow(result_outlier))

result_outlier <- as.data.frame(result_outlier)
# Loop through each row of the table
for (i in 1:nrow(result_outlier)) {
  col_index <- which(result_outlier[i, ] == phacked_outlier_test$zo_test[i])
  if (length(col_index) > 0 && col_index[1] > 1) {
    po_outlier_test[i] <- result_outlier[i, col_index[1] + 1]
    sig_flg_outlier_test[i] <- result_outlier[i, col_index[1] + 2]
  }
}

phacked_outlier_test <- data.table(phacked_outlier_test, po_outlier_test = po_outlier_test, sig_flg_outlier_test = sig_flg_outlier_test)

head(phacked_outlier_test)
head(result_outlier)

# frequency tables for sig_flg 
table(result_outlier$sig_ttest_all_5p)/sum(table(result_outlier$sig_ttest_all_5p))
table(result_outlier$sig_ttest_winsor_5p)/sum(table(result_outlier$sig_ttest_winsor_5p))
table(result_outlier$sig_ttest_rosner_5p)/sum(table(result_outlier$sig_ttest_rosner_5p))
table(result_outlier$sig_ttest_box_5p)/sum(table(result_outlier$sig_ttest_box_5p))
table(phacked_outlier_test$sig_flg_outlier_test)/sum(table(phacked_outlier_test$sig_flg_outlier_test))

result_outlier_teststat <- result_outlier_teststat %>%
  rename(t_test_iqr = t_test_box)

names(result_outlier_teststat)
# Correlation matrix between the t test of the different outlier methods
correlation_matrix_outlier <- cor(result_outlier_teststat)

# Correlation plot
cor_matrix_outlier_ttest <- ggcorrplot(
  correlation_matrix_outlier, 
  hc.order = TRUE, 
  type = "lower",
  lab = TRUE,
  lab_size = 7, 
  colors = brewer.pal(n = 3, name = "PRGn")
) + 
  ggtitle("Correlations of T-Test Statistics by Outlier Handling Method") +  
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 18),  
    plot.title = element_text(size = 19, hjust = 0.5, face = 'bold'),  
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 16),  
    plot.margin = margin(10, 10, 10, 10)  
  ) 

# Create the text box plot
textbox_cor_matrix_outlier_ttest <- ggdraw() +
  draw_text("Short names:\nall = no outlier handling\nwinsor = Winsorizing\nmad = Median Absolute Deviation\niqr = Interquartile Range\nrosner = Rosner Test", 
            x = 0.9, y = -0.8, hjust = 0.9, vjust = -1, size = 15)
cor_matrix_outlier_ttest <- plot_grid(cor_matrix_outlier_ttest, textbox_cor_matrix_outlier_ttest, ncol = 1, rel_heights = c(2, 0.5))

cor_matrix_outlier_ttest

# ggsave(cor_matrix_outlier_ttest, 
       # file = 'plots/cor_matrix/cor_matrix_outlier_ttest.pdf', 
       # width = 10, height = 8, units = "in")


# adjusted p-value: 
m_outlier_test <- 5
min_rho_outlier_test <- min(correlation_matrix_outlier)
max_rho_outlier_test <-find_max_correlation(correlation_matrix_outlier) 

# Apply adjusted_pm 
phacked_outlier_test$pm_minr <- sapply(phacked_outlier_test$zo_test, function(zo) adjusted_pm(zo, m_outlier_test, min_rho_outlier_test))
phacked_outlier_test$pm_maxr <- sapply(phacked_outlier_test$zo_test, function(zo) adjusted_pm(zo, m_outlier_test, max_rho_outlier_test))

phacked_outlier_test <- phacked_outlier_test %>%
  mutate(sig_flg_pm_maxr = ifelse(pm_maxr <= 0.05, 1, 0)
         , sig_flg_change_maxr = ifelse(sig_flg_outlier_test != sig_flg_pm_maxr, 1, 0)
         , diff_po_pm_maxr = pm_maxr - po_outlier_test
         , sig_flg_pm_minr = ifelse(pm_minr <= 0.05, 1, 0)
         , sig_flg_change_minr = ifelse(sig_flg_outlier_test != sig_flg_pm_minr, 1, 0)
         , diff_po_pm_minr = pm_minr - po_outlier_test)

table(phacked_outlier_test$sig_flg_change_maxr)
table(phacked_outlier_test$sig_flg_change_minr)

phacked_outlier_test <- as.data.table(phacked_outlier_test)
head(phacked_outlier_test[phacked_outlier_test$sig_flg_change_maxr ==1], 10)
table(phacked_outlier_test$sig_flg_outlier_test, phacked_outlier_test$sig_flg_pm_maxr) 
table(phacked_outlier_test$sig_flg_outlier_test, phacked_outlier_test$sig_flg_pm_minr) 
plot(density(phacked_outlier_test$diff_po_pm_minr))
plot(density(phacked_outlier_test$diff_po_pm_maxr))

plot(phacked_outlier_test$pm_maxr, phacked_outlier_test$po_outlier_test)
plot(density(result_outlier$z_test))
plot(density(phacked_outlier_test$pm_maxr))


head(phacked_outlier_test)
table(phacked_outlier_test$test_name)/sum(table(phacked_outlier_test$test_name))
table(phacked_outlier_test$test_name, phacked_outlier_test$sig_flg_outlier_test)

#################################################
### outlier handling correlation for linear model
#################################################

# Select the columns with tvalue
result_outlier_lm <- result_outlier %>%
  select(lm_all, lm_box, lm_winsor
         , lm_mad, lm_rosner, lm_cooks)

# Get column names
lm_names <- colnames(result_outlier_lm)
result_outlier_lm <- as.data.frame(result_outlier_lm)

# Initialize empty data table
phacked_outlier_lm <- data.table()

# Loop through each row of the data
for (i in 1:nrow(result_outlier_lm)) {
  # Find the index of the column with the maximum absolute value (vectorized)
  position_max <- which.max(abs(result_outlier_lm[i, ]))
  
  # Extract test name and zo value using vectorization
  lm_name <- lm_names[position_max]
  zo_value <- result_outlier_lm[i, position_max]
  
  # Append the result to the final data table (efficiently)
  phacked_outlier_lm <- rbind(phacked_outlier_lm, data.table(position = position_max,
                                                             lm_name = lm_name
                                                             ,zo_lm = zo_value
  ))
}

head(phacked_outlier_lm)


hist(phacked_outlier_lm$position)
# Reset row names to a simple enumeration
# rownames(phacked_outlier_lm) <- seq_len(nrow(phacked_outlier_lm))

# Initialize the vector to store the previous column values
po_outlier_lm <- rep(NA, nrow(result_outlier))
sig_flg_outlier_lm <- rep(NA, nrow(result_outlier))

result_outlier <- as.data.frame(result_outlier)
# Loop through each row of the table
for (i in 1:nrow(result_outlier)) {
  # Find the column index of the current value in phacked_outlier_lm for current row
  col_index <- which(result_outlier[i, ] == phacked_outlier_lm$zo_lm[i])
  
  # Ensure that the column index is greater than 1
  if (length(col_index) > 0 && col_index[1] > 1) {
    po_outlier_lm[i] <- result_outlier[i, col_index[1] + 1]
    sig_flg_outlier_lm[i] <- result_outlier[i, col_index[1] + 2]
  }
}

# Combine results back into a data.table
phacked_outlier_lm <- data.table(phacked_outlier_lm, po_outlier_lm = po_outlier_lm, sig_flg_outlier_lm = sig_flg_outlier_lm)



# Print the head of the final data.table
head(phacked_outlier_lm)
head(result_outlier)
table(phacked_outlier_lm$lm_name)
table(phacked_outlier_lm$sig_flg_outlier_lm)


# frequency tables for sig_flg 
table(result_outlier$sig_ttest_all_5p)/sum(table(result_outlier$sig_ttest_all_5p))
table(result_outlier$sig_ttest_winsor_5p)/sum(table(result_outlier$sig_ttest_winsor_5p))
table(result_outlier$sig_ttest_rosner_5p)/sum(table(result_outlier$sig_ttest_rosner_5p))
table(result_outlier$sig_ttest_box_5p)/sum(table(result_outlier$sig_ttest_box_5p))
table(phacked_outlier_lm$sig_flg_outlier_lm)/sum(table(phacked_outlier_lm$sig_flg_outlier_lm))



# Correlation matrix between the t test of the different outlier methods
result_outlier_lm <- result_outlier_lm %>%
  rename(lm_iqr = lm_box)

# Calculate the correlation matrix
correlation_matrix_outlier_lm <- cor(result_outlier_lm)

# Create the correlation plot
cor_matrix_outlier_lm <- ggcorrplot(
  correlation_matrix_outlier_lm, 
  hc.order = TRUE, 
  type = "lower",
  lab = TRUE,
  lab_size = 7, 
  colors = brewer.pal(n = 3, name = "PRGn")
) + 
  ggtitle("Correlations of Regression Coefficients by Outlier Handling Method") +  
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 18),  
    plot.title = element_text(size = 20, hjust = 0.5, face = 'bold'),  
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 16),  
    plot.margin = margin(10, 10, 10, 10)  
  ) 

# Create the textbox plot
textbox_cor_matrix_outlier_lm <- ggdraw() +
  draw_text("Short names:\nall = no outlier handling\nwinsor = Winsorizing\nmad = Median Absolute Deviation\niqr = Interquartile Range\nrosner = Rosner Test\ncooks = Cook's Distance", 
            x = 0.9, y = -0.8, hjust = 0.9, vjust = -1, size = 14)


# Combine the correlation plot and the textbox plot
cor_matrix_outlier_lm <- plot_grid(cor_matrix_outlier_lm, textbox_cor_matrix_outlier_lm, ncol = 1, rel_heights = c(2, 0.5))

cor_matrix_outlier_lm

# ggsave(cor_matrix_outlier_lm, 
       # file = 'plots/cor_matrix/cor_matrix_outlier_lm.pdf', 
       # width = 10, height = 8, units = "in")


# adjusted p-value: 
m_outlier_lm <- 5
min_rho_outlier_lm <- min(correlation_matrix_outlier)
max_rho_outlier_lm <-find_max_correlation(correlation_matrix_outlier) 

# Apply adjusted_pm 
phacked_outlier_lm$pm_minr <- sapply(phacked_outlier_lm$zo_lm, function(zo) adjusted_pm(zo, m_outlier_lm, min_rho_outlier_lm))
phacked_outlier_lm$pm_maxr <- sapply(phacked_outlier_lm$zo_lm, function(zo) adjusted_pm(zo, m_outlier_lm, max_rho_outlier_lm))

phacked_outlier_lm <- phacked_outlier_lm %>%
  mutate(sig_flg_pm_maxr = ifelse(pm_maxr <= 0.05, 1, 0)
         , sig_flg_change_maxr = ifelse(sig_flg_outlier_lm != sig_flg_pm_maxr, 1, 0)
         , diff_po_pm_maxr = pm_maxr - po_outlier_lm
         , sig_flg_pm_minr = ifelse(pm_minr <= 0.05, 1, 0)
         , sig_flg_change_minr = ifelse(sig_flg_outlier_lm != sig_flg_pm_minr, 1, 0)
         , diff_po_pm_minr = pm_minr - po_outlier_lm)

table(phacked_outlier_lm$sig_flg_change_maxr)
table(phacked_outlier_lm$sig_flg_change_minr)

phacked_outlier_lm <- as.data.table(phacked_outlier_lm)
head(phacked_outlier_lm[phacked_outlier_lm$sig_flg_change_maxr ==1], 10)
table(phacked_outlier_lm$sig_flg_outlier_lm, phacked_outlier_lm$sig_flg_pm_maxr) 
table(phacked_outlier_lm$sig_flg_outlier_lm, phacked_outlier_lm$sig_flg_pm_minr) 
plot(density(phacked_outlier_lm$diff_po_pm_minr))
plot(density(phacked_outlier_lm$diff_po_pm_maxr))

plot(phacked_outlier_lm$pm_maxr, phacked_outlier_lm$po_outlier_lm)
plot(density(result_outlier$z_lm))
plot(density(phacked_outlier_lm$pm_maxr))


head(phacked_outlier_lm)
table(phacked_outlier_lm$lm_name)/sum(table(phacked_outlier_lm$lm_name))
table(phacked_outlier_lm$lm_name, phacked_outlier_lm$sig_flg_outlier_lm)


