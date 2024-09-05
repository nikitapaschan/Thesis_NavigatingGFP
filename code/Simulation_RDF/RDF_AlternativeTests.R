######################################################
# RDF Alternative hypothesis tests
#####################################################
setwd('/Users/nikitapaschan/Thesis_NavigatingGFP/')

set.seed(34)
n <- 100

# Numerical data, normally distributed
x <- matrix(rnorm(n*10000, 1), n)

# Vector for grouping 
g1 <- sample(c(rep(0, n/2), rep(1, n/2)), n, replace = FALSE)

# Append data set 
data_testing <- as.data.table(cbind(x, g1))

# Descriptive analysis
describeBy(data_testing$V1, data_testing$g1)

# Exclude the grouping column for testing
# data_test <- data_testing[, !("g1")]

# Create an empty data frame to store the results
result_testing <- data.table()

# Iterate through each column and perform tests
for (col in names(data_testing[, !("g1")])) {
  # Z-test
  ztest <- z.test(data_testing[g1==1, get(col)], data_testing[g1==0, get(col)],
                  sigma.x = var(data_testing[g1==1, get(col)]),
                  sigma.y = var(data_testing[g1==0, get(col)]))
  
  # T-test
  ttest <- t.test(data_testing[g1 == 1, get(col)], data_testing[g1 == 0, get(col)])
  
  # Mann-Whitney U test
  formula <- as.formula(paste(col, "~ as.factor(g1)"))
  mwu <- coin::wilcox_test(formula, data = data_testing)
  mwu_test <- statistic(mwu, "standardized")
  pvalue_mwu <- pvalue(mwu)
  
  # Capture results in data set
  result_testing <- rbind(result_testing, 
                          data.table(Variable = col, 
                                     z_test = ztest$statistic, 
                                     pvalue_ztest = ztest$p.value, 
                                     sig_ztest_5p = ifelse(ztest$p.value < 0.05, 1, 0), 
                                     t_test = ttest$statistic, 
                                     pvalue_ttest = ttest$p.value, 
                                     sig_ttest_5p = ifelse(ttest$p.value < 0.05, 1, 0), 
                                     mwu_test = mwu_test, 
                                     pvalue_mwu = pvalue_mwu, 
                                     sig_mwu_5p = ifelse(pvalue_mwu < 0.05, 1, 0)))
}
beep(2)
# Print the result
head(result_testing)
colnames(result_testing)[8] <- 'mwu_test'

## "P-Hacked"/ Min(p) test results
# create new data set
result_testing_teststat <- result_testing %>%
  select(z_test, t_test, mwu_test)
head(result_testing_teststat)
#colnames(result_testing)

# Get column names
test_names <- colnames(result_testing_teststat)

# Initialize empty data frame
result_testing_teststat <- as.data.frame(result_testing_teststat)
phacked_testing <- data.frame()

# Loop through each row of the data
for (i in 1:nrow(result_testing_teststat)) {
  # Find the index of the column with the maximum absolute value (vectorized)
  position_max <- which.max(abs(result_testing_teststat[i, ]))
  
  # Extract test name and zo value using vectorization
  test_name <- test_names[position_max]
  zo_value <- result_testing_teststat[i, position_max]
  
  # Append the result to the final data frame (efficiently)
  phacked_testing <- rbind(phacked_testing, data.frame(position = position_max
                                                       , test_name = test_name
                                                       , zo_test = zo_value)
  )
}

head(phacked_testing)
# Reset row names to a simple enumeration
rownames(phacked_testing) <- seq_len(nrow(phacked_testing))

## Append test statistics to phacked_testing results
# Initialize the vector to store the previous column values
po_testing <- rep(NA, nrow(result_testing))
sig_flg_testing <- rep(NA, nrow(result_testing))
result_testing <- as.data.frame(result_testing)
# Loop through each row of the table
for (i in 1:nrow(result_testing)) {
  col_index <- which(result_testing[i, ] == phacked_testing$zo_test[i])
  if (length(col_index) > 0 && col_index > 1) {
    po_testing[i] <- result_testing[i, col_index + 1]
    sig_flg_testing[i] <- result_testing[i, col_index + 2]
  }
}

# Combine the original vector with the previous column values
phacked_testing <- cbind(phacked_testing, po_testing, sig_flg_testing)


# Print the result
view(phacked_testing)
view(result_testing)

phacked_testing <- as.data.table(phacked_testing)
table(phacked_testing$test_name, phacked_testing$sig_flg_testing) # 179 number of significant mwu maxZ tests 
phacked_testing[position == 3 & sig_flg_testing == 1] 

table(
  # result_testing$sig_ztest_5p,
  result_testing$sig_ttest_5p
  , result_testing$sig_mwu_5p
)

# frequency table: shows slight increase in type 1 error for phacked tests 
table(result_testing$sig_mwu_5p)/sum(table(result_testing$sig_mwu_5p))
table(result_testing$sig_ztest_5p)/sum(table(result_testing$sig_ztest_5p))
table(result_testing$sig_ttest_5p)/sum(table(result_testing$sig_ttest_5p))
table(phacked_testing$sig_flg_testing)/sum(table(phacked_testing$sig_flg_testing))

### correlation matrix 
cor_testing <- cor(result_testing_teststat)
cor_testing


cor_matrix_testing <- ggcorrplot(cor_testing, hc.order = TRUE, type = "lower",
                                 lab = TRUE
                                 , colors = c("#6D9EC1", "white", "#E46726")
)

ggsave('cor_matrix_testing.pdf', width = 14, height = 10, units = "in")


# adjusted p-value: 
m_testing <- 3
min_rho_testing <- min(cor_testing) #NP: or take abs(rho)
max_rho_testing <- find_max_correlation(cor_testing)

# Apply adjusted_pm 
phacked_testing$pm_minr <- sapply(phacked_testing$zo_test, function(zo) adjusted_pm(zo, m_testing, min_rho_testing)) 
phacked_testing$pm_maxr <- sapply(phacked_testing$zo_test, function(zo) adjusted_pm(zo, m_testing, max_rho_testing))

view(phacked_testing)

phacked_testing <- phacked_testing %>%
  mutate(sig_flg_pm_maxr = ifelse(pm_maxr <= 0.05, 1, 0)
         , sig_flg_change = ifelse(sig_flg_testing != sig_flg_pm_maxr, 1, 0)
         , diff_po_pm_maxr = pm_maxr - po_testing
  )
table(phacked_testing$sig_flg_change)
head(phacked_testing[sig_flg_change ==1], 20)
table(phacked_testing$sig_flg_testing, phacked_testing$sig_flg_pm_maxr) # 10 instances where pm is smaller than po with a mean of -0.0005441522 -> almost zero (mostly for very small po)
table(phacked_testing$sig_flg_testing, phacked_testing$test_name)
plot(density(phacked_testing$diff_po_pm_maxr))

plot(phacked_testing$pm_maxr, phacked_testing$po_testing)
plot(density(result_testing$z_test))
plot(density(phacked_testing$pm_maxr))