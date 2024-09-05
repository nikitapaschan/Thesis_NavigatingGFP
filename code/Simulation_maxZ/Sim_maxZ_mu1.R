######################################################
# Simulation maxZ: mu = 1
#####################################################
setwd('/Users/nikitapaschan/Thesis_NavigatingGFP/')
## Parameters
m = c(2, 5, 10, 100, 1000)
rho = c(0, 0.3, 0.5, 0.7, 0.9)
niter = 10

data_phacking_signal_ew1 <- sim_datasets_phacking(m = m, rho = rho, mhu = 1)
# save(data_phacking_signal_ew1, file='datasets/20240627_data_phacking_signal_ew1.RData')

datatable_phacking_ew1 <- as.data.table(data_phacking_signal_ew1$data_all)
unique(datatable_phacking_ew1$df)


# save(datatable_phacking_ew1, 'datasets/20240627_data_phacking_signal_ew1.RData')
# load('datasets/20240627_data_phacking_signal_ew1.RData')


# Use mclapply to process each row in parallel
datatable_phacking_ew1_pm <- do.call(rbind, pblapply(1:nrow(datatable_phacking_ew1_pm), function(i) process_row_mu(datatable_phacking_ew1_pm[i, ]), cl = detectCores() - 1))
head(datatable_phacking_ew1_pm)

datatable_phacking_ew1_pm <- merge(datatable_phacking_ew1_pm, datatable_phacking_ew1_pm, by = c('df', 'zo'))
head(datatable_phacking_ew1_pm)

# save(datatable_phacking_ew1_pm, file = '/Users/nikitapaschan/Master-Thesis/datasets/20240903_datatable_phacking_ew1_pm.Rdata')

# Save the results
save(datatable_phacking_ew1_pm, file='datasets/20240708_datatable_phacking_ew1_pm.Rdata')

# Progress completion sound
beep(2)

# Return the modified dataset
head(datatable_phacking_ew1_pm)

# Merge pm to existing data frame
datatable_phacking_ew1_pm <- merge(datatable_phacking_ew1, datatable_phacking_ew1_pm, by = c('df', 'zo'))
save(datatable_phacking_ew1_pm, file='datasets/20240708_datatable_phacking_ew1_pm.Rdata')

# flags for magnitude of adjustment
datatable_phacking_ew1_pm <- datatable_phacking_ew1_pm %>%
  mutate(diff_pm_po = pm - po
         , diff_pm_pr = pm - pr
         , nonsig_5p_pm = ifelse(sig_5p_zo == 1 & sig_5p_pm == 0, 1, 0)#flag for when pm adjusted po to a non-significant value
         , same_sigflg_pmpr = ifelse(sig_5p_zr == 1 & sig_5p_pm == 1 | sig_5p_zr == 0 & sig_5p_pm == 0, 1, 0)
  )
str(datatable_phacking_ew1_pm)
save(datatable_phacking_ew1_pm, file='datasets/20240708_datatable_phacking_ew1_pm.Rdata')

load('datasets/20240708_datatable_phacking_ew1_pm.Rdata')

# Check data
head(datatable_phacking_ew1_pm)
unique(datatable_phacking_ew1_pm$df)
hist(datatable_phacking_ew1_pm$diff_pm_po) 



######## Descriptive analysis ########
summary_datatable_phacking_ew1_pm <-  datatable_phacking_ew1_pm %>%
  group_by(df) %>%
  summarise(avg_zo = mean(zo)
            , median_zo = median(zo)
            , sd_zo = sd(zo)
            , min_zo = min(zo)
            , max_zo = max(zo)
            #po
            , avg_po = mean(po)
            , median_po = median(po)
            , sd_po = sd(po)
            , min_po = min(po)
            , max_po = max(po)
            #zr
            , avg_zr = mean(zr)
            , median_zr = median(zr)
            , sd_zr = sd(zr)
            , min_zr = min(zr)
            , max_zr = max(zr)
            #pr
            , avg_pr = mean(pr)
            , median_pr = median(pr)
            , sd_pr = sd(pr)
            , min_pr = min(pr)
            , max_pr = max(pr)
            #pm
            , avg_pm = mean(pm)
            , median_pm = median(pm)
            , sd_pm = sd(pm)
            , min_pm = min(pm)
            , max_pm = max(pm)
            , avg_diff_pmpo = mean(diff_pm_po)
            , avg_diff_pmpr = mean(diff_pm_pr)
            ## significance flags 
            , #alpha = 10%
            sum_sig_10p_zo = sum(sig_10p_zo == 1)
            , sum_sig_10p_zr = sum(sig_10p_zr == 1)
            , sum_phacking = sum(p_hacking_10p==1)
            , reverse_sign_z = sum(zo>=0 & zr<=0, zo>=0 & zr<=0) 
            , prob_phacking = (sum(p_hacking_10p == 1) / sum(p_hacking_10p == 1, p_hacking_10p == 0)) 
            , prob_counter_phacking = (sum(sig_10p_zo==0 & sig_10p_zr ==1)/sum(sig_10p_zo == 1, sig_10p_zo == 0))
            , prob_sig_10p_zo = (sum(sig_10p_zo==0 & sig_10p_zr ==0)/sum(sig_10p_zo == 1, sig_10p_zo == 0))
            , prob_notsig_10p_zo = (sum(sig_10p_zo==1 & sig_10p_zr ==1)/sum(sig_10p_zo == 1, sig_10p_zo == 0))
            , type1error_po_10p = (sum(sig_10p_zo == 1) / sum(sig_10p_zo == 1, sig_10p_zo == 0)) 
            , type1error_pr_10p = (sum(sig_10p_zr == 1) / sum(sig_10p_zr == 1, sig_10p_zr == 0))
            #alpha = 5%
            , sum_sig_5p_zo = sum(sig_5p_zo == 1)
            , sum_sig_5p_zr = sum(sig_5p_zr == 1)
            , sum_phacking = sum(p_hacking_5p==1)
            , prob_phacking = (sum(p_hacking_5p == 1) / sum(p_hacking_5p == 1, p_hacking_5p == 0)) 
            , prob_counter_phacking = (sum(sig_5p_zo==0 & sig_5p_zr ==1)/sum(sig_5p_zo == 1, sig_5p_zo == 0))
            , prob_sig_5p_zo = (sum(sig_5p_zo==0 & sig_5p_zr ==0)/sum(sig_5p_zo == 1, sig_5p_zo == 0))
            , prob_notsig_5p_zo = (sum(sig_5p_zo==1 & sig_5p_zr ==1)/sum(sig_5p_zo == 1, sig_5p_zo == 0))
            , type1error_po_5p = (sum(sig_5p_zo == 1) / sum(sig_5p_zo == 1, sig_5p_zo == 0)) 
            , type1error_pr_5p = (sum(sig_5p_zr == 1) / sum(sig_5p_zr == 1, sig_5p_zr == 0))
            #alpha = 1%
            , sum_sig_1p_zo = sum(sig_1p_zo == 1)
            , sum_sig_1p_zr = sum(sig_1p_zr == 1)
            , sum_phacking = sum(p_hacking_1p==1)
            , prob_phacking = (sum(p_hacking_1p == 1) / sum(p_hacking_1p == 1, p_hacking_1p == 0)) 
            , prob_counter_phacking = (sum(sig_1p_zo==0 & sig_1p_zr ==1)/sum(sig_1p_zo == 1, sig_1p_zo == 0))
            , prob_sig_1p_zo = (sum(sig_1p_zo==0 & sig_1p_zr ==0)/sum(sig_1p_zo == 1, sig_1p_zo == 0))
            , prob_notsig_1p_zo = (sum(sig_1p_zo==1 & sig_1p_zr ==1)/sum(sig_1p_zo == 1, sig_1p_zo == 0))
            , type1error_po_1p = (sum(sig_1p_zo == 1) / sum(sig_1p_zo == 1, sig_1p_zo == 0)) 
            , type1error_pr_1p = (sum(sig_1p_zr == 1) / sum(sig_1p_zr == 1, sig_1p_zr == 0))
            #alpha = 0.5%
            , sum_sig_05p_zo = sum(sig_05p_zo == 1)
            , sum_sig_05p_zr = sum(sig_05p_zr == 1)
            , sum_phacking = sum(p_hacking_05p==1)
            , prob_phacking = (sum(p_hacking_05p == 1) / sum(p_hacking_05p == 1, p_hacking_05p == 0)) 
            , prob_counter_phacking = (sum(sig_05p_zo==0 & sig_05p_zr ==1)/sum(sig_05p_zo == 1, sig_05p_zo == 0))
            , prob_sig_05p_zo = (sum(sig_05p_zo==0 & sig_05p_zr ==0)/sum(sig_05p_zo == 1, sig_05p_zo == 0))
            , prob_notsig_05p_zo = (sum(sig_05p_zo==1 & sig_05p_zr ==1)/sum(sig_05p_zo == 1, sig_05p_zo == 0))
            , type1error_po_05p = (sum(sig_05p_zo == 1) / sum(sig_05p_zo == 1, sig_05p_zo == 0)) 
            , type1error_pr_05p = (sum(sig_05p_zr == 1) / sum(sig_05p_zr == 1, sig_05p_zr == 0))
            #parameters
            , param_m = unique(param_m)
            , param_rho = unique(param_rho)
            , param_mhu = unique(param_mhu)
  )

head(summary_datatable_phacking_ew1_pm)
write_csv(summary_datatable_phacking_ew1_pm, file='datasets/20240702_summary_datatable_phacking_ew1_pm.csv')



### Plots ###
names(datatable_phacking_ew1_pm)[names(datatable_phacking_ew1_pm) == "param_m"] <- "m"
names(datatable_phacking_ew1_pm)[names(datatable_phacking_ew1_pm) == "param_rho"] <- "rho"
datatable_phacking_ew1_pm$threshold <- "5% alpha level"

## histogram of zo
hist_ew1_zo <-
  ggplot(datatable_phacking_ew1_pm, aes(x = zo)) +
  geom_histogram(binwidth = 0.3, fill = "blue", color = "black") +
  geom_vline(aes(xintercept = 1.96, color = "5% alpha level"), linetype = "dashed") +
  geom_vline(aes(xintercept = -1.96, color = "5% alpha level"), linetype = "dashed") +
  labs(
    title = expression("Histogram of " *z[o]),
    x = expression(z[o]),
    y = 'Frequency',
    color = "Critical value"
  ) +
  scale_color_manual(
    values = c("5% alpha level" = "red")
  ) +
  facet_grid(rho ~ m, labeller = label_both) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

hist_ew1_zo

# ggsave(
#   filename = "plots/plots_ew1/hist_ew1_zo.pdf",
#   plot = hist_ew1_zo, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )



## histogram of zr
hist_ew1_zr <-
  ggplot(datatable_phacking_ew1_pm, aes(x = zr)) +
  geom_histogram(binwidth = 0.3, fill = "blue", color = "black") +
  
  # Add red dashed lines and map to a dummy variable for legend
  geom_vline(aes(xintercept = 1.96, color = "5% alpha level"), linetype = "dashed") +
  geom_vline(aes(xintercept = -1.96, color = "5% alpha level"), linetype = "dashed") +
  
  # Customize the title and axis labels
  labs(
    title = expression("Histogram of " *z[r]),
    x = expression(z[r]),
    y = 'Frequency',
    color = "Critical value"
  ) +
  
  scale_color_manual(
    values = c("5% alpha level" = "red")
  ) +
  
  facet_grid(rho ~ m, labeller = label_both) +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

hist_ew1_zr

# ggsave(
#   filename = "plots/plots_ew1/hist_ew1_zr.pdf",
#   plot = hist_ew1_zr, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )


# boxplot po

boxplot_rho_po_ew1 <-
  ggplot(datatable_phacking_ew1_pm, aes(x = as.factor(rho), y = po)) +
  geom_boxplot() +
  labs(
    title = expression("Effect of Correlation Parameter on Original Study P-Values (p"[o]*")"),
    x = expression("Correlation parameter " * rho),
    y = expression(p[o])
  ) +
  facet_wrap(~m, labeller = label_both) +
  geom_hline(aes(yintercept = 0.05, color = threshold), linetype = "dashed") +
  scale_color_manual(name = "Legend", values = c("5% alpha level" = "red")) +
  theme_bw() +
  theme(
    legend.position = c(1, 0), 
    legend.justification = c(1, 0), # Adjust legend anchor point
    strip.text = element_text(size= 14, face = "bold"),  
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 16),      
    plot.title = element_text(size = 17, face = "bold"),  
    legend.text = element_text(size = 14),     
    legend.title = element_text(size = 14)     
  )


boxplot_rho_po_ew1

# ggsave(
#   filename = "plots/plots_ew1/boxplot_rho_po_ew1.pdf",
#   plot = boxplot_rho_po_ew1, 
#   width = 10,       
#   height = 8,       
#   units = "in"       
# )


# boxplot pr

boxplot_rho_pr_ew1 <-
  ggplot(datatable_phacking_ew1_pm, aes(x = as.factor(rho), y = pr)) +
  geom_boxplot() +
  labs(
    title = expression("Effect of Correlation Parameter on Replication Study P-Values (p"[r]*")"),
    x = expression("Correlation parameter " * rho),
    y = expression(p[r])
  )+
  facet_wrap(~m, labeller = label_both) +
  geom_hline(aes(yintercept = 0.05, color = threshold), linetype = "dashed") +
  scale_color_manual(name = "Legend", values = c("5% alpha level" = "red")) +
  theme_bw() +
  theme(
    legend.position = c(1, 0), 
    legend.justification = c(1, 0),
    strip.text = element_text(size= 14, face = "bold"),  
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 16),      
    plot.title = element_text(size = 17, face = "bold"),  
    legend.text = element_text(size = 14),     
    legend.title = element_text(size = 14)     
  )

boxplot_rho_pr_ew1

# ggsave(
#   filename = "plots/plots_ew1/boxplot_rho_pr_ew1.pdf",
#   plot = boxplot_rho_pr_ew1, 
#   width = 10,       
#   height = 8,       
#   units = "in"       
# )


### boxplot pm 
boxplot_rho_pm_ew1 <-
  ggplot(datatable_phacking_ew1_pm, aes(x = as.factor(rho), y = pm)) +
  geom_boxplot() +
  labs(
    title = expression("Effect of Correlation Parameter on Adjusted P-Values (p"[m]*")"),
    x = expression("Correlation parameter " * rho),
    y = expression(p[m])
  ) +
  facet_wrap(~m, labeller = label_both) +
  geom_hline(aes(yintercept = 0.05, color = threshold), linetype = "dashed") +
  scale_color_manual(name = "Legend", values = c("5% alpha level" = "red")) +
  theme_bw() +
  theme(
    legend.position = c(1, 0), 
    legend.justification = c(1, 0), 
    strip.text = element_text(size= 14, face = "bold"),  
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 16),      
    plot.title = element_text(size = 17, face = "bold"),  
    legend.text = element_text(size = 14),     
    legend.title = element_text(size = 14)     
  )

boxplot_rho_pm_ew1

# ggsave(
#   filename = "plots/plots_ew1/boxplot_rho_pm_ew1.pdf",
#   plot = boxplot_rho_pm_ew1, 
#   width = 10,       
#   height = 8,       
#   units = "in"       
# )

# Scatter plot po vs pr
phacking_popr <- ggplot(datatable_phacking_ew1_pm, aes(x = po, y = pr, color = factor(p_hacking_5p))) +
  geom_point() +
  geom_vline(aes(xintercept = 0.05, color = threshold), linetype = "dashed", linewidth=0.8) +
  geom_hline(aes(yintercept = 0.05, color = threshold), linetype = "dashed", linewidth=0.8) +
  labs(
    title = expression("Scatter Plot of " *p[o]* " vs " *p[r]),
    x = expression(p[o]),
    y = expression(p[r]),
    color = "Legend"
  ) +
  scale_color_manual(
    values = c("0" = "lightblue2", "1" = "gold", "5% alpha level" = "red"),
    labels = c("1" = "po sig/ pr not sig", "0" = "other", "5% alpha level" = "5% alpha level")  
  ) +  
  facet_grid(rho ~ m, labeller = label_both) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

phacking_popr

# Save the plot
# ggsave(
#   filename = "plots/plots_ew1/phacking_popr.pdf",
#   plot = phacking_popr, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )

## Probabilities of test decisions between po and pr
sum(datatable_phacking_ew1_pm$po < 0.05 & datatable_phacking_ew1_pm$pr < 0.05) / nrow(datatable_phacking_ew1_pm)
sum(datatable_phacking_ew1_pm$po > 0.05 & datatable_phacking_ew1_pm$pr < 0.05) / nrow(datatable_phacking_ew1_pm)
sum(datatable_phacking_ew1_pm$po > 0.05 & datatable_phacking_ew1_pm$pr > 0.05) / nrow(datatable_phacking_ew1_pm)
sum(datatable_phacking_ew1_pm$po < 0.05 & datatable_phacking_ew1_pm$pr > 0.05) / nrow(datatable_phacking_ew1_pm)

plausibility_pmpo_ew1 <- ggplot(datatable_phacking_ew1_pm, aes(x = pm, y = po, color = factor(nonsig_5p_pm))) +
  geom_point() +
  facet_grid(rho ~ m, labeller = label_both) +  
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(
    title = expression("Plausibility Check of Adjusted P-Values (p"[m]*") with Respect to Original P-Values (p"[o]*")"),
    x = expression(p[m]),
    y = expression(p[o]),
    color = "Test decision at alpha 5%"
  ) +
  geom_hline(aes(yintercept = 0.05, color = "5% alpha level"), linetype = "dashed") +  
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +  
  scale_color_manual(
    values = c("1" = "darkorange", "0" = "steelblue1", "5% alpha level" = "red"),  
    labels = c("1" = "Reversed", "0" = "Maintained", "5% alpha level" = "5% alpha level")  
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"), 
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16),  
    plot.title = element_text(size = 17, face = "bold"),  
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1)  
  )


plausibility_pmpo_ew1

# ggsave(
#   filename = "plots/plots_ew1/plausibility_pmpo_ew1.pdf",
#   plot = plausibility_pmpo_ew1, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )

## magnitude of adjustment
magnitude_pmpo_ew1 <- ggplot(datatable_phacking_ew1_pm, aes(x = rho, y = avg_diff_pmpo, color = factor(m))) +
  geom_line(linewidth=2) +
  labs(
    title = "Magnitude of Adjustment",
    x = expression("Correlation Parameter " *rho),
    y =  expression("Average Difference between " *p[m]* " and " *p[m]),
    color = "Nr. Strategies (m)"
  ) +
  scale_x_continuous(breaks = unique(summary_datatable_phacking_ew1_pm$rho)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),  
    axis.title = element_text(size = 18),  
    plot.title = element_text(size = 20, face = "bold"),  
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 16)   
  )

magnitude_pmpo_ew1

# ggsave(
#   filename = "plots/plots_ew1/magnitude_pmpo.pdf",
#   plot = magnitude_pmpo_ew1, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )

####### Power analysis #######
# Set significance level
alpha <- 0.05

# Iterate over unique parameter combinations
for (param_m in unique(datatable_phacking_ew1_pm$param_m)) {
  for (param_rho in unique(datatable_phacking_ew1_pm$param_rho)) {
    subset_data <- datatable_phacking_ew1_pm %>%
      filter(param_m == !!param_m & param_rho == !!param_rho)
    
    # Calculate power for unadjusted
    power_unadjusted <- sum(subset_data$sig_5p_zo) / nrow(subset_data)
    
    # Calculate power for your method (pm)
    rejected_pm <- ifelse(subset_data$pm < alpha, 1, 0)
    power_pm <- sum(rejected_pm) / length(rejected_pm)
    
    # Calculate power for your method (pm)
    rejected_pm <- ifelse(subset_data$pm < alpha, 1, 0)
    power_pm <- sum(rejected_pm) / length(rejected_pm)
    
    # Calculate adjusted p-values for Bonferroni
    po_bonferroni <- pmin(subset_data$po * subset_data$param_m, 1)
    rejected_bonferroni <- ifelse(po_bonferroni < alpha, 1, 0)
    power_bonferroni <- sum(rejected_bonferroni) / length(rejected_bonferroni)
    
    
    # Append results to the power_results data frame
    power_results <- rbind(power_results, data.frame(
      param_m = param_m,
      param_rho = param_rho,
      power_unadjusted = power_unadjusted,
      power_pm = power_pm,
      power_pm = power_pm,
      power_bonferroni = power_bonferroni
    ))
  }
}

head(power_results)

# Plot power curves for different adjustment methods
power_curve_plot <-ggplot(power_results, aes(x = param_rho)) +
  geom_line(aes(y = power_unadjusted, color = "Unadjusted"), linewidth = 1) +
  geom_line(aes(y = power_pm, color = "Distribution-based"), linewidth = 1) +
  geom_line(aes(y = power_bonferroni, color = "Bonferroni"), linewidth = 1) +
  labs(
    title = "Power Curve",
    x = expression("Correlation (" * rho * ")"),
    y = "Power",
    color = "Adjustment Method"
  ) +
  facet_grid(param_m ~ .) +
  scale_x_continuous(
    breaks = unique(power_results$param_rho),  
    labels = unique(power_results$param_rho)   
  ) +
  scale_color_manual(values = c("Unadjusted" = "black", "Distribution-based" = "violetred", "Bonferroni" = "lightseagreen")) +
  theme_bw() +
  theme(
    legend.position = "right",              
    legend.justification = "center",        
    axis.text = element_text(size = 15),     
    axis.title = element_text(size = 16),   
    plot.title = element_text(size = 19, face = "bold"),  
    legend.text = element_text(size = 16),   
    legend.title = element_text(size = 17),  
    strip.text = element_text(size = 17, face = "bold")   
  )


# Print the plot
power_curve_plot

# Save the plot as a PDF file
# ggsave(
#   filename = "plots/plots_ew1/power_curve.pdf",
#   plot = power_curve_plot, 
#   width = 10,        
#   height = 8,        
#   units = "in"       
# )


