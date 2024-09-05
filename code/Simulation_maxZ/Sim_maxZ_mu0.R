######################################################
# Simulation maxZ: mu = 0
#####################################################
setwd('/Users/nikitapaschan/Thesis_NavigatingGFP/')
## Parameters
m = c(2, 5, 10, 100, 1000)
rho = c(0, 0.3, 0.5, 0.7, 0.9)


## Generate data
data_phacking_ew0 <- sim_datasets_phacking(m = m, rho = rho, mhu = 0)
# save(data_phacking_ew0, file='datasets/20240627_data_phacking_ew0.RData')
# load('datasets/20240627_data_phacking_ew0.RData')
head(data_phacking_ew0$data_mvrnorm)
head(data_phacking_ew0$data_all)


## Transform to data table for better handling
datatable_phacking_ew0 <- as.data.table(data_phacking_ew0$data_all)
unique(datatable_phacking_ew0$df)
head(datatable_phacking_ew0)

## load data 
# save(data_phacking_ew0, file='datasets/20240627_data_phacking_ew0.RData')
# load('datasets/20240627_data_phacking_ew0.RData')


# Use mclapply to process each row in parallel
datatable_phacking_ew0_pm_list <- mclapply(
  1:nrow(datatable_phacking_ew0[niter %in% 1:10]),
  function(i) process_row(datatable_phacking_ew0[niter %in% 1:10][i, ]),
  mc.cores = detectCores() / 2
)


datatable_phacking_ew0_pm <- do.call(rbind, pblapply(1:nrow(datatable_phacking_ew0[niter %in% 1:10]), function(i) process_row(datatable_phacking_ew0[niter %in% 1:10][i, ]), cl = detectCores() - 1))



# Combine the list of data frames into one data table
datatable_phacking_ew0_pm <- rbindlist(datatable_phacking_ew0_pm_list)

# Save the results
# save(datatable_phacking_ew0_pm, file='datasets/20240814_datatable_phacking_ew0_pm.Rdata')


# Progress completion sound
beep(2)

# Return the modified dataset
head(datatable_phacking_ew0_pm)

# Merge pm to existing data frame
datatable_phacking_ew0_pm <- merge(datatable_phacking_ew0, datatable_phacking_ew0_pm, by = c('df', 'zo'))
# save(datatable_phacking_ew0_pm, file='datasets/20240714_datatable_phacking_ew0_pm.Rdata')
load('datasets/20240714_datatable_phacking_ew0_pm.Rdata')

# magnitude of adjustment
datatable_phacking_ew0_pm <- datatable_phacking_ew0_pm %>%
  mutate(diff_pm_po = pm - po
         , diff_pm_pr = pm - pr
         , nonsig_5p_pm = ifelse(sig_5p_zo == 1 & sig_5p_pm == 0, 1, 0)#flag for when pm adjusted po to a non-significant value
         , same_sigflg_pmpr = ifelse(sig_5p_zr == 1 & sig_5p_pm == 1 | sig_5p_zr == 0 & sig_5p_pm == 0, 1, 0)
  )

# save(datatable_phacking_ew0_pm, file='datasets/20240814_datatable_phacking_ew0_pm.Rdata')



######## Descriptive analysis ########
summary_datatable_phacking_ew0_pm <-  datatable_phacking_ew0_pm %>%
  group_by(m, rho) %>%
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

head(summary_datatable_phacking_ew0_pm)
# write_csv(summary_datatable_phacking_ew0_pm, file='datasets/20240702_summary_datatable_phacking_ew0_pm.csv')


### Plots ###
names(datatable_phacking_ew0_pm)[names(datatable_phacking_ew0_pm) == "param_m"] <- "m"
names(datatable_phacking_ew0_pm)[names(datatable_phacking_ew0_pm) == "param_rho"] <- "rho"
datatable_phacking_ew0_pm$threshold <- "5% alpha level"

## histogram of zo
hist_ew0_zo <-
  ggplot(datatable_phacking_ew0_pm, aes(x = zo)) +
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

hist_ew0_zo

# ggsave(
#   filename = "plots/plots_ew0/hist_ew0_zo.pdf",
#   plot = hist_ew0_zo, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )



## histogram of zr
hist_ew0_zr <-
  ggplot(datatable_phacking_ew0_pm, aes(x = zr)) +
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

hist_ew0_zr

# ggsave(
#   filename = "plots/plots_ew0/hist_ew0_zr.pdf",
#   plot = hist_ew0_zr, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )


# boxplot po

boxplot_rho_po_ew0 <-
  ggplot(datatable_phacking_ew0_pm, aes(x = as.factor(rho), y = po)) +
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


boxplot_rho_po_ew0

# ggsave(
#   filename = "plots/plots_ew0/boxplot_rho_po_ew0.pdf",
#   plot = boxplot_rho_po_ew0, 
#   width = 10,       
#   height = 8,       
#   units = "in"       
# )


# boxplot pr

boxplot_rho_pr_ew0 <-
  ggplot(datatable_phacking_ew0_pm, aes(x = as.factor(rho), y = pr)) +
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

boxplot_rho_pr_ew0

# ggsave(
#   filename = "plots/plots_ew0/boxplot_rho_pr_ew0.pdf",
#   plot = boxplot_rho_pr_ew0, 
#   width = 10,       
#   height = 8,       
#   units = "in"       
# )


### boxplot pm 
boxplot_rho_pm_ew0 <-
  ggplot(datatable_phacking_ew0_pm, aes(x = as.factor(rho), y = pm)) +
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

boxplot_rho_pm_ew0

# ggsave(
#   filename = "plots/plots_ew0/boxplot_rho_pm_ew0.pdf",
#   plot = boxplot_rho_pm_ew0, 
#   width = 10,       
#   height = 8,       
#   units = "in"       
# )

# Scatter plot po vs pr
phacking_popr <- ggplot(datatable_phacking_ew0_pm, aes(x = po, y = pr, color = factor(p_hacking_5p))) +
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
#   filename = "plots/plots_ew0/phacking_popr.pdf",
#   plot = phacking_popr, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )

## Probabilities of test decisions between po and pr
sum(datatable_phacking_ew0_pm$po < 0.05 & datatable_phacking_ew0_pm$pr < 0.05) / nrow(datatable_phacking_ew0_pm)
sum(datatable_phacking_ew0_pm$po > 0.05 & datatable_phacking_ew0_pm$pr < 0.05) / nrow(datatable_phacking_ew0_pm)
sum(datatable_phacking_ew0_pm$po > 0.05 & datatable_phacking_ew0_pm$pr > 0.05) / nrow(datatable_phacking_ew0_pm)
sum(datatable_phacking_ew0_pm$po < 0.05 & datatable_phacking_ew0_pm$pr > 0.05) / nrow(datatable_phacking_ew0_pm)

plausibility_pmpo_ew0 <- ggplot(datatable_phacking_ew0_pm, aes(x = pm, y = po, color = factor(nonsig_5p_pm))) +
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


plausibility_pmpo_ew0

# ggsave(
#   filename = "plots/plots_ew0/plausibility_pmpo_ew0.pdf",
#   plot = plausibility_pmpo_ew0, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )

## magnitude of adjustment
magnitude_pmpo<- ggplot(datatable_phacking_ew0_pm, aes(x = rho, y = avg_diff_pmpo, color = factor(m))) +
  geom_line(linewidth=2) +
  labs(
    title = "Magnitude of Adjustment",
    x = expression("Correlation Parameter " *rho),
    y =  expression("Average Difference between " *p[m]* " and " *p[m]),
    color = "Nr. Strategies (m)"
  ) +
  scale_x_continuous(breaks = unique(summary_datatable_phacking_ew0_pm$rho)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),  
    axis.title = element_text(size = 18),  
    plot.title = element_text(size = 20, face = "bold"),  
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 16)   
  )

# ggsave(
#   filename = "plots/plots_ew0/magnitude_pmpo.pdf",
#   plot = magnitude_pmpo, 
#   width = 10,        
#   height = 8,        
#   units = "in"      
# )




