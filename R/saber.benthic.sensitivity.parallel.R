library(parallel)
library(pbapply)
library(foreach)
library(doParallel)
library(BayesianTools)

# Define a function to suppress output to a temporary file
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Create cluster
num_cores <- detectCores() - 1  # Use one less than the total number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Load required libraries and source functions on each worker
clusterEvalQ(cl, {
  library(parallel)
  library(pbapply)
  library(foreach)
  library(doParallel)
  library(data.table)
  library(BayesianTools)
  
  source("./R/saber_inversion_final.R")
  source("./R/SABER_forward_fast.R")
  source("./R/solve.objective.inverse_fast.R")
  source("./R/mcmc_bayes_parallel.R")
  
})

# Define the data decomposing function
decompose_data <- function(data, num_groups) {
  split(data, rep(1:num_groups, length.out = nrow(data)))
}

decomposed_data <- decompose_data(obsdata, 11)


set.seed(345)
slope_idx = sample(size = 5, x = seq(1,length(spectral_slope_vec$ad_ag),1))

spectral_slope_vec = spectral_slope_vec[slope_idx,]

rb_idx = sample(size = 2, x = seq(1,length(combinations$V1),1))

combinations = combinations[rb_idx,]



combinations_grid <- expand.grid(i = seq_along(combinations$V1), j = seq_along(spectral_slope_vec$ag))

process_combination <- function(obsdata_row, idx) {
  i <- combinations_grid[idx, "i"]
  j <- combinations_grid[idx, "j"]
  
  # Suppress output to a temporary file
  #quiet({
    interp_rb(wavelength_input = wavelength, bottom_type = as.character(combinations[i,]))
    
    result <- inverse_runBayes(obsdata = as.numeric(obsdata_row), 
                               rrs_type = type_Rrs_below, 
                               max_par = upper.bound, 
                               min_par = lower.bound, 
                               param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", 
                                             "fa3", "pop_sd"),
                               constrain_config = c(FALSE, FALSE, FALSE),
                               bbp_const_val = 0.007, 
                               bgc_const_val = c(5, 1.003, 0.007), 
                               iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                                                  "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
                               qaa_slope = FALSE, 
                               manual_slope = TRUE, 
                               manual_slope_vals = c("s_g" = spectral_slope_vec[j, 4], 
                                                     "s_d" = spectral_slope_vec[j, 3], 
                                                     "gamma" = spectral_slope_vec[j, 2]), 
                               iter_count = 20000, 
                               sampler_mcmc = samplerlist[6], 
                               wavelngth_sim = wavelength, 
                               sa_model = "am03", 
                               hybrid_mode = FALSE, 
                               plot_rrs = FALSE)
  #})
  
  return(result)
}


results_list <- foreach( i = names(decomposed_data),
                        #.combine = rbind, 
                        .verbose = T,
                        .packages = c("data.table", "BayesianTools")) %:% 
                foreach( idx_c = 1:nrow(combinations_grid),
                        #.combine = c, 
                        .verbose = T,
                        .packages = c("data.table", "BayesianTools")) %dopar% {
                          source("./R/SABER_forward_fast.R")
                          source("./R/solve.objective.inverse_fast.R")
                          source("./R/mcmc_bayes_parallel.R")
                          apply(decomposed_data[[i]], 1, process_combination,
                                idx = idx_c)
                        }


# Parallel processing with nested loops using foreach
# results_list <- foreach(i = 1:nrow(obsdata), 
#                         .verbose = T,
#                         #.combine = 'rbind', 
#                         .packages = c('data.table', 'BayesianTools')) %:% 
#   foreach(idx_c = 1:nrow(combinations_grid), 
#           #.combine = 'rbind', 
#           .verbose = T, .packages = c('data.table', 'BayesianTools')) %dopar% {
#             source("./R/SABER_forward_fast.R")
#             source("./R/solve.objective.inverse_fast.R")
#             source("./R/mcmc_bayes_parallel.R")
#             source("./R/saber_inversion_final.R")
#             #apply(obs_group, 1, process_combination, idx = idx)
#             process_combination(obsdata_row = obsdata[i,], idx = idx_c)
#     
#   }


stopCluster(cl)

results_list_bacc = results_list

results_list_df = do.call(c, results_list)
results_list_df = as.data.frame(t(do.call(cbind, results_list_df)))


rrs_names = rownames(obsdata)
rrs_names_count = substr(rrs_names, 5, 10)


# Create a list of list for each rrs observation inverted results for all the slope and benthic ombinations
list_of_50_dfs <- vector("list", length = length(rrs_names))


for (i in seq_along(rrs_names)) {
  rrs_idx = grep(paste0("^Rrs\\.", rrs_names_count[i], ".*"), rownames(results_list_df))
  list_of_50_dfs[[i]] = results_list_df[rrs_idx, ]
}

# Set names of the list for easier access
names(list_of_50_dfs) <- rrs_names


list_of_2380_dfs <- vector("list", length = nrow(combinations_grid))

# Loop over each combination and collect corresponding rows from all rrs_names
for (j in seq_along(list_of_2380_dfs)) {
  list_of_2380_dfs[[j]] <- do.call(rbind, lapply(list_of_50_dfs, function(df) df[j, ]))
}

# Set names of the list for easier access (optional)
names(list_of_2380_dfs) <- apply(combinations_grid, 1, function(x) paste(x, collapse = "_"))


# Initialize a data frame to store the results
results_summary <- data.frame(
  Combination = apply(combinations_grid, 1, function(x) paste(x, collapse = "_")),
  Uncertainty_chl = numeric(nrow(combinations_grid)),
  Bias_chl = numeric(nrow(combinations_grid)),
  R2_chl = numeric(nrow(combinations_grid)),
  Slope_chl = numeric(nrow(combinations_grid)),
  Uncertainty_H = numeric(nrow(combinations_grid)),
  Bias_H = numeric(nrow(combinations_grid)),
  R2_H = numeric(nrow(combinations_grid)),
  Slope_H = numeric(nrow(combinations_grid))
)

# Function to calculate R² and slope
calculate_r2_and_slope <- function(actual, predicted) {
  model <- lm(predicted ~ actual)
  r_squared <- summary(model)$r.squared
  slope <- coef(model)[2]  # Extracting the slope (the coefficient for 'actual')
  return(list(R2 = r_squared, Slope = slope))
}

# Calculate uncertainty, bias, R², and slope for each of the 2380 combinations
for (i in seq_along(list_of_2380_dfs)) {
  
  # Extract the i-th data frame containing predictions
  pred_df <- list_of_2380_dfs[[i]]
  
  # Extract predictions for 'chl' and 'H'
  pred_chl <- pred_df$chl
  pred_H <- pred_df$H
  
  # Extract actual values for 'chl' and 'H' from insitu_val
  actual_chl <- insitu_val$chl
  actual_H <- insitu_val$H
  
  # Calculate Uncertainty (Standard Deviation of predictions)
  uncertainty_chl <- sd(pred_chl)
  uncertainty_H <- sd(pred_H)
  
  # Calculate Bias (Mean of predictions minus actual value)
  bias_chl <- mean(pred_chl) - mean(actual_chl)
  bias_H <- mean(pred_H) - mean(actual_H)
  
  # Calculate R² and Slope for chl
  r2_slope_chl <- calculate_r2_and_slope(actual_chl, pred_chl)
  
  # Calculate R² and Slope for H
  r2_slope_H <- calculate_r2_and_slope(actual_H, pred_H)
  
  # Store the results
  results_summary$Uncertainty_chl[i] <- uncertainty_chl
  results_summary$Bias_chl[i] <- bias_chl
  results_summary$R2_chl[i] <- r2_slope_chl$R2
  results_summary$Slope_chl[i] <- r2_slope_chl$Slope
  results_summary$Uncertainty_H[i] <- uncertainty_H
  results_summary$Bias_H[i] <- bias_H
  results_summary$R2_H[i] <- r2_slope_H$R2
  results_summary$Slope_H[i] <- r2_slope_H$Slope
}

# Print or save the results summary
print(results_summary)

plot(insitu_val$H, list_of_2380_dfs[["24_19"]]$H, xlim = c(0,10), ylim = c(0,10))

# Initialize a data frame to store the results for each observation
optimal_combinations_summary <- data.frame(
  Observation = names(list_of_50_dfs),
  Optimal_Combination = character(length(list_of_50_dfs)),
  Uncertainty_chl = numeric(length(list_of_50_dfs)),
  Bias_chl = numeric(length(list_of_50_dfs)),
  Uncertainty_H = numeric(length(list_of_50_dfs)),
  Bias_H = numeric(length(list_of_50_dfs))
)

# Initialize a list to store detailed metrics for each Rrs observation
detailed_metrics_list <- vector("list", length(list_of_50_dfs))

# Loop through each observation (i.e., each data frame in list_of_50_dfs)
for (obs_idx in seq_along(list_of_50_dfs)) {
  
  pred_df <- list_of_50_dfs[[obs_idx]]
  
  best_combination <- NULL
  min_uncertainty <- Inf
  min_bias <- Inf
  
  # Initialize variables to store the best metrics
  best_uncertainty_chl <- NA
  best_bias_chl <- NA
  best_uncertainty_H <- NA
  best_bias_H <- NA
  
  
  # Initialize a data frame to store metrics for all combinations for the current observation
  detailed_metrics <- data.frame(
    Combination = character(nrow(pred_df)),
    Uncertainty_chl = numeric(nrow(pred_df)),
    Bias_chl = numeric(nrow(pred_df)),
    Uncertainty_H = numeric(nrow(pred_df)),
    Bias_H = numeric(nrow(pred_df))
    
  )
  
  # Loop through each combination of slopes for the current observation
  for (i in 1:nrow(pred_df)) {
    
    pred_chl <- pred_df[i, "chl"]
    pred_H <- pred_df[i, "H"]
    unc_chl = pred_df[i, "sd_chl"]
    unc_H = pred_df[i, "sd_H"]
    
    actual_chl <- insitu_val$chl[obs_idx]
    actual_H <- insitu_val$H[obs_idx]
    
    # Calculate the metrics
    
    uncertainty_chl <- sd(c(actual_chl,pred_chl))
    
    #uncertainty_chl <- unc_chl
    
    bias_chl <- mean(pred_chl) - actual_chl
    
    uncertainty_H <- sd(c(actual_H,pred_H))
    #uncertainty_H <- unc_H
    
    bias_H <- mean(pred_H) - actual_H
    
    
    
    # Store the metrics for this combination
    detailed_metrics$Combination[i] <- paste(combinations_grid[i, ], collapse = "_")
    detailed_metrics$Uncertainty_chl[i] <- uncertainty_chl
    detailed_metrics$Bias_chl[i] <- bias_chl
    detailed_metrics$Uncertainty_H[i] <- uncertainty_H
    detailed_metrics$Bias_H[i] <- bias_H
    
    
    # Update the best combination if the current one is better
    if (uncertainty_chl < min_uncertainty || bias_chl < min_bias) {
      min_uncertainty <- uncertainty_chl
      min_bias <- bias_chl
      
      best_combination <- paste(combinations_grid[i, ], collapse = "_")
      best_uncertainty_chl <- uncertainty_chl
      best_bias_chl <- bias_chl
      best_uncertainty_H <- uncertainty_H
      best_bias_H <- bias_H
    }
  }
  
  # Store the best combination and its metrics for the current observation
  optimal_combinations_summary$Optimal_Combination[obs_idx] <- best_combination
  optimal_combinations_summary$Uncertainty_chl[obs_idx] <- best_uncertainty_chl
  optimal_combinations_summary$Bias_chl[obs_idx] <- best_bias_chl
  optimal_combinations_summary$Uncertainty_H[obs_idx] <- best_uncertainty_H
  optimal_combinations_summary$Bias_H[obs_idx] <- best_bias_H
  
  # Store the detailed metrics for this observation in the list
  detailed_metrics_list[[obs_idx]] <- detailed_metrics
}

# Print the summary of optimal combinations for each observation
print(optimal_combinations_summary)

# Combine the detailed metrics for each observation into a single data frame for easy inspection
all_detailed_metrics <- do.call(rbind, detailed_metrics_list)

# Print the detailed metrics for all observations and combinations
print(all_detailed_metrics)


# Calculate the mean uncertainty and bias for each combination
mean_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, Uncertainty_H, Bias_H) ~ Combination, 
                          data = all_detailed_metrics, 
                          FUN = mean)
min_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, Uncertainty_H, Bias_H) ~ Combination, 
                         data = all_detailed_metrics, 
                         FUN = min)
max_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, Uncertainty_H, Bias_H) ~ Combination, 
                         data = all_detailed_metrics, 
                         FUN = max)

# Rename columns for clarity
colnames(mean_metrics) <- c("Combination", "Mean_Uncertainty_chl", "Mean_Bias_chl", 
                            "Mean_Uncertainty_H", "Mean_Bias_H")

# Print the mean metrics for all combinations
print(mean_metrics)

library(GGally)
scatter_mat <- ggpairs(mean_metrics[, -1]) +
  ggtitle("Scatter Plot Matrix for Mean Uncertainty and Bias")

ggsave(paste0("./outputs/sens_scatter.png"), plot = scatter_mat, scale = 1.7, 
       width = 6, height = 4.5,
       units = "in",dpi = 300)


library(ggplot2)
library(reshape2)

# Melt the data for easy plotting
mean_metrics_melted <- melt(mean_metrics, id.vars = "Combination")

# Line plot
p <- ggplot(mean_metrics_melted, aes(x = as.factor(Combination), y = value, color = variable, 
                                     group = variable)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Trend of Mean Uncertainty and Bias Across Combinations",
       x = "Combination", y = "Value")



ggplot(mean_metrics_melted, aes(x = variable, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of Mean Uncertainty and Bias Metrics",
       x = "Metric", y = "Value")




library(dplyr)
library(ggplot2)
library(reshape2)

# Assuming 'Combination' is in the format "x_y" where x is Benthic Type and y is Spectral Slope
mean_metrics <- mean_metrics[-2381,] %>%
  mutate(Benthic_Type = as.numeric(sub("_.*", "", Combination)),
         Spectral_Slope = as.numeric(sub(".*_", "", Combination)))

# Melt the data for easier plotting
mean_metrics_melted <- melt(mean_metrics, id.vars = c("Benthic_Type", "Spectral_Slope"), 
                            measure.vars = c("Mean_Uncertainty_chl", "Mean_Bias_chl", 
                                             "Mean_Uncertainty_H", "Mean_Bias_H"))

# Plot the Uncertainty for each Benthic Type vs Spectral Slope
unc_benthic_group <- ggplot(mean_metrics_melted %>% filter(grepl("Uncertainty", variable)), 
                            aes(x = Spectral_Slope, y = value, color = as.factor(Benthic_Type), group = Benthic_Type)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y") + 
  labs(title = "Uncertainty for Each Benthic Type vs Spectral Slope",
       x = "Spectral Slope", y = "Uncertainty") +
  theme(legend.position = "right", legend.title = element_blank())

ggsave(paste0("./outputs/sens_benthic_slope_unc.png"), plot = unc_benthic_group, scale = 1.7, 
       width = 6, height = 4,
       units = "in",dpi = 300)

# Plot the Bias for each Benthic Type vs Spectral Slope
bias_benthic_group <- ggplot(mean_metrics_melted %>% filter(grepl("Bias", variable)), 
                             aes(x = Spectral_Slope, y = value, color = as.factor(Benthic_Type), group = Benthic_Type)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y") + 
  labs(title = "Bias for Each Benthic Type vs Spectral Slope",
       x = "Spectral Slope", y = "Bias") +
  theme(legend.position = "right", legend.title = element_blank())

ggsave(paste0("./outputs/sens_benthic_slope_bias.png"), plot = bias_benthic_group, scale = 1.7, 
       width = 6, height = 4,
       units = "in",dpi = 300)







# 
# # Function to process each observation
# process_observation <- function(obsdata) {
#   # Use pblapply to process all combinations for the given obsdata
#   results_list <- pblapply(1:nrow(combinations_grid), function(idx) {
#     process_combination(obsdata, idx)
#   }, cl = cl)
#   
#   # Combine results into a data frame
#   result_df <- as.data.frame(do.call(rbind, results_list))
#   return(result_df)
# }
# 
# # Export required objects to cluster
# clusterExport(cl, c("process_combination", "combinations_grid", "spectral_slope_vec", 
#                     "wavelength", "combinations", "interp_rb", "inverse_runBayes", 
#                     "upper.bound", "lower.bound", "type_Rrs_below", "samplerlist", "quiet"))
# 
# 
# # Create cluster
# num_cores <- detectCores() - 1  # Use one less than the total number of cores
# cl2 <- makeCluster(num_cores)
# registerDoParallel(cl2)
# 
# # Process each observation using foreach with %dopar%
# final_results <- foreach(i = names(rrs_sample_samp), 
#                          .combine = 'rbind', 
#                          .verbose = T,
#                          .packages = c('parallel', 'pbapply', 'BayesianTools', 'data.table'), 
#                          .export = c("process_combination", "combinations_grid", "spectral_slope_vec", 
#                                      "wavelength", "combinations", "interp_rb", "inverse_runBayes", 
#                                      "upper.bound", "lower.bound", "type_Rrs_below", "samplerlist", 
#                                      "quiet")) %dopar% {
#                                        source("./R/SABER_forward_fast.R")
#                                        source("./R/solve.objective.inverse_fast.R")
#                                        source("./R/mcmc_bayes_parallel.R")
#                                        source("./R/saber_inversion_final.R")
#                                        
#                                        apply(rrs_sample_samp[[i]], 1, process_observation)
# }
# 
# # Stop the cluster
# stopCluster(cl)
# 
# # Print the final results
# print(final_results)
