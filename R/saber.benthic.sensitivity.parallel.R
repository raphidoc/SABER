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
results_list <- foreach(i = 1:nrow(obsdata), 
                        .verbose = T,
                        #.combine = 'rbind', 
                        .packages = c('data.table', 'BayesianTools')) %:% 
  foreach(idx_c = 1:nrow(combinations_grid), 
          #.combine = 'rbind', 
          .verbose = T, .packages = c('data.table', 'BayesianTools')) %dopar% {
            source("./R/SABER_forward_fast.R")
            source("./R/solve.objective.inverse_fast.R")
            source("./R/mcmc_bayes_parallel.R")
            source("./R/saber_inversion_final.R")
            #apply(obs_group, 1, process_combination, idx = idx)
            process_combination(obsdata_row = obsdata[i,], idx = idx_c)
    
  }

# Stop the main cluster
stopCluster(cl)

# Print the final results
print(results_list)


# Function to process each observation
process_observation <- function(obsdata) {
  # Use pblapply to process all combinations for the given obsdata
  results_list <- pblapply(1:nrow(combinations_grid), function(idx) {
    process_combination(obsdata, idx)
  }, cl = cl)
  
  # Combine results into a data frame
  result_df <- as.data.frame(do.call(rbind, results_list))
  return(result_df)
}

# Export required objects to cluster
clusterExport(cl, c("process_combination", "combinations_grid", "spectral_slope_vec", 
                    "wavelength", "combinations", "interp_rb", "inverse_runBayes", 
                    "upper.bound", "lower.bound", "type_Rrs_below", "samplerlist", "quiet"))


# Create cluster
num_cores <- detectCores() - 1  # Use one less than the total number of cores
cl2 <- makeCluster(num_cores)
registerDoParallel(cl2)

# Process each observation using foreach with %dopar%
final_results <- foreach(i = names(rrs_sample_samp), 
                         .combine = 'rbind', 
                         .verbose = T,
                         .packages = c('parallel', 'pbapply', 'BayesianTools', 'data.table'), 
                         .export = c("process_combination", "combinations_grid", "spectral_slope_vec", 
                                     "wavelength", "combinations", "interp_rb", "inverse_runBayes", 
                                     "upper.bound", "lower.bound", "type_Rrs_below", "samplerlist", 
                                     "quiet")) %dopar% {
                                       source("./R/SABER_forward_fast.R")
                                       source("./R/solve.objective.inverse_fast.R")
                                       source("./R/mcmc_bayes_parallel.R")
                                       source("./R/saber_inversion_final.R")
                                       
                                       apply(rrs_sample_samp[[i]], 1, process_observation)
}

# Stop the cluster
stopCluster(cl)

# Print the final results
print(final_results)
