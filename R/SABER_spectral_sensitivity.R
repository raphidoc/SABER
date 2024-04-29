# Function :: create multiple partitioned lists for given observation  ----
partition_data_for_multiprocessing <- function(num_groups, data_input){
  
  num_groups = num_groups
  
  library(tidyr)
  #Decompose the concatenated data into lists
  data_split = data_input %>% 
    group_by((row_number()-1) %/% (n()/num_groups)) %>%
    nest %>% pull(data)
  
  names(data_split) = c(paste0((seq(1,num_groups,1))))
  
  return(data_split)
  
}
# Function :: create boxplot showing the error range for each retreival ----
spectral_senstv_boxplot <- function(title_lab = expression(paste("[",italic("chl"), "]")),
                                    input_df = chl_df_sd_long, 
                                    xlab = c("WISE", "MSI", "OLI"),
                                    ymax_error = ceiling(max(adg443_df_sd_long$value) + 0.15*max(adg443_df_sd_long$value))){
  
  df_means <- plyr::ddply(input_df, "variable", 
                          dplyr::summarise, mean_rrs = mean(value, na.rm = T))
  
  g_box <- ggplot(data = input_df, aes(x= reorder((variable), variable)#,as.character(variable)
                                       ))+
    geom_boxplot(aes(y = value, fill = variable),color = "black", alpha = 0.6)+
    geom_point(aes(y = value, shape = variable, fill = variable, color = variable), size = 1.5, 
               alpha = 0.6)+
    scale_x_discrete(name = " ",
                     labels = xlab) +
    
    geom_line(data = df_means, aes(x=as.character(variable), y=as.numeric(mean_rrs),
                                   group=1
    ), linetype = "dashed",
    size=1.3)+
    
    scale_color_viridis(discrete = T, guide = "none")+
    scale_fill_viridis(discrete = T, guide = "none")+
    scale_shape_manual(name = "", values = c(21,22,23), guide = "none")+
    ylim(0,ymax_error)+
    ylab(expression(paste("Standard Deviation (", sigma, ")")))+
    theme_bw()+
    ggtitle(title_lab)+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 15, color = 'black', angle = 0),
          axis.text.y = element_text(size = 15, color = 'black', angle = 0),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.ticks.length = unit(.25, "cm"),
          legend.position=c(0.70, 0.98),
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 20, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5,
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey",
                                          size = 0.5, linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.0,0.5,0.0,0.0), "cm"),
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  return(g_box)
}

# Function :: Convolve the input array to a given RSRF ----
convolve_RSRF <- function(spec_mat, rsrf_mat){
  
  convolve_spec_df = as.data.frame(rho::specconv(spec = spec_mat, rsrf = rsrf_mat))
  
  return(convolve_spec_df)
  
}

# Function :: Convolve the input matrix to a given RSRF ----
convolve_RSRF_vectorized <- function(row, wavelength_obs, rsrf){
  
  obsdata_mat = as.matrix(data.frame("wave" = wavelength_obs, "rrs"= row))
  
  obsdata_convolv_mat = convolve_RSRF(spec_mat = obsdata_mat, rsrf_mat = rsrf)
  obsdata_convolv = obsdata_convolv_mat$rrs[!is.na(obsdata_convolv_mat$rrs)]
  
  return(obsdata_convolv)
  
  
}

# Function :: Interpolate the input matrix to a given wavelength ----
interp_vectorized <- function(row, wavelength_obs, wavelength_out){
  
  obsdata_interp = approx(x= wavelength_obs, y = row, xout = wavelength_out, method = "linear")$y
  
  return(obsdata_interp)
  
  
}

# Function :: Interpolate the bottom reflectance to Rrs wavelengths ----
interp_rb <- function(wavelength_input, bottom_type = c("Mud", "Sand", "Eelgrass")){
  load("./data/EGSL_bottom.RData")
  
  column_numbers <- which(names(rb) %in% bottom_type)
  
  rb = rb[c(1,column_numbers)]
  names(rb) = c("wavelength", "class1", "class2", "class3")
  
  if (!(all(length(rb$wavelength) == length(wavelength_input)) && all(rb$wavelength == wavelength_input))) {
    
    print(paste0("Rb wavelength has " , length(rb$wavelength), " bands; Rrs wavelength has ", length(wavelength_input) ))
    
    rb_interp = data.frame("wavelength" = wavelength_input,
                           "class1" = approx(rb$wavelength, rb$class1, xout = wavelength_input)$y,
                           "class2" = approx(rb$wavelength, rb$class2, xout = wavelength_input)$y,
                           "class3" = approx(rb$wavelength, rb$class3, xout = wavelength_input)$y
    )
    
  }
  rb = rb_interp
  rm(rb_interp)
  
  assign(x = "rb", value = rb, envir = .GlobalEnv)
}



# Create the noise perturbed randomized vector of input variables for Rrs simulation ----
library(sensitivity)
set.seed(123)
vec_size <- seq(1, 500, 1)

param_vec <- data.frame(
  "chl" = sample(seq(from = min(0.5), to = max(30), length.out = length(vec_size))),
  "adg443" = sample(seq(from = min(0.1), to = max(2), length.out = length(vec_size))),
  "bbp555" = sample(seq(from = min(0.002), to = max(0.01), length.out = length(vec_size))),
  "H" = rnorm.trunc(n = vec_size, mean = 5, sd = 4, min = 0.5, max = 12),
  "fa1" = runif(n = length(vec_size), min = 0, max = 1),
  "fa2" = runif(n = length(vec_size), min = 0, max = 1),
  "fa3" = NA  # Initialize fa3 column
)

# Adjust fa1 and fa2 generation to ensure sum of fa1 and fa2 does not exceed 1
param_vec$fa1[param_vec$fa1 + param_vec$fa2 > 1] <- 1 - param_vec$fa2[param_vec$fa1 + param_vec$fa2 > 1]
param_vec$fa2[param_vec$fa1 + param_vec$fa2 > 1] <- 1 - param_vec$fa1[param_vec$fa1 + param_vec$fa2 > 1]

# Generate fa3 based on remaining fraction
param_vec$fa3 <- 1 - param_vec$fa1 - param_vec$fa2

# Verify that sum of fa1, fa2, and fa3 <= 1
sum_check <- with(param_vec, sum(fa1, fa2, fa3))
cat("Sum of fa1, fa2, and fa3:", sum_check, "\n")
head(param_vec)

#param_vec.LUT <- expand.grid(param_vec) 

wavelength = ceiling(seq(from=min(400), to=max(750), 
                         length.out=(80)))

interp_rb(wavelength_input = wavelength)

#Obtain the Rrs vector rom the model param vector
rrs.forward.SABER <- matrix(nrow = length(param_vec$chl),
                            ncol = length(wavelength),0)

rrs.forward.SABER = t(apply(param_vec, 1, function(row){
  
  Saber_forward_fast_sensitivity_test(
    use_true_IOPs = F, use_manual_slope = F, 
    chl = row[1],
    a_dg = row[2],
    bbp.550 = row[3], 
    slope.parametric = F, 
    Rrs_input_for_slope = NA,
    
    z = row[4],
    rb.fraction = row[5:7],
    
    wavelength = wavelength,
    verbose = F
  )
  
}

))
#Convert to a matrix data-frame
colnames(rrs.forward.SABER) = wavelength
rownames(rrs.forward.SABER) = vec_size

rrs.forward.SABER = as.data.frame(rrs.forward.SABER)

rrs_synth = as.data.frame(rrs.forward.SABER)

rrs_synth$id = seq(1, nrow(rrs_synth), 1)

rrs_synth = reshape2::melt(rrs_synth, id.vars = "id")

H_plot_df = data.frame("id" = seq(1,length(param_vec$H), 1), "Depth" = param_vec$H)
H_long = reshape2::melt(H_plot_df, id.vars = "id")

rrs_synth_H = merge(rrs_synth,H_long,by=c("id"))

names(rrs_synth_H) = c("id", "wavelength", "rrs", "shallow_tag", "depth")

xmin <- 400; xmax <- 750;  xstp <- 50
xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0; ymax <- signif(max(rrs_synth_H$rrs) + 0.1*max(rrs_synth_H$rrs), digits = 1) 
ystp <- ymax/4
ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda)^"shallow"[italic("synth")], "[sr"^{-1}, "]"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.70, 0.98)


library(viridis)

g <- ggplot(data = rrs_synth_H) + 
  geom_line(aes(x = as.numeric(as.character(wavelength)), y = rrs, colour = depth, 
                group = id), size = 1.3, show.legend = T, alpha = 0.6)+
  scale_colour_viridis(discrete = F, name = "Depth [m]", direction = -1) +
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax)
              ,expand = FALSE, clip = "on"
  ) +
  #scale_colour_manual(name = "",values = rev(collist))+
  #guides(fill = FALSE)+
  guides( 
    fill = guide_legend(title = "Depth (m)", reverse = T))+
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp))  +
  #scale_y_log10()+
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))  +
  theme_bw()+
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 25, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
        legend.direction = "vertical",
        legend.title = element_text(colour = "black", size = 15, face = "plain"),
        legend.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.0,0.9,0.0,0.0), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

g
ggsave("./outputs/shallow_synthetic_rrs.png", plot = g,scale = 1.7, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)


rrs_shallow_input = data.table(rrs.forward.SABER)

# Convolve to Multispectral wavelengths ----

## LANDSAT-8 OLI ----
l8_srf_meta = read_excel_allsheets("./data/Ball_BA_RSR.v1.2.xlsx")


l8_rsrf = data.frame("wavelength" = seq(400, 2500, 1), "L8_ca" = 0, "L8_blue" = 0,
                     "L8_green" = 0, "L8_red" = 0, "L8_nir" = 0, "L8_cirrus" = 0,  
                     "L8_swir1" = 0, "L8_swir2" = 0)

for (i in 1:length(names(l8_rsrf)[-1])) {
  
  idx <- grep(#paste0("*",paste(l8_srf_meta[[i+1]]$Wavelength, collapse = "|"), "^"), 
    paste("\\b", paste(l8_srf_meta[[i+1]]$Wavelength, collapse = "\\b|\\b"), "\\b", sep = ""),
    l8_rsrf$wavelength)
  
  l8_rsrf[idx,i+1] = l8_srf_meta[[i+1]]$`BA RSR [watts]`
  
}

l8_rsrf[l8_rsrf < 0] = 0
l8_rsrf = as.matrix(l8_rsrf)

rrs_shallow_OLI = apply(rrs_shallow_input, 1, convolve_RSRF_vectorized, 
                        wavelength_obs = wavelength, rsrf = l8_rsrf)

rrs_shallow_OLI = as.data.frame(t(rrs_shallow_OLI))

wavelength_OLI = c(443, 482, 561, 655)

## Sentinel-2 MSI---- 
sentinel_rsrf = read_excel_allsheets("./data/S2-SRF_COPE-GSEG-EOPG-TN-15-0007_3.2.xlsx")
sentinel_rsrf_mat = as.matrix(sentinel_rsrf[[2]])

rrs_shallow_MSI = apply(rrs_shallow_input, 1, convolve_RSRF_vectorized, 
                        wavelength_obs = wavelength, rsrf = sentinel_rsrf_mat)

rrs_shallow_MSI = as.data.frame(t(rrs_shallow_MSI))
wavelength_MSI = c(443, 490, 560, 665, 705, 740)

# COP-S in situ radiometer ----
wavelength_cops = Cops::detection.limit$waves
wavelength_cops = wavelength_cops[wavelength_cops > 400 & wavelength_cops <= 750]

rrs_shallow_cops = apply(rrs_shallow_input, 1, interp_vectorized, wavelength_obs = wavelength,
                         wavelength_out = wavelength_cops)
rrs_shallow_cops = as.data.frame(t(rrs_shallow_cops))


# Spectral plot of Rrs for WISE, MSI and OLI ----
plot(wavelength, as.numeric(rrs_shallow_input[345,]), type = "l", lwd = 2)
lines(wavelength_MSI, as.numeric(rrs_shallow_MSI[345,]), col = "red3", lwd = 2)
lines(wavelength_OLI, as.numeric(rrs_shallow_OLI[345,]), col = "navyblue", lwd = 2)

xmin <- 400; xmax <- 750;  xstp <- 50
xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0; ymax <- 0.0025
ystp <- ymax/5
ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda)^"shallow"[italic("synth")], "[sr"^{-1}, "]"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.70, 0.98)


g <- ggplot()+
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax)
              ,expand = FALSE, clip = "on"
  ) +
  geom_line(data = data.frame("wavelength"= wavelength, "WISE" = as.numeric(rrs_shallow_input[345,])),
            aes(x = wavelength, y = WISE, color = "WISE"), show.legend = T, size = 1.1)+
  
  geom_line(data = data.frame("wavelength"= wavelength_MSI, "MSI" = as.numeric(rrs_shallow_MSI[345,])),
            aes(x = wavelength, y = MSI, color = "MSI"), show.legend = T, size=1.1)+
  geom_point(data = data.frame("wavelength"= wavelength_MSI, "MSI" = as.numeric(rrs_shallow_MSI[345,])),
            aes(x = wavelength, y = MSI, color = "MSI"), show.legend = T, size=2)+
  
  geom_line(data = data.frame("wavelength"= wavelength_OLI, "OLI" = as.numeric(rrs_shallow_OLI[345,])),
             aes(x = wavelength, y = OLI, color = "OLI"), show.legend = T, size=1.1)+
  geom_point(data = data.frame("wavelength"= wavelength_OLI, "OLI" = as.numeric(rrs_shallow_OLI[345,])),
            aes(x = wavelength, y = OLI, color = "OLI"), show.legend = T, size=2)+
  
  scale_color_manual(name = "", values = rev(c("black", "navyblue", "goldenrod2")))+
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp))  +
  #scale_y_log10()+
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))  +
  theme_bw()+
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 25, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
        legend.direction = "vertical",
        legend.title = element_text(colour = "black", size = 15, face = "plain"),
        legend.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.0,0.9,0.0,0.0), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

g
ggsave("./outputs/shallow_rrs_WISE_MSI_OLI.png", plot = g,scale = 1.7, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)

# Run SABER inversion on the simulated Rrs ----

## Partition Rrs with parallel processing ----
rrs_hs_list = partition_data_for_multiprocessing(num_groups = 6, data_input = rrs_shallow_input)
rrs_cops_list = partition_data_for_multiprocessing(num_groups = 6, data_input = rrs_shallow_cops)
rrs_MSI_list = partition_data_for_multiprocessing(num_groups = 6, data_input = rrs_shallow_MSI)
rrs_OLI_list = partition_data_for_multiprocessing(num_groups = 6, data_input = rrs_shallow_OLI)

samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")

## Invert Rrs with parallel processing ----

# run_inverse_mode(rrs_input = as.numeric(rrs_shallow_input[1,]) #+  rnorm(n = length(obsdata_cops), mean = 0, 
#                                                   #sd = sd(obsdata_cops))*0.05
#                                             , 
#                  wavelength_input = wavelength)

## For Hyperspectral Rrs ----
interp_rb(wavelength_input = wavelength)
numCores <- detectCores()
print(paste0("Number of system processor found: ", numCores))

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
print(paste0("The Cluster is created with ", length(cl), " number of processors"))

go_Bayes = T
if (go_Bayes == T) {
  
  start_time = Sys.time()
  
  ### Pure  MCMC inversion
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### UNCONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  # Apply the function to each element in the lists using %dopar%
  output_list_hs <- foreach(i = names(rrs_hs_list)
                         ,
                         #.export = c("par0"),
                         .verbose = T,
                         #.combine = rbind,
                         .packages="data.table"
  ) %dopar% {
    source("./R/SABER_forward_fast.R")
    source("./R/solve.objective.inverse_fast.R")
    source("./R/mcmc_bayes_parallel.R")
    #.GlobalEnv$par0 <- par0
    apply(rrs_hs_list[[i]], 1, inverse_runBayes,
          
          # rrs_type = "shallow", max_par = c(sapply(param_vec[-4], max), 12, 1,1,1,10), 
          # min_par = c(sapply(param_vec, min), 1,1,1,10), constrain_config = c(F, F, F), 
          # qaa_slope = F, manual_slope = F, iter_count = 25000, sampler_mcmc = samplerlist[6], 
          # wavelngth_sim = wavelength, sa_model = "am03", hybrid_mode = F
          
          rrs_type =  "shallow", 
          max_par = c(30,2,0.010, 12, 1,1,1,10), 
          min_par = c(0.5,0.1,0.002, 0.5, 0,0,0,0.00001), 
          param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
          constrain_config = c(F,F,F),
          bbp_const_val = 0.007, 
          bgc_const_val = c(5,1.003, 0.007), 
          
          iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                             "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
          
          qaa_slope = F, 
          manual_slope = F, 
          manual_slope_val = c(0.014, 0.003, 0.5), 
          iter_count = 20000, 
          sampler_mcmc = samplerlist[6], 
          wavelngth_sim = wavelength, 
          sa_model = "am03", 
          hybrid_mode = F, plot_rrs = F
          
    )
  }
  
  end_time = Sys.time(); time_taken <- end_time - start_time
  print(time_taken)
}

stopCluster(cl)


#Save the inversion results as data-frame
fit_results_hs = data.frame(t(do.call(cbind, output_list_hs)))

names(fit_results_hs) = c("chl", "adg443", "bbp555", "H", paste0("fa", (1:rb_count)), "sigma_pop",
                       "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", paste0("sd_fa", (1:rb_count)), 
                       "sd_sigma_pop")

write.csv(x = fit_results_hs, file = "./outputs/shallow_synth_HS_bayes.csv",
          quote = F, col.names = T, sep = ",")

## For COP-S Rrs ----
interp_rb(wavelength_input = wavelength_cops)

numCores <- detectCores()
print(paste0("Number of system processor found: ", numCores))

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
print(paste0("The Cluster is created with ", length(cl), " number of processors"))

go_Bayes = T
if (go_Bayes == T) {
  
  start_time = Sys.time()
  
  ### Pure  MCMC inversion
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### UNCONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  # Apply the function to each element in the lists using %dopar%
  output_list_cops <- foreach(i = names(rrs_cops_list)
                         ,
                         #.export = c("par0"),
                         .verbose = T,
                         #.combine = rbind,
                         .packages="data.table"
  ) %dopar% {
    source("./R/SABER_forward_fast.R")
    source("./R/solve.objective.inverse_fast.R")
    source("./R/mcmc_bayes_parallel.R")
    #.GlobalEnv$par0 <- par0
    apply(rrs_cops_list[[i]], 1, inverse_runBayes,
          
          # rrs_type = "shallow", max_par = c(sapply(param_vec[-4], max), 12, 1,1,1,10), 
          # min_par = c(sapply(param_vec, min), 1,1,1,10), constrain_config = c(F, F, F), 
          # qaa_slope = F, manual_slope = F, iter_count = 25000, sampler_mcmc = samplerlist[6], 
          # wavelngth_sim = wavelength, sa_model = "am03", hybrid_mode = F
          
          rrs_type =  "shallow", 
          max_par = c(30,2,0.010, 12, 1,1,1,10), 
          min_par = c(0.5,0.1,0.002, 0.5, 0,0,0,0.00001),
          param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
          constrain_config = c(F,F,F),
          bbp_const_val = 0.007, 
          bgc_const_val = c(5,1.003, 0.007), 
          
          iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                             "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
          
          qaa_slope = F, 
          manual_slope = F, 
          manual_slope_val = c(0.014, 0.003, 0.5), 
          iter_count = 20000, 
          sampler_mcmc = samplerlist[6], 
          wavelngth_sim = wavelength_cops, 
          sa_model = "am03", 
          hybrid_mode = F, plot_rrs = F
          
    )
  }
  
  end_time = Sys.time(); time_taken <- end_time - start_time
  print(time_taken)
}

stopCluster(cl)


#Save the inversion results as data-frame
fit_results_cops = data.frame(t(do.call(cbind, output_list_cops)))

names(fit_results_cops) = c("chl", "adg443", "bbp555", "H", paste0("fa", (1:rb_count)), "sigma_pop",
                       "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", paste0("sd_fa", (1:rb_count)), 
                       "sd_sigma_pop" )
write.csv(x = fit_results_cops, file = "./outputs/shallow_synth_COPS_bayes.csv",
          quote = F, col.names = T, sep = ",")


## For MSI Rrs ----
interp_rb(wavelength_input = wavelength_MSI)
numCores <- detectCores()
print(paste0("Number of system processor found: ", numCores))

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
print(paste0("The Cluster is created with ", length(cl), " number of processors"))

go_Bayes = T
if (go_Bayes == T) {
  
  start_time = Sys.time()
  
  ### Pure  MCMC inversion
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### UNCONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  # Apply the function to each element in the lists using %dopar%
  output_list_MSI <- foreach(i = names(rrs_MSI_list)
                         ,
                         #.export = c("par0"),
                         .verbose = T,
                         #.combine = rbind,
                         .packages="data.table"
  ) %dopar% {
    source("./R/SABER_forward_fast.R")
    source("./R/solve.objective.inverse_fast.R")
    source("./R/mcmc_bayes_parallel.R")
    #.GlobalEnv$par0 <- par0
    apply(rrs_MSI_list[[i]], 1, inverse_runBayes,
          
          # rrs_type = "shallow", max_par = c(sapply(param_vec[-4], max), 12, 1,1,1,10), 
          # min_par = c(sapply(param_vec, min), 1,1,1,10), constrain_config = c(F, F, F), 
          # qaa_slope = F, manual_slope = F, iter_count = 25000, sampler_mcmc = samplerlist[6], 
          # wavelngth_sim = wavelength, sa_model = "am03", hybrid_mode = F
          
          rrs_type =  "shallow", 
          max_par = c(30,2,0.010, 12, 1,1,1,10), 
          min_par = c(0.5,0.1,0.002, 0.5, 0,0,0,0.00001), 
          param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
          constrain_config = c(F,F,F),
          bbp_const_val = 0.007, 
          bgc_const_val = c(5,1.003, 0.007), 
          
          iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                             "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
          
          qaa_slope = F, 
          manual_slope = F, 
          manual_slope_val = c(0.014, 0.003, 0.5), 
          iter_count = 20000, 
          sampler_mcmc = samplerlist[6], 
          wavelngth_sim = wavelength_MSI, 
          sa_model = "am03", 
          hybrid_mode = F, plot_rrs = F
          
    )
  }
  
  end_time = Sys.time(); time_taken <- end_time - start_time
  print(time_taken)
}

stopCluster(cl)


#Save the inversion results as data-frame
fit_results_MSI = data.frame(t(do.call(cbind, output_list_MSI)))

names(fit_results_MSI) = c("chl", "adg443", "bbp555", "H", paste0("fa", (1:rb_count)), "sigma_pop",
                       "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", paste0("sd_fa", (1:rb_count)), 
                       "sd_sigma_pop" )

write.csv(x = fit_results_MSI, file = "./outputs/shallow_synth_MSI_bayes.csv",
          quote = F, col.names = T, sep = ",")


## For OLI Rrs ----
interp_rb(wavelength_input = wavelength_OLI)
numCores <- detectCores()
print(paste0("Number of system processor found: ", numCores))

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
print(paste0("The Cluster is created with ", length(cl), " number of processors"))

go_Bayes = T
if (go_Bayes == T) {
  
  start_time = Sys.time()
  
  ### Pure  MCMC inversion
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### UNCONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  # Apply the function to each element in the lists using %dopar%
  output_list_OLI <- foreach(i = names(rrs_OLI_list)
                             ,
                             #.export = c("par0"),
                             .verbose = T,
                             #.combine = rbind,
                             .packages="data.table"
  ) %dopar% {
    source("./R/SABER_forward_fast.R")
    source("./R/solve.objective.inverse_fast.R")
    source("./R/mcmc_bayes_parallel.R")
    #.GlobalEnv$par0 <- par0
    apply(rrs_OLI_list[[i]], 1, inverse_runBayes,
          
          # rrs_type = "shallow", max_par = c(sapply(param_vec[-4], max), 12, 1,1,1,10), 
          # min_par = c(sapply(param_vec, min), 1,1,1,10), constrain_config = c(F, F, F), 
          # qaa_slope = F, manual_slope = F, iter_count = 25000, sampler_mcmc = samplerlist[6], 
          # wavelngth_sim = wavelength, sa_model = "am03", hybrid_mode = F
          
          rrs_type =  "shallow", 
          max_par = c(30,2,0.010, 12, 1,1,1,10), 
          min_par = c(0.5,0.1,0.002, 0.5, 0,0,0,0.00001), 
          param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
          constrain_config = c(F,F,F),
          bbp_const_val = 0.007, 
          bgc_const_val = c(5,1.003, 0.007), 
          
          iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                             "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
          
          qaa_slope = F, 
          manual_slope = F, 
          manual_slope_val = c(0.014, 0.003, 0.5), 
          iter_count = 20000, 
          sampler_mcmc = samplerlist[6], 
          wavelngth_sim = wavelength_OLI, 
          sa_model = "am03", 
          hybrid_mode = F, plot_rrs = F
          
    )
  }
  
  end_time = Sys.time(); time_taken <- end_time - start_time
  print(time_taken)
}

stopCluster(cl)


#Save the inversion results as data-frame
fit_results_OLI = data.frame(t(do.call(cbind, output_list_OLI)))

names(fit_results_OLI) = c("chl", "adg443", "bbp555", "H", paste0("fa", (1:rb_count)), "sigma_pop",
                           "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", paste0("sd_fa", (1:rb_count)), 
                           "sd_sigma_pop" )

write.csv(x = fit_results_OLI, file = "./outputs/shallow_synth_OLI_bayes.csv",
          quote = F, col.names = T, sep = ",")

# Create Box-Plots of standard error from inversion ----

# Read inversion results
fit_results_hs = read.csv("./outputs/shallow_synth_HS_bayes.csv", header = T)
fit_results_cops = read.csv("./outputs/shallow_synth_OLI_bayes.csv", header = T)
fit_results_MSI = read.csv("./outputs/shallow_synth_MSI_bayes.csv", header = T)

## chl ----
chl_df = data.frame("chl_real" = param_vec$chl, "chl_hs" = fit_results_hs$chl,
                    "chl_msi" = fit_results_MSI$chl, "chl_cops" = fit_results_cops$chl,
                    
                    "chl_sd_hs" = fit_results_hs$sd_chl,
                    "chl_sd_msi" = fit_results_MSI$sd_chl, "chl_sd_cops" = fit_results_cops$sd_chl)
chl_df_sd = chl_df[(5:7)]
chl_df_sd$id = seq(1,500,1)

chl_df_sd_long = reshape2::melt(chl_df_sd, id.vars = "id")

box_chl = spectral_senstv_boxplot(title_lab = expression(paste("[",italic("chl"), "]")),
                                  input_df = chl_df_sd_long,ymax_error = 10 )
ggsave(paste0("./outputs/sensitivity_box_plot_chl.png"), plot = box_chl,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

## adg443 ----
adg443_df = data.frame("adg443_real" = param_vec$adg443, "adg443_hs" = fit_results_hs$adg443,
                     "adg443_msi" = fit_results_MSI$adg443, "adg443_cops" = fit_results_cops$adg443,
                    
                    "adg443_sd_hs" = fit_results_hs$sd_adg443, 
                    "adg443_sd_msi" = fit_results_MSI$sd_adg443, 
                    "adg443_sd_cops" = fit_results_cops$sd_adg443)

adg443_df_sd = adg443_df[(5:7)]
adg443_df_sd$id = seq(1,500,1)

adg443_df_sd_long = reshape2::melt(adg443_df_sd, id.vars = "id")
box_adg443 = spectral_senstv_boxplot(title_lab = expression(paste(italic("a")["dg"](443))),
                                  input_df = adg443_df_sd_long, ymax_error = 0.25 )

ggsave(paste0("./outputs/sensitivity_box_plot_adg443.png"), plot = box_adg443,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

## bbp555 ----
bbp555_df = data.frame("bbp555_real" = param_vec$bbp555, "bbp555_hs" = fit_results_hs$bbp555,
                     "bbp555_msi" = fit_results_MSI$bbp555, "bbp555_cops" = fit_results_cops$bbp555,
                    
                    "bbp555_sd_hs" = fit_results_hs$sd_bbp555, 
                    "bbp555_sd_msi" = fit_results_MSI$sd_bbp555,
                    "bbp555_sd_cops" = fit_results_cops$sd_bbp555)
bbp555_df_sd = bbp555_df[(5:7)]
bbp555_df_sd$id = seq(1,500,1)

bbp555_df_sd_long = reshape2::melt(bbp555_df_sd, id.vars = "id")
box_bbp555 = spectral_senstv_boxplot(title_lab = expression(paste(italic("b")["bp"](555))),
                                     input_df = bbp555_df_sd_long, ymax_error = 0.003 )

ggsave(paste0("./outputs/sensitivity_box_plot_bbp555.png"), plot = box_bbp555,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

## H ----
H_df = data.frame("H_real" = param_vec$H, "H_hs" = fit_results_hs$H,
                     "H_msi" = fit_results_MSI$H, "H_cops" = fit_results_cops$H,
                    
                    "H_sd_hs" = fit_results_hs$sd_H, 
                    "H_sd_msi" = fit_results_MSI$sd_H,
                  "H_sd_cops" = fit_results_cops$sd_H)
H_df_sd = H_df[(5:7)]
H_df_sd$id = seq(1,500,1)

H_df_sd_long = reshape2::melt(H_df_sd, id.vars = "id")
box_H = spectral_senstv_boxplot(title_lab = expression(paste("[",italic("H"), "]")),
                                  input_df = H_df_sd_long, ymax_error = 2 )
ggsave(paste0("./outputs/sensitivity_box_plot_H.png"), plot = box_H,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)



# library(ggforce)
# 
# g <- ggplot(data = chl_df_sd_long, aes(x=(as.character(variable)), y=value)) +
#   
#   
#   geom_boxplot(aes(fill=(as.character(variable)), 
#                    x=(as.character(variable)), 
#                    y=value, 
#                    group=(as.character(variable))), show.legend = FALSE)+
#   
#   geom_point(aes(y = value, shape = variable, color = variable), size = 1.3)+
#   
#   geom_line(data = df_means, aes(x=as.character(variable), y=as.numeric(mean_rrs),
#                                  group=1
#   ), linetype = "dashed",
#   size=1.3)+
#   
#   labs(x=" ", y=expression(paste("Standard Error"))) +
#   ggtitle(" ")+
#   
#   scale_fill_viridis(discrete = T)+
#   scale_color_viridis(discrete = T)+
#   
#   theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         plot.caption = element_text(face="bold",size=15, hjust = 0, color = 'black'),
#         axis.text.x = element_text(size = 25, color = 'black', angle = 0),
#         axis.text.y = element_text(size = 25, color = 'black', angle = 0),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         axis.ticks.length = unit(.25, "cm"),
#         legend.position=c(0.03, 0.95),
#         legend.direction = "vertical",
#         legend.title = element_text(colour = "black", size = 15, face = "plain"),
#         legend.text = element_text(colour = "black", size = 15, face = "plain"),
#         legend.background = element_rect(fill = NA, size = 0.5,
#                                          linetype = "solid", colour = 0),
#         legend.key = element_blank(),
#         legend.justification = c("left", "top"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour = "grey",
#                                         size = 0.5, linetype = "dotted"),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
# g
# 
# ggsave(paste0("./outputs/seadoo_penetration_depth.png"), plot = g,
#        scale = 1.5, width = 5, height = 4.5, units = "in",dpi = 300)
