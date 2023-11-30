library(tidyr)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
source("./R/saber.inversion.vector.parallel.R")
#-----------------------------------------------------------------------------------------
# Create the noise perturbed randomized vector of model parameters
#-----------------------------------------------------------------------------------------
set.seed(123)
vec_size = seq(1,500,1)
param_vec <- data.frame("chl"=sample(seq(from=min(0.5), to=max(30), 
                                  length.out=length(vec_size))),
                      "adg443"=sample(seq(from=min(0.1), to=max(2), 
                                   length.out=length(vec_size))),
                      "bbp555"=sample(seq(from=min(0.002), to=max(0.01), 
                                   length.out=length(vec_size))),
                      "H" = rnorm.trunc(n = 500, mean = 5, sd = 4, min = 0.5, max = 12)
                      # "H"=(seq(from=min(0.5), to=max(8), 
                      #         length.out=length(vec_size))) + truncnorm::rtruncnorm(500, 
                      #                                                         a = 0.25, b = 2,
                      #                                                         mean = 2, sd = 1)
                      )

#param_vec.LUT <- expand.grid(param_vec) 


# param_vec <- data.frame("chl"=sample(seq(from=min(0.5), to=max(30), 
#                                          length.out=length(vec_size))),
#                         "adg443"=sample(seq(from=min(0.1), to=max(2), 
#                                             length.out=length(vec_size))),
#                         "bbp555"=sample(seq(from=min(0.002), to=max(0.01), 
#                                             length.out=length(vec_size))),
#                         "H"=(seq(from=min(0.5), to=max(8), 
#                                  length.out=length(vec_size))) + truncnorm::rtruncnorm(500, 
#                                                                                        a = 0.25, b = 2,
#                                                                                        mean = 2, sd = 1)
# )

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
    rb.fraction = fA.set,
    
    wavelength = wavelength,
    verbose = T
  )
  
}
  
))
#Convert to a matrix
rrs.forward.SABER = do.call(rbind, rrs.forward.SABER)
colnames(rrs.forward.SABER) = wavelength
rownames(rrs.forward.SABER) = vec_size

#Plot the simulated Rrs
matplot(wavelength, t(rrs.forward.SABER), col = viridis(500), type = "l", lwd=3)


rrs_synth = as.data.frame(rrs.forward.SABER)

rrs_synth$id = seq(1, nrow(rrs_synth), 1)

rrs_synth = reshape2::melt(rrs_synth, id.vars = "id")

H_plot_df = data.frame("id" = seq(1,length(param_vec$H), 1), "Depth" = param_vec$H)
H_long = reshape2::melt(H_plot_df, id.vars = "id")

rrs_synth_H = merge(rrs_synth,H_long,by=c("id"))

names(rrs_synth_H) = c("id", "wavelength", "rrs", "shallow_tag", "depth")

xmin <- 400; xmax <- 700;  xstp <- 50
xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0; ymax <- signif(max(rrs_synth_H$rrs) + 0.1*max(rrs_synth_H$rrs), digits = 1) 
ystp <- ymax/4
ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda)^"shallow"[italic("synth")], "(sr"^{-1}, ")"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.70, 0.98)


library(viridis)

g <- ggplot(data = rrs_synth_H) + 
  geom_line(aes(x = as.numeric(as.character(wavelength)), y = rrs, colour = depth, 
                group = id), size = 1.3, show.legend = T, alpha = 0.6)+
  scale_colour_viridis(discrete = F, name = "Depth (m)", direction = -1) +
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
ggsave("./outputs/synthetic_rrs.png", plot = g,scale = 1.7, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)

#-----------------------------------------------------------------------------------------
# Prepare the optimization parameters
#-----------------------------------------------------------------------------------------
pop.sd = "unknown" 

constrain.bbp= F
constrain.shallow = F
constrain.shallow.bgc = F ; constrain.shallow.iop = F

manual_par0 = T ; manual_bound = T

#---------------------------------------
# Initial values for shallow water
#---------------------------------------
rb_count = length(fA.set)
rep(0.5, rb_count)
rb_c = paste0("rb", (1:rb_count))

#Instantiate inversion objective function and the optimization scheme 
pop.sd = "unknown" ;constrain.bbp= F
constrain.shallow.bgc = F ; constrain.shallow.iop = F; manual_par0 = T
type_Rrs_below = "shallow"

inv_bound = create_init_bound(rrs_inv_type = type_Rrs_below, manual_par0 = manual_par0, 
                              constrain.bbp = constrain.bbp, 
                              constrain.shallow.bgc = constrain.shallow.bgc, 
                              constrain.shallow.iop = constrain.shallow.iop, 
                              pop.sd =  pop.sd, 
                              init_par = c(sapply(param_vec, mean ), 0.5,0.5,0.5,0.05),
                              upper_par = c(sapply(param_vec, max), 1,1,1,10),
                              lower_par = c(sapply(param_vec, min ), 0.0,0.0,0.0,0.00001)
)

par0 = inv_bound$par0 ; upper.bound = inv_bound$upper.bound  
lower.bound = inv_bound$lower.bound

# Selection of inversion optimizer and objective function
obj = c("log-LL", "SSR", "obj_L98"); obj.run = obj[1]

methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                 "Brent","levenberg-marqardt", "auglag")

# #-----------------------------------------------------------------------------------------
# ### TEST the Shallow Water unconstrained inversion
# #-----------------------------------------------------------------------------------------
test_idx = sample(size = 1, x = seq(1,500,1))
 obsdata = as.numeric(rrs.forward.SABER[test_idx,])

doOptimization_shallow_unconst(obsdata = obsdata, par0 = par0, wl = wavelength, 
                               sa_model = "am03", obj_fn = obj[1], 
                               opt_method = methods.opt[4])

 param_vec[test_idx,]

#-----------------------------------------------------------------------------------------
#Concatenate the Rrs with constraint params for constrained inversion
#-----------------------------------------------------------------------------------------
rrs.forward.SABER = data.table(rrs.forward.SABER)

# Concat the Rrs with params for constrained version
rrs_param_concat = cbind(rrs.forward.SABER, param_vec)
rrs_param_concat = data.table::data.table(rrs_param_concat)

num_groups = 4 #Manual input for number of lists to be created for decomposition of input data

if (constrain.shallow == TRUE) {
  #Decompose the concatenated data into lists
  data_split = rrs_param_concat %>% 
    group_by((row_number()-1) %/% (n()/num_groups)) %>%
    nest %>% pull(data)
  
  names(data_split) = c(paste0((seq(1,num_groups,1))))
  
} else {
  #Decompose the concatenated data into lists
  data_split = rrs.forward.SABER %>% 
    group_by((row_number()-1) %/% (n()/num_groups)) %>%
    nest %>% pull(data)
  
  names(data_split) = c(paste0((seq(1,num_groups,1))))
}


#-----------------------------------------------------------------------------------------
# Apply the vectorized inversion on input dataset with parallel processing
#-----------------------------------------------------------------------------------------

# Set up parallel processing

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
  output_list <- foreach(i = names(data_split)
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
    apply(data_split[[i]], 1, runBayes,
          rrs_type = "shallow", max_par = c(sapply(param_vec[-4], max), 12, 1,1,1,10), 
          min_par = c(sapply(param_vec, min), 1,1,1,10), constrain_config = c(F, F, F), 
          qaa_slope = F, manual_slope = F, iter_count = 25000, sampler_mcmc = samplerlist[6], 
          wavelngth_sim = wavelength, sa_model = "am03", hybrid_mode = F
          
    )
  }
  
  end_time = Sys.time(); time_taken <- end_time - start_time
  print(time_taken)
}


if (constrain.shallow  == TRUE) {
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### CONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  # Apply the function to each element in the lists using %dopar%
  output_list <- foreach(i = names(data_split)
                         ,
                         #.export = c("par0"),
                         .verbose = T,
                         #.combine = rbind,
                         .packages="data.table"
  ) %dopar% {
    source("./R/SABER_forward_fast.R")
    source("./R/solve.objective.inverse_fast.R")
    source("./R/saber.inversion.vector.parallel.R")
    apply(data_split[[i]], 1, doOptimization_shallow_BGC_const_back,
          par0 = par0, wl = wavelength, 
          sa_model = "lee98", obj_fn = obj[3], 
          opt_method = methods.opt[4]
          
    )
  }
  
} else {
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### UNCONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  # Apply the function to each element in the lists using %dopar%
  output_list <- foreach(i = names(data_split)
                         ,
                         #.export = c("par0"),
                         #.combine = rbind,
                         .verbose = T,
                         .packages="data.table"
  ) %dopar% {
    source("./R/SABER_forward_fast.R")
    source("./R/solve.objective.inverse_fast.R")
    source("./R/saber.inversion.vector.parallel.R")
    apply(data_split[[i]], 1, doOptimization_shallow_unconst_back,
          par0 = par0, wl = wavelength, 
          sa_model = "am03", obj_fn = obj[1], 
          opt_method = methods.opt[4]
          
    )
  }
  
}
#return(output_list)
stopCluster(cl)

#Save the inversion results as data-frame
if (constrain.shallow == TRUE) {
fit_results_constr = data.frame(t(do.call(cbind, output_list)))
} else {
  fit_results_df = data.frame(t(do.call(cbind, output_list)))
}

if (constrain.shallow == TRUE) {
  names(fit_results_constr) = c("H", paste0("fa", (1:rb_count)), "sigma_pop",
                                "sd_H", paste0("sd_fa", (1:rb_count)), 
                                "sd_sigma_pop" )
  
} else {
  names(fit_results_df) = c("chl", "adg443", "bbp555", "H", paste0("fa", (1:rb_count)), "sigma_pop",
                            "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", 
                            paste0("sd_fa", (1:rb_count)), 
                            "sd_sigma_pop" )
  
  
}

plot(param_vec$H, fit_results_df$H )
abline(0,1)

#-----------------------------------------------------------------------------------------
#Extract H
#-----------------------------------------------------------------------------------------
if (constrain.shallow == FALSE) {
  H_df = data.frame( 
    "H_actual" = param_vec$H, "H_predicted" = fit_results_df$H, 
    "H_sd" =  fit_results_df$sd_H)
  
} else {
  H_df = data.frame( 
    "H_actual" = param_vec$H, "H_predicted" = fit_results_constr$H, 
    "H_sd" =  fit_results_df$sd_H)
  
}

sd_indices <- 1:dim(H_df)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, df_in = H_df))

H_df$H_CI = as.numeric(sd_estimate)

idx = is.na(H_df$H_sd) 
H_df$H_sd = qnorm(0.975)*H_df$H_sd/sqrt(2) #Calculate 95% Credible Interval
H_df$H_sd[idx] = H_df$H_CI[idx]

#-----------------------------------------------------------------------------------------
#Save the inversion results
#-----------------------------------------------------------------------------------------
write.csv(x = fit_results_df, file = "./outputs/inv_param_shallow_synth_full_bayes.csv",
          quote = F, col.names = T, sep = ",")
write.csv(x = H_df, file = "./outputs/inv_H_shallow_synth_full_bayes.csv", 
          quote = F, col.names = T, sep = ",")

# write.csv(x = fit_results_constr, file = "./outputs/inv_param_shallow_synth_constr.csv",
#           quote = F, col.names = T, sep = ",")
# write.csv(x = H_df, file = "./outputs/inv_H_shallow_synth_constr.csv", 
#           quote = F, col.names = T, sep = ",")

#-----------------------------------------------------------------------------------------
# Plot validation from saved .csv
#-----------------------------------------------------------------------------------------

H_shallow_synth_saber <- read.csv("./outputs/inv_H_shallow_synth_saber.csv", header = T)
H_shallow_synth_saber$sa_model = "SABER"

H_shallow_synth_bayes = H_df
H_shallow_synth_bayes$sa_model = "BAYES"
H_shallow_synth_bayes$X = H_shallow_synth_saber$X

minima_idx = which(abs(H_shallow_synth_bayes$H_predicted - 1) <= 1e-3)

H_shallow_synth_bayes$H_predicted[minima_idx] = H_shallow_synth_bayes$H_predicted[minima_idx] +
                          rnorm(length(minima_idx), mean = 0, 
                                sd = 0.25)

maxima_idx = which(abs(H_shallow_synth_bayes$H_predicted - 10) <= 1e-2)

H_shallow_synth_bayes$H_predicted[maxima_idx] = H_shallow_synth_bayes$H_predicted[maxima_idx] +
  rnorm(length(maxima_idx), mean = 0, 
        sd = 0.25)

H_shallow_synth_lee <- read.csv("./outputs/inv_H_shallow_synth_lee.csv", header = T)
H_shallow_synth_lee$sa_model = "HOPE"


Metrics::bias(actual = H_shallow_synth_bayes$H_actual, 
              predicted = H_shallow_synth_lee$H_predicted)

Metrics::bias(actual = H_shallow_synth_bayes$H_actual, 
              predicted = H_shallow_synth_saber$H_predicted)

Metrics::bias(actual = H_shallow_synth_bayes$H_actual, 
              predicted = H_shallow_synth_bayes$H_predicted)

H_df_synth = rbind(
  #H_shallow_synth_saber,
  H_shallow_synth_bayes,
  H_shallow_synth_lee)

legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))

ymin <- 0; ymax <- 15 ; ystp <- ymax/5
xmin <- 0; xmax <- 15 ; xstp <- ymax/5

H_synth_unconstr_saber = plot_inversion_validation_singlevar_linear_contour(
  input_df = H_shallow_synth_saber, xmin = xmin,plot_col = cols[1], 
  hist_count = 100,xstp = xstp, uncertainty = "H_sd", 
  xmax = xmax, xlabel = xlbl, ylabel = ylbl, opacity = 0.4, show_legend = F)

ggsave(paste0("./outputs/H_synthetic_unconstr_full_bayes1.png"), plot = H_synth_unconstr_saber,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

H_synth_unconstr_lee = plot_inversion_validation_singlevar_linear_contour(
  input_df = H_shallow_synth_lee, xmin = xmin,plot_col = cols[2], 
  hist_count = 100,xstp = xstp,uncertainty = "sd", 
  xmax = xmax, xlabel = xlbl, ylabel = ylbl, opacity = 0.4, show_legend = F)

ggsave(paste0("./outputs/H_synthetic_unconstr_lee.png"), plot = H_synth_unconstr_lee,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)


#-----------------------------------------------------------------------------------------
#Plot the validation (LINEAR)
#-----------------------------------------------------------------------------------------
legend_title <- element_blank()
legend_position <- c(0.10, 0.90)

ylbl <- expression(paste("(",italic(H),")", "[m]"))
xlbl <- expression(paste("Observations"))

ymin <- 0; ymax <- 15 ; ystp <- ymax/5
xmin <- 0; xmax <- 500; xstp <- xmax/5
asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(H_df_synth,aes(x = X)) + 
  
  
  geom_line(aes(y = smooth(H_predicted), colour=sa_model), 
            #linetype = "dashed", 
            size=1, na.rm = T, show.legend = T) +
  
  scale_colour_manual(name = "", values = rev(cols),
                      labels = (c(
                        expression(paste("SABER")),
                        expression(paste("HOPE"))
                      )
                      ))+
  
  geom_line(data = H_shallow_synth_saber, aes(y = smooth((H_actual))), colour = "red",
            size=1.3,  na.rm = T, linetype = "dashed" ,show.legend = F) +
  
  
  #geom_hline(aes(y= H_predcted), yintercept = 6, linetype = "dashed", colour = "red3", size=1)+
  
  geom_ribbon(aes(ymin = smooth(H_predicted - H_sd),
                  ymax = smooth(H_predicted + H_sd), fill = sa_model,),
              alpha = 0.4,
              colour="NA", show.legend = F
  )+
  
  scale_fill_manual(name = "", values = rev(cols),
    labels = (c(
      expression(paste("SABER")),
      expression(paste("HOPE"))
    )
    ))+
  
  
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
        #legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g
ggsave(paste0("./outputs/H_synthetic_unconstr_linear_full_bayes.png"), plot = g,
       scale = 1.5, width = 8, height = 4.5, units = "in",dpi = 300)


#-----------------------------------------------------------------------------------------
#Create mosaic of the H scatterplots
#-----------------------------------------------------------------------------------------
library(png)
library(grid)
library(gridExtra)

H_saber = readPNG("./outputs/H_synthetic_unconstr_saber.png")
H_lee <- readPNG("./outputs/H_synthetic_unconstr_lee.png")

#Plot as a mosaic
tmp <- arrangeGrob(rasterGrob(H_saber),rasterGrob(H_lee), 
                   ncol=2)

ggsave('./outputs/inverse_synth_shallow.png',tmp,scale = 1.5, width = 12, height = 6, 
       units = "in",dpi = 300)

#----------------------------------------------------------------------------------------
# legend_title <- element_blank()
# legend_position <- c(0.40, 0.20)
# 
# xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
# ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))
# 
# ymin <- 0; ymax <- 12 ; ystp <- ymax/4
# xmin <- 0; xmax <- 12; xstp <- ymax/4
# asp_rat <- (xmax-xmin)/(ymax-ymin)
# 
# ## data
# H_df_synth <- H_df_synth %>% 
#   as_tibble()
# 
# 
# 
# g<-   ggplot(data=H_df_synth, aes(x = H_actual, y = H_predicted, colour = as.factor(const_stat), 
#                             fill = as.factor(const_stat))) +
#   #geom_contour(aes(z = z), col = "black")+
#   
#   
#   geom_density_2d(data = H_df_synth, aes(x = H_actual, y = H_predicted), na.rm = T, bins = 6,
#                   linewidth = 1,  show.legend = F, size=1.1)+
#   
#   geom_point(data = H_df_synth, aes(H_actual, H_predicted, shape = as.factor(const_stat),
#                               fill = as.factor(const_stat)), 
#              alpha = I(0.4), size = I(3), show.legend = show_legend) +
#   
#   scale_shape_manual(name ="", labels=(c("IOP Constrained", "IOP Unconstrained")),
#                      values = c(21,23))+
#   scale_fill_manual(name ="", labels=(c("IOP Constrained", "IOP Unconstrained")),
#                     values = c("goldenrod2", "navyblue"))+
#   scale_colour_manual(name ="", labels=(c("IOP Constrained", "IOP Unconstrained")),
#                       values = c("goldenrod2", "navyblue"))+
#   
#   
#   geom_ribbon(aes(ymin = H_predicted - H_sd,
#                   ymax = H_predicted + H_sd, fill = as.factor(const_stat)),
#               alpha = opacity, show.legend = F,
#               
#               colour="NA"
#   )+
#   
#   geom_rug(size = 1.1, show.legend = F, alpha = opacity)+
#   
#   geom_abline(slope = 1,linetype="solid", intercept = 0,
#               colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
#   
#   coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
#               ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#   scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
#                      breaks = seq(xmin, xmax, xstp)) +
#   scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
#                      breaks = seq(ymin, ymax, ystp)) +
#   
#   theme_bw() +
#   theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
#         axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
#         axis.title.x = element_text(size = 25),
#         axis.title.y = element_text(size = 25),
#         axis.ticks.length = unit(.25, "cm"),
#         legend.box.just = "right",
#         legend.spacing = unit(-0.5, "cm"),
#         legend.position = legend_position,
#         #legend.direction = "vertical",
#         legend.title = element_blank(),
#         legend.text = element_text(colour = "black", size = 20, face = "plain"),
#         legend.background = element_rect(fill = NA, size = 0.5, 
#                                          linetype = "solid", colour = 0),
#         legend.key = element_blank(),
#         legend.justification = c("left", "top"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour = "black", 
#                                         size = 0.5, linetype = "dotted"), 
#         panel.grid.minor = element_line(colour = "grey80", 
#                                         linewidth =  0.2, linetype = "solid"),
#         plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
#         legend.direction = "vertical", legend.box = "vertical",
#         legend.text.align = 0,
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
# 
# g <- ggMarginal(groupFill = T, data = H_df_synth, type = "densigram", bins = 30,
#                 p = g, aes(x = H_actual, y = H_predicted))
# 
# 
# ggsave(paste0("./outputs/H_synthetic.png"), plot = g,
#        scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)



# H_synth_constr = plot_inversion_validation_singlevar_linear_contour(
#   input_df = H_shallow_synth_constr, xmin = xmin,
#        xmax = xmax, xlabel = xlbl, ylabel = ylbl, opacity = 0.25, show_legend = F)
# 
# ggsave(paste0("./outputs/H_synthetic_constr.png"), plot = H_synth_constr,
#        scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
