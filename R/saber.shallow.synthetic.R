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
                      "H"=(seq(from=min(0.5), to=max(8), 
                              length.out=length(vec_size))) + truncnorm::rtruncnorm(500, 
                                                                              a = 0.25, b = 2,
                                                                              mean = 2, sd = 1)
                      )

#param_vec.LUT <- expand.grid(param_vec) 

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

rrs_synth$id = seq(1, length(rrs_synth$`401`), 1)

rrs_synth = reshape2::melt(rrs_synth, id.vars = "id")

H_plot_df = data.frame("id" = seq(1,length(param_vec$H), 1), "Depth" = param_vec$H)
H_long = reshape2::melt(H_plot_df, id.vars = "id")

rrs_synth_H = merge(rrs_synth,H_long,by=c("id"))

names(rrs_synth_H) = c("id", "wavelength", "rrs", "shallow_tag", "depth")

xmin <- 400; xmax <- 700;  xstp <- 50; xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0; ymax <- 0.015; ystp <- 0.003; ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda),italic("in situ ("), "sr"^{-1}, ")"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.70, 0.98)


library(viridis)

g <- ggplot(data = rrs_synth_H) + 
  geom_line(aes(x = as.numeric(as.character(wavelength)), y = rrs, colour = depth, 
                group = id), size = 1.3, show.legend = T)+
  scale_colour_viridis(discrete = F, name = "Depth (m)") +
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax)
              ,expand = FALSE, clip = "on"
  ) +
  #scale_colour_manual(name = "",values = rev(collist))+
  #guides(fill = FALSE)+
  guides(fill = guide_legend(title = "Depth (m)"))+
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

# BGC CONSTRAINED
if (pop.sd == "unknown" & type_Rrs_below == "shallow" & constrain.shallow.bgc == "TRUE") { 
  par0 = c(z = 2.5, 
           rb_c = rep(0.5, rb_count),
           population.sd = 0.05)
  
  lower.bound <- c(z = 1,
                   rep(0, rb_count),
                   population.sd = 0.00001)
  
  upper.bound <- c(z = 10,
                   rep(1, rb_count),
                   population.sd = 10)
}

#IOP CONSTRAINED
if (pop.sd == "unknown" & type_Rrs_below == "shallow" & constrain.shallow.iop == "TRUE") { 
  par0 = c(z = 2.5, 
           rb_c = rep(0.5, rb_count),
           population.sd = 0.05)
  
  lower.bound <- c(z = 1,
                   rep(0, rb_count),
                   population.sd = 0.00001)
  
  upper.bound <- c(z = 10,
                   rep(1, rb_count),
                   population.sd = 10)
}


# FREE PARAMTER INVERSION 
if (pop.sd == "unknown" & type_Rrs_below == "shallow" & 
    constrain.shallow.bgc == "FALSE" & constrain.shallow.iop == "FALSE") { # <<MANUAL INPUT>>
  par0 = c(apply(param_vec, 2, mean), 
           rb_c = rep(0.5, rb_count),
           population.sd = 0.05)
  
  lower.bound <- c(apply(param_vec, 2, min),
                   rep(0, rb_count),
                   population.sd = 0.00001)
  
  upper.bound <- c(apply(param_vec, 2, max),
                   rep(1, rb_count),
                   population.sd = 10)
}

# Selection of inversion optimizer and objective function
obj = c("log-LL", "SSR", "obj_L98"); obj.run = obj[1]

methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                 "Brent","levenberg-marqardt", "auglag")

# #-----------------------------------------------------------------------------------------
# ### TEST the Shallow Water unconstrained inversion
# #-----------------------------------------------------------------------------------------
# obsdata = as.numeric(rrs.forward.SABER[1,])
# 
# if (type_Rrs_below == "shallow" & constrain.shallow.bgc == FALSE) {
#   inverse_output_shallow_unconstr <- suppressWarnings(solve.objective.inverse.shallow.final.fast(
#     
#     constrain.shallow.iop = F, 
#     
#     unconstrained = T,
#     
#     constrained.bgc = F, 
#     constrain.bgc.value = c(Fit.input$chl, Fit.input$acdom.440+ 
#                               Fit.input$anap.440),
#     
#     
#     initial = as.numeric(par0), 
#     obsdata = obsdata,
#     
#     auto_spectral_slope = F,
#     manual_spectral_slope = F, 
#     
#     manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.003, "gamma"=1),
#     
#     sa.model = "am03", 
#     #obj.fn =obj.run , 
#     obj.fn = obj[1],
#     method.opt = methods.opt[4],
#     
#     lower.b = lower.bound,
#     upper.b = upper.bound, 
#     
#   ))
# }
# param_vec[1,]

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

start_time = Sys.time()
if (constrain.shallow  == TRUE) {
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### CONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  # Apply the function to each element in the lists using %dopar%
  output_list <- foreach(i = names(inv_list)
                         ,
                         #.export = c("par0"),
                         .verbose = T,
                         #.combine = rbind,
                         .packages="data.table"
  ) %dopar% {
    source("./R/SABER_forward_fast.R")
    source("./R/solve.objective.inverse_fast.R")
    source("./R/saber.inversion.vector.parallel.R")
    apply(inv_list[[i]], 1, doOptimization_shallow_BGC_const,
          par0 = par0
          
    )
  }
  
} else {
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### UNCONSTRAINED SHALLOW WATER INVERSION #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  # Apply the function to each element in the lists using %dopar%
  output_list <- foreach(i = names(inv_list)
                         ,
                         #.export = c("par0"),
                         #.combine = rbind,
                         .verbose = T,
                         .packages="data.table"
  ) %dopar% {
    source("./R/SABER_forward_fast.R")
    source("./R/solve.objective.inverse_fast.R")
    source("./R/saber.inversion.vector.parallel.R")
    apply(inv_list[[i]], 1, doOptimization_shallow_unconst,
          par0 = par0
          
    )
  }
  
}
end.time = Sys.time(); time_taken <- end.time - start.time
cat(paste0("\033[0;33m","Time Taken for the inversion: ",time_taken," secs. \033[0m","\n"))
return(output_list)
stopCluster(cl)

#Save the inversion results as data-frame
fit_results_constr = data.frame(t(do.call(cbind, output_list)))

if (constrain.shallow == TRUE) {
  names(fit_results_constr) = c("H", paste0("fa", (1:rb_count)), "sigma_pop",
                                "sd_H", paste0("sd_fa", (1:rb_count)), 
                                "sd_sigma_pop" )
  
} else {
  names(fit_results_df) = c("chl", "adg443", "bbp555", "H", paste0("fa", (1:rb_count)), "sigma_pop",
                            "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", paste0("sd_fa", (1:rb_count)), 
                            "sd_sigma_pop" )
  
  
}

# startloc <- 1
# indices <- startloc:dim(rrs.forward.SABER)[1]
# 
# # Apply the optimization function and create a dataframe of obtained params
# fit_results <- t(sapply(indices,  doOptimization))  
# fit_results_df <- as.data.frame(fit_results)

plot(param_vec$H, fit_results_df$H )
abline(0,1)
#-----------------------------------------------------------------------------------------
#Extract H
#-----------------------------------------------------------------------------------------
H_df = data.frame( 
                  "H_actual" = param_vec$H, "H_predicted" = fit_results_df$H, 
                  "H_sd" =  fit_results_df$sd_H)

H_df = data.frame( 
  "H_actual" = param_vec$H, "H_predicted" = fit_results_constr$H, "H_sd" =  fit_results_df$sd_H)

#Calculate Confidence Interval
cal_sd <- function(i){
  errorb <- qnorm(0.975)*sd(H_df[i,-3])/
    sqrt(length(H_df[i,-3]))
}

sd_indices <- 1:dim(H_df)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd)) 

H_df$H_CI = as.numeric(sd_estimate)

idx = is.na(H_df$H_sd) 
H_df$H_sd = qnorm(0.975)*H_df$H_sd/sqrt(2) #Calculate 95% Credible Interval
H_df$H_sd[idx] = H_df$H_CI[idx]

#-----------------------------------------------------------------------------------------
#Plot the validation
#-----------------------------------------------------------------------------------------
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))

ymin <- 0; ymax <- 12 ; ystp <- ymax/4
xmin <- 0; xmax <- 12 ; xstp <- ymax/4
asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(H_df,aes(x = H_actual)) + 
  
  stat_density_2d(aes(y = H_predicted, fill = ..level..),
                  geom = "polygon", colour="gray", show.legend = F, alpha = 0.5, size=1.1)+

  scale_fill_distiller(palette = 2, direction = 1)+
  
  geom_point(aes(y = H_predicted), shape=21, fill="goldenrod2", 
             size=3.0, na.rm = T, show.legend = F) +
  
  # geom_errorbar( aes(ymin = H_predicted - (2.5*H_sd),
  #                    ymax = H_predicted + (2.5*H_sd)),
  #                color="black",
  #                
  #                width = 0.1, 
  #                size=1.3
  #                
  # )+
  geom_ribbon(aes(ymin = H_predicted - H_sd,
                  ymax = H_predicted + H_sd),
              alpha = 0.5,
              fill = "navyblue",
              colour="NA"
  )+
  geom_abline(slope = 1,linetype="solid", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  
  # geom_smooth(size=1.5,level = 0.95,show.legend = F,linetype = "dashed",
  #             color="black", data = H_df,
  #             se= T, method = "lm", aes(x=H_actual, y=H_predicted))+
  
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
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
ggsave(paste0("./outputs/H_synthetic_constr.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#-----------------------------------------------------------------------------------------
#Save the inversion results
#-----------------------------------------------------------------------------------------
write.csv(x = fit_results_df, file = "./outputs/inv_param_shallow_synth.csv",
          quote = F, col.names = T, sep = ",")
write.csv(x = H_df, file = "./outputs/inv_H_shallow_synth.csv", 
          quote = F, col.names = T, sep = ",")

write.csv(x = fit_results_constr, file = "./outputs/inv_param_shallow_synth_constr.csv",
          quote = F, col.names = T, sep = ",")
write.csv(x = H_df, file = "./outputs/inv_H_shallow_synth_constr.csv", 
          quote = F, col.names = T, sep = ",")

#-----------------------------------------------------------------------------------------
# Plot validation from saved .csv
#-----------------------------------------------------------------------------------------

H_shallow_synth_unconstr <- read.csv("./outputs/inv_H_shallow_synth.csv", header = T)
H_shallow_synth_unconstr$const_stat = "unconstrained"

H_shallow_synth_constr <- read.csv("./outputs/inv_H_shallow_synth_constr.csv", header = T)
H_shallow_synth_constr$const_stat = "constrained"

H_df_synth = rbind(H_shallow_synth_constr, H_shallow_synth_unconstr)


legend_title <- element_blank()
legend_position <- c(0.40, 0.20)

xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))

ymin <- 0; ymax <- 12 ; ystp <- ymax/4
xmin <- 0; xmax <- 12; xstp <- ymax/4
asp_rat <- (xmax-xmin)/(ymax-ymin)

## data
H_df_synth <- H_df_synth %>% 
  as_tibble()



g<-   ggplot(data=H_df_synth, aes(x = H_actual, y = H_predicted, colour = as.factor(const_stat), 
                            fill = as.factor(const_stat))) +
  #geom_contour(aes(z = z), col = "black")+
  
  
  geom_density_2d(data = H_df_synth, aes(x = H_actual, y = H_predicted), na.rm = T, bins = 6,
                  linewidth = 1,  show.legend = F, size=1.1)+
  
  geom_point(data = H_df_synth, aes(H_actual, H_predicted, shape = as.factor(const_stat),
                              fill = as.factor(const_stat)), 
             alpha = I(0.4), size = I(3), show.legend = show_legend) +
  
  scale_shape_manual(name ="", labels=(c("IOP Constrained", "IOP Unconstrained")),
                     values = c(21,23))+
  scale_fill_manual(name ="", labels=(c("IOP Constrained", "IOP Unconstrained")),
                    values = c("goldenrod2", "navyblue"))+
  scale_colour_manual(name ="", labels=(c("IOP Constrained", "IOP Unconstrained")),
                      values = c("goldenrod2", "navyblue"))+
  
  
  geom_ribbon(aes(ymin = H_predicted - H_sd,
                  ymax = H_predicted + H_sd, fill = as.factor(const_stat)),
              alpha = opacity, show.legend = F,
              
              colour="NA"
  )+
  
  geom_rug(size = 1.1, show.legend = F, alpha = opacity)+
  
  geom_abline(slope = 1,linetype="solid", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  
  theme_bw() +
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
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

g <- ggMarginal(groupFill = T, data = H_df_synth, type = "densigram", bins = 30,
                p = g, aes(x = H_actual, y = H_predicted))


ggsave(paste0("./outputs/H_synthetic.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)



H_synth_constr = plot_inversion_validation_singlevar_linear_contour(
  input_df = H_shallow_synth_constr, xmin = xmin,
       xmax = xmax, xlabel = xlbl, ylabel = ylbl, opacity = 0.25, show_legend = F)

ggsave(paste0("./outputs/H_synthetic_constr.png"), plot = H_synth_constr,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

H_synth_unconstr = plot_inversion_validation_singlevar_linear_contour(
  input_df = H_shallow_synth_unconstr, xmin = xmin,
  xmax = xmax, xlabel = xlbl, ylabel = ylbl, opacity = 0.25, show_legend = F)

ggsave(paste0("./outputs/H_synthetic_unconstr.png"), plot = H_synth_unconstr,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#-----------------------------------------------------------------------------------------
#Create mosaic of the H scatterplots
#-----------------------------------------------------------------------------------------
library(png)
library(grid)
library(gridExtra)

H_unconst = readPNG("./outputs/H_synthetic_unconstr.png")
H_const <- readPNG("./outputs/H_synthetic_constr.png")

#Plot as a mosaic
tmp <- arrangeGrob(rasterGrob(H_const),rasterGrob(H_unconst), 
                   ncol=2)

ggsave('./outputs/inverse_shallow.png',tmp,scale = 1.5, width = 12, height = 6, 
       units = "in",dpi = 300)
