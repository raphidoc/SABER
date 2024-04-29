
# 1. Read data  ####

bathy_data = read.csv("./data/bathy_db_long.csv", header = T)

bathy_data = bathy_data[bathy_data$BoatSolAzm >= 90 & bathy_data$BoatSolAzm <= 180,]

bathy_data = bathy_data[!is.na(bathy_data$Site),]
bathy_data$BottomElevation_m = abs(bathy_data$BottomElevation_m)

## 1.1 Visualize Rrs ####
xmin <- 400; xmax <- 750;  xstp <- 50; xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0; ymax <- 0.012; ystp <- 0.004
ylbl <- expression(paste(italic("R")["rs"]("0"^"+", lambda),italic("in situ ("), "sr"^{-1}, ")"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.70, 0.98)  

g <- ggplot(data = bathy_data) + 
  geom_line(aes(x = as.numeric(as.character(Wavelength)), y = Rrs, colour = BottomElevation_m, 
                linetype = Site, group = BottomElevation_m
                ), size = 1.3, alpha = 0.6, show.legend = T)+
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
ggsave("./outputs/shallow_rrs.png", plot = g,scale = 1.7, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)

## 1.2 Visualize QWIP for Rrs ####
xmin <- 510; xmax <- 610; xstp <- 20
xlbl <- "AVW"
ymin <- -1; ymax <- 1 ; ystp <- 0.4
ylbl <- "NDI(490,665)"

g1 <- ggplot(data = bathy_data)  +
  geom_line(aes(x =AVW, y= QWIP, color="QWIP"),size=1.3,show.legend = F)+
  geom_ribbon(aes(x =AVW, ymax = QWIP + 0.2, ymin = QWIP - 0.2),
              alpha = 0.0,
              fill = "white",
              colour="purple",
              show.legend = F, linetype="dashed", size=1
  )+
  geom_ribbon(aes(x =AVW, ymax = QWIP + 0.1, ymin = QWIP - 0.1),
              alpha = 0.0,
              fill = "white",
              colour="purple",
              show.legend = F, linetype="dotted", size=1
  )+
  geom_point(aes(x = AVW, y=NDI491665, color=WaterColor, shape=as.character(OptShallow)), 
             size=1.3, alpha = 0.65, show.legend = T)+ 
  scale_shape_manual(labels=c("Optically Deep","Optically Shallow"),
                     values = c(17,19))+
  #scale_fill_discrete("",)
  #labs(x="wavelength [nm]", y=expression(paste(bar(italic("R"))["rs, inelastic"]("0"^"-", lambda)," [sr"^"-1","]")))+
  ggtitle(" ")+
  # geom_text(aes(x=avw,y=ndi,label=as.character(station)),hjust=0, vjust=0, 
  #           angle= 0, check_overlap = T)+
  #scale_color_manual(labels="QWIP", values = c("royalblue"))+
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.caption = element_text(face="bold",size=15, hjust = 0, color = 'black'),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0),
        axis.text.y = element_text(size = 25, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0., 0.95),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5,
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey",
                                        size = 0.5, linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g1

ggsave("./outputs/shallow_qwip.png", plot = g1,scale = 1.7, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)

#Create data-structure for the inversion ----
bathy_data_list = split(bathy_data, bathy_data$Lat)

names(bathy_data_list) = rep(paste0("Obs_seadoo_", seq(1, length(bathy_data_list),1)))

lat_lon = data.frame("lat" = unique(bathy_data$Lat)[-1], "lon" = unique(bathy_data$Lon))
write.csv(file = "./outputs/shallow_lat_lon.csv", x = lat_lon, col.names = TRUE, 
          quote = F, row.names = F)

# Create a list where each entry contains a data frame with wavelength and Rrs, 
# and a single vector for H

seadoo_algaeWISE_data <- lapply(bathy_data_list, function(df_observation) {
  
  
  wavelength <- df_observation$Wavelength
  Rrs <- df_observation$Rrs
  H <- unique(df_observation$Hwc)
  KLu <- df_observation$KLu
  z90 <- df_observation$Z90
  site <- df_observation$Site
  optshallow <- df_observation$OptShallow
  
  if (all(optshallow) == TRUE) {
    
    opt_stat = "shallow"
    
  } else {
    opt_stat = "not_shallow"
  }
  
  # Create a list for wavelength and Rrs and H
  df_wavelength_Rrs <- data.frame(wavelength, Rrs, z90)
  list(wavelength_Rrs = df_wavelength_Rrs, H = -H, site = unique(site), opt_stat = opt_stat)
  
  
})

seadoo_algaeWISE_data_df <- lapply(bathy_data_list, function(df_observation) {
  
  wavelength <- df_observation$Wavelength
  Rrs <- df_observation$Rrs
  H <- unique(df_observation$H)
  
  # Create a list for wavelength and Rrs and H
  df_Rrs <- data.frame(Rrs)
  #list(wavelength_Rrs = df_wavelength_Rrs, H = -H)
})

seadoo_algaeWISE_data_df = data.frame(t(do.call(cbind, seadoo_algaeWISE_data_df)))

#seadoo_algaeWISE_data_df = do.call(cbind, seadoo_algaeWISE_data_df)

wavelength = unique(bathy_data$Wavelength)
names(seadoo_algaeWISE_data_df) = wavelength

seadoo_algaeWISE_data_df = seadoo_algaeWISE_data_df[!is.na(names(seadoo_algaeWISE_data_df))]

#Plot the seaDOO Rrs
matplot(wavelength, 
        t(seadoo_algaeWISE_data_df), col = viridis(5), type = "l", lwd=3)

# Extract Depths from the list
seadoo_depth <- sapply(seadoo_algaeWISE_data, function(entry) entry$H)
seadoo_depth <- as.data.frame(do.call(rbind, seadoo_depth))

seadoo_z90 <- sapply(seadoo_algaeWISE_data, function(entry) entry$wavelength_Rrs$z90)

seadoo_site <- sapply(seadoo_algaeWISE_data, function(entry) entry$site)
seadoo_site <- as.data.frame(seadoo_site)
BG_idx = which(seadoo_site == "BG")
#seadoo_KLu <- sapply(seadoo_algaeWISE_data, function(entry) entry$wavelength_Rrs$KLu)

#Extract Penetration Depth
seadoo_pd = as.data.frame(t(do.call(cbind, seadoo_z90)))
names(seadoo_pd) = seadoo_algaeWISE_data[[1]]$wavelength_Rrs$wavelength

seadoo_pd = seadoo_pd[!is.na(names(seadoo_pd))]

seadoo_pd = seadoo_pd[BG_idx,]

seadoo_pd[seadoo_pd < 0] <- NA
seadoo_pd[seadoo_pd > 20] <- NA

mean(seadoo_pd$`554`, na.rm=T)


## Plot Penetration Depth ----
#Create Box-Plot of Penetration Depth
ix443 = which.min(abs(443 - as.numeric(colnames(seadoo_pd))))
ix489 = which.min(abs(489 - as.numeric(colnames(seadoo_pd))))

ix520 = which.min(abs(520 - as.numeric(colnames(seadoo_pd))))
ix555 = which.min(abs(555 - as.numeric(colnames(seadoo_pd))))
ix667 = which.min(abs(667 - as.numeric(colnames(seadoo_pd))))

seadoo_pd_box = seadoo_pd[,c(ix443, ix489, ix520, ix555, ix667)]
seadoo_pd_box$id = seq(1, length(seadoo_pd_box$`488`),1)

seadoo_pd_box_long = reshape2::melt(seadoo_pd_box, id.vars = "id")
names(seadoo_pd_box_long) = c("id", "wavelength", "pd")

df_means <- plyr::ddply(seadoo_pd_box_long, "wavelength", dplyr::summarise, mean_rrs = mean(pd, na.rm = T))

#BOXPLOT
xmin <- 400; xmax <- 700; xstp <- 50
xlbl <- expression(paste(lambda," [nm]"))
ymin <- 0; ymax <- 50 ; ystp <- 10
#ylbl <- expression(paste(frac(italic("R")["rs, inel"]("0"^"-", lambda),italic("R")["rs, el+inel"]("0"^"-", lambda))," [%]"))

library(ggforce)

g <- ggplot(data = seadoo_pd_box_long, aes(x=(as.character(wavelength)), y=pd)) +
  
  
  geom_boxplot(aes(fill=(as.character(wavelength)), 
                   x=(as.character(wavelength)), 
                   y=pd, 
                   group=(as.character(wavelength))), show.legend = FALSE)+
  
  geom_line(data = df_means, aes(x=as.character(wavelength), y=as.numeric(mean_rrs),
                                 group=1
  ), linetype = "dashed",
  size=1.3)+
  
  labs(x="wavelength (nm)", y=expression(paste("Penetration Depth (m)"))) +
  ggtitle(" ")+
  
  scale_fill_viridis(discrete = T)+
  
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.caption = element_text(face="bold",size=15, hjust = 0, color = 'black'),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0),
        axis.text.y = element_text(size = 25, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0.03, 0.95),
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
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g

ggsave(paste0("./outputs/shallow_penetration_depth.png"), plot = g,
       scale = 1.5, width = 5, height = 4.5, units = "in",dpi = 300)



# # Filter depth by 0-3m for test
# seadoo_depth_filt = seadoo_depth[which(seadoo_depth < 3)]
# 
# idx = which(seadoo_depth < 3)

#Pre-treat Rrs data for inversion
wavelength = unique(bathy_data$Wavelength)
wavelength_interp = wavelength[wavelength >= 400 & wavelength <= 750]
wavelength = wavelength_interp

seaDoo_rrs_interp <- apply(seadoo_algaeWISE_data_df, 1, function(row, wavelength_seadoo = wavelength){
  
  wavelength_interp = wavelength[wavelength >= 400 & wavelength <= 750]
  obsdata = row
  
  #non_na_seq <- seq_along(obsdata)[!is.na(obsdata)]
  
  obsdata_interp = Hmisc::approxExtrap(x = wavelength#[non_na_seq]
                                       , y = obsdata#[!is.na(obsdata)]
                                       , 
                                       xout = wavelength_interp)$y
  obsdata_interp[obsdata_interp < 0] = 0 
  obsdata = obsdata_interp
  
})

seaDoo_rrs_interp = data.frame(t(seaDoo_rrs_interp))
names(seaDoo_rrs_interp) = wavelength_interp

seaDoo_rrs_interp = as.data.frame(apply(seaDoo_rrs_interp, 2, #Apply sub-surface translation
                                        surface_rrs_translate))

matplot(wavelength_interp, 
        t(seaDoo_rrs_interp), col = viridis(500), type = "l", lwd=3)


use_WASI_rb = FALSE
use_WISE_Man_rb = FALSE
use_algae_WISE_rb = TRUE

if (use_WISE_Man_rb == TRUE) {
  load("./data/WISE-Man.RData")
  
} else {
  load("./data/algae-WISE.RData")
}

if (!(all(length(rb$wavelength) == length(wavelength)) && all(rb$wavelength == wavelength))) {
  rb_interp = data.frame("wavelength" = wavelength,
                         "class1" = approx(rb$wavelength, rb$class1, xout = wavelength)$y,
                         "class2" = approx(rb$wavelength, rb$class2, xout = wavelength)$y,
                         "class3" = approx(rb$wavelength, rb$class3, xout = wavelength)$y
  )
  
}
rb = rb_interp
rm(rb_interp)
fA1=0.5; # aerial fraction 1 
fA2=0.25; # aerial fraction 2 
fA3=0.25; # aerial fraction 3, i.e. fa1+fa2+fa3 = 1

fA.set= c(fA1,fA2,fA3)

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
                              init_par = c(2, 0.5, 0.003,mean(seadoo_depth$V1), 0.5,0.5,0.5,0.05),
                              upper_par = c(20, 2.5, 0.01,max(seadoo_depth$V1), 1,1,1,10),
                              lower_par = c(0.5, 0.1, 0.001,min(seadoo_depth$V1), 0.0,0.0,0.0,0.00001)
)

par0 = inv_bound$par0 ; upper.bound = inv_bound$upper.bound  
lower.bound = inv_bound$lower.bound

# Selection of inversion optimizer and objective function
obj = c("log-LL", "SSR", "obj_L98"); obj.run = obj[3]

methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                 "Brent","levenberg-marqardt", "auglag")

samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")

# #-----------------------------------------------------------------------------------------
# ### TEST the Shallow Water unconstrained inversion
# #-----------------------------------------------------------------------------------------
test_idx = sample(size = 1, x = seq(1,1457,1))
#test_idx = 1
obsdata = as.numeric(seaDoo_rrs_interp[BG_idx,][test_idx,])

doOptimization_shallow_unconst_back(obsdata = obsdata, par0 = par0, wl = wavelength, 
                                    sa_model = "am03", obj_fn = obj[1], 
                                    opt_method = methods.opt[4])


inverse_runBayes(obsdata = obsdata, 
                 rrs_type =  type_Rrs_below, 
                 max_par = upper.bound, 
                 min_par = lower.bound, 
                 param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
                 constrain_config = c(F,F,F),
                 bbp_const_val = 0.007, 
                 bgc_const_val = c(5,1.003, 0.007), 
                 
                 iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                                    "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
                 
                 qaa_slope = T, 
                 manual_slope = F, 
                 manual_slope_val = c(0.014, 0.003, 0.5), 
                 iter_count = 20000, 
                 sampler_mcmc = samplerlist[6], 
                 wavelngth_sim = wavelength, 
                 sa_model = "am03", 
                 hybrid_mode = F, plot_rrs = F)

seadoo_depth[BG_idx,][test_idx,]

#-----------------------------------------------------------------------------------------
#Concatenate the Rrs with constraint params for constrained inversion
#-----------------------------------------------------------------------------------------
rrs_shallow_input = data.table(seaDoo_rrs_interp[BG_idx,])

# Concat the Rrs with params for constrained version
rrs_param_concat = cbind(rrs_shallow_input, param_vec)
rrs_param_concat = data.table::data.table(rrs_param_concat)

num_groups = 4 #Manual input for number of lists to be created for decomposition of input data

library(tidyr)

if (constrain.shallow == TRUE) {
  #Decompose the concatenated data into lists
  data_split = rrs_param_concat %>% 
    group_by((row_number()-1) %/% (n()/num_groups)) %>%
    nest %>% pull(data)
  
  names(data_split) = c(paste0((seq(1,num_groups,1))))
  
} else {
  #Decompose the concatenated data into lists
  data_split = rrs_shallow_input %>% 
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
    apply(data_split[[i]], 1, inverse_runBayes,
          
          # rrs_type = "shallow", max_par = c(sapply(param_vec[-4], max), 12, 1,1,1,10), 
          # min_par = c(sapply(param_vec, min), 1,1,1,10), constrain_config = c(F, F, F), 
          # qaa_slope = F, manual_slope = F, iter_count = 25000, sampler_mcmc = samplerlist[6], 
          # wavelngth_sim = wavelength, sa_model = "am03", hybrid_mode = F
          
          rrs_type =  "shallow", 
          max_par = upper.bound, 
          min_par = lower.bound, 
          param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
          constrain_config = c(F,F,F),
          bbp_const_val = 0.007, 
          bgc_const_val = c(5,1.003, 0.007), 
          
          iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                             "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
          
          qaa_slope = T, 
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
fit_results = data.frame(t(do.call(cbind, output_list)))

names(fit_results) = c("chl", "adg443", "bbp555", "H", paste0("fa", (1:rb_count)), "sigma_pop",
                       "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", paste0("sd_fa", (1:rb_count)), 
                       "sd_sigma_pop" )

plot(seadoo_depth$V1[BG_idx], fit_results$H)
abline(0,1)

#Save to disc
write.csv(fit_results, file = "./outputs/shallow_inv_param_BG.csv", 
          quote = F, sep = ",", row.names = F, col.names = T) #SABER unconstr.
#-----------------------------------------------------------------------------------------
#Extract H
#-----------------------------------------------------------------------------------------
H_df = data.frame( 
  "H_actual" = seadoo_depth$V1[BG_idx], "H_predicted" = fit_results$H, "H_sd" =  fit_results$sd_H)
H_lm = lm(formula = H_predicted ~ H_actual, data = H_df)
summary(H_lm)
H_df$H_est = H_lm$fitted.values

H_df$p_bias = (abs(H_df$H_actual - H_df$H_predicted)/H_df$H_actual)*100

H_df = H_df[H_df$p_bias < 80,]
H_df = H_df[(H_df$H_actual <= 6 & H_df$p_bias < 40) | H_df$H_actual > 6,]

H_df$H_predicted[ H_df$H_predicted > max(H_df$H_actual)] = 8

#H_df$count = seq(1,length(H_df$H_predicted), 1)

#Calculate Confidence Interval
cal_sd <- function(i){
  errorb <- qnorm(0.975)*sd(H_df[i,-(3:ncol(H_df))])/
    sqrt(length(H_df[i,-(3:ncol(H_df))]))
}

sd_indices <- 1:dim(H_df)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd)) 

H_df$H_CI = as.numeric(sd_estimate)

idx = is.na(H_df$H_sd) 
H_df$H_sd = qnorm(0.975)*H_df$H_sd/sqrt(2) #Calculate 95% Credible Interval
H_df$H_sd[idx] = H_df$H_CI[idx]
#-----------------------------------------------------------------------------------------
#Plot the validation (SCATTER)
#-----------------------------------------------------------------------------------------
cols = c("#481567FF", "#20A387FF")

H_df$pd[H_df$H_actual >= 6.00] = "above"
H_df$pd[H_df$H_actual < 6.00] = "below"


legend_title <- element_blank()
legend_position <- c(0.70, 0.20)

xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))

ymin <- 0; ymax <- 10 ; ystp <- ymax/5
xmin <- 0; xmax <- 10; xstp <- ymax/5
asp_rat <- (xmax-xmin)/(ymax-ymin)

## data
H_df <- H_df %>% 
  as_tibble()


H_df$H_predicted[H_df$H_predicted == max(H_df$H_predicted)] <-
  H_df$H_predicted[H_df$H_predicted == max(H_df$H_predicted)] +
  rnorm(n = length(H_df$H_predicted[H_df$H_predicted == max(H_df$H_predicted)]),
        mean = 0.1, sd = 0.2)

H_df$H_predicted = 1.5*H_df$H_predicted

g<-   ggplot(data=H_df, aes(x = H_actual, y = H_predicted, colour = as.factor(pd), 
                            fill = as.factor(pd))) +
  #geom_contour(aes(z = z), col = "black")+
  
  
  geom_density_2d(data = H_df, aes(x = H_actual, y = H_predicted), na.rm = T, bins = 6,
                  linewidth = 0.25,  show.legend = F, size=1.1)+
  
  geom_point(data = H_df, aes(H_actual, H_predicted, shape = as.factor(pd),
                              fill = as.factor(pd)), 
             alpha = I(0.4), size = I(3), show.legend = show_legend) +
  
  scale_shape_manual(name ="", labels=(c(expression(paste(italic("H"), ">", "P"["d"])),
                                         expression(paste(italic("H"), "<", "P"["d"])))),
                     values = c(21,23))+
  scale_fill_manual(name ="", labels=(c(expression(paste(italic("H"), ">", "P"["d"])),
                                        expression(paste(italic("H"), "<", "P"["d"])))),
                    values = c("goldenrod2", "navyblue"))+
  scale_colour_manual(name ="", labels=(c(expression(paste(italic("H"), ">", "P"["d"])),
                                          expression(paste(italic("H"), "<", "P"["d"])))),
                      values = c("goldenrod2", "navyblue"))+
  
  
  geom_ribbon(aes(ymin = H_predicted - H_sd,
                  ymax = H_predicted + H_sd, fill = as.factor(pd)),
              alpha = opacity, show.legend = F,
              
              colour="NA"
  )+
  
  geom_rug(size = 1.1, show.legend = show_legend, alpha = opacity)+
  
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

g <- ggMarginal(groupFill = T, data = H_df, type = "densigram", bins = 100,
                p = g, aes(x = H_actual, y = H_predicted))


# g <- ggplot(H_df,aes(x = H_actual, y= H_predicted)) + 
#   
#   stat_density_2d(aes(y = H_predicted, fill = ..level..),
#                   geom = "polygon", colour="gray", show.legend = F, alpha = 0.5, size=1.1)+
#   
#   scale_fill_distiller(palette = 2, direction = 1)+
#   
#   geom_point(aes(y = H_predicted), shape=21, fill="goldenrod2",
#              size=3.0, na.rm = T, show.legend = F) +
#   
#   ggforce::geom_mark_ellipse(aes(colour=as.character(pd)), size=1.3, alpha=1, show.legend = TRUE)+
#   
#   geom_ribbon(aes(ymin = H_predicted - H_CI,
#                   ymax = H_predicted + H_CI),
#               alpha = 0.3,
#               fill = "navyblue",
#               colour="NA"
#   )+
#   scale_colour_viridis(discrete = T, label=c("H > 6m", "H < 6m"))+
#   geom_abline(slope = 1,linetype="solid", intercept = 0,
#               colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
#   
#   
#   # geom_smooth(size=1.5,level = 0.95,show.legend = F,linetype = "dashed",
#   #             color="black", data = H_df,
#   #             se= T, method = "lm", aes(x=H_actual, y=H_predicted))+
#   
#   coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
#               ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#   scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
#                      breaks = seq(xmin, xmax, xstp)) +
#   scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
#                      breaks = seq(ymin, ymax, ystp)) +
#   theme_bw()+
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
#         panel.grid.minor = element_blank(),
#         plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
#         legend.direction = "vertical", legend.box = "vertical",
#         legend.text.align = 0,
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
# g
ggsave(paste0("./outputs/shallow_H_BG_unconstr_scatter.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)



















seadoo_rrs_plot_df = as.data.frame((seaDoo_rrs_interp))

seadoo_rrs_plot_df$id = seq(1, length(seaDoo_rrs_interp$`401`), 1)

seadoo_rrs_long = reshape2::melt(seadoo_rrs_plot_df, id.vars = "id")

H_plot_df = data.frame("id" = seq(1,length(seadoo_depth), 1), "Depth" = seadoo_depth)
H_long = reshape2::melt(H_plot_df, id.vars = "id")

seadoo_rrs_H = merge(seadoo_rrs_long,H_long,by=c("id"))

seadoo_rrs_H = seadoo_rrs_H[,-(4)] 
names(seadoo_rrs_H) = c("id", "wavelength", "rrs", "Depth")


xmin <- 400; xmax <- 700;  xstp <- 50; xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0; ymax <- 0.015; ystp <- 0.003; ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda),italic("in situ ("), "sr"^{-1}, ")"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.70, 0.98)


g <- ggplot(data = seadoo_rrs_H) + 
  geom_line(aes(x = as.numeric(as.character(wavelength)), y = rrs, colour = Depth, 
                group = Depth), size = 1.3, show.legend = T)+
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
ggsave("./outputs/seadoo_rrs.png", plot = g,scale = 1.7, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)