library(dplyr)
library(tidyr)
library(pbapply)
library(parallel)
library(doParallel)
library(foreach)

# Function :: Interpolate the bottom reflectance to Rrs wavelengths ----
interp_rb <- function(wavelength_input, bottom_type = c("Mud_2019", "Sand_2019", "Eelgrass_2019")){
  load("./data/EGSL_all.RData")
  
  column_numbers <- which(names(rb) %in% bottom_type)
  
  rb = rb[c(1,column_numbers)]
  
  print(paste0("Classnames reassigned as: Class1:", colnames(rb)[2], ", Class2:", colnames(rb)[3],
               ", Class3:", colnames(rb)[4]))
  
  names(rb) = c("wavelength", "class1", "class2", "class3")
  
  
  if (!(all(length(rb$wavelength) == length(wavelength_input)) && all(rb$wavelength == wavelength_input))) {
    
    print(paste0("Rb wavelength has " , length(rb$wavelength), " bands; Rrs wavelength has ", length(wavelength_input) ))
    
    rb_interp = data.frame("wavelength" = wavelength_input,
                           "class1" = approx(rb$wavelength, rb$class1, xout = wavelength_input)$y,
                           "class2" = approx(rb$wavelength, rb$class2, xout = wavelength_input)$y,
                           "class3" = approx(rb$wavelength, rb$class3, xout = wavelength_input)$y
    )
    
  }
  else {
    print("The desired and existing wavelength are same")
    rb_interp = rb
  }
  rb = rb_interp
  rm(rb_interp)
  
  assign(x = "rb", value = rb, envir = .GlobalEnv)
}

# Read shallow water database  ----

site_input = "BG"
year = 2022

if (year == 2022) {
  
  bathy_data = read.csv("./data/bathy_db_long.csv", header = T)
  
}

if (year == 2023) {
  
  bathy_data = read.csv("./data/algae_wise_2023.csv", header = T)
  
}


if (site_input == "BG") {
  
  bathy_data = bathy_data[bathy_data$BoatSolAzm >= 90 & bathy_data$BoatSolAzm <= 180,]
  
} 

bathy_data = bathy_data[!is.na(bathy_data$Site),]
bathy_data = bathy_data[!is.na(bathy_data$OptShallow),]
bathy_data$BottomElevation_m = abs(bathy_data$BottomElevation_m)

bathy_data <- bathy_data %>%
  mutate(DateOnly = format(ymd_hms(DateTime), "%d_%m_%y"))

ggplot(data = bathy_data[bathy_data$Site == "BG",]#[BG_HS_idx$...1,]
) +
  geom_line(aes(x = as.numeric(as.character(Wavelength)), y = BRI, colour = Hwc,
                linetype = Site, group = Hwc
  ), size = 1.3, alpha = 0.6, show.legend = T)+
  scale_colour_viridis(discrete = F, name = "Depth (m)")+
  ylim(c(0,1))

# Visualize Rrs ----
xmin <- 400; xmax <- 750;  xstp <- 50; xlbl <- expression(paste("Wavelength (", lambda, ") [nm]"))
ymin <- 0; ymax <- 0.010; ystp <- 0.002
ylbl <- expression(paste(italic("R")["rs"]("0"^"+", lambda),italic("in situ ("), "sr"^{-1}, ")"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.70, 0.98)

g <- ggplot(data = bathy_data) +
  geom_line(aes(x = as.numeric(as.character(Wavelength)), y = Rrs, colour = Hwc,
                linetype = DateOnly,group = Hwc
  ), size = 1.0, alpha = 0.6, show.legend = T)+
  scale_colour_viridis(discrete = F, name = "Depth (m)") +
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax)
              ,expand = FALSE, clip = "on"
  ) +
  #scale_colour_manual(name = "",values = rev(collist))+
  guides(linetype = guide_legend(title = ""))+
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
ggsave(paste0("./outputs/shallow_rrs_", year,".png"), plot = g, scale = 1.7, 
       width = 4.5, height = 4.5,
       units = "in",dpi = 300)

# Visualize QWIP for Rrs ----
xmin <- 510; xmax <- 570; xstp <- 20
xlbl <- "AVW"
ymin <- -1; ymax <- 1 ; ystp <- 0.5
ylbl <- "NDI(490,665)"

g1 <- ggplot(data = bathy_data[bathy_data$Site == "BG",])  +
  geom_line(aes(x =AVW, y= QWIP),color = "navyblue", size=1.3,show.legend = F)+
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
  geom_point(aes(x = AVW, y=NDI491665, color= as.numeric(Hwc), shape=as.character(DateOnly)
  ),
  size=1.8, alpha = 0.60, show.legend = T)+
  
  geom_hline(yintercept = 0, col = "grey", linetype = "dashed", size = 1.3)+
  
  scale_color_viridis(discrete = F, name = "Depth (m)") +
  
  # scale_shape_manual(labels=c("Optically Deep","Optically Shallow"),
  #                    values = c(17,19))+
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
        legend.position=c(0.6, 0.95),
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

ggsave(paste0("./outputs/shallow_qwip_", year,".png"), 
       plot = g1,scale = 1.7, width = 4.5, height = 4.5,
       units = "in",dpi = 300)

#Create data-structure for the inversion ----
bathy_data_list = split(bathy_data, bathy_data$Lat)

names(bathy_data_list) = rep(paste0("Obs_shallow_", seq(1, length(bathy_data_list),1)))

lat_lon = data.frame("lat" = unique(bathy_data$Lat)[-1], 
                     "lon" = unique(bathy_data$Lon))
write.csv(file = paste0("./outputs/shallow_lat_lon-",site_input ,"_",year,".csv"), 
          x = lat_lon, col.names = TRUE,
          quote = F, row.names = F)

## Create a list where each entry contains a data frame with wavelength and Rrs and a single vector for H ----
seadoo_algaeWISE_data <- lapply(bathy_data_list, function(df_observation) {
  
  
  wavelength <- df_observation$Wavelength
  Rrs <- df_observation$Rrs
  H <- unique(df_observation$Hwc)
  KLu <- df_observation$KLu
  z90 <- df_observation$Z90
  site <- df_observation$Site
  optshallow <- df_observation$OptShallow
  
  benthic_veg_frac <- df_observation$PercentCoverage
  plant_height <- df_observation$PlantHeight_m
  
  chl = df_observation$Chl
  vsf_700 = df_observation$VSF700
  
  
  if (all(optshallow) == TRUE) {
    
    opt_stat = "shallow"
    
  } else {
    opt_stat = "not_shallow"
  }
  
  # Create a list for wavelength and Rrs and H
  df_wavelength_Rrs <- data.frame(wavelength, Rrs, z90)
  list(wavelength_Rrs = df_wavelength_Rrs, H = -H, site = unique(site), opt_stat = opt_stat,
       chl = unique(chl), vsf700 = unique(vsf_700),
       benthic_veg_frac = unique(benthic_veg_frac), 
       plant_height = unique(plant_height))
  
  
})

## Convert the list into data-frame  ----
seadoo_algaeWISE_data_df <- lapply(bathy_data_list, function(df_observation) {
  
  wavelength <- df_observation$Wavelength
  Rrs <- df_observation$Rrs
  H <- unique(df_observation$H)
  
  # Create a list for wavelength and Rrs and H
  df_Rrs <- data.frame(Rrs)
  #list(wavelength_Rrs = df_wavelength_Rrs, H = -H)
})

seadoo_algaeWISE_data_df = data.frame(t(do.call(cbind, seadoo_algaeWISE_data_df)))

filtered_list <- lapply(seadoo_algaeWISE_data_df, function(df) {
  if (ncol(df) == 149) {
    return(df)
  } else {
    return(NULL)
  }
})

filtered_list <- Filter(Negate(is.null), filtered_list)


seadoo_algaeWISE_data_df <- data.frame(t(do.call(cbind, filtered_list)))


#seadoo_algaeWISE_data_df = do.call(cbind, seadoo_algaeWISE_data_df)

wavelength = unique(bathy_data$Wavelength)
names(seadoo_algaeWISE_data_df) = wavelength

seadoo_algaeWISE_data_df = seadoo_algaeWISE_data_df[!is.na(names(seadoo_algaeWISE_data_df))]

#Plot the Rrs
matplot(wavelength, 
        t(seadoo_algaeWISE_data_df), col = viridis(5), type = "l", lwd=3)

# Extract Depths from the list
seadoo_depth <- sapply(seadoo_algaeWISE_data, function(entry) entry$H)

if(year == 2022){
  seadoo_depth <- as.data.frame(do.call(rbind, seadoo_depth))
  
} else {
  seadoo_depth <- as.data.frame(seadoo_depth)
  
}


seadoo_z90 <- sapply(seadoo_algaeWISE_data, function(entry) entry$wavelength_Rrs$z90)

seadoo_site <- sapply(seadoo_algaeWISE_data, function(entry) entry$site)
seadoo_site <- as.data.frame(seadoo_site)


seadoo_chl <- sapply(seadoo_algaeWISE_data, function(entry) entry$chl)
if(year == 2022){
  seadoo_chl <- as.data.frame(do.call(rbind, seadoo_chl))
  
} else {
  seadoo_chl <- as.data.frame(seadoo_chl)
  
}

seadoo_vsf700 <- sapply(seadoo_algaeWISE_data, function(entry) entry$vsf700)
if(year == 2022){
  seadoo_vsf700 <- as.data.frame(do.call(rbind, seadoo_vsf700))
  
} else {
  seadoo_vsf700 <- as.data.frame(seadoo_vsf700)
  
}


seadoo_bottom_veg_cover <- sapply(seadoo_algaeWISE_data, function(entry) entry$benthic_veg_frac)
if(year == 2022){
  seadoo_bottom_veg_cover <- as.data.frame(do.call(rbind, seadoo_bottom_veg_cover))
  
} else {
  seadoo_bottom_veg_cover <- as.data.frame(seadoo_bottom_veg_cover)
  
}


seadoo_plant_height <- sapply(seadoo_algaeWISE_data, function(entry) entry$plant_height)

if(year == 2022){
  seadoo_plant_height <- as.data.frame(do.call(rbind, seadoo_plant_height))
  
} else {
  seadoo_plant_height <- as.data.frame(seadoo_plant_height)
  
}


BG_idx = which(seadoo_site == site_input)

# Apply SABER on the shallow water dataset of Rrs---- 

## Pre-treat Rrs data for inversion ----
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
        t(seaDoo_rrs_interp[BG_idx,]), col = viridis(5), type = "l", lwd=3)

## Prepare the optimization parameters ----
pop.sd = "unknown" 

constrain.bbp= F
constrain.shallow = F
constrain.shallow.bgc = F ; constrain.shallow.iop = F

manual_par0 = T ; manual_bound = T


## Set-up initial values and run parameters for shallow water inversion ----

rb_count = length(fA.set)
rep(0.5, rb_count)
rb_c = paste0("rb", (1:rb_count))

pop.sd = "unknown" ;constrain.bbp= F
constrain.shallow.bgc = F ; constrain.shallow.iop = F; manual_par0 = T
type_Rrs_below = "shallow"

inv_bound = create_init_bound(rrs_inv_type = type_Rrs_below, manual_par0 = manual_par0, 
                              constrain.bbp = constrain.bbp, 
                              constrain.shallow.bgc = constrain.shallow.bgc, 
                              constrain.shallow.iop = constrain.shallow.iop, 
                              pop.sd =  pop.sd, 
                              init_par = c(#mean(seadoo_chl$seadoo_chl[BG_idx], na.rm = T),
                                          mean(seadoo_chl$V1[BG_idx], na.rm = T), 
                                           0.25, 
                                           0.003,
                                           mean(seadoo_depth$V1[BG_idx], na.rm = T), 
                                           #mean(seadoo_depth$seadoo_depth, na.rm = T),
                                           
                                           0.5,0.5,0.5,0.05),
                              
                              upper_par = c(#max(seadoo_chl$seadoo_chl[BG_idx], na.rm = T),
                                            max(seadoo_chl$V1[BG_idx], na.rm = T), 
                                            0.5, 
                                            0.005,
                                            max(seadoo_depth$V1[BG_idx], na.rm = T), 
                                            #max(seadoo_depth$seadoo_depth, na.rm = T),
                                            1,1,1,10),
                              
                              lower_par = c(min(seadoo_chl$V1[BG_idx], na.rm = T),
                                            #min(seadoo_chl$seadoo_chl[BG_idx], na.rm = T), 
                                            0.05, 
                                            0.001,
                                            min(seadoo_depth$V1[BG_idx], na.rm = T),
                                            #min(seadoo_depth$seadoo_depth, na.rm = T),
                                            0.0,0.0,0.0,0.00001)
)

par0 = inv_bound$par0 ; upper.bound = inv_bound$upper.bound  
lower.bound = inv_bound$lower.bound

# Selection of inversion optimizer and objective function
obj = c("log-LL", "SSR", "obj_L98"); obj.run = obj[3]

methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                 "Brent","levenberg-marqardt", "auglag")

samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")


na_idx = is.na(seadoo_chl$V1)

# TEST the Shallow Water unconstrained inversion
set.seed(345)
test_idx = sample(size = 50, x = seq(1,length(BG_idx),1))

if(inel_cor == "OFF") {
  #obsdata = as.numeric(seaDoo_rrs_interp[BG_idx,][test_idx,])
  obsdata = (seaDoo_rrs_interp[BG_idx,][test_idx,])
} else {
  
  obsdata = as.numeric(seaDoo_rrs_interp_inel_corrected[test_idx,])
}

insitu_val = data.frame("chl" = seadoo_chl$V1[BG_idx][test_idx],
                        "H" = seadoo_depth$V1[BG_idx][test_idx],
                        "vsf700" = seadoo_vsf700$V1[BG_idx][test_idx]
                        )

Rb_set = read.csv("./data/EGSL_Rb.csv", header = T)
benthic_classes_2019 <- colnames(Rb_set)[grep("2019", colnames(Rb_set))]
benthic_classes_2022 <- colnames(Rb_set)[grep("2022", colnames(Rb_set))]
benthic_classes_2023 <- colnames(Rb_set)[grep("2023", colnames(Rb_set))]

combinations <- combn(benthic_classes_2023, 3, simplify = FALSE)

combinations <- lapply(combinations, function(x) sort(x))
combinations <- as.data.frame(do.call(rbind, combinations), stringsAsFactors = FALSE)

#combinations = as.data.frame(do.call(rbind, combinations))


param_vec_grad = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0, 
                       "H" = 0, "fa1" = 0, "fa2" = 0, "fa3" = 0, "pop_sd" = 0)

param_vec_bayes = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0, 
                            "H" = 0, "fa1" = 0, "fa2" = 0, "fa3" = 0, "pop_sd" = 0,
                            "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
                            "sd_H" = 0, "sd_fa1" = 0, "sd_fa2" = 0, "sd_fa3" = 0, "sd_pop_sd" = 0) 

# spectral_slope_vec = expand.grid(data.frame("ag" = seq(0.001,0.020, 0.001*2),
#                                             "ad" = seq(0.0005,0.005, 0.000225*2)[-10],
#                                             "gamma" = seq(0.1,1,0.045*2)[-10]
#                                             ))


ad_ag_values <- seq(0.0015, 0.025, length.out = 10)
gamma_values <- seq(0.1, 1, length.out = 10)
spectral_slope_vec <- expand.grid(ad_ag = ad_ag_values, gamma = gamma_values)

set.seed(123)  
spectral_slope_vec$ad <- runif(n = nrow(spectral_slope_vec), min = 0.0005, max = 0.005)
spectral_slope_vec$ag <- spectral_slope_vec$ad_ag - spectral_slope_vec$ad

spectral_slope_vec <- spectral_slope_vec %>%
  filter(ad != ag & ag >= 0.001 & ag <= 0.02 & ad + ag == ad_ag)

#spectral_slope_vec <- spectral_slope_vec[1:100, c("ag", "ad", "gamma")]

param_vec_bayes <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(param_vec_bayes) <- c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd", 
                               "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", "sd_fa1", "sd_fa2", 
                               "sd_fa3", "sd_pop_sd")

num_cores <- detectCores() - 1  
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Create a grid of all combinations of i and j
combinations_grid <- expand.grid(i = seq_along(combinations$V1), j = seq_along(spectral_slope_vec$ag))

# Define a function to run your process for a given combination of i and j
process_combination <- function(idx) {
  i <- combinations_grid[idx, "i"]
  j <- combinations_grid[idx, "j"]
  
  # Suppress output from the functions
  capture.output({
    interp_rb(wavelength_input = wavelength, bottom_type = as.character(combinations[i,]))
    
    result <- inverse_runBayes(obsdata = as.numeric(obsdata[10,]
    ), 
    rrs_type = type_Rrs_below, 
    max_par = upper.bound, 
    min_par = lower.bound, 
    param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", 
                  "pop_sd"),
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
  }, file = NULL) 
  
  return(result)
}

results_list <- foreach(idx = 1:nrow(combinations_grid), 
                        .combine = rbind, 
                        .verbose = T,
                        .packages = c("data.table", "BayesianTools")) %dopar% {
                          source("./R/SABER_forward_fast.R")
                          source("./R/solve.objective.inverse_fast.R")
                          source("./R/mcmc_bayes_parallel.R")
                          process_combination(idx)
}

param_vec_bayes <- as.data.frame(results_list)
stopCluster(cl)



library(parallel)
library(pbapply)
library(foreach)
library(doParallel)

# Create cluster
num_cores <- detectCores() - 1  # Use one less than the total number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Load required libraries on each worker
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

rrs_sample = partition_data_for_multiprocessing(num_groups = 8, data_input = obsdata)

# Function to process combinations
process_combination <- function(obsdata, idx) {
  i <- combinations_grid[idx, "i"]
  j <- combinations_grid[idx, "j"]
  
  # Suppress output from the functions
  capture.output({
    interp_rb(wavelength_input = wavelength, bottom_type = as.character(combinations[i,]))
    
    result <- inverse_runBayes(obsdata = as.numeric(obsdata), 
                               rrs_type = type_Rrs_below, 
                               max_par = upper.bound, 
                               min_par = lower.bound, 
                               param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
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
  }, file = NULL) 
  
  return(result)
}

# Create a grid of all combinations of i and j
combinations_grid <- expand.grid(i = seq_along(combinations$V1), j = seq_along(spectral_slope_vec$ag))

# Function to process each observation
process_observation <- function(obsdata) {
  results_list <- pblapply(1:nrow(combinations_grid), function(idx) {
    process_combination(obsdata, idx)
  }, cl = cl)
  
  # Combine results into a data frame
  result_df <- as.data.frame(do.call(rbind, results_list))
  return(result_df)
}

# results_list <- pblapply(1:nrow(combinations_grid), function(idx) {
#   process_combination(idx)
# }, cl = cl) 
# 
# 
# param_vec_bayes <- as.data.frame(do.call(rbind, results_list))

# Export required objects to cluster
clusterExport(cl, c("process_combination", "combinations_grid", 
                    "spectral_slope_vec", "create_init_bound",
                    "wavelength", "combinations", "interp_rb", "inverse_runBayes", 
                    "upper.bound", "lower.bound", "type_Rrs_below", "samplerlist", 
                    "fA.set", "obsdata", "type_case_water"))

obs_test = process_observation(obsdata = as.numeric(obsdata[10,]))

final_results <- foreach(obsdata = rrs_sample_samp, 
                         #.combine = 'rbind', 
                         .verbose = T,
                         .packages = c('BayesianTools','parallel', 'pbapply')) %dopar% {
                           source("./R/SABER_forward_fast.R")
                           source("./R/solve.objective.inverse_fast.R")
                           source("./R/mcmc_bayes_parallel.R")
                           process_observation(obsdata)
}

# Stop the cluster
stopCluster(cl)

# Print the final results
print(final_results)


























# # Create a grid of all combinations of i and j
# combinations_grid <- expand.grid(i = seq_along(combinations$V1), j = seq_along(spectral_slope_vec$ag))
# 
# # Define a function to run your process for a given combination of i and j
# process_combination <- function(idx) {
#   i <- combinations_grid[idx, "i"]
#   j <- combinations_grid[idx, "j"]
#   
#   # Suppress output from the functions
#   capture.output({
#     interp_rb(wavelength_input = wavelength, bottom_type = as.character(combinations[i,]))
#     
#     result <- inverse_runBayes(obsdata = as.numeric(obsdata[10,]
#                                                     ), 
#                                rrs_type = type_Rrs_below, 
#                                max_par = upper.bound, 
#                                min_par = lower.bound, 
#                                param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", 
#                                              "pop_sd"),
#                                constrain_config = c(FALSE, FALSE, FALSE),
#                                bbp_const_val = 0.007, 
#                                bgc_const_val = c(5, 1.003, 0.007), 
#                                iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
#                                                   "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
#                                qaa_slope = FALSE, 
#                                manual_slope = TRUE, 
#                                manual_slope_vals = c("s_g" = spectral_slope_vec[j, 4], 
#                                                      "s_d" = spectral_slope_vec[j, 3], 
#                                                      "gamma" = spectral_slope_vec[j, 2]), 
#                                iter_count = 20000, 
#                                sampler_mcmc = samplerlist[6], 
#                                wavelngth_sim = wavelength, 
#                                sa_model = "am03", 
#                                hybrid_mode = FALSE, 
#                                plot_rrs = FALSE)
#   }, file = NULL) 
#   
#   return(result)
# }
# 
# 
# results_list <- pblapply(1:nrow(combinations_grid), function(idx) {
#   process_combination(idx)
# }, cl = cl) 
# 
# 
# param_vec_bayes <- as.data.frame(do.call(rbind, results_list))


seadoo_chl$seadoo_chl[test_idx][10]
seadoo_depth$seadoo_depth[test_idx][10]


rrs_sample = partition_data_for_multiprocessing(num_groups = 6, data_input = obsdata)


# for (i in 1:length(combinations$V1)) {
#   
#   interp_rb(wavelength_input = wavelength, bottom_type = as.character(combinations[i,]))
#   
#   # param_vec_grad[i,] = doOptimization_shallow_unconst(obsdata = as.numeric(obsdata[1,]), 
#   #                                                     init_par = par0, 
#   #                                                 max_par = upper.bound, 
#   #                                                 min_par = lower.bound,
#   #                                QAA_mode = F,
#   #                                 wavelength_sim = wavelength, qaa_prefit = F, 
#   #                                qaa_slope = F, manual_slope = F, sa_model = "am03", 
#   #                                obj_fn = obj[1], opt_method = methods.opt[8] )
#   
#   # doOptimization_shallow_unconst_back(obsdata = obsdata, par0 = par0, wl = wavelength, 
#   #                                     sa_model = "am03", obj_fn = obj[1], 
#   #                                     opt_method = methods.opt[8])
#   for(j in 1:length(spectral_slope_vec$ag))
#   param_vec_bayes[i,j] = inverse_runBayes(obsdata = as.numeric(obsdata[10,]), 
#                    rrs_type =  type_Rrs_below, 
#                    max_par = upper.bound, 
#                    min_par = lower.bound, 
#                    param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
#                    constrain_config = c(F,F,F),
#                    bbp_const_val = 0.007, 
#                    bgc_const_val = c(5,1.003, 0.007), 
#                    
#                    iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
#                                       "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
#                    
#                    qaa_slope = F, 
#                    manual_slope = T, 
#                    manual_slope_vals  = c("s_g"=spectral_slope_vec[j,1], 
#                                           "s_d"=spectral_slope_vec[j,2], 
#                                           "gamma"=spectral_slope_vec[j,3]), 
#                    iter_count = 20000, 
#                    sampler_mcmc = samplerlist[6], 
#                    wavelngth_sim = wavelength, 
#                    sa_model = "am03", 
#                    hybrid_mode = F, plot_rrs = F)
#   
#   
#   
# }
