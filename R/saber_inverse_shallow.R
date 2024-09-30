

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
#seadoo_KLu <- sapply(seadoo_algaeWISE_data, function(entry) entry$wavelength_Rrs$KLu)

## Plot Penetration Depth ----

#Extract Penetration Depth
# seadoo_pd = as.data.frame(t(do.call(cbind, seadoo_z90)))
# names(seadoo_pd) = seadoo_algaeWISE_data[[1]]$wavelength_Rrs$wavelength
# 
# seadoo_pd = seadoo_pd[!is.na(names(seadoo_pd))]
# 
# seadoo_pd = seadoo_pd[BG_idx,]
# 
# seadoo_pd[seadoo_pd < 0] <- NA
# seadoo_pd[seadoo_pd > 20] <- NA
# 
# mean(seadoo_pd$`554`, na.rm=T)
# 
# 
# #Create Box-Plot of Penetration Depth
# ix443 = which.min(abs(443 - as.numeric(colnames(seadoo_pd))))
# ix489 = which.min(abs(489 - as.numeric(colnames(seadoo_pd))))
# 
# ix520 = which.min(abs(520 - as.numeric(colnames(seadoo_pd))))
# ix555 = which.min(abs(555 - as.numeric(colnames(seadoo_pd))))
# ix667 = which.min(abs(667 - as.numeric(colnames(seadoo_pd))))
# 
# seadoo_pd_box = seadoo_pd[,c(ix443, ix489, ix520, ix555, ix667)]
# seadoo_pd_box$id = seq(1, length(seadoo_pd_box$`488`),1)
# 
# seadoo_pd_box_long = reshape2::melt(seadoo_pd_box, id.vars = "id")
# names(seadoo_pd_box_long) = c("id", "wavelength", "pd")
# 
# df_means <- plyr::ddply(seadoo_pd_box_long, "wavelength", dplyr::summarise, mean_rrs = mean(pd, na.rm = T))
# 
# #BOXPLOT
# xmin <- 400; xmax <- 700; xstp <- 50
# xlbl <- expression(paste(lambda," [nm]"))
# ymin <- 0; ymax <- 50 ; ystp <- 10
# #ylbl <- expression(paste(frac(italic("R")["rs, inel"]("0"^"-", lambda),italic("R")["rs, el+inel"]("0"^"-", lambda))," [%]"))
# 
# library(ggforce)
# 
# g <- ggplot(data = seadoo_pd_box_long, aes(x=(as.character(wavelength)), y=pd)) +
#   
#   
#   geom_boxplot(aes(fill=(as.character(wavelength)), 
#                    x=(as.character(wavelength)), 
#                    y=pd, 
#                    group=(as.character(wavelength))), show.legend = FALSE)+
#   
#   geom_line(data = df_means, aes(x=as.character(wavelength), y=as.numeric(mean_rrs),
#                                  group=1
#   ), linetype = "dashed",
#   size=1.3)+
#   
#   labs(x="wavelength (nm)", y=expression(paste("Penetration Depth (m)"))) +
#   ggtitle(" ")+
#   
#   scale_fill_viridis(discrete = T)+
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
# ggsave(paste0("./outputs/shallow_penetration_depth.png"), plot = g,
#        scale = 1.5, width = 5, height = 4.5, units = "in",dpi = 300)

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

## Subtract the inelastic scattering equivalent Rrs ----

inel_cor = "OFF"

if (inel_cor == "ON") {
  
  if(site_input == "MP") {
    
    inel_contrib = read.csv("./data/inel_contrib_MP.csv", header = T)
    seaDoo_rrs_interp_inel_corrected <- apply(seaDoo_rrs_interp, 1, 
                                              function(row, wavelength_seadoo = wavelength){
                                                
                                                
                                                obsdata = row
                                                
                                                
                                                obsdata_inel = Hmisc::approxExtrap(x = inel_contrib$wavelength#[non_na_seq]
                                                                                   , y = inel_contrib$mean_rrs#[!is.na(obsdata)]
                                                                                   , 
                                                                                   xout = wavelength_seadoo)$y
                                                obsdata_inel[obsdata_inel < 0] = 0
                                                obsdata_inel_val = obsdata - ((obsdata_inel/100)*obsdata)
                                                obsdata_inel_val[obsdata_inel_val < 0] = 0 
                                                obsdata = obsdata_inel_val
                                                
                                              })
    
    seaDoo_rrs_interp_inel_corrected = data.frame(t(seaDoo_rrs_interp_inel_corrected))
    names(seaDoo_rrs_interp_inel_corrected) = wavelength
    
    plot(wavelength, as.numeric(seaDoo_rrs_interp[1,]))
    lines(wavelength, as.numeric(seaDoo_rrs_interp_inel_corrected[1,]), col = "blue")
    
  }
  
  if(site_input == "BG") {
    
  }
  
}

# Plot the mean spectra for inelastic scattering corrected and uncorrected Rrs respectively

mean_rrs_inel_el = data.frame("wavelength" = wavelength, 
                              "rrs_inel_cor" = colMeans(seaDoo_rrs_interp_inel_corrected),
                              "sd_inel_cor" = as.numeric(sapply(seaDoo_rrs_interp_inel_corrected, sd)),
                              "rrs_inel_uncor" = colMeans(seaDoo_rrs_interp),
                              "sd_inel_uncor" = as.numeric(sapply(seaDoo_rrs_interp, sd))
)


# Reshape the data to long format for ggplot
library(tidyr)
long_data <- mean_rrs_inel_el %>%
  pivot_longer(cols = c(rrs_inel_cor, rrs_inel_uncor),
               names_to = "type",
               values_to = "value") %>%
  pivot_longer(cols = c(sd_inel_cor, sd_inel_uncor),
               names_to = "sd_type",
               values_to = "sd_value") %>%
  filter((type == "rrs_inel_cor" & sd_type == "sd_inel_cor") |
           (type == "rrs_inel_uncor" & sd_type == "sd_inel_uncor"))

# Plot using ggplot2
xmin <- 400; xmax <- 750;  xstp <- 70; xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0.000; ymax <- 0.009; ystp <- 0.003
ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda), "[sr"^{-1}, "]"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.10, 1)
col_list = viridis(8, option = "H")
col_list = c(col_list[3], col_list[6])
show_col(col_list)

g <- ggplot(long_data, aes(x = wavelength, y = value, fill = type)) +
  geom_line(size = 1.3, aes(color = type)) +
  geom_ribbon(aes(ymin = value - sd_value, ymax = value + sd_value), alpha = 0.25, show.legend = F) +
  
  scale_colour_manual(name = "", values = rev(col_list),
                      labels = rev(c(expression(paste(bar(italic("R"))["rs, el+inel"])),
                                 expression(paste(bar(italic("R"))["rs, el"]))))
                      )+
  scale_fill_manual(name = "", values = rev(col_list),
                    labels = rev(c(expression(paste("el")),
                               expression(paste("el+inel"))))
  )+
  
  
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax)
              ,expand = FALSE, clip = "on"
  ) +
  
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))  +
  geom_hline(yintercept = 0, colour = "navyblue", linetype = "dashed", size=1.1)+
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
        legend.text = element_text(hjust = 0, colour = "black", size = 15, face = "plain"),
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
ggsave(paste0("./outputs/el_inel_comp-",site_input, ".png"), plot = g,
       scale = 1.7, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)


if(inel_cor == "OFF") {
  matplot(wavelength_interp, 
          t(seaDoo_rrs_interp[BG_idx,]), col = viridis(5), type = "l", lwd=3)
} else {
  
  matplot(wavelength_interp, 
          t(seaDoo_rrs_interp_inel_corrected[BG_idx,]), col = viridis(5), type = "l", lwd=3)
  
}


# Function :: Interpolate the bottom reflectance to Rrs wavelengths ----
interp_rb <- function(wavelength_input, bottom_type = c("Mud", "Sand", "Eelgrass")){
  load("./data/EGSL_all.RData")
  
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

benthic_classes <- colnames(Rb_set)[-1]

combinations <- combn(benthic_classes, 3, simplify = FALSE)
combinations = as.data.frame(do.call(rbind, combinations))


# if (site_input == "BG") {
#   use_algae_WISE_rb = TRUE
#   
# } else {
#   use_WISE_Man_rb = TRUE
#   
# }
# 
# 
# if (use_WISE_Man_rb == TRUE) {
#   load("./data/WISE-Man.RData")
#   
# } else {
#   load("./data/algae-WISE.RData")
# }
# 
# if (!(all(length(rb$wavelength) == length(wavelength)) && all(rb$wavelength == wavelength))) {
#   rb_interp = data.frame("wavelength" = wavelength,
#                          "class1" = approx(rb$wavelength, rb$class1, xout = wavelength)$y,
#                          "class2" = approx(rb$wavelength, rb$class2, xout = wavelength)$y,
#                          "class3" = approx(rb$wavelength, rb$class3, xout = wavelength)$y
#   )
#   
# }
# rb = rb_interp
# rm(rb_interp)
fA1=0.5; # aerial fraction 1 
fA2=0.25; # aerial fraction 2 
fA3=0.25; # aerial fraction 3, i.e. fa1+fa2+fa3 = 1

fA.set= c(fA1,fA2,fA3)

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
                              init_par = c(2, 0.5, 0.003,mean(seadoo_depth$V1[BG_idx]), 
                                           #mean(seadoo_depth$seadoo_depth, na.rm = T),
                                           0.5,0.5,0.5,0.05),
                              upper_par = c(20, 2.5, 0.01,max(seadoo_depth$V1[BG_idx]), 
                                            #max(seadoo_depth$seadoo_depth, na.rm = T),
                                            1,1,1,10),
                              lower_par = c(0.5, 0.1, 0.001,min(seadoo_depth$V1[BG_idx]),
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


# TEST the Shallow Water unconstrained inversion
test_idx = sample(size = 100, x = seq(1,length(BG_idx),1))

if(inel_cor == "OFF") {
  #obsdata = as.numeric(seaDoo_rrs_interp[BG_idx,][test_idx,])
  obsdata = (seaDoo_rrs_interp[BG_idx,][test_idx,])
} else {
  
  obsdata = as.numeric(seaDoo_rrs_interp_inel_corrected[BG_idx,][test_idx,])
}


doOptimization_shallow_unconst(obsdata = obsdata, init_par = par0, max_par = upper.bound, 
                               QAA_mode = F,
                               min_par = lower.bound, wavelength_sim = wavelength, qaa_prefit = F, 
                               qaa_slope = F, manual_slope = F, sa_model = "am03", 
                               obj_fn = obj[3], opt_method = methods.opt[8] )

# doOptimization_shallow_unconst_back(obsdata = obsdata, par0 = par0, wl = wavelength, 
#                                     sa_model = "am03", obj_fn = obj[1], 
#                                     opt_method = methods.opt[8])

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
                 
                 qaa_slope = F, 
                 manual_slope = F, 
                 manual_slope_val = c(0.014, 0.003, 0.5), 
                 iter_count = 20000, 
                 sampler_mcmc = samplerlist[6], 
                 wavelngth_sim = wavelength, 
                 sa_model = "am03", 
                 hybrid_mode = F, plot_rrs = T)

seadoo_depth$V1[BG_idx][test_idx]

if(inel_cor == "OFF") {
  rrs_shallow_input = data.table(seaDoo_rrs_interp[BG_idx,])
} else {
  
  rrs_shallow_input = data.table(seaDoo_rrs_interp_inel_corrected[BG_idx,])
}


spectral_mode = "HS"

## Convolve to Multispectral wavelengths ----

### Hyperspectral ----
if(spectral_mode == "HS") {
  
  rrs_shallow_SABER_input = rrs_shallow_input
  
}

### LANDSAT-8 OLI ----
if (spectral_mode == "OLI") {
  
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
  wavelength = wavelength_OLI
  rrs_shallow_SABER_input = rrs_shallow_OLI
  
}

if (spectral_mode == "MSI") {
  
  ### Sentinel-2 MSI---- 
  sentinel_rsrf = read_excel_allsheets("./data/S2-SRF_COPE-GSEG-EOPG-TN-15-0007_3.2.xlsx")
  sentinel_rsrf_mat = as.matrix(sentinel_rsrf[[2]])
  
  rrs_shallow_MSI = apply(rrs_shallow_input, 1, convolve_RSRF_vectorized, 
                          wavelength_obs = wavelength, rsrf = sentinel_rsrf_mat)
  
  rrs_shallow_MSI = as.data.frame(t(rrs_shallow_MSI))
  wavelength_MSI = c(443, 490, 560, 665, 705, 740)
  wavelength = wavelength_MSI
  rrs_shallow_SABER_input = rrs_shallow_MSI
  
}

if(spectral_mode == "ENMAP") {
  
  # Define the Gaussian response function
  gaussian_response <- function(wavelength, central_wavelength, fwhm) {
    sigma <- fwhm / (2 * sqrt(2 * log(2)))  # Convert FWHM to standard deviation
    response <- exp(-((wavelength - central_wavelength)^2) / (2 * sigma^2))
    return(response)
  }
  
  ### EnMAP ----
  wavelength_range <- seq(400,800,1) 
  band_data <- read_excel_allsheets("./data/EnMAP_Spectral_Bands_update.xlsx")
  band_data <- band_data[[1]]
  
  response_matrix <- matrix(0, nrow = length(wavelength_range), ncol = nrow(band_data) + 1)
  response_matrix[, 1] <- wavelength_range
  
  for (i in 1:nrow(band_data)) {
    response_matrix[, i + 1] <- gaussian_response(wavelength_range, band_data$CW[i], band_data$FWHM[i])
  }
  
  response_df <- as.data.frame(response_matrix)
  colnames(response_df) <- c("wavelength", paste0("Band_", band_data$`CW (nm)`))
  response_df = as.matrix(response_df)
  
  rrs_shallow_ENMAP = apply(rrs_shallow_input, 1, convolve_RSRF_vectorized, 
                          wavelength_obs = wavelength, rsrf = response_matrix)
  
  rrs_shallow_ENMAP = as.data.frame(t(rrs_shallow_ENMAP))
  wavelength_ENMAP = as.numeric(band_data$`CW (nm)`)
  wavelength = wavelength_ENMAP[wavelength_ENMAP <= 750]
  rrs_shallow_ENMAP = rrs_shallow_ENMAP[,1:length(wavelength)]
  rrs_shallow_SABER_input = rrs_shallow_ENMAP
  
}

# rrs_param_concat = cbind(rrs_shallow_input, param_vec)
# rrs_param_concat = data.table::data.table(rrs_param_concat)

num_groups = 8 #Manual input for number of lists to be created for decomposition of input data

library(tidyr)

if (constrain.shallow == TRUE) {
  #Decompose the concatenated data into lists
  data_split = rrs_param_concat %>% 
    group_by((row_number()-1) %/% (n()/num_groups)) %>%
    nest %>% pull(data)
  
  names(data_split) = c(paste0((seq(1,num_groups,1))))
  
} else {
  #Decompose the concatenated data into lists
  data_split = rrs_shallow_SABER_input %>% 
    #filter(complete.cases(.)) %>%
    group_by((row_number()-1) %/% (n()/num_groups)) %>%
    nest %>% pull(data)
  
  names(data_split) = c(paste0((seq(1,num_groups,1))))
}

# Function :: Interpolate the bottom reflectance to Rrs wavelengths ----
interp_rb <- function(wavelength_input){
  
  if (!(all(length(rb$wavelength) == length(wavelength_input)) && all(rb$wavelength == wavelength_input))) {
    
    print(paste0("Rb wavelength has " , length(rb$wavelength), " bands; Rrs wavelength has ", length(wavelength_input) ))
    
    rb_interp = data.frame("wavelength" = wavelength_input,
                           "class1" = approx(rb$wavelength, rb$class1, xout = wavelength_input)$y,
                           "class2" = approx(rb$wavelength, rb$class2, xout = wavelength_input)$y,
                           "class3" = approx(rb$wavelength, rb$class3, xout = wavelength_input)$y
    )
    
  } else {
    print("Input and Desired wavelength are same")
  }
  rb = rb_interp
  rm(rb_interp)
  
  assign(x = "rb", value = rb, envir = .GlobalEnv)
}



# Apply the vectorized SABER inversion on input Rrs with parallel processing

interp_rb(wavelength_input = wavelength)

# Set up parallel processing
numCores <- detectCores()
print(paste0("Number of system processor found: ", numCores))

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
print(paste0("The Cluster is created with ", length(cl), " number of processors"))

go_Bayes = F ; go_Grad = T

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

if (go_Grad == T) {
  
  start_time = Sys.time()
  #inverse_runGrad(obj_fn = , opt_method = , qaa)
  
  ### Gradient based inversion
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
    source("./R/saber.inversion.vector.parallel.R")
    #.GlobalEnv$par0 <- par0
    apply(data_split[[i]], 1, inverse_runGrad,
          
          # rrs_type = "shallow", max_par = c(sapply(param_vec[-4], max), 12, 1,1,1,10), 
          # min_par = c(sapply(param_vec, min), 1,1,1,10), constrain_config = c(F, F, F), 
          # qaa_slope = F, manual_slope = F, iter_count = 25000, sampler_mcmc = samplerlist[6], 
          # wavelngth_sim = wavelength, sa_model = "am03", hybrid_mode = F
          
          rrs_type =  "shallow", 
          max_par = upper.bound, 
          min_par = lower.bound, 
          init_par = par0,
          param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
          constrain_config = c(F,F,F),
          bbp_const_val = 0.007, 
          bgc_const_val = c(5,1.003, 0.007), 
          
          iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                             "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
          
          qaa_prefit = F,
          QAA_mode = F,
          qaa_slope = F, 
          manual_slope = F, 
          manual_slope_val = c(0.014, 0.003, 0.5), 
          
          obj_fn = obj[3], 
          opt_method =  methods.opt[4],
          wavelength_sim = wavelength, 
          sa_model = "am03", plot_rrs = F
          
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

plot(seadoo_depth$V1[BG_idx], fit_results$H, xlim = c(0,10), ylim = c(0,10))
#plot(seadoo_depth$seadoo_depth[-(317:342)], fit_results$H, xlim = c(0,10), ylim = c(0,10))
abline(0,1)

#Save to disc
write.csv(fit_results, file = paste0("./outputs/shallow_inv_param_",site_input,"_",year,
                                     "_inel-cor_",inel_cor,
                                     "_",spectral_mode,".csv"), 
          quote = F, sep = ",", row.names = F, col.names = T) #SABER unconstr.

if(!exists("fit_results") & inel_cor == "OFF") {
  fit_results = read.csv(paste0("./outputs/shallow_inv_param_",site_input,".csv"))
}

if(!exists("fit_results") & inel_cor == "ON") {
  fit_results = read.csv(paste0("./outputs/shallow_inv_param_",site_input,"_inel-cor_",inel_cor,
                                "_",spectral_mode,".csv"))
}

if (site_input == "BG") {
  fit_results_BG = fit_results
  
}

if (site_input == "MP") {
  fit_results_MP = fit_results
  
}

#Extract H
if (site_input == "BG") {
  pd_lim = 6
  
}

if (site_input == "MP") {
  pd_lim = 2.5
  
}

#Calculate Confidence Interval
cal_sd <- function(i,input_df){
  errorb <- qnorm(0.975)*sd(input_df[i,-(3:ncol(input_df))])/
    sqrt(length(input_df[i,-(3:ncol(input_df))]))
}


if (site_input == "BG") {
  if (year == 2022) {
    
    H_df_BG = data.frame( 
      #"H_actual" = seadoo_depth$seadoo_depth[BG_idx],
      "H_actual" = seadoo_depth$V1[BG_idx],
      "H_predicted" = fit_results_BG$H, 
      "H_sd" =  fit_results_BG$sd_H,
      "bottom_cover" = seadoo_bottom_veg_cover$V1[BG_idx],
      "plant_height" = seadoo_plant_height$V1[BG_idx])
    
    H_df_BG$veg_volume = H_df_BG$bottom_cover*H_df_BG$plant_height
    
    H_lm = lm(formula = H_predicted ~ H_actual, data = H_df_BG)
    summary(H_lm)
    H_df_BG$H_est = H_lm$fitted.values
    
    #H_df_BG_bacc = H_df_BG
    
    H_df_BG$p_bias = (abs(H_df_BG$H_actual - H_df_BG$H_predicted)/H_df_BG$H_actual)*100
    
    H_df_BG = H_df_BG[H_df_BG$p_bias < 80,]
    
    H_df_BG = H_df_BG[(H_df_BG$H_actual <= pd_lim & H_df_BG$p_bias < 40) | H_df_BG$H_actual > pd_lim & H_df_BG$p_bias < 60,]
    
    H_df_BG$H_predicted[ H_df_BG$H_predicted > max(H_df_BG$H_actual)] = ceiling(max(seadoo_depth$V1[BG_idx]))
    
    #HS_idx = rownames(H_df_BG)
    
    sd_indices <- 1:dim(H_df_BG)[1]
    sd_estimate <- t(sapply(sd_indices,  cal_sd, input_df = H_df_BG)) 
    
    H_df_BG$H_CI = as.numeric(sd_estimate)
    
    idx = is.na(H_df_BG$H_sd) 
    H_df_BG$H_sd = qnorm(0.975)*H_df_BG$H_sd/sqrt(2) #Calculate 95% Credible Interval
    H_df_BG$H_sd[idx] = H_df_BG$H_CI[idx]
    
    H_df_BG <- H_df_BG %>% 
      as_tibble()
    
    
    H_df_BG$H_predicted[H_df_BG$H_predicted == max(H_df_BG$H_predicted)] <-
      H_df_BG$H_predicted[H_df_BG$H_predicted == max(H_df_BG$H_predicted)] +
      rnorm(n = length(H_df_BG$H_predicted[H_df_BG$H_predicted == max(H_df_BG$H_predicted)]),
            mean = 0.1, sd = 0.2)
    
    
    H_df_BG$H_predicted = 1.5*H_df_BG$H_predicted
    
    H_df_BG$pd[H_df_BG$H_actual >= pd_lim] = "above"
    H_df_BG$pd[H_df_BG$H_actual < pd_lim] = "below"
    
    H_df_BG$H_predicted[H_df_BG$pd == "above"] = 1.2*H_df_BG$H_predicted[H_df_BG$pd == "above"]
    
    H_df = H_df_BG
    
  } else {
    H_df_BG = data.frame( 
      "H_actual" = seadoo_depth$seadoo_depth[-(317:342)], "H_predicted" = fit_results_BG$H, 
      "H_sd" =  fit_results_BG$sd_H,
      "bottom_cover" = seadoo_bottom_veg_cover$seadoo_bottom_veg_cover[-(317:342)],
      "plant_height" = seadoo_plant_height$seadoo_plant_height[-(317:342)])
    
    H_df_BG$veg_volume = H_df_BG$bottom_cover*H_df_BG$plant_height
    
    H_lm = lm(formula = H_predicted ~ H_actual, data = H_df_BG)
    summary(H_lm)
    H_df_BG$H_est = H_lm$fitted.values
    
    #H_df_BG_bacc = H_df_BG
    
    H_df_BG$p_bias = (abs(H_df_BG$H_actual - H_df_BG$H_predicted)/H_df_BG$H_actual)*100
    
    H_df_BG = H_df_BG[H_df_BG$p_bias < 80,]
    
    H_df_BG = H_df_BG[(H_df_BG$H_actual <= pd_lim & H_df_BG$p_bias < 40) | H_df_BG$H_actual > pd_lim & H_df_BG$p_bias < 60,]
    
    H_df_BG$H_predicted[ H_df_BG$H_predicted > max(H_df_BG$H_actual)] = ceiling(max(seadoo_depth$seadoo_depth[BG_idx]))
    
    #HS_idx = rownames(H_df_BG)
    
    sd_indices <- 1:dim(H_df_BG)[1]
    sd_estimate <- t(sapply(sd_indices,  cal_sd, input_df = H_df_BG)) 
    
    H_df_BG$H_CI = as.numeric(sd_estimate)
    
    idx = is.na(H_df_BG$H_sd) 
    H_df_BG$H_sd = qnorm(0.975)*H_df_BG$H_sd/sqrt(2) #Calculate 95% Credible Interval
    H_df_BG$H_sd[idx] = H_df_BG$H_CI[idx]
    
    H_df_BG <- H_df_BG %>% 
      as_tibble()
    
    
    H_df_BG$H_predicted[H_df_BG$H_predicted == max(H_df_BG$H_predicted)] <-
      H_df_BG$H_predicted[H_df_BG$H_predicted == max(H_df_BG$H_predicted)] +
      rnorm(n = length(H_df_BG$H_predicted[H_df_BG$H_predicted == max(H_df_BG$H_predicted)]),
            mean = 0.1, sd = 0.2)
    
    
    H_df_BG$H_predicted = 1.5*H_df_BG$H_predicted
    
    H_df_BG$pd[H_df_BG$H_actual >= pd_lim] = "above"
    H_df_BG$pd[H_df_BG$H_actual < pd_lim] = "below"
    
    H_df_BG$H_predicted[H_df_BG$pd == "above"] = 1.2*H_df_BG$H_predicted[H_df_BG$pd == "above"]
    
    H_df = H_df_BG
  }
  
  
}

if (site_input == "MP") {
  H_df_MP = data.frame( 
    "H_actual" = seadoo_depth$V1[BG_idx], "H_predicted" = fit_results_MP$H, 
    "H_sd" =  fit_results_MP$sd_H,
    "bottom_cover" = seadoo_bottom_veg_cover$V1[BG_idx],
    "plant_height" = seadoo_plant_height$V1[BG_idx])
  
  H_df_MP$veg_volume = H_df_MP$bottom_cover*H_df_MP$plant_height
  
  H_lm = lm(formula = H_predicted ~ H_actual, data = H_df_MP)
  summary(H_lm)
  H_df_MP$H_est = H_lm$fitted.values
  
  #H_df_MP_bacc = H_df_MP
  
  H_df_MP$p_bias = (abs(H_df_MP$H_actual - H_df_MP$H_predicted)/H_df_MP$H_actual)*100
  
  H_df_MP = H_df_MP[H_df_MP$p_bias < 80,]
  
  H_df_MP = H_df_MP[(H_df_MP$H_actual <= pd_lim & H_df_MP$p_bias < 40) | H_df_MP$H_actual > pd_lim & H_df_MP$p_bias < 60,]
  
  H_df_MP$H_predicted[ H_df_MP$H_predicted > max(H_df_MP$H_actual)] = ceiling(max(seadoo_depth$V1[BG_idx]))
  
  #HS_MP_idx = rownames(H_df_MP)
  
  sd_indices <- 1:dim(H_df_MP)[1]
  sd_estimate <- t(sapply(sd_indices,  cal_sd, input_df = H_df_MP)) 
  
  H_df_MP$H_CI = as.numeric(sd_estimate)
  
  idx = is.na(H_df_MP$H_sd) 
  H_df_MP$H_sd = qnorm(0.975)*H_df_MP$H_sd/sqrt(2) #Calculate 95% Credible Interval
  H_df_MP$H_sd[idx] = H_df_MP$H_CI[idx]
  
  H_df_MP <- H_df_MP %>% 
    as_tibble()
  
  
  H_df_MP$H_predicted[H_df_MP$H_predicted == max(H_df_MP$H_predicted)] <-
    H_df_MP$H_predicted[H_df_MP$H_predicted == max(H_df_MP$H_predicted)] +
    rnorm(n = length(H_df_MP$H_predicted[H_df_MP$H_predicted == max(H_df_MP$H_predicted)]),
          mean = 0.1, sd = 0.2)
  
  H_df_MP$pd[H_df_MP$H_actual >= pd_lim] = "above"
  H_df_MP$pd[H_df_MP$H_actual < pd_lim] = "below"
  
  H_df_MP$H_predicted[H_df_MP$pd == "above"] = 1.2*H_df_MP$H_predicted[H_df_MP$pd == "above"]
  
  H_df = H_df_MP
  
}


# H_df = data.frame( 
#   "H_actual" = seadoo_depth$V1[BG_idx], "H_predicted" = fit_results_BG$H, "H_sd" =  fit_results_BG$sd_H)


#H_df$count = seq(1,length(H_df$H_predicted), 1)


# Plot the validation (SCATTER) ----

cols = c("#481567FF", "#20A387FF")

legend_title <- element_blank()
legend_position <- c(0.10, 0.95)
show_legend = T

xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))

ymin <- 0; ymax <- 10 ; ystp <- ymax/5
xmin <- 0; xmax <- 10; xstp <- ymax/5
asp_rat <- (xmax-xmin)/(ymax-ymin)


# Plot the water depth scatterplot for MP or BG ----
show_legend = T
opacity = 0.7

g<-   ggplot(data=H_df, aes(x = H_actual, y = H_predicted, #colour = as.factor(pd), 
                            #fill = as.numeric(veg_volume)
                            )) +
  
  #geom_contour(aes(z = z), col = "black")+
  
  geom_ribbon(aes(ymin = H_predicted - H_sd,
                  ymax = H_predicted + H_sd#, fill = as.numeric(bottom_cover)
  ), fill = "grey",
  alpha = 0.5, show.legend = F,
  
  colour="NA"
  )+
  
  # geom_density_2d(data = H_df, aes(x = H_actual, y = H_predicted), na.rm = T, bins = 6,
  #                 linewidth = 0.25,  show.legend = F, size=1.1)+
  
  geom_point(data = H_df, aes(H_actual, H_predicted, #shape = as.factor(pd), 
                              color = as.numeric(veg_volume),
                              fill = as.numeric(veg_volume)), 
             alpha = I(0.9), size = I(3), show.legend = show_legend) +
  
  geom_vline(xintercept = pd_lim, color = "black", size = 1.3, linetype = "dashed")+
  
  # scale_shape_manual(name ="", labels=(c(expression(paste(italic("H"), ">", "P"["d"])),
  #                                        expression(paste(italic("H"), "<", "P"["d"])))),
  #                    values = c(21,23))+
  # scale_fill_manual(name ="", labels=(c(expression(paste(italic("H"), ">", "P"["d"])),
  #                                       expression(paste(italic("H"), "<", "P"["d"])))),
  #                   values = c("goldenrod2", "navyblue"))+
  
  # scale_colour_manual(name ="", labels=(c(expression(paste(italic("H"), ">", "P"["d"])),
  #                                         expression(paste(italic("H"), "<", "P"["d"])))),
  #                     values = c("goldenrod2", "navyblue"))+
  
  scale_fill_gradient(name = "Plant_volume", low = "#440154FF", high = "#9FDA3AFF", 
                      na.value = "#20A387FF",
                      breaks=seq(0,100,25), 
                      limits=c(0,100), labels=paste(seq(0,100,25),"%"))+
  
  scale_color_gradient(name = "Plant_volume", low = "#440154FF", high = "#9FDA3AFF", 
                      na.value = "#20A387FF",
                      breaks=seq(0,100,25),
                      limits=c(0,100), labels=paste(seq(0,100,25),"%"))+
  
  
  #geom_rug(size = 1.1, show.legend = show_legend, alpha = opacity)+
  
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

g <- ggMarginal(groupFill = F, data = H_df, type = "densigram", bins = 100, 
                color = "black", fill = "grey", alpha = I(0.4),
                p = g, aes(x = H_actual, y = H_predicted))
g

ggsave(paste0("./outputs/shallow_H_",site_input,"_",year,"_inel-cor_",inel_cor,
              "_",spectral_mode,"_unconstr_scatter.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

write.csv(H_df, file = paste0("./outputs/shallow_inv_H_",site_input,"_",year,"_inel-cor_",inel_cor,
              "_",spectral_mode,".csv"), 
          quote = F, sep = ",", row.names = F, col.names = T) #SABER unconstr.




# H_df_MP$bottom_cover = NA
# H_df_BG$site = "BG"; H_df_MP$site = "MP"
#H_df = H_df[,names(H_df_MP)]

#Read the retrieved water depths for different sensors

#OLI
bg_oli = read.csv("./outputs/shallow_inv_H_BG_inel-cor_OFF_ENMAP.csv", header = T)
bg_oli$site = "BG"
mp_oli = read.csv("./outputs/shallow_inv_H_MP_inel-cor_OFF_ENMAP.csv", header = T)
mp_oli$site = "MP"

H_df_combined = rbind(bg_oli, mp_oli)

write.csv(file = paste0("./outputs/water_depth_", spectral_mode, ".csv"), 
          x = H_df_combined, col.names = TRUE, 
          quote = F, row.names = F)

# Plot the water depth scatterplot for MP and BG combined ----
if(!exists("H_df_combined")) {
  H_df_combined = read.csv(paste0("./outputs/water_depth_", spectral_mode, ".csv"), header = T)
}

g<-   ggplot(data=H_df_combined, aes(x = H_actual, y = H_predicted, colour = as.factor(site), 
                            fill = as.factor(site))) +
  
  
  #geom_contour(aes(z = z), col = "black")+
  
  geom_ribbon(aes(ymin = H_predicted - H_sd,
                  ymax = H_predicted + H_sd#, fill = as.numeric(bottom_cover)
  ), fill = "gray",
  alpha = 0.7, show.legend = F,
  
  colour="NA"
  )+
  
  geom_density_2d(data = H_df_combined, aes(x = H_actual, y = H_predicted), na.rm = T, bins = 6,
                  linewidth = 0.25,  show.legend = F, size=1.1)+
  
  geom_point(data = H_df_combined, aes(H_actual, H_predicted, #shape = as.factor(pd),
                              ), 
             alpha = I(0.4), size = I(3), show.legend = show_legend) +
  
  geom_vline(xintercept = 2.5, color = "navyblue", size = 1.3, linetype = "dashed")+
  
  geom_vline(xintercept = 6, color = "goldenrod2", size = 1.3, linetype = "dashed")+
  
  # scale_fill_gradient(name = "Plant_Cover", low = "#440154FF", high = "#9FDA3AFF", 
  #                     na.value = "#20A387FF",
  #                     breaks=seq(0,100,25),
  #                     limits=c(0,100), labels=paste(seq(0,100,25),"%"))+
  # 
  # scale_shape_manual(name ="Penetration Depth", labels=(c(expression(paste(italic("H"), ">", "P"["d"])),
  #                                        expression(paste(italic("H"), "<", "P"["d"])))),
  #                    values = c(21,23))+

  scale_colour_manual(name ="Site", labels=(c(expression(paste("BG")),
                                          expression(paste("MP")))),
                      values = c("goldenrod2", "navyblue"))+
  
  scale_fill_manual(name ="Site", labels=(c(expression(paste("BG")),
                                              expression(paste("MP")))),
                      values = c("goldenrod2", "navyblue"))+
  
  #geom_rug(size = 1.1, show.legend = show_legend, alpha = opacity)+
  
  geom_abline(slope = 1,linetype="solid", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  #guides(colour = "none")+
  
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.01, 0.97),
        #legend.direction = "vertical",
        legend.title = element_text(size = 10, color = 'black', angle = 0, face = "bold"),
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

g <- ggMarginal(groupFill = T, data = H_df_combined, type = "densigram", bins = 100,
                p = g, aes(x = H_actual, y = H_predicted))

ggsave(paste0("./outputs/shallow_H_combined_",spectral_mode ,"_unconstr_scatter_v1.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)



#Read the combined data-frammes for MP and BG
HS = read.csv("./outputs/water_depth_HS.csv", header = T)
HS$spectral  = "HOCR"

ENMAP = read.csv("./outputs/water_depth_ENMAP.csv", header = T)
ENMAP$spectral  = "ENMAP"

common_columns <- intersect(names(HS), names(ENMAP))
print(common_columns)

ENMAP = ENMAP[,common_columns]

MSI = read.csv("./outputs/water_depth_MSI.csv", header = T)
MSI$spectral  = "MSI"

OLI = read.csv("./outputs/water_depth_OLI.csv", header = T)
OLI$spectral  = "OLI"


spectral_H_df = rbind(HS, ENMAP, MSI, OLI)

  df_means <- plyr::ddply(spectral_H_df, "spectral", 
                          dplyr::summarise, mean_rrs = mean(H_sd, na.rm = T))
  
  g_box <- ggplot(data = spectral_H_df, aes(x= as.character(spectral)
                                              #,reorder((spectral), spectral)
                                            )
  )+
    geom_boxplot(aes(y = H_sd, fill = spectral),color = "black", alpha = 0.6)+
    geom_point(aes(y = H_CI, shape = spectral, fill = spectral, color = spectral), size = 1.5, 
               alpha = 0.6)+
    scale_x_discrete(name = " ",
                     labels =  c("HOCR", "ENMAP", "MSI","OLI")) +
    
    geom_line(data = df_means, aes(x=as.character(spectral), y=as.numeric(mean_rrs),
                                   group=1
    ), linetype = "dashed",
    size=1.3)+
    
    scale_color_viridis(discrete = T, guide = "none")+
    scale_fill_viridis(discrete = T, guide = "none")+
    scale_shape_manual(name = "", values = c(21,22,23,19), guide = "none")+
    ylim(0,10)+
    ylab(expression(paste("Standard Deviation (", sigma, ")")))+
    theme_bw()+
    ggtitle(expression(paste("[",italic("H"), "]")))+
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



g_box

ggsave(paste0("./outputs/sensitivity_box_plot_H_insitu.png"), plot = g_box,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

# Create Ternary plots for the aerial fraction of Rb ----
lat_lon = read.csv("./outputs/shallow_lat_lon-BG.csv", header = T)

if (year == 2022) {
  lat_lon = lat_lon[BG_idx,]
} else {
  lat_lon = lat_lon[-(317:341),]
}

shallow_param_df = read.csv("./outputs/shallow_inv_param_BG.csv", header = T)

shallow_param_df = read.csv("./outputs/shallow_inv_param_BG_2023_inel-cor_OFF_HS.csv", header = T)

shallow_param_df = read.csv("./outputs/shallow_inv_param_BG_2022_inel-cor_OFF_HS_grad.csv",
                            header = T)
shallow_param_df = shallow_param_df[BG_HS_idx$...1,]

shallow_param_df$fa1_scale = shallow_param_df$fa1/max(shallow_param_df$fa1) #-  
  # rnorm.trunc(n = nrow(shallow_param_df), mean = 0.5, sd = 2,
  #                                   min = 0, max = 0.2)

shallow_param_df$fa1_scale[shallow_param_df$fa1_scale > 1] = 0.999

shallow_param_df$fa2_scale = shallow_param_df$fa2/max(shallow_param_df$fa2)

shallow_param_df$fa3_scale = shallow_param_df$fa3/max(shallow_param_df$fa3) #+ 
  # rnorm.trunc(n = nrow(shallow_param_df), mean = 0.1, sd = 2,
  #             min = 0, max = 0.2)
shallow_param_df$fa3_scale[shallow_param_df$fa3_scale > 1] = 0.999

# shallow_param_df$fa1_scale = runif(n = length(shallow_param_df$fa1_scale),
#                                     min = 0.3, max = 1)
# shallow_param_df$fa2_scale = runif(n = length(shallow_param_df$fa1_scale),
#                                    min = 0.05, max = 0.3)
# shallow_param_df$fa3_scale = runif(n = length(shallow_param_df$fa1_scale),
#                                    min = 0.1, max = 0.4)
# 
# # Adjust fa1_scale_scale and fa2_scale_scale generation to ensure sum of fa1_scale_scale and fa2_scale_scale does not exceed 1
# shallow_param_df$fa1_scale[shallow_param_df$fa1_scale + shallow_param_df$fa2_scale_scale > 1] <- 1 - shallow_param_df$fa2_scale[shallow_param_df$fa1_scale + shallow_param_df$fa2_scale > 1]
# shallow_param_df$fa2_scale[shallow_param_df$fa1_scale + shallow_param_df$fa2_scale > 1] <- 1 - shallow_param_df$fa1_scale[shallow_param_df$fa1_scale + shallow_param_df$fa2_scale > 1]
# 
# # Generate fa3 based on remaining fraction
# shallow_param_df$fa3_scale <- 1 - shallow_param_df$fa1_scale - shallow_param_df$fa2_scale

#rnorm.trunc(n = 500, mean = 5, sd = 4, min = 0.5, max = 12)

shallow_param_df$veg_volume = H_df_BG$veg_volume
shallow_param_df$veg_frac = H_df_BG$bottom_cover


summary(shallow_param_df$fa1_scale)
summary(shallow_param_df$fa2_scale)
summary(shallow_param_df$fa3_scale)

shallow_param_df = cbind(lat_lon[BG_HS_idx$...1,],shallow_param_df)
# 
# shallow_param_df$veg_volume[(shallow_param_df$fa3_scale > 0.15 | shallow_param_df$fa1_scale > 0.35) & 
#                               (shallow_param_df$fa2_scale < 0.01)] = 
#   #shallow_param_df$veg_volume[shallow_param_df$fa3_scale > 0.15 | shallow_param_df$fa1_scale > 0.35] + 
#   rnorm.trunc(n = nrow(shallow_param_df[(shallow_param_df$fa3_scale > 0.15 | shallow_param_df$fa1_scale > 0.35) & 
#                                           (shallow_param_df$fa2_scale < 0.01),]), 
#               mean = 20, sd = 5, 
#               min = 40, max = 70)
# 
# shallow_param_df$veg_volume[(shallow_param_df$fa3_scale > 0.2 | shallow_param_df$fa1_scale > 0.6) & 
#                               (shallow_param_df$fa2_scale < 0.01)] = 
#   #shallow_param_df$veg_volume[shallow_param_df$fa3_scale > 0.15 | shallow_param_df$fa1_scale > 0.6] + 
#   rnorm.trunc(n = nrow(shallow_param_df[(shallow_param_df$fa3_scale > 0.2 | shallow_param_df$fa1_scale > 0.6) & 
#                                           (shallow_param_df$fa2_scale < 0.01),]),
#               mean = 75, sd = 10, 
#               min = 70, max = 99)
# 
# shallow_param_df$veg_volume[shallow_param_df$veg_volume > 100] = 99

# shallow_param_df$veg_volume[shallow_param_df$fa3_scale > 0.05] = 
#   shallow_param_df$veg_volume[shallow_param_df$fa3_scale > 0.05] + 
#   rnorm.trunc(n = nrow(shallow_param_df[shallow_param_df$fa3_scale > 0.05,]), 
#               mean = 20, sd = 5, 
#               min = 10, max = 40)
# shallow_param_df$fa3_scale[shallow_param_df$fa3_scale > 1] = 0.999


library(ggtern)

g_tern <- ggtern(data=shallow_param_df#[HS_idx_BG,]
                 ,
                 aes(x=fa1_scale,y=fa2_scale,
                                             z=fa3_scale),aes(x,y,z)) + 
  
  
  geom_point(aes(fill=veg_volume),size=4,shape=21, alpha = 0.7) + 
  
  # stat_density_tern(geom="polygon",
  #                   aes(fill=..level..,
  #                       #weight=veg_volume,
  #                       #alpha=abs(..level..)
  #                   ),
  #                   alpha = 0.3,
  #                   na.rm = TRUE,
  #                   bdl = 0.02) + 
  geom_density_tern(aes(#weight=veg_volume,
    color=..level..),
    #n=200,
    #                 binwidth=100,
    bdl = 0.02) +
  
  scale_fill_gradient(name = "Plant volume (%)", low = "#440154FF", high = "#9FDA3AFF", 
                      na.value = "#20A387FF",
                      
                      breaks=seq(0,100,25),
                      limits=c(0,100), labels=paste(seq(0,100,25))
                      
                      # breaks=seq(0,5,2),
                      # limits=c(0,5), labels=paste(seq(0,5,2),"(m)")
                      )+
  scale_color_gradient(name = "Plant volume (%)", low = "#440154FF", high = "#9FDA3AFF", 
                       na.value = "#20A387FF",
                       
                       breaks=seq(0,100,25),
                       limits=c(0,100), labels=paste(seq(0,100,25))
                       
                       # breaks=seq(0,5,2),
                       # limits=c(0,5), labels=paste(seq(0,5,2),"(m)")
  )+
  
  #coord_fixed()+
  
  theme_rgbw() + 
  theme(
    axis.text.x = element_text(size = 15, color = 'black', angle = 0),
    axis.text.y = element_text(size = 15, color = 'black', angle = 0),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.ticks.length = unit(.25, "cm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.key.width = unit(2, "cm"),  # Adjust the width of the legend key
    legend.key.height = unit(0.5, "cm")  # Adjust the height of the legend key
  ) +
  
  guides(fill = guide_colorbar(order=1),
         alpha= guide_legend(order=2),
         color="none") + 
  labs(fill = "Vegetation Volume", x = "Lithom.", y = "Rock", "z" = "Sacha")

g_tern

ggsave(paste0("./outputs/rb_ternary_",year,".png"), plot = g_tern,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)



# Create a spatial plot for the SABER retrieved benthic reflectance ----
library(ggplot2)
library(sf)
library(raster)
library(viridis)

shallow_param_df = shallow_param_df[-length(shallow_param_df$lat),] 

set.seed(100)

rand <- shallow_param_df[,19:21]

rand <- shallow_param_df[,21:23]

comb_col = (rand$fa1_scale + rand$fa2_scale + rand$fa3_scale)/3

rand <- rand / apply(rand, 1, sum) # Make sure the numbers sum to one in each row

rand <- cbind(shallow_param_df[,1:2], rand)
names(rand) = c("lat", "lon", "fa1_scale", "fa2_scale", "fa3_scale")


# Create an RGB color based on fa1_scale, fa2_scale, and fa3_scale
rand <- rand %>%
  mutate(color = rgb(fa1_scale, fa2_scale, fa3_scale))

# Convert RGB to Viridis
rand$viridis_color <- viridis_pal()(100)[as.numeric(cut(col2rgb(rand$color)[1,], breaks=100))]

get_viridis_color <- function(x, opt="D", begin=0, end=1, reverse=FALSE) {
  x <- x * (end - begin) + begin
  cmap <- viridisLite::viridis.map[viridisLite::viridis.map$opt == opt,]
  if (reverse) cmap <- cmap[rev(seq_len(nrow(cmap))),]
  map_rgbs <- grDevices::rgb(cmap$R, cmap$G, cmap$B)
  ramp <- grDevices::colorRamp(map_rgbs, space="Lab", interpolate="spline")
  out_rgbs <- ramp(x) / 255
  grDevices::rgb(out_rgbs[,1], out_rgbs[,2], out_rgbs[,3])
}

rand$color_v = get_viridis_color(x = rand$fa1_scale)


# Plot using ggplot2
ggplot(rand[rand$lon > -68,]) +
  geom_tile(aes(x = lon, y = lat, fill = viridis_color),width=0.001,height=0.001, alpha = 0.6) +
  #geom_point(aes(x = lon, y = lat, fill = color, color = color), alpha = 0.6, show.legend = T) +
  
  scale_fill_identity() +
  #scale_color_identity()+
  labs(title = "Spatial Plot of FA Scales",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Convert data frame to sf object for each FA scale
fa1_sf <- st_as_sf(rand, coords = c("lon", "lat"), crs = 4326)
fa2_sf <- fa1_sf
fa3_sf <- fa1_sf

# Assign the correct data for each FA scale
fa1_sf <- fa1_sf %>% select(fa1_scale)
fa2_sf <- fa2_sf %>% select(fa2_scale)
fa3_sf <- fa3_sf %>% select(fa3_scale)

# Write each sf object to a shapefile
st_write(fa1_sf, "fa1_scale.shp")
st_write(fa2_sf, "fa2_scale.shp")
st_write(fa3_sf, "fa3_scale.shp")

d <- data.frame(a = rand[,3],
                b = rand[,4],
                c = rand[,5],
                #expand.grid(x=shallow_param_df$lat, y=shallow_param_df$lon)
                x = rand$lat, y = rand$lon
                )

g_spat <- ggplot(d[d$y > -64,], aes(x=x, y=y))+
  theme_bw() + #For a clearer background
  
  # geom_sf(data = world, fill = "antiquewhite") +
  # coord_sf(xlim = c(48, 50), ylim = c(-70, -65), expand = FALSE) +
  
  geom_tile(alpha=d$a[d$y > -64], fill="brown",width=0.001,height=0.001, show.legend = T)+
  geom_tile(alpha=d$b[d$y > -64], fill="goldenrod",width=0.001,height=0.001, show.legend = T)+
  geom_tile(alpha=d$c[d$y > -64], fill="seagreen",width=0.001,height=0.001, show.legend = T)+
  theme(panel.grid.major=element_blank())+
  labs(title = "",
       x = "Latitude(o)", y = "Longtitude(o)") 


g_spat <- ggplot(d[d$y > -68,], aes(x=x, y=y))+
  theme_bw() + #For a clearer background
  
  # geom_sf(data = world, fill = "antiquewhite") +
  # coord_sf(xlim = c(48, 50), ylim = c(-70, -65), expand = FALSE) +
  
  geom_tile(alpha=d$a[d$y > -68], fill="brown",width=0.001,height=0.001, show.legend = T)+
  geom_tile(alpha=d$b[d$y > -68], fill="goldenrod",width=0.001,height=0.001, show.legend = T)+
  geom_tile(alpha=d$c[d$y > -68], fill="seagreen",width=0.001,height=0.001, show.legend = T)+
  theme(panel.grid.major=element_blank())+
  labs(title = "",
       x = "Latitude(o)", y = "Longtitude(o)") 
g_spat
ggsave(paste0("./outputs/rb_spatial-",site_input,"_",year,".png"), plot = g_spat,
       scale = 1.5, width = 6.5, height = 4.5, units = "in",dpi = 300)


# Load necessary libraries
library(ggplot2)
library(raster)
library(viridis)
library(dplyr)

# Create SpatialPointsDataFrame
coordinates(shallow_param_df) <- ~lon + lat

# Define a raster template
r <- raster(extent(shallow_param_df), res = 0.0001)

# Create raster layers for each fa scale
r_fa1 <- rasterize(shallow_param_df, r, field = "fa1_scale", fun = mean, na.rm = TRUE)
r_fa2 <- rasterize(shallow_param_df, r, field = "fa2_scale", fun = mean, na.rm = TRUE)
r_fa3 <- rasterize(shallow_param_df, r, field = "fa3_scale", fun = mean, na.rm = TRUE)

# Stack the rasters
r_stack <- stack(r_fa1, r_fa2, r_fa3)
names(r_stack) <- c("fa1_scale", "fa2_scale", "fa3_scale")

# Convert raster stack to data frame for ggplot
r_df <- as.data.frame(r_stack, xy = TRUE, na.rm = TRUE)

# Create a composite color based on the three FA scales
r_df$composite_color <- with(r_df, fa1_scale + fa2_scale + fa3_scale)

# Plotting using ggplot2 with composite color
ggplot(d) +
  geom_tile(aes(x = x, y = y, fill = b), alpha = 0.8) +
  scale_fill_viridis(name = "Composite FA Scale", option = "D") +
  labs(title = "Spatial Plot of FA Scales",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")









































# seadoo_rrs_plot_df = as.data.frame((seaDoo_rrs_interp))
# 
# seadoo_rrs_plot_df$id = seq(1, length(seaDoo_rrs_interp$`401`), 1)
# 
# seadoo_rrs_long = reshape2::melt(seadoo_rrs_plot_df, id.vars = "id")
# 
# H_plot_df = data.frame("id" = seq(1,length(seadoo_depth), 1), "Depth" = seadoo_depth)
# H_long = reshape2::melt(H_plot_df, id.vars = "id")
# 
# seadoo_rrs_H = merge(seadoo_rrs_long,H_long,by=c("id"))
# 
# seadoo_rrs_H = seadoo_rrs_H[,-(4)] 
# names(seadoo_rrs_H) = c("id", "wavelength", "rrs", "Depth")
# 
# 
# xmin <- 400; xmax <- 700;  xstp <- 50; xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
# ymin <- 0; ymax <- 0.015; ystp <- 0.003; ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda),italic("in situ ("), "sr"^{-1}, ")"))
# asp_rat <-  (xmax-xmin)/(ymax-ymin)
# legend_title <- element_blank()
# legend_position <- c(0.70, 0.98)
# 
# 
# g <- ggplot(data = seadoo_rrs_H) + 
#   geom_line(aes(x = as.numeric(as.character(wavelength)), y = rrs, colour = Depth, 
#                 group = Depth), size = 1.3, show.legend = T)+
#   scale_colour_viridis(discrete = F, name = "Depth (m)") +
#   coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
#               ylim = c(ymin, ymax)
#               ,expand = FALSE, clip = "on"
#   ) +
#   #scale_colour_manual(name = "",values = rev(collist))+
#   #guides(fill = FALSE)+
#   guides(fill = guide_legend(title = "Depth (m)"))+
#   scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
#                      breaks = seq(xmin, xmax, xstp))  +
#   #scale_y_log10()+
#   scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
#                      breaks = seq(ymin, ymax, ystp))  +
#   theme_bw()+
#   theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
#         axis.text.x = element_text(size = 25, color = 'black', angle = 0), 
#         axis.text.y = element_text(size = 25, color = 'black', angle = 0), 
#         axis.title.x = element_text(size = 25),
#         axis.title.y = element_text(size = 25),
#         axis.ticks.length = unit(.25, "cm"),
#         legend.box.just = "right",
#         legend.spacing = unit(-0.5, "cm"),
#         legend.position = legend_position,
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
#         plot.margin = unit(c(0.0,0.9,0.0,0.0), "cm"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
# 
# g
# ggsave("./outputs/seadoo_rrs.png", plot = g,scale = 1.7, width = 4.5, height = 4.5, 
#        units = "in",dpi = 300)