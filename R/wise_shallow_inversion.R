# Read, restructure, interpolate Rrs data ----

## WISE-Man 2019 COPS data ----
wise_cops = read.csv("./data/Rb_spectral/copsdata.csv", header = T)
wise_cops_wide = read.csv("./data/Rb_spectral/COPS_db.csv", header = T)

kildir_iop = read.csv("./data/Kildir_IOP.Process_Log.csv", header = T)
saucier_iop = read.csv("./data/saucierQC.csv", header = T)

#Input the depth statistics and simulation notes from WISE-Man IOP data
wise_depth_stat = read.csv("./data/Depth_wise.csv", header = T)
wise_depth_stat$correct_statname = gsub(wise_depth_stat$Station, 
                                        pattern = "\\_", replacement = "-") 
wise_depth_stat <- wise_depth_stat %>% dplyr::select(correct_statname, everything())

#Select the QC passed stations for WISE-Man
wise_depth_stat_shallow = wise_depth_stat[wise_depth_stat$Shallow == "T",]

kildir_iop_qc = kildir_iop[kildir_iop$ASPH == "Y" & kildir_iop$HS6 == "Y",]
saucier_iop_qc = saucier_iop[saucier_iop$station.keep == "TRUE",]

# qc_statnames_bind = rbind(list(kildir_iop_qc$StationID), 
#                           list(saucier_iop_qc$Station))

qc_statnames_bind = rbind(list(kildir_iop$StationID), 
                          list(saucier_iop$Station))

qc_statnames_bind = data.frame("qcstat" = as.character(unlist(qc_statnames_bind)))


shallow_qc_keep_idx = vector()

for (i in 1:length(qc_statnames_bind$qcstat)) {
  tryCatch({
    shallow_qc_keep_idx[i] = grep(x =  wise_depth_stat_shallow$correct_statname, 
                                  pattern = qc_statnames_bind$qcstat[i])
    
    if (is.na(shallow_qc_keep_idx[i])) {
      cat(paste0("\033[0;31m","The station ", wise_depth_stat_shallow$correct_statname[i]," is not shallow, coerce to NA","\033[0m","\n"))
      shallow_qc_keep_idx[i] = NA
    } }, 
    error = function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

shallow_qc_stations_grep = qc_statnames_bind$qcstat[!is.na(shallow_qc_keep_idx)]
shallow_qc_stations_in = qc_statnames_bind[qc_statnames_bind$qcstat %in% 
                                             wise_depth_stat_shallow$correct_statname,]
#shallow_qc_stations = union(shallow_qc_stations_grep, shallow_qc_stations_in)

shallow_qc_stations =  shallow_qc_stations_in

shallow_qc_stations = gsub(shallow_qc_stations, 
                                        pattern = "\\-", replacement = ".")

#shallow_qc_stations = shallow_qc_stations[-which(shallow_qc_stations == "MAN.F0")]

wise_cops_shallow =  wise_cops[,colnames(wise_cops) %in%  
                                 # wise_depth_stat_shallow$correct_statname
                                 shallow_qc_stations
                               ] 

wise_cops_wide_shallow = wise_cops_wide[,colnames(wise_cops_wide) %in% 
                                          colnames(wise_cops_shallow)]

rownames(wise_cops_wide_shallow) = wise_cops_wide[,1]

wise_cops_shallow = wise_cops_shallow[(1:351),]

rm(qc_statnames_bind)

## WISE-Man 2019 PSR data ----
wise_psr = read.csv("./data/Rrs_psr_long_wiseman.csv", header = T)

wise_psr = split(wise_psr, wise_psr$station)

## PSR multicast data ----
multicast_idx <- sapply(wise_psr, function(df) {
  nrow(df) > 781 && (nrow(df) %% 781 == 0)
})

multicast_idx <- which(multicast_idx)

wise_psr_multicast = wise_psr[multicast_idx]

wise_psr_multicast_df <- lapply(wise_psr_multicast, function(df_observation) {
  multicast_df_split <- split(df_observation, df_observation$station_alt)
  
  lapply(multicast_df_split, function(df_obs){
    
    wavelength <- df_obs$Wavelength
    Rrs <- df_obs$Rrs
    station = df_obs$station_alt
    df_Rrs <- data.frame(Rrs)
    #do.call(rbind, multicast_df_split)
  })
  
})

extract_rrs <- function(data_list) {
  
  all_rrs <- data.frame()
  

  for(df in data_list) {
   
    df_d = data.frame(t(do.call(cbind, df)))
    rownames(df_d) = names(df)
    colnames(df_d) = unique(wise_psr_multicast$`EXP-WISE`$wavelength)
    all_rrs = rbind(all_rrs, df_d)
    
  }
  
  return(all_rrs)
}

wise_psr_multicast_df = extract_rrs(wise_psr_multicast_df)

## PSR unicast data ----
wise_psr_unicast_df <- lapply(wise_psr[-multicast_idx], function(df_observation) {
  wavelength <- df_observation$Wavelength
  Rrs <- df_observation$Rrs
  station = df_observation$station
  df_Rrs <- data.frame(Rrs)
})

## combined PSR data ----
wise_psr_unicast_df = data.frame(t(do.call(cbind, wise_psr_unicast_df)))

rownames(wise_psr_unicast_df)  = names(wise_psr[-multicast_idx])
colnames(wise_psr_unicast_df) = wise_psr[[1]]$wavelength


wise_psr_rrs_df = rbind(wise_psr_unicast_df, wise_psr_multicast_df)

wise_psr_rrs_df = wise_psr_rrs_df[colnames(wise_psr_rrs_df) %in% seq(400,750,1)]
wise_psr_rrs_df = data.frame(t(wise_psr_rrs_df))


# 
# wise_psr_shallow =  wise_psr_rrs_df[,colnames(wise_psr_rrs_df) %in%  
#                                               # wise_depth_stat_shallow$correct_statname
#                                       bgc_wiseman$correct_statname] 
# 
wise_psr_shallow =  wise_psr_rrs_df[, 1:47][,colnames(wise_psr_rrs_df)[1:47] %in%
                                 # wise_depth_stat_shallow$correct_statname
                                 shallow_qc_stations]

wise_psr_shallow = cbind(wise_psr_shallow, wise_psr_rrs_df[, -(1:47)])

rm(wise_psr_multicast_df, wise_psr_unicast_df)

#Retrieve PSR station depths
bgc_wiseman = read.csv("./data/biogeochemistry_wiseman.csv", header = T)
bgc_wiseman$correct_statname = gsub(bgc_wiseman$station,
                                    pattern = "\\-", replacement = ".")

bgc_wiseman_shallow = bgc_wiseman[bgc_wiseman$correct_statname %in% colnames(wise_psr_shallow),]
bgc_wiseman_shallow = rbind(bgc_wiseman_shallow, bgc_wiseman[nrow(bgc_wiseman),])
bgc_wiseman_shallow = bgc_wiseman_shallow[bgc_wiseman_shallow$depth < 1, ]

wise_psr_shallow_known_H = wise_psr_shallow[colnames(wise_psr_shallow) %in% 
                                              bgc_wiseman_shallow$correct_statname]

wise_psr_shallow_known_H = cbind(wise_psr_shallow_known_H, 
                                 rowMeans(wise_psr_rrs_df[,48:53]))

colnames(wise_psr_shallow_known_H)[10] = "EXP.WISE"

wise_psr_shallow_known_H = wise_psr_shallow_known_H[,order(colnames(wise_psr_shallow_known_H))]
bgc_wiseman_shallow = bgc_wiseman_shallow[order(bgc_wiseman_shallow$correct_statname),]



#Apply sub-surface translation
wise_cops_shallow = as.data.frame(apply(wise_cops_shallow, 2, 
                                        surface_rrs_translate))

wise_psr_shallow = as.data.frame(apply(#wise_psr_shallow, 
                                      wise_psr_shallow_known_H,
                                       2, 
                                        surface_rrs_translate))

## Exclude the inelastic scattering contribution from Rrs ----
subtract_wise_inel <- function(input_df){
  
  inel_contrib = read.csv("./data/inel_contrib_MP.csv", header = T)
  
  
  output_df <- apply(input_df, 2, 
                                          function(row, 
                                                   wavelength_seadoo = as.numeric(wise_cops$wave)[1:351]){
                                            
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
  
  output_df = data.frame((output_df))
  return(output_df)
  
}

### COPS ----
wise_cops_shallow_incor = subtract_wise_inel(input_df = wise_cops_shallow)

### PSR ----
wise_psr_shallow_incor = subtract_wise_inel(input_df = wise_psr_shallow)


wise_rrs_df = cbind(wise_cops_shallow, wise_psr_shallow)

wise_rrs_incor_df = cbind(wise_cops_shallow_incor, wise_psr_shallow_incor)


# Plot the mean spectra for inelastic scattering corrected and uncorrected Rrs respectively

mean_rrs_inel_el = data.frame("wavelength" = seq(400,750,1), 
                              "rrs_inel_cor" = rowMeans(wise_psr_shallow_incor),
                              "sd_inel_cor" = as.numeric(apply(wise_psr_shallow_incor, 1,sd)),
                              "rrs_inel_uncor" = rowMeans(wise_psr_shallow),
                              "sd_inel_uncor" = as.numeric(apply(wise_psr_shallow, 1, sd))
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
ggsave(paste0("./outputs/wise_el_inel_comp_PSR.png"), plot = g,
       scale = 1.25, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)

inel_contrib = data.frame("wavelength" = seq(400,750,1), "inel_ratio" = inel_contrib)
inel_contrib$inel_ratio = 100*inel_contrib$inel_ratio


legend_title <- element_blank()
legend_position <- c(0.10, 0.95)
show_legend = T

xlbl <- expression(paste("Wavelength (", lambda, ") [nm]"))
ylbl <- expression(paste(italic("R")["rs,inel"]("0"^"-", lambda), "[%]"))

ymin <- 0; ymax <- 100 ; ystp <- ymax/5
xmin <- 400; xmax <- 750; xstp <- 70
asp_rat <- (xmax-xmin)/(ymax-ymin)

g_inel = ggplot(data = inel_contrib, aes(x = wavelength, y = inel_ratio))+
  geom_line(size = 1.3, color = "purple3")+
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
        legend.position = c(0.01,1),
        #legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 10, face = "plain"),
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
        legend.direction = "horizontal", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
ggsave(paste0("./outputs/wise_2019_inel_contrib.png"), plot = g_inel,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)



## Retrieve BGC parameters for given stations ----
### COPS ----
bgc_val_vec = data.frame()
shallow_statist = gsub(colnames(wise_cops_shallow_incor), 
                       pattern = "\\.", replacement = "-")
shallow_statist[c(which(shallow_statist == "MAN-F03-5"),
                  which(shallow_statist == "MAN-R11-5"))] = c("MAN-F03.5", "MAN-R11.5")

for (i in 1:length(shallow_statist)) {
  
  temp_insitu_param = get_in_situ_params(station_name = shallow_statist[i])
  temp_bgc_val = data.frame("chl" = temp_insitu_param$chl_invivo,
                            "adg443" = temp_insitu_param[[1]]$a_cdom + 
                              temp_insitu_param[[1]]$a_nap,
                            "bbp555" = temp_insitu_param$bbp555)
  bgc_val_vec = rbind(bgc_val_vec, temp_bgc_val)
  
}

### PSR ----

# Retrieve the stations common between COPS and PSR for comparative study----
common_columns <- intersect(colnames(wise_cops_shallow_incor), colnames(wise_psr_shallow_incor))

wise_cops_shallow_common_incor <- wise_cops_shallow_incor[colnames(wise_cops_shallow_incor) 
                                              %in% common_columns]
wise_cops_wide_shallow_common_incor <- wise_cops_wide_shallow[colnames(wise_cops_wide_shallow) 
                                                        %in% common_columns]
wise_psr_shallow_common_incor <- wise_psr_shallow_incor[colnames(wise_psr_shallow_incor) %in% 
                                                    common_columns]

## Plot the COPS and PSR common Rrs together ----
color_var <- seq(400, 750, length.out = nrow(wise_cops_shallow_common_incor)) 

wise_cops_shallow_common_incor_long <- wise_cops_shallow_common_incor %>%
  mutate(Wavelength = color_var) %>%
  pivot_longer(cols = -Wavelength, names_to = "variable", values_to = "COPS")

wise_psr_shallow_common_incor_long <- wise_psr_shallow_common_incor %>%
  mutate(Wavelength = color_var) %>%
  pivot_longer(cols = -Wavelength, names_to = "variable", values_to = "PSR")

combined_df <- wise_cops_shallow_common_incor_long %>%
  left_join(wise_psr_shallow_common_incor_long, by = c("Wavelength", "variable"))

combined_df_long <- combined_df %>%
  pivot_longer(cols = c(COPS, 
                        PSR), names_to = "source", values_to = "value")

rm(wise_cops_shallow_common_incor_long, wise_psr_shallow_common_incor_long,
   combined_df)

g_rrs = ggplot(combined_df_long, aes(x = Wavelength, y = value, color = source)) +
  geom_line(show.legend = T, size = 1.3) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_y_continuous(limits = c(0, 0.010), breaks = seq(0, 0.010, by = 0.002)) +
  scale_x_continuous(limits = c(400, 750), breaks = seq(400, 750, by = 70)) +
  labs(x = "Wavelength (nm)", y = expression(paste(italic("R")["rs"]("0"^"+", lambda),italic("in situ ["), "sr"^{-1}, "]")),
       color = " ") +
  scale_color_manual(values = c("purple3", "green3"), labels = c("COPS", "PSR"))+
  theme_bw()+
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0),
        axis.text.y = element_text(size = 25, color = 'black', angle = 0),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.42,0.9),
        legend.direction = "horizontal",
        legend.title = element_text(colour = "black", size = 15, face = "plain"),
        legend.text = element_text(colour = "black", size = 18, face = "bold"),
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

ggsave(paste0("./outputs/rrs_shallow_cops-psr_comp.png"), plot = g_rrs, scale = 1.7, 
       width = 8, height = 4.5,
       units = "in",dpi = 300)

## Prepare the optimization parameters ----

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
                              init_par = c(colMeans(bgc_val_vec),
                                           mean(as.numeric(wise_cops_wide_shallow[4,]), 
                                                na.rm = T),
                                           0.5,0.5,0.5,0.1
                                           ),
                              
                              upper_par = c( apply(bgc_val_vec, 2, max),
                                            max(as.numeric(wise_cops_wide_shallow[4,]), 
                                                 na.rm = T),
                                            1,1,1,10
                              ),
                              
                              lower_par = c(apply(bgc_val_vec, 2, min),
                                            min(as.numeric(wise_cops_wide_shallow[4,]), 
                                                 na.rm = T),
                                            0,0,0,0.001
                              )
)

par0 = inv_bound$par0 ; upper.bound = inv_bound$upper.bound  
lower.bound = inv_bound$lower.bound

# Selection of inversion optimizer and objective function
obj = c("log-LL", "SSR", "obj_L98"); obj.run = obj[3]

methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                 "Brent","levenberg-marqardt", "auglag")

samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")


inv_res_df_bayes = data.frame(#"chl" = 0, "adg443" = 0, "bbp555" = 0,
                           "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
                           #"sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
                           "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
                           "sd_pop_sd" = 0)
   
inv_res_df_bayes_unconst = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
     "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
     "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
     "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
     "sd_pop_sd" = 0)
   
inv_res_df_grad_unconst = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
     "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
     "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
     "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
     "sd_pop_sd" = 0)

# inv_res_df_bayes_cops = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
#   "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
#   "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
#   "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
#   "sd_pop_sd" = 0)
# 
# inv_res_df_bayes_psr = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
#   "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
#   "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
#   "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
#   "sd_pop_sd" = 0)


inverse_rrs_input = as.data.frame(cbind(wise_cops_shallow_incor, wise_psr_shallow_incor))

for(j in 1:ncol(inverse_rrs_input)) {
  
  # param_grad = doOptimization_shallow_BGC_const(obsdata = 
  #                                                 as.numeric(inverse_rrs_input[,j]),
  #                                               init_par = par0[-(1:3)],
  #                                               max_par = upper.bound[-(1:3)],
  #                                               min_par = lower.bound[-(1:3)],
  #                                               QAA_mode = F, bgc_const_val = as.numeric(bgc_val_vec[j,]),
  #                                               wavelength_sim = seq(400,750,1) #as.numeric(wise_cops$wave[1:351])
  #                                               , 
  #                                               qaa_prefit = F,
  #                                               qaa_slope = F, manual_slope = F, sa_model = "am03",
  #                                               obj_fn = obj[1], opt_method = methods.opt[4] )
  # 
  # param_bayes = inverse_runBayes(obsdata = 
  #                                         as.numeric(inverse_rrs_input[,j]),
  #                                       rrs_type =  type_Rrs_below,
  #                                       max_par = upper.bound,
  #                                       min_par = lower.bound,
  #                                       param_lab = c(#"chl", "adg443", "bbp555", 
  #                                                     "H", "fa1", "fa2", "fa3", 
  #                                                     "pop_sd"),
  #                                       
  #                                       constrain_config = c(F,T,F),
  #                                       bbp_const_val = 0.007,
  #                                       bgc_const_val = as.numeric(bgc_val_vec[j,]),
  #                                       
  #             iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
  #                                "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
  #                                       
  #                                       qaa_slope = F,
  #                                       manual_slope = T,
  #                                       manual_slope_vals  = c("s_g"=0.014,
  #                                                              "s_d"=0.003,
  #                                                              "gamma"=0.5),
  #                                       iter_count = 20000,
  #                                       sampler_mcmc = samplerlist[6],
  #                                       wavelngth_sim = seq(400,750,1) #as.numeric(wise_cops$wave[1:351])
  #             , 
  #                                       sa_model = "am03",
  #                                       hybrid_mode = F, plot_rrs = F)
  
  
  param_grad_unconst = doOptimization_shallow_unconst(obsdata =
                                           as.numeric(inverse_rrs_input[,j]),
                                           init_par = c(3.5,1.3,0.017,4.5,0.6,0.2,
                                                        0.2,0.05),
                                           max_par = c(6.5, 2.5, 0.01, 8, 1, 1, 
                                                       1, 10),
                                           min_par = c(1.5, 1, 0.003, 2, 0.1, 0.1, 
                                                       0.1, 0.01),
                                                QAA_mode = F, 
                                                bgc_const_val = as.numeric(bgc_val_vec[j,]),
                                                wavelength_sim = seq(400,750,1),
                                                qaa_prefit = F,
                                                qaa_slope = F, manual_slope = T, 
                                                manual_slope_vals = c("s_g"=0.014, "s_d"=0.003, 
                                                                      "gamma"=0.5),
                                                sa_model = "am03",
                                                obj_fn = obj[1], opt_method = methods.opt[4] )

  param_bayes_unconst = inverse_runBayes(obsdata = 
                                   as.numeric(inverse_rrs_input[,j]),
                                 rrs_type =  type_Rrs_below,
                                 max_par = c(6.5, 2.5, 0.01, 8, 1, 1, 1, 10),
                                 min_par = c(1.5, 1, 0.003, 2, 0.1, 0.1, 0.1, 0.01),
                                 param_lab = c("chl", "adg443", "bbp555", 
                                   "H", "fa1", "fa2", "fa3", 
                                   "pop_sd"),
                                 
                                 constrain_config = c(F,F,F),
                                 bbp_const_val = 0.007,
                                 bgc_const_val = as.numeric(bgc_val_vec[j,]),
                                 
                                 iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                                                    "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
                                 
                                 qaa_slope = F,
                                 manual_slope = T,
                                 manual_slope_vals  = c("s_g"=0.014,
                                                        "s_d"=0.003,
                                                        "gamma"=0.5),
                                 iter_count = 25000,
                                 sampler_mcmc = samplerlist[6],
                                 wavelngth_sim = seq(400,750,1), 
                                 sa_model = "am03",
                                 hybrid_mode = F, plot_rrs = T)
  
  # param_bayes_unconst_psr = inverse_runBayes(obsdata = 
  #                                               as.numeric(wise_psr_shallow_common_incor[,j]),
  #                                             rrs_type =  type_Rrs_below,
  #                                             max_par = upper.bound,
  #                                             min_par = lower.bound,
  #                                             param_lab = c("chl", "adg443", "bbp555", 
  #                                                           "H", "fa1", "fa2", "fa3", 
  #                                                           "pop_sd"),
  #                                             
  #                                             constrain_config = c(F,F,F),
  #                                             bbp_const_val = 0.007,
  #                                             bgc_const_val = as.numeric(bgc_val_vec[j,]),
  #                                             
  #                                             iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
  #                                                                "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
  #                                             
  #                                             qaa_slope = F,
  #                                             manual_slope = T,
  #                                             manual_slope_vals  = c("s_g"=0.014,
  #                                                                    "s_d"=0.003,
  #                                                                    "gamma"=0.5),
  #                                             iter_count = 20000,
  #                                             sampler_mcmc = samplerlist[6],
  #                                             wavelngth_sim = seq(400,750,1) #as.numeric(wise_cops$wave[1:351])
  #                                             , 
  #                                             sa_model = "am03",
  #                                             hybrid_mode = F, plot_rrs = T)
  
  inv_res_df_grad_unconst = rbind(inv_res_df_grad_unconst, param_grad_unconst)
  inv_res_df_bayes_unconst = rbind(inv_res_df_bayes_unconst, param_bayes_unconst)
  #inv_res_df_bayes_unconst = rbind(inv_res_df_bayes_unconst, param_bayes_unconst)
  
  #inv_res_df_bayes_cops = rbind(inv_res_df_bayes_cops, param_bayes_unconst_cops)
  #inv_res_df_bayes_psr = rbind(inv_res_df_bayes_psr, param_bayes_unconst_psr)
  
}

inv_res_df_bayes_unconst = inv_res_df_bayes_unconst[-1,]
inv_res_df_bayes_unconst$sensor = c(rep("COPS", 16), rep("PSR", 10))

inv_res_df_grad_unconst = inv_res_df_grad_unconst[-1,]
inv_res_df_grad_unconst$sensor = c(rep("COPS", 16), rep("PSR", 10))


inv_res_df_bayes_unconst$H_actual = c(as.numeric(wise_cops_wide_shallow[4,]), 
                                      as.numeric(bgc_wiseman_shallow$water_depth))

inv_res_df_grad_unconst$H_actual = c(as.numeric(wise_cops_wide_shallow[4,]), 
                                      as.numeric(bgc_wiseman_shallow$water_depth))

bgc_wiseman = bgc_wiseman[bgc_wiseman$depth == 0 ,]
shallow_chl = data.frame("station" = bgc_wiseman$correct_statname[bgc_wiseman$correct_statname %in% 
                                    inv_res_df_grad_unconst$station][-18],
                         "chl" = bgc_wiseman$Chl[bgc_wiseman$correct_statname %in% 
                                     inv_res_df_grad_unconst$station][-18])
                         
inv_res_df_bayes_unconst = read.csv("./outputs/inv_param_bayes_wise_shallow.csv", header = T)

#Calculate Confidence Interval
cal_sd <- function(i,input_df, col_no ){
  errorb <- qnorm(0.975)*sd(input_df[i,col_no])/
    sqrt(length(input_df[i,col_no]))
}

sd_indices <- 1:dim(inv_res_df_bayes_unconst)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, input_df = inv_res_df_bayes_unconst, 
                        col_no = c(4,18))) 

inv_res_df_bayes_unconst$H_CI = as.numeric(sd_estimate)
inv_res_df_bayes_unconst$station = c(colnames(wise_cops_shallow_incor),
                                    colnames(wise_psr_shallow_incor))

inv_res_df_bayes_unconst$station_type <- ifelse(duplicated(inv_res_df_bayes_unconst$station) | 
                                                 duplicated(inv_res_df_bayes_unconst$station, 
                                                            fromLast = TRUE), 
                                               "duplicate", 
                                               paste(inv_res_df_bayes_unconst$sensor, "unique", sep = "_"))


inv_res_df_bayes_unconst$statname_new <- ifelse(inv_res_df_bayes_unconst$station_type == 
                                                 "duplicate", 
                                                inv_res_df_bayes_unconst$station, NA)

sd_estimate <- t(sapply(sd_indices,  cal_sd, input_df = inv_res_df_grad_unconst, 
                        col_no = c(4,18))) 

inv_res_df_grad_unconst$H_CI = as.numeric(sd_estimate)
inv_res_df_grad_unconst$station = c(colnames(wise_cops_shallow_incor),
                                                           colnames(wise_psr_shallow_incor))


inv_res_df_grad_unconst$station_type <- ifelse(duplicated(inv_res_df_grad_unconst$station) | 
                                                 duplicated(inv_res_df_grad_unconst$station, 
                                                            fromLast = TRUE), 
                          "duplicate", 
                          paste(inv_res_df_grad_unconst$sensor, "unique", sep = "_"))


inv_res_df_grad_unconst$statname_new <- ifelse(inv_res_df_grad_unconst$station_type == 
                                                 "duplicate", 
                                               inv_res_df_grad_unconst$station, NA)

inv_res_df_grad_unconst$chl_actual = inv_res_df_bayes_unconst$chl_actual
# inversion_results = data.frame("station" = c(colnames(wise_cops_shallow_incor),
#                                              colnames(wise_psr_shallow_incor)),
#                                # "chl_actual" = rep(bgc_val_vec$chl,2), 
#                                # "chl_predicted" = inv_res_df_bayes$chl,
#                                # 
#                                # "adg443_actual" = rep(bgc_val_vec$adg443, 2),
#                                # "adg443_predicted" = inv_res_df_bayes$adg443,
#                                # 
#                                # "bbp555_actual" = rep(bgc_val_vec$bbp555,2), 
#                                # "bbp555_predicted" = inv_res_df_bayes$bbp555,
#                                
#                                "H_actual" =  wise_psr,
#                                "H_predicted" = inv_res_df_bayes$H
#                                
#                                )


#inversion_results_merged = cbind(inversion_results,inv_res_df_bayes[5:17])


plot_inversion_validation_multivar_linear_contour <- function(
    
                      input_df = inv_res_df_grad_unconst, 
                      input_x = "H_actual", input_y = "H",
                      uncertainty = "H_CI", xmin = 0, xmax = 12, 
                     xlabel = expression(paste(italic("H"["actual"]), "[m]")),
                     ylabel = expression(paste(italic("H"["predicted"]), "[m]")), 
                     opacity = 0.8, plot_col = "purple4", xstp = 2, ystp = 2, 
                     show_legend = T, hist_count = 30){
    ## data
    d <- input_df %>% 
      as_tibble() 
    
    ymin = xmin; ymax = xmax
    
    asp_rat <- (xmax-xmin)/(ymax-ymin)
    
    g<-   ggplot(data=d, aes(x = .data[[input_x]], y = .data[[input_y]])) +
      
      geom_density_2d(data = d, aes(x = .data[[input_x]], y = .data[[input_y]]),na.rm = T, 
                      bins = 6,
                      linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
      
      geom_point(data = d, aes(.data[[input_x]], .data[[input_y]],
                               shape=sensor, fill=(station_type), color = station_type),
                 alpha = I(opacity), size = I(5), show.legend = show_legend) +
      
      geom_text(data = d, size = 2.5,  
                aes(label=as.character(statname_new)),hjust=0, vjust=2, angle= 45,
                check_overlap = T)+
      
      # geom_ribbon(data = d, aes(fill = sensor,
      #                           x = .data[[input_x]], y = .data[[input_y]], 
      #                           ymin = (.data[[input_y]] - .data[[uncertainty]]),
      #                           ymax = (.data[[input_y]] + .data[[uncertainty]])),
      #             
      #             alpha = 0.3, show.legend = F,
      #             colour="NA"
      # )+
      
      scale_shape_manual(name = " ",
                         labels=(c(expression(paste("COPS")),
                                   expression(paste("PSR")))),
                         values = c(21,23))+
      
      # scale_color_manual(name = " ", 
      #                    labels=(c(expression(paste("COPS")),expression(paste("PSR")))),
      #                    values = c("black", "red"))+
      
      scale_color_manual(name = " ",
                         labels=(c("COPS_unique", "Duplicate", "PSR_unique")),
                         values = c("purple3", "cyan3", "green3"))+
      
      scale_fill_manual(name = " ",
                        labels=(c("COPS_unique", "Duplicate", "PSR_unique")),
                        values = c("purple3", "cyan3", "green3"))+
      
      geom_abline(slope = 1,linetype="dashed", intercept = 0,
                  colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
      #geom_vline(xintercept = 2.54, color = "orange3", linetype = "dashed", size = 1.2)+
      
      geom_smooth(#aes(color = sensor),
        size=1,level = 0.95,show.legend = F,linetype = "solid",
        color="red4",
        se= T, method = "lm")+
      
      coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                  ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
      scale_x_continuous(name = xlabel, limits = c(xmin, xmax), 
                         breaks = seq(xmin, xmax, xstp)) +
      scale_y_continuous(name = ylabel, limits = c(ymin, ymax), 
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
            legend.position = c(0.01,1),
            #legend.direction = "vertical",
            legend.title = element_blank(),
            legend.text = element_text(colour = "black", size = 10, face = "plain"),
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
            legend.direction = "horizontal", legend.box = "vertical",
            legend.text.align = 0,
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
    
    g <- ggMarginal(groupFill = F, data = d, type = "densigram", bins = hist_count, color = "grey",
                    p = g, aes(x = .data[[input_x]], y = .data[[input_y]]))
    return(g)
  }
  
g1 = plot_inversion_validation_multivar_linear_contour(
  input_df = inv_res_df_bayes_unconst, 
                                                   input_x = "H_actual", input_y = "H",
                                                   uncertainty = "H_CI", 
  xmin = 0, xmax = 12, 
                                      
  xlabel = expression(paste(italic("H"["actual"]), "[m]")),
                                      
  ylabel = expression(paste(italic("H"["predicted"]), "[m]")), 
                                      
  opacity = 0.8, plot_col = "purple4", xstp = 2, ystp = 2, 
                                      
  show_legend = T, hist_count = 30)

ggsave(paste0("./outputs/wise_2019_shallow_H_cops-psr_bayes.png"), plot = g1, scale = 1.25, 
       width = 4.5, height = 4.5,
       units = "in",dpi = 300)

g2 = plot_inversion_validation_multivar_linear_contour(input_df = inv_res_df_bayes_unconst, 
                                                   input_x = "chl_actual", input_y = "chl",
                                                   uncertainty = "sd_chl", xmin = 0, xmax = 10, 
                              xlabel = expression(paste(italic("[chl]"["actual"]), "mg/m^3")),
                              ylabel = expression(paste(italic("[chl]"["predicted"]), "mg/m^3")), 
                                                   opacity = 0.8, plot_col = "green4", xstp = 2, 
                              ystp = 2, 
                                                   show_legend = F, hist_count = 20)


ggsave(paste0("./outputs/wise_2019_shallow_chl_cops-psr_grad.png"), plot = g2, scale = 1.25, 
       width = 4.5, height = 4.5,
       units = "in",dpi = 300)

g3 = plot_inversion_validation_singlevar_linear_contour(input_df = inversion_results_merged, 
                                                   input_x = "adg443_actual", 
                                                   input_y = "adg443_predicted",
                                                   uncertainty = "sd_H", xmin = 0, xmax = 4, 
                               xlabel = expression(paste(italic("a"["dg,actual"](443)), "m"^{-1})),
                              ylabel = expression(paste(italic("a"["dg,predicted"](443)), "m"^{-1})), 
                                             opacity = 0.8, plot_col = "goldenrod4", 
                              xstp = 1, ystp = 1, 
                                                   show_legend = F, hist_count = 30)
ggsave(paste0("./outputs/wise_2019_shallow_adg443_comp.png"), plot = g, scale = 1.25, 
       width = 4.5, height = 4.5,
       units = "in",dpi = 300)

g4 = plot_inversion_validation_singlevar_linear_contour(input_df = inversion_results_merged, 
                                                   input_x = "bbp555_actual", 
                                                   input_y = "bbp555_predicted",
                                                   uncertainty = "sd_H", xmin = 0, xmax = 0.02, 
                                  xlabel = expression(paste(italic("b"["bp,actual"](555)), "m"^{-1})),
                                ylabel = expression(paste(italic("b"["bp,predicted"](555)), "m"^{-1})), 
                                   opacity = 0.8, plot_col = "blue4", xstp = 0.01, ystp = 0.01, 
                                                   show_legend = F, hist_count = 30)

ggsave(paste0("./outputs/wise_2019_shallow_bbp555_comp.png"), plot = g, scale = 1.25, 
       width = 4.5, height = 4.5,
       units = "in",dpi = 300)


wise_h_hocr = data.frame("H_actual" = water_depth_HS$H_actual[water_depth_HS$site == "MP"],
                         "H_predicted" = water_depth_HS$H_predicted[water_depth_HS$site == "MP"],
                         "H_sd" = water_depth_HS$H_sd[water_depth_HS$site == "MP"],
                         "H_CI" = water_depth_HS$H_CI[water_depth_HS$site == "MP"],
                         "sensor" = "HOCR")
wise_h_cops = data.frame("H_actual" = inv_res_df_bayes_unconst$H_actual[inv_res_df_bayes_unconst$sensor 
                                                                        == "COPS"],
                         "H_predicted" = inv_res_df_bayes_unconst$H[inv_res_df_bayes_unconst$sensor 
                                                                           == "COPS"],
                         "H_sd" = inv_res_df_bayes_unconst$sd_H[inv_res_df_bayes_unconst$sensor 
                                                                    == "COPS"],
                         "H_CI" = inv_res_df_bayes_unconst$H_CI[inv_res_df_bayes_unconst$sensor 
                                                                    == "COPS"],
                         "sensor" = "COPS")

wise_h_psr = data.frame("H_actual" = inv_res_df_bayes_unconst$H_actual[inv_res_df_bayes_unconst$sensor 
                                                                        == "PSR"],
                         "H_predicted" = inv_res_df_bayes_unconst$H[inv_res_df_bayes_unconst$sensor 
                                                                    == "PSR"],
                         "H_sd" = inv_res_df_bayes_unconst$sd_H[inv_res_df_bayes_unconst$sensor 
                                                                == "PSR"],
                         "H_CI" = inv_res_df_bayes_unconst$H_CI[inv_res_df_bayes_unconst$sensor 
                                                                == "PSR"],
                         "sensor" = "PSR")

wise_h_combined = rbind(wise_h_hocr, wise_h_cops, wise_h_psr)

wise_h_combined$sensor <- factor(wise_h_combined$sensor, levels = c("HOCR", "COPS", "PSR"))

legend_title <- element_blank()
legend_position <- c(0.10, 0.95)
show_legend = T

xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))

ymin <- 0; ymax <- 12 ; ystp <- ymax/4
xmin <- 0; xmax <- 12; xstp <- ymax/4
asp_rat <- (xmax-xmin)/(ymax-ymin)

g<-   ggplot(data=wise_h_combined, aes(x = H_actual, y = H_predicted, colour = as.factor(sensor), 
                                       shape = as.factor(sensor),
                                     fill = as.factor(sensor))) +
  
  
  #geom_contour(aes(z = z), col = "black")+
  
  geom_ribbon(aes(ymin = H_predicted - H_CI,
                  ymax = H_predicted + H_CI#, fill = as.numeric(bottom_cover)
  ), fill = "gray",
  alpha = 0.7, show.legend = F,
  
  colour="NA"
  )+
  
  geom_density_2d(data = wise_h_combined, aes(x = H_actual, y = H_predicted), na.rm = T, bins = 6,
                  linewidth = 0.25,  show.legend = F, size=1.1)+
  
  geom_point(data = wise_h_combined, aes(H_actual, H_predicted, #shape = as.factor(pd),
  ), 
  alpha = I(0.6), size = I(4), show.legend = T) +
  
  geom_vline(xintercept = 2.5, color = "navyblue", size = 1.3, linetype = "dashed")+
  
  
  # scale_fill_gradient(name = "Plant_Cover", low = "#440154FF", high = "#9FDA3AFF", 
  #                     na.value = "#20A387FF",
  #                     breaks=seq(0,100,25),
  #                     limits=c(0,100), labels=paste(seq(0,100,25),"%"))+
  # 
  # scale_shape_manual(name ="Penetration Depth", labels=(c(expression(paste(italic("H"), ">", "P"["d"])),
  #                                        expression(paste(italic("H"), "<", "P"["d"])))),
  #                    values = c(21,23))+
  
  scale_colour_viridis(name =" ", discrete = T, labels=(c("HOCR", "COPS", "PSR")),
                      )+
  
  scale_fill_viridis(name =" ", discrete = T, labels=(c("HOCR", "COPS", "PSR")),
                    )+
  scale_shape_manual(name = " ", values = c(21,22,23), labels= c("HOCR", "COPS", "PSR") ) +
  
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
        legend.position = c(0.65, 0.30),
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

g <- ggMarginal(groupFill = F, data = H_df_combined, type = "densigram", bins = 60, color = "grey",
                p = g, aes(x = H_actual, y = H_predicted))

ggsave(paste0("./outputs/shallow_H_MP_unconstr_scatter.png"), plot = g,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)


library(dplyr)
library(tidyr)
library(pbapply)
library(parallel)
library(doParallel)
library(foreach)

# Function :: Interpolate the bottom reflectance to Rrs wavelengths ----
interp_rb <- function(wavelength_input, bottom_type = c("Mud_2019", "Sand_2019", "Eelgrass_2019")){
  load("./data/EGSL_all.RData")
  
  #column_numbers <- which(names(rb) %in% bottom_type)
  
  column_numbers <- match(bottom_type, names(rb))
  
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


Rb_set = read.csv("./data/EGSL_Rb.csv", header = T)
benthic_classes_2019 <- colnames(Rb_set)[grep("2019", colnames(Rb_set))]
benthic_classes_2022 <- colnames(Rb_set)[grep("2022", colnames(Rb_set))]
benthic_classes_2023 <- colnames(Rb_set)[grep("2023", colnames(Rb_set))]

combinations <- combn(benthic_classes_2023, 3, simplify = FALSE)

combinations <- lapply(combinations, function(x) sort(x))
combinations <- as.data.frame(do.call(rbind, combinations), stringsAsFactors = FALSE)


sdg_values <- seq(0.0015, 0.025, length.out = 10)
gamma_values <- seq(0.1, 1, length.out = 10)
sd_values <- seq(0.0005, 0.005, length.out = 10)
sg_values <- sdg_values - sd_values


spectral_slope_vec <- data.frame("gamma" = gamma_values,
                                 "sd" = sd_values,
                                 "sg" = sg_values)


param_vec_bayes <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(param_vec_bayes) <- c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd", 
                               "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", "sd_fa1", "sd_fa2", 
                               "sd_fa3", "sd_pop_sd")


# Define the process_combination function (as before)
process_combination <- function(obsdata, idx) {
  # Suppress output from the functions
  capture.output({
    interp_rb(wavelength_input = seq(400, 750, 1), 
              bottom_type = c("Sand_2019", "Mud_2019", "Eelgrass_2019"))
    
    result <- inverse_runBayes(obsdata = as.numeric(obsdata), 
                               rrs_type = type_Rrs_below, 
                               max_par = c(6.5, 2.5, 0.01, 8, 1, 1, 1, 10),
                               min_par = c(1.5, 1, 0.003, 2, 0.1, 0.1, 0.1, 0.01), 
                               param_lab = c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd"),
                               constrain_config = c(FALSE, FALSE, FALSE),
                               bbp_const_val = 0.007, 
                               bgc_const_val = c(5, 1.003, 0.007), 
                               iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                                                  "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
                               qaa_slope = FALSE, 
                               manual_slope = TRUE, 
                               manual_slope_vals = c("s_g" = spectral_slope_vec[idx, 3], 
                                                     "s_d" = spectral_slope_vec[idx, 2], 
                                                     "gamma" = spectral_slope_vec[idx, 1]), 
                               iter_count = 20000, 
                               sampler_mcmc = samplerlist[6], 
                               wavelngth_sim = seq(400, 750, 1), 
                               sa_model = "am03", 
                               hybrid_mode = FALSE, 
                               plot_rrs = FALSE)
  }, file = NULL) 
  
  return(result)
}


results_list <- list()


for (idx in 1:nrow(spectral_slope_vec)) {
  

  fit_results <- apply(inverse_rrs_input[,1:2], 2, process_combination, idx = idx)
  

  results_list[[idx]] <- fit_results
}





