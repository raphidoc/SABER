#==================================================================================================
# main_Rb_spectral.R performs spectral retrieval of bottom reflectance (Rb) based on Radiative 
# transfer in optically complex case-II waters. There are two modes to retrieve Rb:

# 1. Model the Remote-Sensing Reflectance (Rrs) in optically complex waters from user-given [chl],
# CDOM and NAP absorption at 440nm wavelength and depth and subtract the model value from observed 
# Rrs using algebric reparametrization of AM03. 

# 2. Model the Remote-Sensing Reflectance (Rrs) in optically complex waters from user-given spectral
# IOP and depth and subtract the model value using algebric reparametrization of AM03.

#The function to calculate Rb is referred as "Saber_retrieve_rb_wise".
#==================================================================================================

###======================================================================
# 1. Retrieve IOPs, AOPs and BGC variables for the station under process
###======================================================================

#Input the WISE-Man IOP QC sheets for both sets of IOP cages
kildir_iop = read.csv("/home/musk0001/Kildir_IOP.Process_Log.csv", header = T)
saucier_iop = read.csv("/home/musk0001/saucierQC.csv", header = T)

#Input the depth statistics and simulation notes from WISE-Man IOP data
wise_depth_stat = read.csv("/home/musk0001/Depth_wise.csv", header = T)
wise_depth_stat$correct_statname = gsub(wise_depth_stat$Station, pattern = "\\_", replacement = "-") 
wise_depth_stat <- wise_depth_stat %>% dplyr::select(correct_statname, everything())

#Select the QC passed stations for WISE-Man
wise_depth_stat_shallow = wise_depth_stat[wise_depth_stat$Shallow == "T",]
kildir_iop_qc = kildir_iop[kildir_iop$ASPH == "Y" & kildir_iop$HS6 == "Y",]
saucier_iop_qc = saucier_iop[saucier_iop$station.keep == "TRUE",]

qc_statnames_bind = rbind(list(kildir_iop_qc$StationID), 
                                                       list(saucier_iop_qc$Station))

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

shallow_qc_stations = qc_statnames_bind$qcstat[!is.na(shallow_qc_keep_idx)]
shallow_qc_stations = shallow_qc_stations[-(3:4)] #To get rid of stations with no Rrs

#=============================================================
## 2. BATCH Rb retrieval for WISE-Man in situ
#=============================================================
Saber_generate_rb_batch(staionlist = shallow_qc_stations)

#=======================================================================================================
#=======================================================================================================
# Saber_generate_rb_batch <- function(staionlist, save_plot=T, save_data=T) {
#   
#   shallow_qc_stations = staionlist
#   
#   for (i in 1:length(shallow_qc_stations)) {
#     stat_temp = shallow_qc_stations[i]
#     cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
#     cat(paste0("\033[0;36m","##################### Rb retrieval start for Station: ", stat_temp, " #########################","\033[0m","\n"))
#     cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
#     
#     iop_aop_rb = suppressWarnings(read.surface.IOPs.wise(station.args = stat_temp,save_on_disc = T, 
#                                                          verbose = F))
#     
#     bgc_params_rb = get_in_situ_params(station_name = stat_temp, use_bb_nup = T)
#     
#     IOP_files = list.files("./data/Rb_spectral/surface_iops/", full.names = T)
#     
#     idx_a_rb = grep(IOP_files, pattern = paste0("abs_surf_",stat_temp, ".csv$"))
#     idx_bb_rb = grep(IOP_files, pattern = paste0("bb_surf_",stat_temp, ".csv$"))
#     
#     ###======================================================================
#     # 4. Retrieve spectral Rb by providing IOPs and Depth
#     ###======================================================================
#     if (type_Rrs_below_rb == "shallow") {
#       #Generate Rb using forward model parameters
#       if (insitu.present == FALSE) {
#         
#         Rrs_bottom_gen = Saber_generate_rb(lambda = wavelength, fA = fA.set)
#         print(paste0("Bottom composition is linear mixture of sand = ", fA.set[1], " , sediment = ",
#                      fA.set[2], " , Chara contraria = ", fA.set[3],
#                      " , Potamogeton perfoliatus = ", fA.set[4],
#                      " , Potamogeton pectinatus = ", fA.set[5]))
#         
#         print("in situ Rrs is not supplied, user simulated Rrs is used for bottom reflectance retrieval")
#         
#         if (type_Rrs_water == "below_surface") {
#           
#           obs_rrs =  rrs.forward.am.param.conc.true_iop
#         } else {
#           obs_rrs = surface_rrs_translate(Rrs = rrs.forward.am.param.conc.true_iop)
#         }
#         
#       } else {
#         print("Historic data of the bottom composition is unknown for the station")
#         print("in situ Rrs is supplied for bottom reflectance retrieval")
#         obs_rrs = surface_rrs_translate(Rrs = iop_aop_rb$Rrs_0p)
#       }
#       
#       cat(paste0("\033[0;36m","=======================PARAMETRIC RETIREVAL=================================================","\033[0m","\n"))
#       #Retrieve spectral Rb from forward SABER generated Rrs with BGC parametrization
#       Rrs_bottom_est_bgc = Saber_retrieve_rb_wise(use_true_IOPs = F, 
#                                                   a_non_water_path = IOP_files[idx_a_rb],
#                                                   bb_non_water_path = IOP_files[idx_bb_rb],
#                                                   
#                                                   chl = bgc_params_rb$chl_invivo, 
#                                                   acdom440 =bgc_params_rb$abs_invivo$a_cdom, 
#                                                   anap440 = bgc_params_rb$abs_invivo$a_nap,
#                                                   bbp.550 = bgc_params_rb$bbp555,
#                                                   
#                                                   slope.parametric = F, dg_composite = F,
#                                                   
#                                                   #z = 2,
#                                                   z = iop_aop_rb$zB_COPS,
#                                                   
#                                                   #obs_rrs = rrs.forward.am.param.conc.true_iop,
#                                                   obs_rrs =  obs_rrs
#       )
#       
#       cat(paste0("\033[0;36m","=======================IOP based RETIREVAL=================================================","\033[0m","\n"))
#       #Retrieve spectral Rb from forward SABER generated Rrs with True IOPs
#       Rrs_bottom_est_iop = Saber_retrieve_rb_wise(use_true_IOPs = T, 
#                                                   a_non_water_path = IOP_files[idx_a_rb],
#                                                   bb_non_water_path = IOP_files[idx_bb_rb],
#                                                   
#                                                   chl = bgc_params_rb$chl_invivo, 
#                                                   acdom440 =bgc_params_rb$abs_invivo$a_cdom, 
#                                                   anap440 = bgc_params_rb$abs_invivo$a_nap,
#                                                   bbp.550 = bgc_params_rb$bbp555,
#                                                   
#                                                   slope.parametric = F, dg_composite = F,
#                                                   
#                                                   #z = 2,
#                                                   z = iop_aop_rb$zB_COPS,
#                                                   
#                                                   #obs_rrs = rrs.forward.am.param.conc.true_iop,
#                                                   obs_rrs = obs_rrs
#       )
#       
#       cat(paste0("\033[0;36m","=======================SAVING OUTPUT=================================================","\033[0m","\n"))
#       #Compare the actual vs retrieved Rb as a function of wavelength (#Change the plot to ggplot)
#       if (insitu.present == TRUE) {
#         bottom_ref_df = data.frame("wave"= as.numeric(wavelength),
#                                    "rb_est_bgc" = as.numeric(Rrs_bottom_est_bgc$rb),
#                                    "rb_est_iop" = as.numeric(Rrs_bottom_est_iop$rb),
#                                    "rrs_shallow" = obs_rrs) 
#         
#         xmin = min(bottom_ref_df$wave); xmax= max(bottom_ref_df$wave); xstp=100
#         ymin= 0; ymax=1.25*max(bottom_ref_df$rb_est_iop[bottom_ref_df$wave < "650"], na.rm = T)
#         ystp= signif(ymax/5, digits = 2)
#         asp_rat <- (xmax-xmin)/(ymax-ymin)
#         #browser()
#         g <- ggplot()  + geom_line(data = bottom_ref_df,aes(y=rrs_shallow,color="xx1",x = wave),
#                                    size=1.3,show.legend = TRUE) +
#           geom_line(data = bottom_ref_df,aes(y=rb_est_bgc,x = wave,color="xx3"), 
#                     size=1.3,show.legend = TRUE)+
#           geom_line(data = bottom_ref_df,aes(y=rb_est_iop,x = wave,color="xx4"), 
#                     size=1.3,show.legend = TRUE)+
#           
#           scale_colour_manual(labels = c(expression(paste(italic("R")["rs,observed"])),
#                                          expression(paste(italic("R")["B, model_BOPT"])),
#                                          expression(paste(italic("R")["B, model_IOP"]))
#           ), 
#           values = c("purple","navyblue", "red")) +
#           ggtitle(paste0(stat_temp)) +
#           scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
#                              limits = c(xmin, xmax), 
#                              breaks = seq(xmin, xmax, xstp))  +
#           scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")[", sr^-1,"]")) , limits = c(ymin, ymax),
#                              breaks = seq(ymin, ymax, ystp))+
#           coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
#                       ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#           theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#                 axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
#                 axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
#                 axis.title.x = element_text(size = 25),
#                 axis.title.y = element_text(size = 25),
#                 axis.ticks.length = unit(.25, "cm"),
#                 legend.position=c(0.6,1),
#                 legend.direction = "vertical",
#                 legend.title = element_blank(),
#                 legend.text = element_text(colour = "black", size = 20, face = "plain"),
#                 legend.background = element_rect(fill = NA, size = 0.5, 
#                                                  linetype = "solid", colour = 0),
#                 legend.key = element_blank(),
#                 legend.justification = c("left", "top"),
#                 panel.background = element_blank(),
#                 panel.grid.major = element_line(colour = "grey", 
#                                                 size = 0.5, linetype = "dotted"), 
#                 panel.grid.minor = element_blank(),
#                 #legend.spacing.y = unit(2.0, 'cm'),
#                 plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
#                 legend.text.align = 0,
#                 panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
#         
#         
#       } else {
#         bottom_ref_df = data.frame("wave"= as.numeric(wavelength), 
#                                    "rb_actual" = as.numeric(Rrs_bottom_gen$rb),
#                                    "rb_est_bgc" = as.numeric(Rrs_bottom_est_bgc$rb),
#                                    "rb_est_iop" = as.numeric(Rrs_bottom_est_iop$rb),
#                                    "rrs_shallow" =obs_rrs) 
#         
#         xmin = min(bottom_ref_df$wave); xmax= max(bottom_ref_df$wave); xstp=100
#         ymin= 0; ymax=1.25*max(bottom_ref_df$rb_actual[bottom_ref_df$wave < "650"], na.rm = T)
#         ystp= signif(ymax/5, digits = 2)
#         asp_rat <- (xmax-xmin)/(ymax-ymin)
#         
#         g <- ggplot()  + geom_line(data = bottom_ref_df,aes(y=rrs_shallow,color="xx1",x = wave),
#                                    size=1.3,show.legend = TRUE) +
#           geom_line(data = bottom_ref_df,aes(y=rb_est_bgc,x = wave,color="xx3"), 
#                     size=1.3,show.legend = TRUE)+
#           geom_point(data = bottom_ref_df,aes(y=rb_actual,x = wave,color="xx2"), 
#                      size=3,show.legend = TRUE)+
#           geom_line(data = bottom_ref_df,aes(y=rb_est_iop,x = wave,color="xx4"), 
#                     size=1.3,show.legend = TRUE)+
#           
#           scale_colour_manual(labels = c(expression(paste(italic("R")["rs,observed"])),
#                                          expression(paste(italic("R")["B, actual"])),
#                                          expression(paste(italic("R")["B, model_BOPT"])),
#                                          expression(paste(italic("R")["B, model_IOP"]))
#           ), 
#           values = c("purple","black","navyblue", "red")) +
#           ggtitle(paste0(stat_temp)) +
#           scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
#                              limits = c(xmin, xmax), 
#                              breaks = seq(xmin, xmax, xstp))  +
#           scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")[", sr^-1,"]")) , limits = c(ymin, ymax),
#                              breaks = seq(ymin, ymax, ystp))+ 
#           coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
#                       ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#           theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#                 axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
#                 axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
#                 axis.title.x = element_text(size = 25),
#                 axis.title.y = element_text(size = 25),
#                 axis.ticks.length = unit(.25, "cm"),
#                 legend.position=c(0.6,1),
#                 legend.direction = "vertical",
#                 legend.title = element_blank(),
#                 legend.text = element_text(colour = "black", size = 20, face = "plain"),
#                 legend.background = element_rect(fill = NA, size = 0.5, 
#                                                  linetype = "solid", colour = 0),
#                 legend.key = element_blank(),
#                 legend.justification = c("left", "top"),
#                 panel.background = element_blank(),
#                 panel.grid.major = element_line(colour = "grey", 
#                                                 size = 0.5, linetype = "dotted"), 
#                 panel.grid.minor = element_blank(),
#                 #legend.spacing.y = unit(2.0, 'cm'),
#                 plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
#                 legend.text.align = 0,
#                 panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
#       }
#       
#       g 
#       cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
#       cat(paste0("\033[0;32m","##################### Rb retrieval end for Station: ", stat_temp, " #########################","\033[0m","\n"))
#       cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
#       
#       #Save the bottom ref to disk
#       if (save_data == TRUE) {
#         write.csv(file = paste0("./Outputs/Bottom_ref/Rb_spectral_data_",stat_temp,".csv"), 
#                   x = bottom_ref_df, quote = F, sep = ",", col.names = T, row.names = F)
#       }
#       
#       
#       if (save_plot == TRUE) {
#         ggsave(filename = paste0("./Outputs/Bottom_ref/Rb_spectral_station_",stat_temp,".png"), plot = g,
#                scale = 1.5, width = 4.5, height = 4.5, units = "in", dpi = 300)
#       }
#     } else {
#       print("Rb retrieval is not applicable for deep water column")
#     }
#     
#   }
# }
# 
# stationID = "MAN-F16"
# iop_aop_rb = suppressWarnings(read.surface.IOPs.wise(station.args = stationID,save_on_disc = T))
# 
# bgc_params_rb = get_in_situ_params(station_name = stationID, use_bb_nup = T)
# 
# IOP_files = list.files("./data/Rb_spectral/surface_iops/", full.names = T)
# 
# idx_a_rb = grep(IOP_files, pattern = paste0("abs_surf_",stationID))
# idx_bb_rb = grep(IOP_files, pattern = paste0("bb_surf_",stationID))
# 
# ###======================================================================
# # 4. Retrieve spectral Rb by providing IOPs and Depth
# ###======================================================================
# if (type_Rrs_below_rb == "shallow") {
#   #Generate Rb using forward model parameters
#   if (insitu.present == FALSE) {
#     
#     Rrs_bottom_gen = Saber_generate_rb(lambda = wavelength, fA = fA.set)
#     print(paste0("Bottom composition is linear mixture of sand = ", fA.set[1], " , sediment = ",
#                  fA.set[2], " , Chara contraria = ", fA.set[3],
#                  " , Potamogeton perfoliatus = ", fA.set[4],
#                  " , Potamogeton pectinatus = ", fA.set[5]))
#     
#   } else {
#     print("Historic data of the bottom composition is unknown for the station")
#   }
#   
#   
#   #Retrieve spectral Rb from forward SABER generated Rrs with BGC parametrization
#   Rrs_bottom_est_bgc = Saber_retrieve_rb_wise(use_true_IOPs = F, 
#                        a_non_water_path = IOP_files[idx_a_rb],
#                        bb_non_water_path = IOP_files[idx_bb_rb],
#                                               
#                        chl = bgc_params_rb$chl_invivo, 
#                        acdom440 =bgc_params_rb$abs_invivo$a_cdom, 
#                        anap440 = bgc_params_rb$abs_invivo$a_nap,
#                        bbp.550 = bgc_params_rb$bbp555,
#                                               
#                       slope.parametric = F, dg_composite = F,
#                                               
#                       #z = 2,
#                       z = iop_aop_rb$zB_COPS,
#                                               
#                       #obs_rrs = rrs.forward.am.param.conc.true_iop,
#                       obs_rrs =  surface_rrs_translate(Rrs = iop_aop_rb$Rrs_0p)
#   )
#   
#   
#   #Retrieve spectral Rb from forward SABER generated Rrs with True IOPs
#   Rrs_bottom_est_iop = Saber_retrieve_rb_wise(use_true_IOPs = T, 
#                                               a_non_water_path = IOP_files[idx_a_rb],
#                                               bb_non_water_path = IOP_files[idx_bb_rb],
#                                               
#                        chl = bgc_params_rb$chl_invivo, 
#                        acdom440 =bgc_params_rb$abs_invivo$a_cdom, 
#                        anap440 = bgc_params_rb$abs_invivo$a_nap,
#                        bbp.550 = bgc_params_rb$bbp555,
#                                               
#                        slope.parametric = F, dg_composite = F,
#                                               
#                        #z = 2,
#                        z = iop_aop_rb$zB_COPS,
#                                               
#                        #obs_rrs = rrs.forward.am.param.conc.true_iop,
#                        obs_rrs =  surface_rrs_translate(Rrs = iop_aop_rb$Rrs_0p)
#   )
#   
#   #Compare the actual vs retrieved Rb as a function of wavelength (#Change the plot to ggplot)
#   if (insitu.present == FALSE) {
#     bottom_ref_df = data.frame("wave"= as.numeric(wavelength), 
#                                "rb_actual" = as.numeric(Rrs_bottom_gen$rb),
#                                "rb_est_bgc" = as.numeric(Rrs_bottom_est_bgc$rb),
#                                "rb_est_iop" = as.numeric(Rrs_bottom_est_iop$rb),
#                                "rrs_shallow" =surface_rrs_translate(Rrs = iop_aop_rb$Rrs_0p)) 
#     
#     xmin = min(bottom_ref_df$wave); xmax= max(bottom_ref_df$wave); xstp=100
#     ymin= 0; ymax=1.25*max(bottom_ref_df$rb_actual[bottom_ref_df$wave < "650"])
#     ystp= signif(ymax/5, digits = 2)
#     asp_rat <- (xmax-xmin)/(ymax-ymin)
#     
#     g <- ggplot()  + geom_line(data = bottom_ref_df,aes(y=rrs_shallow,color="xx1",x = wave),
#                                size=1.3,show.legend = TRUE) +
#       geom_line(data = bottom_ref_df,aes(y=rb_est_bgc,x = wave,color="xx3"), 
#                 size=1.3,show.legend = TRUE)+
#       geom_point(data = bottom_ref_df,aes(y=rb_actual,x = wave,color="xx2"), 
#                  size=3,show.legend = TRUE)+
#       geom_line(data = bottom_ref_df,aes(y=rb_est_iop,x = wave,color="xx4"), 
#                 size=1.3,show.legend = TRUE)+
#       
#       scale_colour_manual(labels = c(expression(paste(italic("R")["rs,observed"])),
#                                      expression(paste(italic("R")["B, actual"])),
#                                      expression(paste(italic("R")["B, model_BOPT"])),
#                                      expression(paste(italic("R")["B, model_IOP"]))
#       ), 
#       values = c("purple","black","navyblue", "red")) +
#       #ggtitle(paste0(stationlist[i])) +
#       scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
#                          limits = c(xmin, xmax), 
#                          breaks = seq(xmin, xmax, xstp))  +
#       scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")[", sr^-1,"]")) , limits = c(ymin, ymax),
#                          breaks = seq(ymin, ymax, ystp))+ 
#       coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
#                   ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#       theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#             axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
#             axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
#             axis.title.x = element_text(size = 25),
#             axis.title.y = element_text(size = 25),
#             axis.ticks.length = unit(.25, "cm"),
#             legend.position=c(0.6,1),
#             legend.direction = "vertical",
#             legend.title = element_blank(),
#             legend.text = element_text(colour = "black", size = 20, face = "plain"),
#             legend.background = element_rect(fill = NA, size = 0.5, 
#                                              linetype = "solid", colour = 0),
#             legend.key = element_blank(),
#             legend.justification = c("left", "top"),
#             panel.background = element_blank(),
#             panel.grid.major = element_line(colour = "grey", 
#                                             size = 0.5, linetype = "dotted"), 
#             panel.grid.minor = element_blank(),
#             #legend.spacing.y = unit(2.0, 'cm'),
#             plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
#             legend.text.align = 0,
#             panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
#     
#   } else {
#     bottom_ref_df = data.frame("wave"= as.numeric(wavelength),
#                                "rb_est_bgc" = as.numeric(Rrs_bottom_est_bgc$rb),
#                                "rb_est_iop" = as.numeric(Rrs_bottom_est_iop$rb),
#                                "rrs_shallow" = rrs.forward.am.param.conc.true_iop) 
#     
#     xmin = min(bottom_ref_df$wave); xmax= max(bottom_ref_df$wave); xstp=100
#     ymin= 0; ymax=1.25*max(bottom_ref_df$rb_est_iop[bottom_ref_df$wave < "650"])
#     ystp= signif(ymax/5, digits = 2)
#     asp_rat <- (xmax-xmin)/(ymax-ymin)
#     
#     g <- ggplot()  + geom_line(data = bottom_ref_df,aes(y=rrs_shallow,color="xx1",x = wave),
#                                size=1.3,show.legend = TRUE) +
#       geom_line(data = bottom_ref_df,aes(y=rb_est_bgc,x = wave,color="xx3"), 
#                 size=1.3,show.legend = TRUE)+
#       geom_line(data = bottom_ref_df,aes(y=rb_est_iop,x = wave,color="xx4"), 
#                 size=1.3,show.legend = TRUE)+
#       
#       scale_colour_manual(labels = c(expression(paste(italic("R")["rs,observed"])),
#                                      expression(paste(italic("R")["B, model_BOPT"])),
#                                      expression(paste(italic("R")["B, model_IOP"]))
#       ), 
#       values = c("purple","navyblue", "red")) +
#       #ggtitle(paste0(stationlist[i])) +
#       scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
#                          limits = c(xmin, xmax), 
#                          breaks = seq(xmin, xmax, xstp))  +
#       scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")[", sr^-1,"]")) , limits = c(ymin, ymax),
#                          breaks = seq(ymin, ymax, ystp))+ 
#       coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
#                   ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#       theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#             axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
#             axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
#             axis.title.x = element_text(size = 25),
#             axis.title.y = element_text(size = 25),
#             axis.ticks.length = unit(.25, "cm"),
#             legend.position=c(0.6,1),
#             legend.direction = "vertical",
#             legend.title = element_blank(),
#             legend.text = element_text(colour = "black", size = 20, face = "plain"),
#             legend.background = element_rect(fill = NA, size = 0.5, 
#                                              linetype = "solid", colour = 0),
#             legend.key = element_blank(),
#             legend.justification = c("left", "top"),
#             panel.background = element_blank(),
#             panel.grid.major = element_line(colour = "grey", 
#                                             size = 0.5, linetype = "dotted"), 
#             panel.grid.minor = element_blank(),
#             #legend.spacing.y = unit(2.0, 'cm'),
#             plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
#             legend.text.align = 0,
#             panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
#   }
#   
#   g 
#   #Save the bottom ref to disk
#   write.csv(file = paste0("./Outputs/Bottom_ref/Rb_spectral_data_",stat_temp,".csv"), 
#             x = bottom_ref_df, quote = F, sep = ",", col.names = T, row.names = F)
#   
#   if (plot == "TRUE") {
#     ggsave(filename = paste0("./Outputs/Bottom_ref/Rb_spectral_station_",stationID,".png"), plot = g,
#            scale = 1.5, width = 4.5, height = 4.5, units = "in", dpi = 300)
#   }
# } else {
#   print("Rb retrieval is not applicable for deep water column")
# }

# #1.1 Read Kildir absorption data
# absdata_kildir = read.csv("./data/Rb_spectral/absorption_final/At-w.All_Kildir_surf-QC.csv",
#                           header = T)
# 
# abs_surf_depth_kildir = colnames(absdata_kildir) #extract depth from where acquisition started
# 
# colnames(absdata_kildir) =  t(absdata_kildir[1,])
# absdata_kildir = absdata_kildir[-1,] #Set correct column-names
# 
# abs_surf_depth_kildir_df = data.frame("station" = colnames(absdata_kildir)[-1],
#                                       "surface_depth" = substr(abs_surf_depth_kildir[-1],2,4),
#                                       "boat" = "kildir") #Create surface depth data-frame
# 
# #1.2 Read Saucier absorption data
# absdata_saucier = read.csv("./data/Rb_spectral/absorption_final/At-w.All_Saucier_surf-QC_new.csv", 
#                            header = T)
# abs_surf_depth_saucier = colnames(absdata_saucier) #extract depth from where acquisition started
# 
# colnames(absdata_saucier) =  t(absdata_saucier[1,])
# absdata_saucier = absdata_saucier[-1,] #Set correct column-names
# 
# abs_surf_depth_saucier_df = data.frame("station" = colnames(absdata_saucier)[-1],
#                                       "surface_depth" = substr(abs_surf_depth_saucier[-1],2,4),
#                                       "boat" = "saucier") #Create surface depth data-frame
# 
# #1.3 Create unified long data-format for absorption
# absdata_kildir_long = reshape2::melt(absdata_kildir, id.vars = "wavelength") #kildir
# colnames(absdata_kildir_long) = c("wave", "station", "at_w")
# absdata_kildir_long$boat = "kildir"
# 
# absdata_saucier_long = reshape2::melt(absdata_saucier, id.vars = "wavelength") #saucier
# colnames(absdata_saucier_long) = c("wave", "station", "at_w")
# absdata_saucier_long$boat = "saucier"
# 
# absdata_all_surf_long = rbind(absdata_kildir_long, absdata_saucier_long) #Final at_w
# absdata_all_surf_depth = rbind(abs_surf_depth_kildir_df, abs_surf_depth_saucier_df) #Final depth
# 
# rm(absdata_kildir, abs_surf_depth_kildir, abs_surf_depth_kildir_df,
#    absdata_saucier, abs_surf_depth_saucier, abs_surf_depth_saucier_df,
#    absdata_kildir_long, absdata_saucier_long)
# 
# ###======================================================================
# # 2. READ BACKSCATTER DATA
# ###======================================================================
# 
# #2.1 Read Kildir backscatter data
# bbdata_kildir = read.csv("./data/Rb_spectral/backscatter_final/Bbp.All_Kildir_surf_new.csv",
#                          header = T)
# bb_surf_depth_kildir = colnames(bbdata_kildir) #extract depth from where acquisition started
# 
# colnames(bbdata_kildir) =  t(bbdata_kildir[1,])
# bbdata_kildir = bbdata_kildir[-1,] #Set correct column-names
# 
# bb_surf_depth_kildir_df = data.frame("station" = colnames(bbdata_kildir)[-1],
#                                       "surface_depth" = substr(bb_surf_depth_kildir[-1],2,4),
#                                       "boat" = "kildir") #Create surface depth data-frame
# 
# #2.2 Read Saucier backscatter data
# bbdata_saucier = read.csv("./data/Rb_spectral/backscatter_final/Bbp.All_Saucier_surf-QC_new.csv",
#                          header = T)
# bb_surf_depth_saucier = colnames(bbdata_saucier) #extract depth from where acquisition started
# colnames(bbdata_saucier) =  t(bbdata_saucier[1,])
# bbdata_saucier = bbdata_saucier[-1,] #Set correct column-names
# 
# bb_surf_depth_saucier_df = data.frame("station" = colnames(bbdata_saucier)[-1],
#                                      "surface_depth" = substr(bb_surf_depth_saucier[-1],2,4),
#                                      "boat" = "saucier") #Create surface depth data-frame
# 
# #2.3 Create unified long data-format for backscatter
# bbdata_kildir_long = reshape2::melt(bbdata_kildir, id.vars = "wavelength") #kildir
# colnames(bbdata_kildir_long) = c("wave", "station", "bbp")
# bbdata_kildir_long$boat = "kildir"
# 
# bbdata_saucier_long = reshape2::melt(bbdata_saucier, id.vars = "wavelength") #saucier
# colnames(bbdata_saucier_long) = c("wave", "station", "bbp")
# bbdata_saucier_long$boat = "saucier"
# 
# bbdata_all_surf_long = rbind(bbdata_kildir_long, bbdata_saucier_long) #Final at_w
# bbdata_all_surf_depth = rbind(bb_surf_depth_kildir_df, bb_surf_depth_saucier_df) #Final depth
# 
# rm(bbdata_kildir, bb_surf_depth_kildir, bb_surf_depth_kildir_df,
#    bbdata_saucier, bb_surf_depth_saucier, bb_surf_depth_saucier_df,
#    bbdata_kildir_long, bbdata_saucier_long)
# 
# ###======================================================================
# # 3. READ Rrs DATA
# ###======================================================================
# 
# #3.1 Read COP-S actual database and hyperspectral database
# cops_db =  read.csv("./data/Rb_spectral/COPS_db.csv",header = T)
# 
# cops_hs =  read.csv("./data/Rb_spectral/copsdata.csv",header = T)
# 
# correct_statnames = gsub( names(cops_hs)[-1], pattern = "\\.", replacement = "-") #replace "." 
#                                                                                   #by "-"
# colnames(cops_hs) = c("wave", correct_statnames) #assign correct station-names
# colnames(cops_db) = c("station", correct_statnames)
# 
# #3.2 Retrieve COP-S acquired depths (CHECK b4 USE)
# cops_db_flip = as.data.frame(t(cops_db)) 
# colnames(cops_db_flip) = cops_db$stationID
# cops_db_flip = cops_db_flip[-1,]
# colnames(cops_db_flip) = cops_db$station
# 
# #3.3 Extract shallow stations as per COP-S db 
# #[!!May not be accurate as it uses COP-S cast depth!!]
# cops_db_shallow_idx = !is.na(cops_db_flip$Bottom.depth)
# cops_db_shallow_idx = which(cops_db_shallow_idx == TRUE)
# 
# cops_db_shallow = cops_db_flip[cops_db_shallow_idx,]
# cops_db_shallow$stationID = rownames(cops_db_shallow)
# cops_db_hs_shallow = cops_hs[,c(1,cops_db_shallow_idx+1)]
# 
# 
# ###======================================================================
# # 4. READ desired IOPs and Rrs DATA
# ###======================================================================
# 
# #4.1 Select the station-name
# stationID = "OUT-R01"
# 
# #4.2 absorption
# absdata_sample = absdata_all_surf_long[absdata_all_surf_long$station == 
#                                          paste0("Station",stationID),]
# 
# # #abs_data_idx = grep(stationID, colnames(absdata_kildir))
# # abs_surf_temp = data.frame("wave" = absdata_kildir$wavelength, 
# #                            "a"= (absdata_kildir[abs_data_idx]))
# # colnames(abs_surf_temp) = c('wave','a')
# 
# write.csv(file = paste0("./data/Rb_spectral/surface_iops/abs_surf_", stationID, ".csv"), 
#           x = absdata_sample, row.names = F, quote = F)
# 
# #4.3 Backscatter
# bbdata_sample = bbdata_all_surf_long[bbdata_all_surf_long$station == 
#                                        paste0("Station",stationID),]
# 
# # bb_data_idx = grep(stationID, colnames(bbdata_kildir))
# # 
# # bb_surf_temp = data.frame("wave" = bbdata_kildir$wavelength, 
# #                            "bb"= (bbdata_kildir[bb_data_idx]))
# # colnames(bb_surf_temp) = c('wave','bb')
# 
# write.csv(file = paste0("./data/Rb_spectral/surface_iops/abs_surf_", stationID, ".csv"), 
#           x = bbdata_sample, row.names = F, quote = F)
# 
# 
# #4.4 Rrs
# rrs_data_idx = grep(stationID, correct_statnames)
# rrs_hs_temp = as.matrix(cops_hs[(rrs_data_idx+1)])
# 
# rrs_hs_temp_interp = approx(x=cops_hs$wave, y = rrs_hs_temp,
#                             xout = wavelength, method = "linear" )$y
# 
# #4.5 Depth from COP-S
# 
# cops_db_correct_statnames = rownames(cops_db_flip)
# 
# depth_idx = grep(stationID, cops_db_correct_statnames)
# depth_input = as.numeric(as.character(cops_db_flip$Bottom.depth[depth_idx]))
