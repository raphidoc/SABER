
read.surface.IOPs.wise <- function(absdata_kildir = "./data/Rb_spectral/absorption_final/At-w.All_Kildir_surf-QC.csv",
                                   absdata_saucier = "./data/Rb_spectral/absorption_final/At-w.All_Saucier_surf-QC_new.csv",
                                   bbdata_kildir = "./data/Rb_spectral/backscatter_final/Bbp.All_Kildir_surf_new.csv",
                                   bbdata_saucier = "./data/Rb_spectral/backscatter_final/Bbp.All_Saucier_surf-QC_new.csv",
                                   station.args = "OUT-F18",
                                   verbose = T,
                                   save_on_disc = F){
  ###======================================================================
  # 1. READ ABSORPTION DATA
  ###======================================================================
  #1.1 Read Kildir absorption data
  absdata_kildir = read.csv(absdata_kildir,
                            header = T)
  
  abs_surf_depth_kildir = colnames(absdata_kildir) #extract depth from where acquisition started
  
  colnames(absdata_kildir) =  t(absdata_kildir[1,])
  absdata_kildir = absdata_kildir[-1,] #Set correct column-names
  
  abs_surf_depth_kildir_df = data.frame("station" = colnames(absdata_kildir)[-1],
                                        "surface_depth" = substr(abs_surf_depth_kildir[-1],2,4),
                                        "boat" = "kildir") #Create surface depth data-frame
  
  #1.2 Read Saucier absorption data
  absdata_saucier = read.csv(absdata_saucier, 
                             header = T)
  abs_surf_depth_saucier = colnames(absdata_saucier) #extract depth from where acquisition started
  
  colnames(absdata_saucier) =  t(absdata_saucier[1,])
  absdata_saucier = absdata_saucier[-1,] #Set correct column-names
  
  abs_surf_depth_saucier_df = data.frame("station" = colnames(absdata_saucier)[-1],
                                         "surface_depth" = substr(abs_surf_depth_saucier[-1],2,4),
                                         "boat" = "saucier") #Create surface depth data-frame
  if (verbose == TRUE) {
    print("Absorption for all stations loaded")
  }
  
  
  #1.3 Create unified long data-format for absorption
  absdata_kildir_long = reshape2::melt(absdata_kildir, id.vars = "wavelength") #kildir
  colnames(absdata_kildir_long) = c("wave", "station", "at_w")
  absdata_kildir_long$boat = "kildir"
  
  absdata_saucier_long = reshape2::melt(absdata_saucier, id.vars = "wavelength") #saucier
  colnames(absdata_saucier_long) = c("wave", "station", "at_w")
  absdata_saucier_long$boat = "saucier"
  
  absdata_all_surf_long = rbind(absdata_kildir_long, absdata_saucier_long) #Final at_w
  absdata_all_surf_depth = rbind(abs_surf_depth_kildir_df, abs_surf_depth_saucier_df) #Final depth
  
  rm(absdata_kildir, abs_surf_depth_kildir, abs_surf_depth_kildir_df,
     absdata_saucier, abs_surf_depth_saucier, abs_surf_depth_saucier_df,
     absdata_kildir_long, absdata_saucier_long)
  
  if (verbose == TRUE) {
    print("Absorption in long format generated")
  }
 
  
  ###======================================================================
  # 2. READ BACKSCATTER DATA
  ###======================================================================
  
  #2.1 Read Kildir backscatter data
  bbdata_kildir = read.csv(bbdata_kildir,
                           header = T)
  bb_surf_depth_kildir = colnames(bbdata_kildir) #extract depth from where acquisition started
  
  colnames(bbdata_kildir) =  t(bbdata_kildir[1,])
  bbdata_kildir = bbdata_kildir[-1,] #Set correct column-names
  
  bb_surf_depth_kildir_df = data.frame("station" = colnames(bbdata_kildir)[-1],
                                       "surface_depth" = substr(bb_surf_depth_kildir[-1],2,4),
                                       "boat" = "kildir") #Create surface depth data-frame
  
  #2.2 Read Saucier backscatter data
  bbdata_saucier = read.csv(bbdata_saucier,
                            header = T)
  bb_surf_depth_saucier = colnames(bbdata_saucier) #extract depth from where acquisition started
  colnames(bbdata_saucier) =  t(bbdata_saucier[1,])
  bbdata_saucier = bbdata_saucier[-1,] #Set correct column-names
  
  bb_surf_depth_saucier_df = data.frame("station" = colnames(bbdata_saucier)[-1],
                                        "surface_depth" = substr(bb_surf_depth_saucier[-1],2,4),
                                        "boat" = "saucier") #Create surface depth data-frame
  
  if (verbose == TRUE) {
    print("Backscatter for all stations loaded")
  }
  
  
  #2.3 Create unified long data-format for backscatter
  bbdata_kildir_long = reshape2::melt(bbdata_kildir, id.vars = "wavelength") #kildir
  colnames(bbdata_kildir_long) = c("wave", "station", "bbp")
  bbdata_kildir_long$boat = "kildir"
  
  bbdata_saucier_long = reshape2::melt(bbdata_saucier, id.vars = "wavelength") #saucier
  colnames(bbdata_saucier_long) = c("wave", "station", "bbp")
  bbdata_saucier_long$boat = "saucier"
  
  bbdata_all_surf_long = rbind(bbdata_kildir_long, bbdata_saucier_long) #Final at_w
  bbdata_all_surf_depth = rbind(bb_surf_depth_kildir_df, bb_surf_depth_saucier_df) #Final depth
  
  rm(bbdata_kildir, bb_surf_depth_kildir, bb_surf_depth_kildir_df,
     bbdata_saucier, bb_surf_depth_saucier, bb_surf_depth_saucier_df,
     bbdata_kildir_long, bbdata_saucier_long)
  
  if (verbose == TRUE) {
    print("Backscatter in long format generated")
  }
  
  
  ###======================================================================
  # 3. READ Rrs DATA
  ###======================================================================
  
  #3.1 Read COP-S actual database and hyperspectral database
  cops_db =  read.csv("./data/Rb_spectral/COPS_db.csv",header = T)
  
  cops_hs =  read.csv("./data/Rb_spectral/copsdata.csv",header = T)
  
  correct_statnames = gsub( names(cops_hs)[-1], pattern = "\\.", replacement = "-") #replace "." 
  #by "-"
  colnames(cops_hs) = c("wave", correct_statnames) #assign correct station-names
  colnames(cops_db) = c("station", correct_statnames)
  
  #3.2 Retrieve COP-S acquired depths (CHECK b4 USE)
  cops_db_flip = as.data.frame(t(cops_db)) 
  colnames(cops_db_flip) = cops_db$stationID
  cops_db_flip = cops_db_flip[-1,]
  colnames(cops_db_flip) = cops_db$station
  
  #3.3 Extract shallow stations as per COP-S db 
  #[!!May not be accurate as it uses COP-S cast depth!!]
  cops_db_shallow_idx = !is.na(cops_db_flip$Bottom.depth)
  cops_db_shallow_idx = which(cops_db_shallow_idx == TRUE)
  
  cops_db_shallow = cops_db_flip[cops_db_shallow_idx,]
  cops_db_shallow$stationID = rownames(cops_db_shallow)
  cops_db_hs_shallow = cops_hs[,c(1,cops_db_shallow_idx+1)]
  
  if (verbose == TRUE) {
    print("Rrs in long format generated")
  }
  
  
  ###======================================================================
  # 4. READ desired IOPs and Rrs DATA
  ###======================================================================
  
  #4.1 Select the station-name
  stationID = station.args
  
  #4.2 absorption
  absdata_sample = absdata_all_surf_long[absdata_all_surf_long$station == 
                                           paste0("Station",stationID),]
  #browser()
  if (rlang::is_empty(absdata_sample$wave) == TRUE) {
    print(paste0("!!!!Absorption for station ",stationID, " is not found in absorption database!!!!"))
    print(paste0("!!!!ABSORPTION VECTOR WILL BE EMPTY; INTERPOLATION WILL FAIL!!!!"))
  }
  
  # #abs_data_idx = grep(stationID, colnames(absdata_kildir))
  # abs_surf_temp = data.frame("wave" = absdata_kildir$wavelength, 
  #                            "a"= (absdata_kildir[abs_data_idx]))
  # colnames(abs_surf_temp) = c('wave','a')
  absdata_sample$at_w[absdata_sample$at_w < 0] = 0
  write.csv(file = paste0("./data/Rb_spectral/surface_iops/abs_surf_", stationID, ".csv"), 
            x = absdata_sample, row.names = F, quote = F)
  
  if (verbose == TRUE) {
    print(paste0("Absorption for station ", stationID, " is exported"))
  }
  
  
  #4.3 Backscatter
  bbdata_sample = bbdata_all_surf_long[bbdata_all_surf_long$station == 
                                         paste0("Station",stationID),]
  
  if (rlang::is_empty(bbdata_sample$wave) == TRUE) {
    print(paste0("!!!!backscatter for station ",stationID, " is not found in backscatter database!!!!"))
    print(paste0("!!!!BACKSCATTER VECTOR WILL BE EMPTY; INTERPOLATION WILL FAIL!!!!"))
  }
  
  # bb_data_idx = grep(stationID, colnames(bbdata_kildir))
  # 
  # bb_surf_temp = data.frame("wave" = bbdata_kildir$wavelength, 
  #                            "bb"= (bbdata_kildir[bb_data_idx]))
  # colnames(bb_surf_temp) = c('wave','bb')
  
  write.csv(file = paste0("./data/Rb_spectral/surface_iops/bb_surf_", stationID, ".csv"), 
            x = bbdata_sample, row.names = F, quote = F)
  
  if (verbose == TRUE) {
    print(paste0("Backscatter for station ", stationID, " is exported"))
  }
  
  
  #4.4 Rrs
  rrs_data_idx = grep(paste0("^",stationID,"$"), correct_statnames)
  
  if (rlang::is_empty(rrs_data_idx) == TRUE) {
    print(paste0("!!!!Rrs for station ",stationID, " is not found in Rrs database!!!!"))
    print(paste0("!!!!ERROR INCOMING!!!!"))
    
  }
  rrs_hs_temp = as.matrix(cops_hs[(rrs_data_idx+1)])
  
  #browser()
  rrs_hs_temp_interp = approx(x=cops_hs$wave, y = rrs_hs_temp,
                              xout = wavelength, method = "linear" )$y
  
  if (verbose == TRUE) {
    print(paste0("Rrs for station ", stationID, " is exported"))
  }
  
  
  #4.5 Depth from COP-S
  
  cops_db_correct_statnames = rownames(cops_db_flip)
  
  depth_idx = grep(stationID, cops_db_correct_statnames)
  
  if (rlang::is_empty(depth_idx) == TRUE) {
    print(paste0("!!!Depth for station ",stationID, " is not found in Depth database!!!"))
    print(paste0("!!!!ERROR INCOMING!!!!"))
    
  }
  
  depth_input = as.numeric(as.character(cops_db_flip$Bottom.depth[depth_idx]))
  
  if (verbose == TRUE) {
    print(paste0("Depth for station ", stationID, " is exported"))
  }
  
  
  return(list("a_data" = absdata_sample, "bb_data" = bbdata_sample, "zB_COPS" = depth_input,
              "Rrs_0p" = rrs_hs_temp_interp))
  
}

#=======================================================================================
get_in_situ_params <- function(station_name = "OUT-F18", use_bb_nup = TRUE) {
  
  statname = station_name
  #Load WISE-Man BGC data
  bgc_wiseman <- read.csv("./data/biogeochemistry_wiseman.csv",header = T)
  
  a_d <- read.csv("./data/ad_long_wiseman.csv", header = T)
  a_p <- read.csv("./data/ap_long_wiseman.csv", header = T)
  a_g <- read.csv("./data/ag_long_wiseman.csv", header = T)
  a_g <- subset(a_g,wavelength >= 290)
  
  #Create data structure
  a_t <- a_p$ap + a_g$ag
  a_t_long <- data.frame("station"= a_d$station,
                         "depth"=a_d$depth,"wavelength"=a_d$wavelength, 
                         "a_phi"=a_p$ap - a_d$ad,
                         "a_cdom"=a_g$ag, "a_nap"=a_d$ad,
                         "a_total"=a_t)
  
  rm(a_t)
  
  invivo_abs = a_t_long[a_t_long$station == statname & a_t_long$depth == 0 & 
                          a_t_long$wavelength == "443",]
  
  invivo_chl = bgc_wiseman$Chl[bgc_wiseman$station == statname & bgc_wiseman$depth < 1]
  
  #Load relevant backscatter
  if (use_bb_nup == TRUE ) {
    
    #invisible(capture.output(test_IOP_func = suppressWarnings(read.surface.IOPs.wise(station.args = statname))))  
    test_IOP_func = suppressWarnings(read.surface.IOPs.wise(station.args = statname, verbose = F))
    bbdata = test_IOP_func$bb_data
    
    #### Compute bbp and bb spectral slope
    if (length(bbdata$wave) == 6) {
      sensor_wl = c(394, 420, 470, 532, 620, 700)
      x = 555/sensor_wl
      print("HS-6 VSF is used")
    }
    if (length(bbdata$wave) == 9) {
      sensor_wl = c(412, 440, 488, 510, 532, 595, 650, 676, 715)
      x = 555/sensor_wl
      print("BB-9 VSF is used")
    }
    
    #nz=length(IOP.fitted.down$Depth)
    nz=1 #As we only are interested for surface bb
    bbP555.down =rep(0,nz)
    nuP.down    =rep(0,nz)
    
    for (i in 1:nz) {
      #y = IOP.fitted.down$HS6$bbP[i,]
      y = as.numeric(as.matrix(bbdata$bbp)[,i])
      y[y < 0] = NA
      y[y > 0.1] = NA
      if (any(is.na(y))) {
        bbP555.down[i]=NA
        nuP.down[i]=NA
      } else {
        model = nls(y~b*x^z, start = list(b = y[3], z = 1),data=data.frame(x,y), 
                    control=list(maxiter=100, warnOnly=T))
        bbP555.down[i]=coef(model)[1]
        nuP.down[i]=coef(model)[2]
      }
    }
  } else {
    print("WARNING!!! ONLY USE WHEN QC CHECKED DEPTH PROFILE IS AVAILABLE")
    
    bbpath = list.files(path = "./data/Rb_spectral/backscatter/", 
                        pattern = "2019*", full.names = T, )
    idx = grep(statname, bbpath)
    if (is.na(idx)) {
      print("The backscatter csv file is not available on disc, use generate.IOP.DB.funcmode.R function to generate compatiable bb file-structure")
    }
    bbdata = read.csv(bbpath[idx], header = T)
    depth_bbdata = as.numeric(substr(colnames(bbdata)[-1],start = 2, stop = 20))
    bbdata_trim = bbdata[,-1]
    
    if (nrow(bbdata) == 6) {
      bbdata$wl = c(394, 420, 470, 532, 620, 700)
    }
    if (nrow(bbdata) == 9) {
      bbdata$wl = c(412, 440, 488, 510, 532, 595, 650, 676, 715)
    }
    #extrapolate bb for z=0 depth and desired waves
    bb = matrix(nrow=length(depth_bbdata), ncol=length(wavelength))
    
    for (i in 1:length(depth_bbdata)) {
      s     <- smooth.spline(bbdata$wl, bbdata_trim[,i])
      bb[i,] <- predict(s,wavelength)$y
    }
    
    bb = as.data.frame(t(bb))  
    colnames(bb) = depth_bbdata
    bb$wavelength = wavelength
    
    bb <- bb %>% dplyr::select(wavelength, everything())
    optimal_depth = which.min(abs(0 - depth_bbdata))
    optimal_wave = which.min(abs(550 - bb$wavelength))
    bbP555.down = bb[optimal_wave,optimal_depth+1]
  }
  
  return(list("abs_invivo"= invivo_abs, "chl_invivo" = invivo_chl, "bbp555"=bbP555.down))
  
}
