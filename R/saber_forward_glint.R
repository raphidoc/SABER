

saber_forward_surface <- function( calc_glint = T,
                                   rrs_subsurface, #Rrs values obtained for subsurface (0-) from
                                   # Saber_forward_fast.R
                                   
                                # Sun-Sensor Geometry above surface
                                view_above = 15, #viewing angle of sensor above water surface
                                sun_above = 30, #sun angle of sensor above water surface
                                
                                # Water surface roughness
                                q_surf = pi, #[sr]
    
                                # Atmospheric conditions
                                # Irradiance intensities [1/sr]
                                g_dd=0.05, g_dsr=0, g_dsa=0,
                                
                                # Intensities of light sources 
                                f_dd= 1, f_ds= 1,
                                
                                # Angstrom exponent
                                alpha = 1.317,
                                
                                # Atmospheric pressure 
                                P = 1013.25, # [mbar]
                                
                                # Relative Humidity
                                RH = 0.60,
                                
                                # Scale height for ozone
                                Hoz = 0.300, # [cm]
                                
                                # Scale height of the precipitate water in the atmosphere
                                WV= 2.500, # [cm])
                                
                                lambda = wavelength,
                                
                                plot_rrs = T #TRUE if want to plot the modeled Rrs
  ) {
  #--------------------------------------------------------------------------
  ## Remote sensing reflectance  (L_sky/Ed) above the surface
  #--------------------------------------------------------------------------
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### Calculation of Sky reflectance : START #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  
  
  # Extraterrestrial solar irradiance [mW/m^2 nm]
  E0 <-  read.table("./data/input-spectra/E0.txt", header = F)
  E0.wavelength <-  E0$V1
  E0.ett <-  E0$V2
  E0 <- rep(length(lambda), 0)
  
  E0 <-  Hmisc::approxExtrap(x =E0.wavelength, y =E0.ett, xout = lambda, method = "linear")$y
  
  # Oxygen absorption [1/cm]
  absO2 <- read.table("./data/input-spectra/absO2.A", header = F)
  absO2.wavelength <-  absO2$V1
  absO2.oxy <-  absO2$V2
  abs_O2 <- rep(length(lambda), 0)
  
  abs_O2 <-  Hmisc::approxExtrap(x =absO2.wavelength, y =absO2.oxy, 
                                 xout = lambda, method = "linear")$y
  
  # Ozone absorption [1/cm]
  absO3 <- read.table("./data/input-spectra/absO3.A", header = F)
  absO3.wavelength <-  absO3$V1
  absO3.oxy <-  absO3$V2
  abs_O3 <- rep(length(lambda), 0)
  
  abs_O3 <-  Hmisc::approxExtrap(x =absO3.wavelength, y =absO3.oxy, 
                                 xout = lambda, method = "linear")$y
  
  # Water vapour absorption [1/cm]
  absWV <- read.table("./data/input-spectra/absWV.A", header = F)
  absWV.wavelength <-  absWV$V1
  absWV.wv <-  absWV$V2
  abs_WV <- rep(length(lambda), 0)
  
  abs_WV <-  Hmisc::approxExtrap(x =absWV.wavelength, y =absWV.wv, 
                                 xout = lambda, method = "linear")$y
  
  # angles from deg to rad
  geometry_above = snell_law(view = view_above, sun = sun_above)
  
  view_rad=geometry_above$view_w # rad
  sun_rad=geometry_above$sun_w   # rad
  rho_L_above = geometry_above$rho_L #fresnel reflectance
  
  cat(paste0("\033[0;34m**************************************************************************\033[0m","\n"))
  cat(paste0("\033[0;34m Sun-Sensor Geometry ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Input Viewing Zenith: ",view_above," degrees ===> ", signif(view_rad, digits = 2) , " radians \033[0m","\n"))
  cat(paste0("\033[0;32m","Input Solar Zenith: ",sun_above," degrees ===> ", signif(sun_rad, digits = 2) , " radians \033[0m","\n"))
  cat(paste0("\033[0;34m**************************************************************************\033[0m","\n"))
  
  
  # Downwelling Irradiance [mW/m^2 nm]
  M <- 1/(cos(sun_rad)+(0.50572*((90+ 6.079975-sun_rad)^(-1.253))))
  M1 <-  (M*P)/1013.25
  Moz <-  1.0035/ (((cos(sun_rad)^2)+0.007)^0.5)
  
  # Air mass type
  if (type_case_water == "1") { 
    AM <- 1
  } else{
    
    AM <- 1
  }
  
  omega_a <- ((-0.0032*AM) + 0.972)*exp(RH*3.06*(10^-4))
  
  Ha <- 1; V <- 15# [km] aerosol scale height and horizontal visibility
  beta <- 3.91*(Ha/V)
  tau_a <-  beta*((lambda/550)^(-alpha))
  
  Tr <-   exp(-M1/((115.6406*(lambda^(4)))-(1.335*(lambda^(2)))))
  Taa <-  exp(-(1-omega_a)*tau_a*M)
  Tas <-  exp(-omega_a*tau_a*M)
  Toz <-  exp(-abs_O3*Hoz*Moz)
  To <-   exp((-1.41*abs_O2*M1)/((1+(118.3*abs_O2*M1))^0.45))
  Twv <-  exp((-0.2385*abs_WV*WV*M)/((1+(20.07*abs_WV*WV*M))^0.45))
  
  
  B3 <- 0.82 - (0.1417*alpha)
  B1 <-  B3*(1.459 +(B3*(0.1595+(0.4129*B3))))
  B2 <-  B3*(0.0783 +(B3*(-0.3824-(0.5874*B3))))
  
  Fa <-  1-(0.5*exp((B1+(B2*cos(sun_rad)))*cos(sun_rad)))
  
  Edd <-  E0*Tr*Taa*Tas*Toz*To*Twv*cos(sun_rad)
  Edsr <- (1/2)*E0*(1-(Tr^(0.95)))*Taa*Tas*Toz*To*Twv*cos(sun_rad)
  Edsa <-  E0*(Tr^(1/2))*Taa*(1-Tas)*Toz*To*Twv*cos(sun_rad)*Fa
  Eds <-  Edsr+Edsa# diffuse downwelling irradiance (sum of Rayleigh and aerosol)
  
  Ed <-  (f_dd*Edd) + (f_ds*Eds)# downwelling irradiance
  
  # Sky Radiance
  Ls <-  (g_dd*Edd) + (g_dsr*Edsr) + (g_dsa*Edsa)# [mW/sr m^2 nm]
  
  # Remote sensing reflectance above the surface [1/sr]
  
  Rrs_above <-  rho_L_above*(Ls/Ed)
  
  cat(paste0("\033[0;34m Sky Glint equivalent Rrs calculated ::::\033[0m","\n"))
  
  #--------------------------------------------------------------------------
  ## Final Remote sensing reflectance computation
  #--------------------------------------------------------------------------
  sigma <- 0.03; nW <- 1.33; rho_U <- 0.54; Q <- q_surf#[sr]
  
  rrs_surf = (((1-sigma)*(1-rho_L)/(nW^2))* (rrs_subsurface/(1-rho_U*Q*rrs_subsurface)))
  
  cat(paste0("\033[0;32m Sub-surface Rrs (0-) is translated to surface (0+)  ::::\033[0m","\n"))

  if (calc_glint == TRUE) {
    
    Rrs <-  rrs_surf+ Rrs_above
    
    cat(paste0("\033[0;32m Above surface (0^+) Rrs with glint is calculated  ::::\033[0m","\n"))
    
    
  } else {
    
    Rrs <-  rrs_surf
    
    cat(paste0("\033[0;32m Above surface (0^+) Rrs without glint is calculated  ::::\033[0m","\n"))
  }
  
  rrs_df = data.frame("wave" = lambda, "rrs_subsurf" = rrs_subsurface, "rrs_glint" = Rrs_above,
                      "rrs_surf" = Rrs)
  
  if (plot_rrs == TRUE) {
    
    xmin = min(rrs_df$wave); xmax = max(rrs_df$wave); xstp = 100
    ymin = 0 
    
    #Rrs above and under
    ymax = signif(max(rrs_df$rrs_subsurf) + 0.2*max(rrs_df$rrs_subsurf), digits = 2)
    
    ystp = signif(ymax/4, digits = 2)
    
    asp_rat <- (xmax-xmin)/(ymax-ymin)
    
    g_rrs<- ggplot(rrs_df, aes(x=wave)) +
      #geom_line(aes(y=rrs_glint, col="xx1"), lwd = 1.3) +
      geom_line(aes(y=rrs_subsurf, col="xx2"), lwd = 1.3) +
      geom_line(aes(y=rrs_surf, col = "xx3"), lwd = 1.3) +
      
      scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,
                                                ")[", sr^-1,"]")) , 
                         limits = c(ymin, ymax),
                         breaks = seq(ymin, ymax, ystp))+ 
      
      scale_colour_viridis(discrete = T,
                           labels = c(#expression(paste(italic("R"),{}[rs],"(glint)")),
                             expression(paste(italic("R"),{}[rs],"(",0^"-",")")),
                             expression(paste(italic("R"),{}[rs],"(",0^"+",")"))
                                                 ) 
      ) +
      
      scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
                         limits = c(xmin, xmax), 
                         breaks = seq(xmin, xmax, xstp))  +
      
      coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                  ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
      
      theme_bw()+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
            axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
            axis.title.x = element_text(size = 25),
            axis.title.y = element_text(size = 25),
            axis.ticks.length = unit(.25, "cm"),
            legend.position=c(0.62, 0.905),
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
            #legend.spacing.y = unit(2.0, 'cm'),
            plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
            legend.text.align = 0,
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
    
    ###### GLINT
    ymin = signif(min(rrs_df$rrs_glint), digits = 5)
    ymax = signif(max(rrs_df$rrs_glint) + 0.2*max(rrs_df$rrs_glint), digits = 5)
    
    ystp = signif(ymax/4, digits = 5)
    
    asp_rat <- (xmax-xmin)/(ymax-ymin)
    
    g_glint<- ggplot(rrs_df, aes(x=wave)) +
      geom_line(aes(y=rrs_glint, col="xx1"), lwd = 1.3) +
      
      scale_y_continuous(name =" " 
                         # , 
                         # limits = c(ymin, ymax),
                         # breaks = seq(ymin, ymax, ystp)
                         )+ 
      
      scale_colour_viridis(discrete = T,
                           labels = c(expression(paste(italic("R"),{}[rs],"(glint)"))
                             
                           ) 
      ) +
      
      scale_x_continuous(name = expression(paste(lambda)), 
                         limits = c(xmin, xmax), 
                         breaks = seq(xmin, xmax, xstp))  +
      
      # coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
      #             ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
      
      theme_bw()+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
            axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
            axis.title.x = element_text(size = 25),
            axis.title.y = element_text(size = 25),
            axis.ticks.length = unit(.25, "cm"),
            legend.position=c(0.62, 0.905),
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
            #legend.spacing.y = unit(2.0, 'cm'),
            plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
            legend.text.align = 0,
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
      
    print(g_rrs + g_glint)
    
  }
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### Calculation of Sky reflectance : FINISH #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  
  return(rrs_df)
}
