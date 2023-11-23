
#Function to trnaslate above water Sun Zenith angle to sub-surface Zenith angle                            
sunzen_below = function(sun_zen_aove=45){
  sun_zen_below_rad = asin(sin(sun_zen_aove*(pi/180))/1.333)
  sun_zen_below_deg = sun_zen_below_rad*(180/pi)
  return(sun_zen_below_deg)
}

calc_Ed_analytic <- function(sunzen_Ed = 60, lat_Ed = 49, lon_Ed = -68,
                             date_time_Ed = "2019-08-18 20:50 GMT",
                             
                             wavelength=seq(400,800,10), plot_ed = T){
  
  print("SICF: NOTE: Gregg & Carder 1990 is used to estimate Ed")
  
  Cops::GreggCarder.data()
  
  library(lubridate)
  library(dplyr)
  
  x = as_datetime(as.POSIXct(date_time_Ed))
  
  jday_no = yday(x) #Get Jullien Day no.
  time_no = format(x, "%T") #Get time in UTC
  
  time_dec = sapply(strsplit(as.character(time_no),":"), #Convert time in decimal format
                    function(x) {
                      x <- as.numeric(x)
                      x[1]+x[2]/60
                    }
  )
  
  if (sunzen_Ed < 0) {
    print("SICF: Sun Zenith was not provided, calculated from Geometry")
    sunzen_Ed = Cops::GreggCarder.sunang(rad = 180/pi, iday = jday_no, 
                                         xlon = lon_Ed, ylat = lat_Ed, hr = time_dec)
  }
  
  tryCatch({
    #Calculate the Ed following Gregg & Carder 1990
    test_Ed = GreggCarder.f.modified(the = sunzen_Ed, 
                                     lam.sel = lambda, hr = time_dec,
                                     jday = jday_no, rlon = lon_Ed, rlat = lat_Ed, debug = T)
    
    if (all(is.na(test_Ed))) {
      cat(paste0("Sun is below horizon with given sun-earth geometry, reseting to default"))
      
      sunzen_Ed = -999; lat_Ed = 49; lon_Ed = -68;
      date_time_Ed = "2019-08-18 20:50 GMT"
      
      
      library(lubridate)
      library(dplyr)
      
      x = as_datetime(as.POSIXct(date_time_Ed))
      
      jday_no = yday(x) #Get Jullien Day no.
      time_no = format(x, "%T") #Get time in UTC
      
      time_dec = sapply(strsplit(as.character(time_no),":"), #Convert time in decimal format
                        function(x) {
                          x <- as.numeric(x)
                          x[1]+x[2]/60
                        }
      )
      
      test_Ed = GreggCarder.f.modified(the = sunzen_Ed, 
                                       lam.sel = lambda, hr = time_dec,
                                       jday = jday_no, rlon = lon_Ed, rlat = lat_Ed, debug = F)
    } }, 
    error = function(e){cat("ERROR in Irradiance model :",conditionMessage(e), "\n")}
  )
  
  
  Ed0.0p = test_Ed$Ed #Total Ed at 0+
  Ed0_dir.0p = test_Ed$Edir #Direct Ed at 0+
  Ed0_dif.0p = test_Ed$Edif #Diffused Ed at 0+
  
  rhoF <- GreggCarder.sfcrfl(rad = 180.0/pi, theta=sunzen_Ed, ws=5) #Calculate fresnel reflectance
  sun_view_op = snell_law(view = 0 , sun = sunzen_Ed) #calculate sun angle 0+ in radian
  
  #Translate Ed0+ to Ed0-
  Ed0.0m = (Ed0_dir.0p *(1 - rhoF$rod)) + # direct
    (Ed0_dif.0p * (1 - rhoF$ros))   # diffuse
  
  #Convert it to scalar irradiance
  sun_view_om = snell_law(view = 0 , sun = sunzen_below(sunzen_Ed))
  Ed0m_0 = (Ed0.0m) / cos(sun_view_om$sun_w)
  
  #Plot the simulated Ed profiles
  if (plot_ed == TRUE) {
    plot(lambda, Ed0.0p, col="black", pch=19, ylab = "Ed", ylim=c(0,1.25*max(Ed0m_0)))
    lines(lambda, Ed0.0m, col="red", lwd=2.5, lty="dashed")
    lines(lambda, Ed0m_0, col="navyblue", lwd=2.5)
    
    legend(x = "bottom",          # Position
           legend = c("Ed0+","Ed0-", "E0"),  # Legend texts
           lty = c(1, 2, 1),           # Line types
           col = c("black", "red", "navyblue"), # Line colors
           lwd = 2, bty = "n")
    
  }
  return(list("sunzen_above_w" = sun_view_op$sun_w, "sunzen_below_w" = sun_view_om$sun_w,
              "Ed_above_w" = Ed0.0p, "Ed_below_w" = Ed0.0m, "E0_below_w" = Ed0m_0))
}


wise_datetime = "2019-08-18 20:03:56Z"

# Convert string to POSIXct object
utc_date_time <- as.POSIXct(wise_datetime, format = "%Y-%m-%d %H:%M:%S")

# Format the POSIXct object as a string with the desired format
wise_datetime <- format(utc_date_time, "%Y-%m-%d %H:%M GMT")

wise_Lon=-68.32577;  wise_Lat=49.32577

wise_Ed = calc_Ed_analytic(sunzen_Ed = -999, lat_Ed = wise_Lat, lon_Ed = wise_Lon, 
                           date_time_Ed = wise_datetime, plot_ed = T)




