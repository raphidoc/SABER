#--------------------------------------------------------------------------------
##Test the concentration-scaled model to obtain a_chl with non-scaled model
#--------------------------------------------------------------------------------
lambda = seq(400,800,10)
#Create a random sample of [chl]
C_ph_vec = rweibull(n=500, scale = 4, shape=2)
#hist(C_ph_vec)

abs_ph_actual_vec = matrix(0, ncol = length(wavelength), nrow = length(C_ph_vec))
abs_ph_norm_vec =  matrix(0, ncol = length(wavelength), nrow = length(C_ph_vec))

#Obtain concentration-scaled abs_chl
for (i in 1:length(C_ph_vec)) {
  
  aph_440 <- 0.06*(C_ph_vec[i])^0.65# [mg/m^3]
  
  abs_ph <- rep(0,length(lambda))
  
  # Compute the value of plankton absorption as function of the concentration and wavelength
  for (j in 1:length(lambda)){
    
    abs_ph[j] <- (a0[j] + a1[j]*log(aph_440))*aph_440
    
  }
  abs_ph[abs_ph < 0] <- 0
  
  abs_ph_actual_vec[i,] = abs_ph
}



#Obtain concentration non-scaled abs_chl
for (i in 1:length(C_ph_vec)) {
  
  aph_440 <- 1# [mg/m^3]
  abs_ph_norm <- rep(0,length(lambda))
  
  # Compute the value of plankton absorption as function of the concentration and wavelength
  for (j in 1:length(lambda)){
    
    abs_ph_norm[j] <- (a0[j] + a1[j]*log(aph_440))*aph_440
    
  }
  
  C_ph_norm=0.06*(C_ph_vec[i])^0.65
  abs_ph = C_ph_norm*abs_ph_norm
  abs_ph[abs_ph < 0] <- 0
  abs_ph_norm_vec[i,] = abs_ph
}

#Test the scatter
plot(abs_ph_actual_vec, abs_ph_norm_vec, pch=19, 
     xlab="Actual a_ph", ylab="Predicted a_ph",
     col="navyblue")
abline(0,1, col="red", lwd=3)

##--------------------------------------------------
abs_chl_actual_df = data.frame(abs_ph_actual_vec)

colnames(abs_chl_actual_df) = wavelength
abs_chl_actual_df = as.data.frame(t(abs_chl_actual_df))

abs_chl_actual_df$wave = wavelength
abs_chl_actual_df <- abs_chl_actual_df %>% dplyr::select(wave, everything())

##--------------------------------------------------

abs_chl_norm_df = data.frame(abs_ph_norm_vec)

colnames(abs_chl_norm_df) = wavelength
abs_chl_norm_df = as.data.frame(t(abs_chl_norm_df))

abs_chl_norm_df$wave = wavelength
abs_chl_norm_df <- abs_chl_norm_df %>% dplyr::select(wave, everything())
##--------------------------------------------------

abs_chl_actual_df_long = reshape2::melt(abs_chl_actual_df, id.vars = "wave")
abs_chl_actual_df_long$type = "actual"

abs_chl_norm_df_long = reshape2::melt(abs_chl_norm_df, id.vars = "wave")
abs_chl_norm_df_long$type = "normalized"

abs_chl_long = rbind(abs_chl_actual_df_long, abs_chl_norm_df_long)

plot(abs_chl_long$value[abs_chl_long$type == "actual"],
     abs_chl_long$value[abs_chl_long$type == "normalized"])

#Starting values from linear coefficients from log-transformed data for both Y and X
yx.lm = lm(log(abs_chl_long$value[abs_chl_long$type == "normalized"]) ~ log(abs_chl_long$value[abs_chl_long$type == "actual"]))
summary(yx.lm)

ydata = abs_chl_long$value[abs_chl_long$type == "normalized"]
xdata = abs_chl_long$value[abs_chl_long$type == "actual"]
test.df = data.frame("x"=xdata, "y"=ydata)

#Do the GWLS NLS with lm obtained starting values
m <- nls( y ~ a*(x^b), test.df, 
         start = list(a = exp(coef(yx.lm)[1]), b = coef(yx.lm)[2])) # power formula: y = a*x^b
summary(m)




m[[1]]$formula()
m[[1]]$getPars()


#TesT the FIT
plot(test.df$y ~ test.df$x)
curve(predict(m, newdata = data.frame(x = x)), add = TRUE, col="red")

#Do the Plot
test.df$y_est = m[[1]]$getPars()[[1]] * test.df$x^m[[1]]$getPars()[[2]]
plot(test.df$x, test.df$y_est)

abs_ph_est_vec = matrix(0, ncol = length(wavelength), nrow = length(C_ph_vec))

for (i in 1:length(C_ph_vec)) {
  abs_ph_est_vec[i,] = m[[1]]$getPars()[[1]] * abs_ph_actual_vec[i,]^m[[1]]$getPars()[[2]]
}


plot(wavelength, abs_ph_actual_vec[1,])
lines(wavelength, abs_ph_norm_vec[1,], col="blue", lwd=3)
lines(wavelength, abs_ph_est_vec[1,], col="goldenrod2", lwd=3)

errorb <- qnorm(0.975)*sd(acdom_data$acdom412_est)/sqrt(length(acdom_data$acdom412_est))


# #Create the stat_function to draw the curve (!!!!!!NOT FINISHED!!!!!!)
# 
# powerFun <- function(x) m[[1]]$getPars()[[1]]*x^m[[1]]$getPars()[[2]]
# 
# legend_title <- element_blank()
# legend_position <- c(0.65, 0.50)
# 
# xmin <- -0.5; xmax <- 15; xstp <- xmax/5
# ymin <- -0.5; ymax <- 6; ystp <- ymax/5
# 
# xlabel <- expression(paste(italic("R")["rs"]("0"^"+", frac(490,667)),italic(" in situ"), " (sr"^-1," )"))
# ylabel <- expression(paste(italic("a")["CDOM"]("412"),italic("in vivo"), " (m"^-1," )"))
# 
# 
# xlbl <- xlabel; ylbl <-ylabel
# #ymin <- 0; ymax <- 0.010 ; ystp <- ymax/5
# 
# asp_rat <- (xmax-xmin)/(ymax-ymin)
# 
# g <- ggplot(acdom_data,aes(x = band_ratio)) + 
#   
#   geom_point(aes(y = acdom412, 
#                  #shape=factor(cruise), 
#                  fill=factor(cruise), 
#                  color = factor(cruise)), 
#              size=3.5, na.rm = FALSE, show.legend = T) +
#   
#   
#   
#   geom_ribbon(aes(ymin = acdom412_est - errorb,
#                   ymax = acdom412_est + errorb),
#               alpha = 0.6,
#               fill = "grey70",
#               colour="NA"
#   )+
#   
#   stat_function(fun = powerFun, size=1.3, col="#4e0707")+
#   geom_hline(yintercept =0, linetype=2, size=1.1, col="grey")+
#   geom_vline(xintercept =0, linetype=2, size=1.1, col="grey")+
#   #scale_shape_manual(name = "Cruise",labels = c("HB2018", "NE2006"), values = c(21, 24)) +
#   scale_color_manual(name = "Cruise",labels = c("HB2018", "NE2006"), values = c("navyblue", "goldenrod2")) +
#   scale_fill_manual(name = "Cruise",labels = c("HB2018", "NE2006"), values = c("navyblue", "goldenrod2")) +
#   labs(shape="cruise")+
#   
#   coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
#               ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#   scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
#                      breaks = seq(xmin, xmax, xstp)) +
#   scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
#                      breaks = seq(ymin, ymax, ystp)) +
#   guides(shape = guide_legend(title = NULL, order = 1, 
#                               keyheight = unit(0.45, 'inch')), 
#          colors = guide_legend(title = NULL, order = 2, label.vjust = 1, 
#                                keyheight = unit(0.2, 'inch'))) +
#   
#   theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         axis.text.x = element_text(size = 15, color = 'black', angle = 0), 
#         axis.text.y = element_text(size = 15, color = 'black', angle = 0), 
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         axis.ticks.length = unit(.25, "cm"),
#         legend.box.just = "right",
#         legend.spacing = unit(-0.5, "cm"),
#         legend.position = legend_position,
#         legend.direction = "vertical",
#         legend.title = element_text(colour = "black", size = 18, face = "plain"),
#         legend.text = element_text(colour = "black", size = 15, face = "plain"),
#         legend.background = element_rect(fill = NA, size = 0.5, 
#                                          linetype = "solid", colour = 0),
#         legend.key = element_blank(),
#         legend.justification = c("left", "top"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour = "grey", 
#                                         size = 0.5, linetype = "dotted"), 
#         panel.grid.minor = element_blank(),
#         plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
# g



#--------------------------------------------------------------------------------
##Test the concentration-scaled model to obtain a_CDOM(440) with non-scaled model
#--------------------------------------------------------------------------------
slope.parametric = T

abs_CDM_440 <-  1 # [1/m], CDOM abs. coeff. at 440 [nm] normalized

abs_CDM_norm_vec <- rep(0,length(lambda)); abs_CDM_norm <- rep(0,length(lambda))

for (i in 1:length(lambda)){
  
  abs_CDM_norm[i]  <- abs_CDM_440*exp(-S_CDM*(lambda[i] - 440))
  
}
abs_CDM_norm_vec = Ga_CDOM*abs_CDM_norm

plot(wavelength, abs_CDM_norm_vec)

Ga_CDOM <- 1# [m^2/mg]
Oa_CDOM <- 0

if (slope.parametric == TRUE) {
  #parametric formula to retrieve spectral slope of CDOM + NAP
  S_CDM = 0.015 + (0.002/(0.6 + (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))])))
  print(paste0("The spectral slope for CDOM + NAP is calculated as:", S_CDM))
  
} else {
  S_CDM <- 0.014 #<< USER INPUT >>
}

abs_CDM_440 <-  (Ga_CDOM*base_CDM)+Oa_CDOM# [1/m], CDOM abs. coeff. at 440 [nm]

abs_CDM <- rep(0,length(lambda))

for (i in 1:length(lambda)){
  
  abs_CDM[i]  <- abs_CDM_440*exp(-S_CDM*(lambda[i] - 440))
  
}
lines(wavelength, abs_CDM)