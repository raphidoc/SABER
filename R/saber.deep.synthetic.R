plot_inversion_validation <- function(input_df, xmin, xmax, xlabel, ylabel){
  
  xmin = xmin; xmax = xmax; xstp <- xmax/4
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g <- ggplot(input_df,aes(x = actual)) + 
    
    stat_density_2d(aes(y = predicted, fill = ..level..),
                    geom = "polygon", colour="gray", show.legend = F, alpha = 0.5, size=1.1)+
    
    scale_fill_distiller(palette = 2, direction = 1)+
    
    geom_point(aes(y = predicted), shape=21, fill="goldenrod2", 
               size=3.0, na.rm = T, show.legend = F) +
    
    geom_ribbon(aes(ymin = predicted - CI,
                    ymax = predicted + CI),
                alpha = 0.3,
                fill = "navyblue",
                colour="NA"
    )+
    geom_abline(slope = 1,linetype="solid", intercept = 0,
                colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
    
    # geom_smooth(size=1.5,level = 0.95,show.legend = F,linetype = "dashed",
    #             color="black", data = H_df,
    #             se= T, method = "lm", aes(x=H_actual, y=H_predicted))+
    
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
    
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    xlab(xlbl)+
    ylab(ylbl)+
    annotation_logticks()+
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
  
  return(g)
}


#---------------------------------------------
### NOMAD parameter validation
#---------------------------------------------
# fit_param_saber = read.csv("./outputs/inv_param_NOMAD_qaa_initial_new.csv", header = T)
# fit_param_saber = read.csv("./outputs/inv_param_NOMAD.csv", header = T)

fit_param_qaa = read.csv("./outputs/QAA_param_NOMAD.csv", header = T)

fit_param_saber = read.csv("./outputs/inv_param_NOMAD_new.csv", header = T)
iop_param_qc = iop_nomad_qc[-c(nomad_qaa_idx),]

bbp_idx = (which(is.na(iop_param_qc$bb555) == FALSE))

#===========================================
# 1. BACKSCATTER
#===========================================
inv_bbp = data.frame( 
  "actual" = iop_param_qc$bb555[bbp_idx], "predicted" = fit_param_saber$bbp550[bbp_idx], 
  "sd" =  fit_param_saber$bbp550.sd[bbp_idx])
#Calculate Confidence Interval
cal_sd <- function(df_in,i){
  errorb <- qnorm(0.975)*sd(df_in[i,-3])/
    sqrt(length(df_in[i,-3]))
}

sd_indices <- 1:dim(inv_bbp)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, df_in = inv_bbp)) 

inv_bbp$CI = as.numeric(sd_estimate)

idx = is.na(inv_bbp$sd) 
inv_bbp$sd = qnorm(0.975)*(inv_bbp$sd*10^3)/sqrt(2) #Calculate 95% Credible Interval

inv_bbp$sd[idx] = inv_bbp$CI[idx]
inv_bbp = inv_bbp[-(1:2),]

#--------------------------------------
#Plot bbp(555) validation
#--------------------------------------
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

xmin <- 10^-3; xmax <- 10^-1.5; xstp <- xmax/4
xlbl <- expression(paste(italic("b")["bp"](555),italic("actual"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("b")["bp"](555),italic("predicted"), " [", "m"^-1, "]"))

bbp_ioccg = plot_inversion_validation(input_df = inv_bbp, xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl)

ggsave(paste0("./outputs/bbp_nomad.png"), plot = bbp_ioccg,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#===========================================
# 2. chl
#===========================================
inv_chl = data.frame( 
  "actual" = iop_param_qc$chl, "predicted" = fit_param_saber$chl, 
  "sd" =  fit_param_saber$chl.sd)

#Calculate Confidence Interval
cal_sd <- function(df_in,i){
  errorb <- qnorm(0.975)*sd(df_in[i,-3])/
    sqrt(length(df_in[i,-3]))
}

sd_indices <- 1:dim(inv_chl)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, df_in = inv_chl)) 

inv_chl$CI = as.numeric(sd_estimate)

idx = is.na(inv_chl$sd) 
inv_chl$sd = qnorm(0.975)*(inv_chl$sd)/sqrt(2) #Calculate 95% Credible Interval

inv_chl$sd[idx] = inv_chl$CI[idx]
inv_chl$sd_min = inv_chl$predicted - inv_chl$sd
inv_chl$sd_max = inv_chl$predicted + inv_chl$sd

inv_chl$CI_min = inv_chl$predicted - inv_chl$CI
inv_chl$CI_max = inv_chl$predicted + inv_chl$CI

#--------------------------------------
#Plot [chl] validation
#--------------------------------------
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

xmin <- 10^-2.5; xmax <- 10^2.5; xstp <- xmax/4
xlbl <- expression(paste("[",italic("chl"),"]",italic("actual"), " [","mg"^1,"m"^-3,"]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste("[",italic("chl"),"]",italic("predicted"), " [","mg"^1,"m"^-3,"]"))

chl_ioccg = plot_inversion_validation(input_df = inv_chl, xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl)

ggsave(paste0("./outputs/chl_nomad.png"), plot = chl_ioccg,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#===========================================
# 3. adg(443)
#===========================================
inv_adg = data.frame( 
  "actual" = iop_param_qc$adg443, "predicted" = fit_param_saber$adg.440, 
  "sd" =  fit_param_saber$adg440.sd)

#Calculate Confidence Interval
cal_sd <- function(df_in,i){
  errorb <- qnorm(0.975)*sd(df_in[i,-3])/
    sqrt(length(df_in[i,-3]))
}

sd_indices <- 1:dim(inv_adg)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, df_in = inv_adg)) 

inv_adg$CI = as.numeric(sd_estimate)

idx = is.na(inv_adg$sd) 
inv_adg$sd = qnorm(0.975)*(inv_adg$sd)/sqrt(2) #Calculate 95% Credible Interval

inv_adg$sd[idx] = inv_adg$CI[idx]
inv_adg$sd_min = inv_adg$predicted - inv_adg$sd
inv_adg$sd_max = inv_adg$predicted + inv_adg$sd

inv_adg$CI_min = inv_adg$predicted - inv_adg$CI
inv_adg$CI_max = inv_adg$predicted + inv_adg$CI

#--------------------------------------
#Plot adg(443) validation
#--------------------------------------
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

xmin <- 10^-2.5; xmax <- 10^1; xstp <- xmax/4
xlbl <- expression(paste(italic("a")["dg"](443),italic("actual"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("a")["dg"](443),italic("predicted"), " [", "m"^-1, "]"))

adg_ioccg = plot_inversion_validation(input_df = inv_adg, xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl)

ggsave(paste0("./outputs/adg_nomad.png"), plot = adg_ioccg,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)


#-------------------------------------------------
# IOCCG parameter validation
#-------------------------------------------------
fit_param_saber = read.csv("./outputs/inv_param_IOCCG_new.csv", header = T)

iop_param_qc = HL.deep.iop
iop_param_qc$adg443 = iop_param_qc$acdom440 + iop_param_qc$anap440

#===========================================
# 1. BACKSCATTER
#===========================================
inv_bbp = data.frame( 
  "actual" = iop_param_qc$bbp550, "predicted" = fit_param_saber$bbp550, 
  "sd" =  fit_param_saber$chl.sd)

#Calculate Confidence Interval
cal_sd <- function(df_in,i){
  errorb <- qnorm(0.975)*sd(df_in[i,-3])/
    sqrt(length(df_in[i,-3]))
}

sd_indices <- 1:dim(inv_bbp)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, df_in = inv_bbp)) 

inv_bbp$CI = as.numeric(sd_estimate)

idx = is.na(inv_bbp$sd) 
inv_bbp$sd = qnorm(0.975)*(inv_bbp$sd*10^3)/sqrt(2) #Calculate 95% Credible Interval

inv_bbp$sd[idx] = inv_bbp$CI[idx]


#--------------------------------------
#Plot bbp(555) validation
#--------------------------------------
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

xmin <- 10^-3.5; xmax <- 10^-0.5; xstp <- xmax/4
xlbl <- expression(paste(italic("b")["bp"](555),italic("actual"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("b")["bp"](555),italic("predicted"), " [", "m"^-1, "]"))

bbp_ioccg = plot_inversion_validation(input_df = inv_bbp, xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl)

ggsave(paste0("./outputs/bbp_ioccg.png"), plot = bbp_ioccg,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#===========================================
# 2. chl
#===========================================
inv_chl = data.frame( 
  "actual" = iop_param_qc$chl, "predicted" = fit_param_saber$chl, 
  "sd" =  fit_param_saber$chl.sd)

#Calculate Confidence Interval
cal_sd <- function(df_in,i){
  errorb <- qnorm(0.975)*sd(df_in[i,-3])/
    sqrt(length(df_in[i,-3]))
}

sd_indices <- 1:dim(inv_chl)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, df_in = inv_chl)) 

inv_chl$CI = as.numeric(sd_estimate)

idx = is.na(inv_chl$sd) 
inv_chl$sd = qnorm(0.975)*(inv_chl$sd)/sqrt(2) #Calculate 95% Credible Interval

inv_chl$sd[idx] = inv_chl$CI[idx]
inv_chl$sd_min = inv_chl$predicted - inv_chl$sd
inv_chl$sd_max = inv_chl$predicted + inv_chl$sd

inv_chl$CI_min = inv_chl$predicted - inv_chl$CI
inv_chl$CI_max = inv_chl$predicted + inv_chl$CI

#--------------------------------------
#Plot [chl] validation
#--------------------------------------
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

xmin <- 10^-2; xmax <- 10^2; xstp <- xmax/4
xlbl <- expression(paste("[",italic("chl"),"]",italic("actual"), " [","mg"^1,"m"^-3,"]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste("[",italic("chl"),"]",italic("predicted"), " [","mg"^1,"m"^-3,"]"))

chl_ioccg = plot_inversion_validation(input_df = inv_chl, xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl)

ggsave(paste0("./outputs/chl_ioccg.png"), plot = chl_ioccg,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#===========================================
# 3. adg(443)
#===========================================
inv_adg = data.frame( 
  "actual" = iop_param_qc$adg443, "predicted" = fit_param_saber$adg.440, 
  "sd" =  fit_param_saber$acdom.440.sd + fit_param_saber$anap440.sd)

#Calculate Confidence Interval
cal_sd <- function(df_in,i){
  errorb <- qnorm(0.975)*sd(df_in[i,-3])/
    sqrt(length(df_in[i,-3]))
}

sd_indices <- 1:dim(inv_adg)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, df_in = inv_adg)) 

inv_adg$CI = as.numeric(sd_estimate)

idx = is.na(inv_adg$sd) 
inv_adg$sd = qnorm(0.975)*(inv_adg$sd)/sqrt(2) #Calculate 95% Credible Interval

inv_adg$sd[idx] = inv_adg$CI[idx]
inv_adg$sd_min = inv_adg$predicted - inv_adg$sd
inv_adg$sd_max = inv_adg$predicted + inv_adg$sd

inv_adg$CI_min = inv_adg$predicted - inv_adg$CI
inv_adg$CI_max = inv_adg$predicted + inv_adg$CI

#--------------------------------------
#Plot adg(443) validation
#--------------------------------------
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

xmin <- 10^-2.5; xmax <- 10^1; xstp <- xmax/4
xlbl <- expression(paste(italic("a")["dg"](443),italic("actual"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("a")["dg"](443),italic("predicted"), " [", "m"^-1, "]"))

adg_ioccg = plot_inversion_validation(input_df = inv_adg, xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl)

ggsave(paste0("./outputs/adg_ioccg.png"), plot = adg_ioccg,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#-----------------------------------------------------------------------------------------
#Create mosaic of the deep parameter scatterplots
#-----------------------------------------------------------------------------------------
library(png)
library(grid)
library(gridExtra)

chl_ioccg = readPNG("./outputs/chl_ioccg.png")
chl_nomad <- readPNG("./outputs/chl_nomad.png")

adg_ioccg = readPNG("./outputs/adg_ioccg.png")
adg_nomad <- readPNG("./outputs/adg_nomad.png")

bbp_ioccg = readPNG("./outputs/bbp_ioccg.png")
bbp_nomad <- readPNG("./outputs/bbp_nomad.png")

#Plot together
grid.arrange()
tmp <- arrangeGrob(rasterGrob(chl_ioccg),rasterGrob(chl_nomad), 
                   rasterGrob(adg_ioccg), rasterGrob(adg_nomad), 
                   rasterGrob(bbp_ioccg), rasterGrob(bbp_nomad),
                   ncol=2)

ggsave('./outputs/inverse_deep.png',tmp,scale = 1.5, width = 8, height = 12, 
       units = "in",dpi = 300)

