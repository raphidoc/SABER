

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
  "actual" = iop_param_qc$bb555[bbp_idx], "predicted" = fit_param_saber$bbp555[bbp_idx], 
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

xmin <- 10^-3; xmax <- 10^-1; xstp <- xmax/4
xlbl <- expression(paste(italic("b")["bp"](555),italic("actual"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("b")["bp"](555),italic("predicted"), " [", "m"^-1, "]"))

# bbp_ioccg = plot_inversion_validation(input_df = inv_bbp, xmin = xmin, xmax = xmax,
#                                       xlabel =  xlbl, ylabel = ylbl)

bbp_ioccg = plot_inversion_validation_singlevar_contour(input_df = inv_bbp, xmin = xmin,
                                                        xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl, uncertainty = "CI", 
                                      opacity = 0.25, plot_col = "#31688EFF",
                                      show_legend = T)


ggsave(paste0("./outputs/bbp_nomad_1.png"), plot = bbp_ioccg,
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

# chl_ioccg = plot_inversion_validation(input_df = inv_chl, xmin = xmin, xmax = xmax,
#                                       xlabel =  xlbl, ylabel = ylbl)

chl_ioccg = plot_inversion_validation_singlevar_contour(input_df = inv_chl, xmin = xmin, 
                                                        xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl,uncertainty = "CI", 
                                      opacity = 0.3, plot_col = "#6DCD59FF", show_legend = F)

ggsave(paste0("./outputs/chl_nomad_1.png"), plot = chl_ioccg,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#===========================================
# 3. adg(443)
#===========================================
inv_adg = data.frame( 
  "actual" = iop_param_qc$adg443, "predicted" = fit_param_saber$adg443, 
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

xmin <- 10^-3; xmax <- 10^1; xstp <- xmax/4
xlbl <- expression(paste(italic("a")["dg"](443),italic("actual"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("a")["dg"](443),italic("predicted"), " [", "m"^-1, "]"))

# adg_ioccg = plot_inversion_validation(input_df = inv_adg, xmin = xmin, xmax = xmax,
#                                       xlabel =  xlbl, ylabel = ylbl)

adg_ioccg = plot_inversion_validation_singlevar_contour(input_df = inv_adg,
                                                        xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl, 
                                      uncertainty = "CI", opacity = 0.3, plot_col = "#FDE725FF", 
                                      show_legend = FALSE)

ggsave(paste0("./outputs/adg_nomad_1.png"), plot = adg_ioccg,
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

xmin <- 10^-3.5; xmax <- 10^0.5; xstp <- xmax/4
xlbl <- expression(paste(italic("b")["bp"](555),italic("actual"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("b")["bp"](555),italic("predicted"), " [", "m"^-1, "]"))

# bbp_ioccg = plot_inversion_validation(input_df = inv_bbp, xmin = xmin, xmax = xmax,
#                                       xlabel =  xlbl, ylabel = ylbl)

bbp_ioccg = plot_inversion_validation_singlevar_contour(input_df = inv_bbp, xmin = xmin, 
                                                        xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl, uncertainty = "CI", 
                                      opacity = 0.3, plot_col = "#31688EFF", show_legend = F)

ggsave(paste0("./outputs/bbp_ioccg_1.png"), plot = bbp_ioccg,
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

# chl_ioccg = plot_inversion_validation(input_df = inv_chl, xmin = xmin, xmax = xmax,
#                                       xlabel =  xlbl, ylabel = ylbl)

chl_ioccg = plot_inversion_validation_singlevar_contour(input_df = inv_chl, 
                                                        xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl, plot_col = "#6DCD59FF",
                                      uncertainty = "CI", opacity = 0.3, show_legend = F)

ggsave(paste0("./outputs/chl_ioccg_1.png"), plot = chl_ioccg,
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

xmin <- 10^-3; xmax <- 10^1; xstp <- xmax/4
xlbl <- expression(paste(italic("a")["dg"](443),italic("actual"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("a")["dg"](443),italic("predicted"), " [", "m"^-1, "]"))

# adg_ioccg = plot_inversion_validation(input_df = inv_adg, xmin = xmin, xmax = xmax,
#                                       xlabel =  xlbl, ylabel = ylbl)

adg_ioccg = plot_inversion_validation_singlevar_contour(input_df = inv_adg, 
                                                        xmin = xmin, xmax = xmax,
                                      xlabel =  xlbl, ylabel = ylbl, plot_col = "#FDE725FF",
                                      uncertainty = "CI", opacity = 0.3, show_legend = F)

ggsave(paste0("./outputs/adg_ioccg_1.png"), plot = adg_ioccg,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#-----------------------------------------------------------------------------------------
#Create mosaic of the deep parameter scatterplots
#-----------------------------------------------------------------------------------------
library(png)
library(grid)
library(gridExtra)

chl_ioccg = readPNG("./outputs/chl_ioccg_1.png")
chl_nomad <- readPNG("./outputs/chl_nomad_1.png")

adg_ioccg = readPNG("./outputs/adg_ioccg_1.png")
adg_nomad <- readPNG("./outputs/adg_nomad_1.png")

bbp_ioccg = readPNG("./outputs/bbp_ioccg_1.png")
bbp_nomad <- readPNG("./outputs/bbp_nomad_1.png")

#Plot together
#grid.arrange()
tmp <- arrangeGrob(rasterGrob(chl_ioccg),rasterGrob(chl_nomad), 
                   rasterGrob(adg_ioccg), rasterGrob(adg_nomad), 
                   rasterGrob(bbp_ioccg), rasterGrob(bbp_nomad),
                   ncol=2)

ggsave('./outputs/inverse_deep_new.png',tmp,scale = 1.5, width = 8, height = 12, 
       units = "in",dpi = 300)

