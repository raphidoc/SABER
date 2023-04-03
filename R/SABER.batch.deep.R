#=================================================================================================================
# Test the S.A.B.E.R. inversion for deep water in batch mode for synthetic & in situ data-sets. The synthetic data
# used here is IOCCG 2003 data and the in situ data used here consists of GLORIA and NOMAD dataset.

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#=================================================================================================================
rm(list=ls(all=TRUE))

require(dplyr)
require(readxl)
require(stats4)
require(MASS)
require(dglm)
library(fitdistrplus)
require(Riops)
require(Cops)
library(lmodel2)

setwd("/home/musk0001/R_inverse_wasi")

#=====================================================
# Function for goodness of fit for inversion retrieved
# parameters for QAA and SABER
#=====================================================
goodness_of_fit <- function(actual, predicted, bbp.exist = F) {
  
  #BIAS
  bias_chl = Metrics::bias(actual = actual$chl, #chl
                               predicted = predicted$chl)
  
  bias_adg443 = Metrics::bias(actual = actual$adg443, #adg443
                              predicted = predicted$adg443)
  
  
  #%-bias
  p_bias_chl = Metrics::percent_bias(actual = actual$chl, #chl
                                     predicted = predicted$chl)
  
  p_bias_adg443 = Metrics::percent_bias(actual = actual$adg443, #adg443
                                        predicted = predicted$adg443)
  
  #RMSE
  rmse_chl = Metrics::rmse(actual = actual$chl, #chl
                           predicted = predicted$chl)
  
  rmse_adg443 = Metrics::rmse(actual = actual$adg443, #adg443
                              predicted = predicted$adg443)
  
  #R^2
  yx.lmodel2_chl <- lmodel2(log10(predicted$chl) ~ log10(actual$chl)) #chl
  
  r_square_chl = yx.lmodel2_chl$rsquare
  slope_chl = yx.lmodel2_chl$regression.results$Slope[2]
  
  yx.lmodel2_adg443 <- lmodel2(log10(predicted$adg443) ~ log10(actual$adg443)) #adg443
  
  r_square_adg443 = yx.lmodel2_adg443$rsquare
  slope_adg443 = yx.lmodel2_adg443$regression.results$Slope[2]
  
  if (bbp.exist == T) {
    bias_bbp555 = Metrics::bias(actual = actual$bbp555, #bbp555
                                predicted = predicted$bbp555)
    
    p_bias_bbp555 = Metrics::percent_bias(actual = actual$bbp555, #bbp555
                                       predicted = predicted$bbp555)
    
    rmse_bbp555 = Metrics::rmse(actual = actual$bbp555, #bbp555
                             predicted = predicted$bbp555)
    
    yx.lmodel2_bbp555 <- lmodel2(log10(predicted$bbp555) ~ log10(actual$bbp555)) #bbp555
    
    r_square_bbp555 = yx.lmodel2_bbp555$rsquare
    slope_bbp555 = yx.lmodel2_bbp555$regression.results$Slope[2]
    
  }
  
  if (bbp.exist == T) {
    bias = c(bias_chl, bias_adg443, bias_bbp555)
    p_bias = c(p_bias_chl, p_bias_adg443, p_bias_bbp555)
    rmse = c(rmse_chl, rmse_adg443, rmse_bbp555)
    r_square = c(r_square_chl, r_square_adg443, r_square_bbp555)
    slope = c(slope_chl, slope_adg443, slope_bbp555)
    errorstat = data.frame(bias, p_bias, rmse, r_square, slope)
    errorstat$parameter = c("chl", "adg443", "bbp555")
    
  } else {
    bias = c(bias_chl, bias_adg443)
    p_bias = c(p_bias_chl, p_bias_adg443)
    rmse = c(rmse_chl, rmse_adg443)
    r_square = c(r_square_chl, r_square_adg443)
    slope = c(slope_chl, slope_adg443)
    errorstat = data.frame(bias, p_bias, rmse, r_square, slope)
    errorstat$parameter = c("chl", "adg443")
    
  }
  
  return(errorstat)
  
}

# ====================================================
# 1. For IOCCG 2003 dataset
# ====================================================
#Load IOP and AOP data

#water IOP
insitu.data <- read_excel_allsheets(paste0(getwd(),"/data/IOP_AOP_Sun60.xls"))
water.data <- insitu.data$Basics[6:46, 1:3]                           
names(water.data) <- c("wavelength", "a_w", "bb_w") 

#Store spectral slopes
a_g_vec = insitu.data$a_g$Sg ; a_g_HL_mean = mean(a_g_vec)
a_d_vec = insitu.data$a_d$Sd; a_d_HL_mean = mean(a_d_vec) ; a_dg_HL_mean = a_g_HL_mean + a_d_HL_mean


HL.deep.iop = read_IOCCG_data()

#sub-surface (0^-) Rrs
rrs.HL <- as.matrix(insitu.data$r_rs); rrs.HL.wl <- as.numeric(names(insitu.data$r_rs))

#---------------------------------------------
### 1.1 Inversion of IOCCG dataset using SABER
#---------------------------------------------
#Data frame to store inversion results
Fit.optimized.ssobj.batch <- data.frame("chl"=0, 
                                  "adg.440"=0,
                                  "bbp550"=0,
                                  "chl.sd"=0, "adg440.sd"=0, "bbp550.sd"=0)
for (j in 1:dim(rrs.HL)[1]) {
  
  #Initial values from pre-fit
  
  par0 = c(chl = mean(Hl.deep.iop$chl), adg440 = mean(Hl.deep.iop$acdom440 + Hl.deep.iop$anap440), 
             bbp550 = 0.025, "pop_sd" = 0.13129354/100)
  
  upper.bound <- c(30,10,0.5,0.01)  
  lower.bound <- c(0.01, 0.01, 0.001, 0.0001)
  
  ##Do the Optimization
  obj = "log-LL"
  methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")
  
  inverse_output <- suppressWarnings(solve.objective.inverse.deep.final(
                                initial = par0, 
                                bbp.constrain  = F,
                                bbp.550 = HL.deep.iop$bbp550[j],
                                
                                obsdata = rrs.HL[j,],
                                
                                sa.model = "am03", 
                                
                                obj.fn =obj.run ,
                                obj.fn ="obj_L98",
                                
                                auto_spectral_slope = F,
                                manual_spectral_slope = T, 
                                
                                manual_spectral_slope_vals = c("s_g"=a_g_vec[i], "s_d"=a_d_vec[i], "gamma"=1),
                                
                                method.opt = methods.opt[4],
                                lower.b = lower.bound,
                                upper.b = upper.bound, 
                                batch = FALSE, pop.sd = FALSE))
  
  Fit.optimized.ssobj.batch <- rbind(Fit.optimized.ssobj.batch, c(inverse_output[[1]]$estimates[1:3],
                                                      inverse_output[[1]]$`sd(+/-)`[1:3]))
  cat(paste0("\033[0;43m","############### INVERSION FINISHED for ", j, " no. spectra ################","\033[0m","\n"))
}

#Export the SABER outputs to disc
Fit.optimized.ssobj.batch <- Fit.optimized.ssobj.batch[-1,]

write.csv(file = "./Outputs/SABER_param_IOCCG.csv", x=Fit.optimized.ssobj.batch, sep = ",", quote = F,
          row.names = F, col.names = T)

#---------------------------------------------
### 1.2 Inversion of IOCCG dataset using QAA
#---------------------------------------------
#Create data-frame to store QAA outputs
Fit.optimized.qaa.ioccg <- data.frame("chl"=0, 
                                      "adg.440"=0,
                                      "bbp550"=0
)

negative_status_vec_qaa = vector()


#Run QAA
for (j in 1:dim(rrs.HL)[1]) {
  
  rrs_inverse_input = (as.numeric(rrs.HL[j,]))
  
  ix405 = which.min(abs(405 - rrs.HL.wl))
  ix411 = which.min(abs(411 - rrs.HL.wl))
  ix670 = which.min(abs(670 - rrs.HL.wl))
  ix683 = which.min(abs(683 - rrs.HL.wl))
  interp_ignore = c(ix405, ix411, ix670, ix683)
  
  
  
  if (any(rrs_inverse_input[-interp_ignore] < 0)) {
    print("!!!!!! WARNING: The Rrs observation has negative value between 420-665nm, discarding data")
    negative_status_vec_qaa[j] = "neg_rrs_green"
    j = j+1
    
    next
    
  } else {
    if (any(rrs_inverse_input[interp_ignore] < 0)) {
      print("!!!!!! WARNING: The Rrs observation has negative value outside 420-665nm, Set to Zero")
      
      negidx = rrs_inverse_input[interp_ignore] < 0
      negloc = which(negidx == TRUE)
      negative_status_vec_qaa[negloc] = 0
      
      if (any(rrs_inverse_input[interp_ignore[1:2]] < 0)) {
        negative_status_vec_qaa[j] = "neg_rrs_blue"
      }
      if (any(rrs_inverse_input[interp_ignore[3:4]] < 0)) {
        negative_status_vec_qaa[j] = "neg_rrs_red"
      }
      if (all(rrs_inverse_input[interp_ignore] < 0)) {
        negative_status_vec_qaa[j] = "neg_rrs_blue-red"
      }
      
    }
    
    if (all(rrs_inverse_input > 0)) {
      print(paste0(j, " no. Rrs is good observation"))
      negative_status_vec_qaa[j] = "no_neg"
    }
    
  }
  
  qaa_output = QAA.v5(waves = wavelength, Rrs = rrs_inverse_input)
  ix443 = which.min(abs(443 - wavelength))
  a_ph_443 = qaa_output$a_phi[ix443]
  
  chl_est = (a_ph_443/0.06)^(1/0.65)
  
  par0 = c(chl = chl_est, adg440 =qaa_output$a_dg_443, 
           bbp550 = qaa_output$b_bp_555)
  
  if (any(par0 == Inf) | any(is.na(par0))) {
    print("QAA failed, default mean values wil be used")
    par0 = c(chl = mean(HL.deep.iop$chl), adg440 = mean(HL.deep.iop$acdom440 + HL.deep.iop$anap440), 
             bbp550 = mean(HL.deep.iop$bbp550, na.rm = T))
    
    
  } 
  
  Fit.optimized.qaa.ioccg <- rbind(Fit.optimized.qaa.ioccg, par0)
  cat(paste0("\033[0;43m","############### QAA FINISHED for ", j, " no. spectra ################","\033[0m","\n"))
  
}

#Export the QAA outputs to disc
Fit.optimized.qaa.ioccg <- Fit.optimized.qaa.ioccg[-1,]
write.csv(file = "./Outputs/QAA_param_IOCCG.csv", x=Fit.optimized.qaa.ioccg, sep = ",", quote = F,
          row.names = F, col.names = T)


#---------------------------------------------
### 1.3 #Validate inversion parameters
#---------------------------------------------
fit_param_saber = read.csv("./Outputs/SABER_param_IOCCG.csv", header = T)
fit_param_qaa = read.csv("./Outputs/QAA_param_IOCCG.csv", header = T)
iop_param_qc = HL.deep.iop
iop_param_qc$adg443 = iop_param_qc$acdom440 + iop_param_qc$anap440

#====================
# 1.3.1 Plot [chl]
#====================
chl_saber = data.frame("id" = seq(1, length(fit_param_saber$chl),1),
                       "chl" = fit_param_saber$chl)
chl_saber_long = reshape2::melt(data = chl_saber, id.vars = "id")

chl_qaa = data.frame("id" = seq(1, length(fit_param_qaa$chl),1),
                     "chl" = fit_param_qaa$chl)
chl_qaa_long = reshape2::melt(data = chl_qaa, id.vars = "id")

chl_IOCCG = data.frame("id" = seq(1, length(iop_param_qc$chl),1),
                       "chl" = iop_param_qc$chl)
chl_IOCCG_long = reshape2::melt(data = chl_IOCCG, id.vars = "id")

zprime <- merge(chl_saber_long,chl_qaa_long,by=c("id", "variable"))
zprime <- merge(zprime, chl_IOCCG_long, by=c("id", "variable"))

names(zprime) <- c("id", "param", "SABER", "QAA", "IOCCG")

rm(chl_saber, chl_saber_long, chl_qaa, chl_qaa_long, chl_IOCCG, chl_IOCCG_long)

# Plot the validation
legend_title <- element_blank()
legend_position <- c(0.05, 0.85)

xmin <- 10^-3.5; xmax <- 10^3.5; xstp <- xmax/4
xlbl <- expression(paste(italic("[chl]")["IOCCG"],italic("in vivo"), " [","mg"^1,"m"^-3,"]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("[chl]")["MODEL"],italic("estimated"), " [","mg"^1,"m"^-3,"]"))

asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(zprime,aes(x = IOCCG)) +
  geom_abline(slope = 1,intercept = 0,colour="black", na.rm = FALSE, show.legend = FALSE) +
  
  # For varying shapes for SABER and QAA
  geom_point(aes(y = QAA, shape="xx2",fill = "yy2"),
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  geom_point(aes(y = SABER, shape="xx1", fill = "yy1"), 
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  
  scale_shape_manual(name = "Model Type",
                     labels=(c(expression(paste("SABER")),expression(paste("QAAv5")))),
                     values = c(21,23))+
  
  # scale_color_manual(labels=rev(c(expression(paste("SABER")),expression(paste("QAAv5")))),
  #                    values = c("navyblue","goldenrod2"))+
  
  scale_fill_manual(name = "Model Type",
                    labels=c(expression(paste("SABER")),expression(paste("QAAv5"))),
                    values = rev(c("navyblue","goldenrod2")))+
  
  
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
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 0),
        axis.text.y = element_text(size = 15, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
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
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g

ggsave(paste0("./Outputs/IOCCG_chl_valid_v1.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)


#====================
# 1.3.2 Plot [adg443]
#====================
adg.440_saber = data.frame("id" = seq(1, length(fit_param_saber$adg.440),1),
                           "adg.440" = fit_param_saber$adg.440)
adg.440_saber_long = reshape2::melt(data = adg.440_saber, id.vars = "id")

adg.440_qaa = data.frame("id" = seq(1, length(fit_param_qaa$adg.440),1),
                         "adg.440" = fit_param_qaa$adg.440)
adg.440_qaa_long = reshape2::melt(data = adg.440_qaa, id.vars = "id")

adg.440_nomad = data.frame("id" = seq(1, length(iop_param_qc$adg443),1),
                           "adg.440" = iop_param_qc$adg443)
adg.440_nomad_long = reshape2::melt(data = adg.440_nomad, id.vars = "id")

zprime <- merge(adg.440_saber_long,adg.440_qaa_long,by=c("id", "variable"))
zprime <- merge(zprime, adg.440_nomad_long, by=c("id", "variable"))

names(zprime) <- c("id", "param", "SABER", "QAA", "IOCCG")

rm(adg.440_saber, adg.440_saber_long, adg.440_qaa, adg.440_qaa_long, adg.440_nomad, adg.440_nomad_long)

# Plot the validation
legend_title <- element_blank()
legend_position <- c(0.05, 0.85)

xmin <- 10^-4; xmax <- 10^2; xstp <- xmax/4
xlbl <- expression(paste(italic("a")["dg"](443),italic("IOCCG"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("a")["dg"](443),italic("MODEL"), " [", "m"^-1, "]"))

asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(zprime,aes(x = IOCCG)) +
  geom_abline(slope = 1,intercept = 0,colour="black", na.rm = FALSE, show.legend = FALSE) +
  
  # For varying shapes for SABER and QAA
  geom_point(aes(y = QAA, shape="xx2",fill = "yy2"),
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  geom_point(aes(y = SABER, shape="xx1", fill = "yy1"), 
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  
  scale_shape_manual(name = "Model Type",
                     labels=(c(expression(paste("SABER")),expression(paste("QAAv5")))),
                     values = c(21,23))+
  
  # scale_color_manual(labels=rev(c(expression(paste("SABER")),expression(paste("QAAv5")))),
  #                    values = c("navyblue","goldenrod2"))+
  
  scale_fill_manual(name = "Model Type",
                    labels=c(expression(paste("SABER")),expression(paste("QAAv5"))),
                    values = rev(c("navyblue","goldenrod2")))+
  
  
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
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 0),
        axis.text.y = element_text(size = 15, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
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
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g

ggsave(paste0("./Outputs/IOCCG_adg443_valid_v1.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#====================
# 1.3.3 Plot [bbp555]
#====================
bbp.555_saber = data.frame("id" = seq(1, length(fit_param_saber$bbp550),1),
                           "bbp.555" = fit_param_saber$bbp550)
bbp.555_saber_long = reshape2::melt(data = bbp.555_saber, id.vars = "id")

bbp.555_qaa = data.frame("id" = seq(1, length(fit_param_qaa$bbp550),1),
                         "bbp.555" = fit_param_qaa$bbp550)
bbp.555_qaa_long = reshape2::melt(data = bbp.555_qaa, id.vars = "id")

bbp.555_ioccg = data.frame("id" = seq(1, length(HL.deep.iop$bbp550),1),
                           "bbp.555" = HL.deep.iop$bbp550)
bbp.555_ioccg_long = reshape2::melt(data = bbp.555_ioccg, id.vars = "id")

zprime <- merge(bbp.555_saber_long,bbp.555_qaa_long,by=c("id", "variable"))
zprime <- merge(zprime, bbp.555_ioccg_long, by=c("id", "variable"))

names(zprime) <- c("id", "param", "SABER", "QAA", "IOCCG")

rm(bbp.555_saber, bbp.555_saber_long, bbp.555_qaa, bbp.555_qaa_long, bbp.555_ioccg, bbp.555_ioccg_long )

# Plot the validation
legend_title <- element_blank()
legend_position <- c(0.05, 0.85)

xmin <- 10^-4; xmax <- 10^0; xstp <- xmax/4
xlbl <- expression(paste(italic("b")["bp"](555),italic("IOCCG"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("b")["bp"](555),italic("MODEL"), " [", "m"^-1, "]"))

asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(zprime,aes(x = IOCCG)) +
  geom_abline(slope = 1,intercept = 0,colour="black", na.rm = FALSE, show.legend = FALSE) +
  
  # For varying shapes for SABER and QAA
  geom_point(aes(y = QAA, shape="xx2",fill = "yy2"),
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  geom_point(aes(y = SABER, shape="xx1", fill = "yy1"), 
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  
  scale_shape_manual(name = "Model Type",
                     labels=(c(expression(paste("SABER")),expression(paste("QAAv5")))),
                     values = c(21,23))+
  
  # scale_color_manual(labels=rev(c(expression(paste("SABER")),expression(paste("QAAv5")))),
  #                    values = c("navyblue","goldenrod2"))+
  
  scale_fill_manual(name = "Model Type",
                    labels=c(expression(paste("SABER")),expression(paste("QAAv5"))),
                    values = rev(c("navyblue","goldenrod2")))+
  
  
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
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 0),
        axis.text.y = element_text(size = 15, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
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
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g

ggsave(paste0("./Outputs/IOCCG_bbp555_valid_v1.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#----------------------------------------------------
# 2.4 Statistics on the goodness of Fit for inversion 
# retrieved parameters
#----------------------------------------------------
names(fit_param_saber) = c("chl", "adg443", "bbp555", "chl_sd", "adg443_sd", "bbp555_sd")
names(fit_param_qaa) = c("chl", "adg443", "bbp555")
names(iop_param_qc) = c("chl", "adg443", "bbp555")

errorstat_IOCCG_SABER = goodness_of_fit(actual = iop_param_qc, predicted = fit_param_saber, bbp.exist = T)
errorstat_IOCCG_SABER <- errorstat_IOCCG_SABER %>% dplyr::select(parameter, everything())

errorstat_IOCCG_QAA = goodness_of_fit(actual = iop_param_qc, predicted = fit_param_qaa, bbp.exist = T)
errorstat_IOCCG_QAA <- errorstat_IOCCG_QAA %>% dplyr::select(parameter, everything())

# ====================================================
# 2. For NOMAD dataset
# ====================================================
nomad = read.csv("./data/nomad_dataset_simplified.csv", header = T)

#Load NOMAD Rrs
es_idx = grep("es", colnames(nomad))
colnames(nomad)[es_idx]

lw_idx = grep("lw", colnames(nomad))
colnames(nomad)[lw_idx]

rrs_nomad = data.frame("nomad_id"=nomad$id, "id" = seq(1, length(nomad$year),1), 
                       "Rrs"= nomad[lw_idx]/nomad[es_idx])

Rrs_idx = grep("Rrs", colnames(rrs_nomad))

vars  <- colnames(rrs_nomad)[Rrs_idx]
nomad_wave = as.numeric(gsub("\\D", "", vars))

rrs_nomad_qc = rrs_nomad[vars]  %>% naniar::replace_with_na_all(condition = ~.x == 1)

ix405 = which.min(abs(405 - nomad_wave))
ix411 = which.min(abs(411 - nomad_wave))
ix670 = which.min(abs(670 - nomad_wave))
ix683 = which.min(abs(683 - nomad_wave))
interp_ignore = c(ix405, ix411, ix670, ix683)


rrs_nomad_qc_interp = matrix(data = 0, nrow = nrow(rrs_nomad_qc), ncol = ncol(rrs_nomad_qc))
for (i in 1:nrow(rrs_nomad_qc)) {
  
  rrs_nomad_qc_interp[i,] = Hmisc::approxExtrap(x=  nomad_wave ,as.numeric(rrs_nomad_qc[i,]), 
                                   xout = nomad_wave,
                                   method = "linear", na.rm = T)$y
}

colnames(rrs_nomad_qc_interp) = paste0("Rrs_",nomad_wave)

rrs_nomad_qc_interp_df = data.frame(rrs_nomad_qc_interp)

rrs_nomad_qc_interp_df$nomad_id = rrs_nomad$nomad_id
rrs_nomad_qc_interp_df$id = rrs_nomad$id


rm(rrs_nomad, rrs_nomad_qc, rrs_nomad_qc_interp)


#Load NOMAD IOPs
ag443_idx = grep("ag443", colnames(nomad))
colnames(nomad)[ag443_idx]

ad443_idx = grep("ad443", colnames(nomad))
colnames(nomad)[ad443_idx]

chl_idx = grep("chl", colnames(nomad))
colnames(nomad)[chl_idx][1]

bb555_idx = grep("bb555", colnames(nomad))
colnames(nomad)[bb555_idx]

iop_nomad = data.frame("nomad_id"=nomad$id, "id" = seq(1, length(nomad$year),1), 
                       "chl" = nomad[chl_idx][1], 
                       "ag443"= nomad[ag443_idx], "ad443"= nomad[ad443_idx],
                       "bbp555"=nomad[bb555_idx])
iop_nomad_qc = iop_nomad  %>% naniar::replace_with_na_all(condition = ~.x == -999)

iop_nomad_qc = iop_nomad_qc[complete.cases(iop_nomad_qc[-length(iop_nomad_qc)]),]
rrs_nomad_qc_interp_df_iop_idx = rrs_nomad_qc_interp_df[rrs_nomad_qc_interp_df$id %in% iop_nomad_qc$id,]

rrs_nomad_final = rrs_nomad_qc_interp_df_iop_idx[-(21:22)]
iop_nomad_qc$adg443  = iop_nomad_qc$ag443 + iop_nomad_qc$ad443


#---------------------------------------------
### 2.1 Inversion of NOMAD dataset using SABER
#---------------------------------------------
#Create data-frame to store SABER outputs
Fit.optimized.ssobj.nomad <- data.frame("chl"=0, 
                                        "adg.440"=0,
                                        "bbp550"=0,
                                        "chl.sd"=0, "adg440.sd"=0, "bbp550.sd"=0)
negative_status_vec = vector()

#Run Conditions
qaa_prefit = TRUE
manual_par = FALSE
wavelength = nomad_wave

for (j in 1:dim(rrs_nomad_final)[1]) {
  
  ##Do the optimization 

  obj = c("log-LL",  "SSR" ,    "obj_L98")
  methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")
  
  rrs_inverse_input = surface_rrs_translate(as.numeric(rrs_nomad_final[j,]))
  
  ix405 = which.min(abs(405 - nomad_wave))
  ix411 = which.min(abs(411 - nomad_wave))
  ix670 = which.min(abs(670 - nomad_wave))
  ix683 = which.min(abs(683 - nomad_wave))
  interp_ignore = c(ix405, ix411, ix670, ix683)
  
  
  
  if (any(rrs_inverse_input[-interp_ignore] < 0)) {
    print("!!!!!! WARNING: The Rrs observation has negative value between 420-665nm, discarding data")
    negative_status_vec[j] = "neg_rrs_green"
    j = j+1
     
    next
    
  } else {
    if (any(rrs_inverse_input[interp_ignore] < 0)) {
      print("!!!!!! WARNING: The Rrs observation has negative value outside 420-665nm, Set to Zero")
      
      negidx = rrs_inverse_input[interp_ignore] < 0
      negloc = which(negidx == TRUE)
      rrs_inverse_input[negloc] = 0
      
      if (any(rrs_inverse_input[interp_ignore[1:2]] < 0)) {
        negative_status_vec[j] = "neg_rrs_blue"
      }
      if (any(rrs_inverse_input[interp_ignore[3:4]] < 0)) {
        negative_status_vec[j] = "neg_rrs_red"
      }
      if (all(rrs_inverse_input[interp_ignore] < 0)) {
        negative_status_vec[j] = "neg_rrs_blue-red"
      }
      
    }
    
    if (all(rrs_inverse_input > 0)) {
      negative_status_vec[j] = "no_neg"
    }
    
  }
  
  if (qaa_prefit == TRUE) {
    qaa_output = QAA.v5(waves = wavelength, Rrs = rrs_inverse_input)
    par0 = c(chl = qaa_output$chl, adg440 =qaa_output$a_dg_443, 
             bbp550 = qaa_output$b_bp_555, "pop_sd" = 0.001)
    
    if (any(par0 == Inf) | any(is.na(par0))) {
      print("Prefit failed, default mean values wil be used")
      par0 = c(chl = mean(iop_nomad_qc$chl), adg440 = mean(iop_nomad_qc$ag443 + iop_nomad_qc$ad443), 
               bbp550 = mean(iop_nomad_qc$bb555, na.rm = T), "pop_sd" = 0.001)
      
      upper.bound <- c(60,10,0.1,0.01)  
      lower.bound <- c(0.05, 0.01, 0.001, 0.0001)
    } else {
      print(paste0("prefit values are computed as chl: ", signif(par0[1], digits = 3), " adg443: ", signif(par0[2], digits = 3), " bbp555: ",signif(par0[3], digits = 3), " pop_sd: ",signif(par0[4], digits = 3)))
      
      upper.bound <- c(par0[-length(par0)] + 0.15*par0[-length(par0)],
                       1)  
      lower.bound <- c(par0[-length(par0)] - 0.15*par0[-length(par0)], 
                       0.0001)
    }
    
  } else {
    if (manual_par == T) {
      par0 = c(chl = mean(iop_nomad_qc$chl), adg440 = mean(iop_nomad_qc$ag443 + iop_nomad_qc$ad443), 
               bbp550 = mean(iop_nomad_qc$bb555, na.rm = T), "pop_sd" = 0.13129354/100)
      
      upper.bound <- c(50,10,0.5,0.01)  
      lower.bound <- c(0.05, 0.01, 0.001, 0.0001)
  }
  
    
  }
  
  #Optimize
  inverse_output <- suppressWarnings(solve.objective.inverse.deep.final(
                            initial = par0, 
                            bbp.constrain  = F,
                            
                            #bbp.550 = ,
                            
                            #obsdata = as.numeric(rrs_nomad_final[j,]),
                            obsdata = rrs_inverse_input,
                            
                            sa.model = "am03", 
                            
                            #obj.fn =obj.run ,
                            obj.fn =obj[1],
                            
                            auto_spectral_slope = T,
                            manual_spectral_slope = F, 
                            
                            manual_spectral_slope_vals = c("s_g"=a_g_vec[i], "s_d"=a_d_vec[i], "gamma"=1),
                            
                            method.opt = methods.opt[4],
                            lower.b = lower.bound,
                            upper.b = upper.bound, 
                            batch = FALSE, pop.sd = FALSE))
  
  Fit.optimized.ssobj.nomad <- rbind(Fit.optimized.ssobj.nomad, c(inverse_output[[1]]$estimates[1:3],
                                                                  inverse_output[[1]]$`sd(+/-)`[1:3]))
  cat(paste0("\033[0;43m","############### INVERSION FINISHED for ", j, " no. spectra ################","\033[0m","\n"))
}

######## Save SABER output to the disc
Fit.optimized.ssobj.nomad <- Fit.optimized.ssobj.nomad[-1,]

nomad_inverse_idx = which(negative_status_vec == "neg_rrs_green") 

length(iop_nomad_qc$nomad_id[-nomad_inverse_idx])

write.csv(file = "./Outputs/inv_param_NOMAD_qaa_initial_new.csv", x=Fit.optimized.ssobj.nomad, sep = ",", quote = F,
          row.names = F, col.names = T)

#---------------------------------------------
### 2.2 Inversion of NOMAD dataset using QAA
#---------------------------------------------
#Create data-frame to store QAA outputs
Fit.optimized.qaa.nomad <- data.frame("chl"=0, 
                                        "adg.440"=0,
                                        "bbp550"=0
                                        )
negative_status_vec_qaa = vector()

#Run QAA
for (j in 1:dim(rrs_nomad_final)[1]) {
  
  rrs_inverse_input = surface_rrs_translate(as.numeric(rrs_nomad_final[j,]))
  
  ix405 = which.min(abs(405 - nomad_wave))
  ix411 = which.min(abs(411 - nomad_wave))
  ix670 = which.min(abs(670 - nomad_wave))
  ix683 = which.min(abs(683 - nomad_wave))
  interp_ignore = c(ix405, ix411, ix670, ix683)
  
  
  
  if (any(rrs_inverse_input[-interp_ignore] < 0)) {
    print("!!!!!! WARNING: The Rrs observation has negative value between 420-665nm, discarding data")
    negative_status_vec_qaa[j] = "neg_rrs_green"
    j = j+1
    
    next
    
  } else {
    if (any(rrs_inverse_input[interp_ignore] < 0)) {
      print("!!!!!! WARNING: The Rrs observation has negative value outside 420-665nm, Set to Zero")
      
      negidx = rrs_inverse_input[interp_ignore] < 0
      negloc = which(negidx == TRUE)
      negative_status_vec_qaa[negloc] = 0
      
      if (any(rrs_inverse_input[interp_ignore[1:2]] < 0)) {
        negative_status_vec_qaa[j] = "neg_rrs_blue"
      }
      if (any(rrs_inverse_input[interp_ignore[3:4]] < 0)) {
        negative_status_vec_qaa[j] = "neg_rrs_red"
      }
      if (all(rrs_inverse_input[interp_ignore] < 0)) {
        negative_status_vec_qaa[j] = "neg_rrs_blue-red"
      }
      
    }
    
    if (all(rrs_inverse_input > 0)) {
      negative_status_vec_qaa[j] = "no_neg"
    }
    
  }
  
    qaa_output = QAA.v5(waves = wavelength, Rrs = rrs_inverse_input)
    ix443 = which.min(abs(443 - wavelength))
    a_ph_443 = qaa_output$a_phi[ix443]
    
    chl_est = (a_ph_443/0.06)^(1/0.65)
    
    par0 = c(chl = chl_est, adg440 =qaa_output$a_dg_443, 
             bbp550 = qaa_output$b_bp_555)
    
    if (any(par0 == Inf) | any(is.na(par0))) {
      print("QAA failed, default mean values wil be used")
      par0 = c(chl = mean(iop_nomad_qc$chl), adg440 = mean(iop_nomad_qc$ag443 + iop_nomad_qc$ad443), 
               bbp550 = mean(iop_nomad_qc$bb555, na.rm = T))
      
      
    } 
    
    Fit.optimized.qaa.nomad <- rbind(Fit.optimized.qaa.nomad, par0)
    cat(paste0("\033[0;43m","############### QAA FINISHED for ", j, " no. spectra ################","\033[0m","\n"))
    
}

######## Save QAA output to the disc
Fit.optimized.qaa.nomad <- Fit.optimized.qaa.nomad[-1,]

nomad_qaa_idx = which(negative_status_vec_qaa == "neg_rrs_green") 

length(iop_nomad_qc$nomad_id[-nomad_qaa_idx])

write.csv(file = "./Outputs/QAA_param_NOMAD.csv", x=Fit.optimized.qaa.nomad, sep = ",", quote = F,
          row.names = F, col.names = T)
  

#---------------------------------------------
### 2.3 #Validate inversion parameters
#---------------------------------------------

fit_param_saber = read.csv("./Outputs/inv_param_NOMAD_qaa_initial_new.csv", header = T)
fit_param_qaa = read.csv("./Outputs/QAA_param_NOMAD.csv", header = T)
iop_param_qc = iop_nomad_qc[-c(nomad_qaa_idx),]

#====================
# 2.3.1 Plot [chl]
#====================
chl_saber = data.frame("id" = seq(1, length(fit_param_saber$chl),1),
                       "chl" = fit_param_saber$chl)
chl_saber_long = reshape2::melt(data = chl_saber, id.vars = "id")

chl_qaa = data.frame("id" = seq(1, length(fit_param_qaa$chl),1),
                       "chl" = fit_param_qaa$chl)
chl_qaa_long = reshape2::melt(data = chl_qaa, id.vars = "id")

chl_nomad = data.frame("id" = seq(1, length(iop_param_qc$chl),1),
                       "chl" = iop_param_qc$chl)
chl_nomad_long = reshape2::melt(data = chl_nomad, id.vars = "id")

zprime <- merge(chl_saber_long,chl_qaa_long,by=c("id", "variable"))
zprime <- merge(zprime, chl_nomad_long, by=c("id", "variable"))

names(zprime) <- c("id", "param", "SABER", "QAA", "NOMAD")


# Plot the validation
legend_title <- element_blank()
legend_position <- c(0.05, 0.85)

xmin <- 10^-3.5; xmax <- 10^3.5; xstp <- xmax/4
xlbl <- expression(paste(italic("[chl]")["NOMAD"],italic("in vivo"), " [","mg"^1,"m"^-3,"]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("[chl]")["MODEL"],italic("estimated"), " [","mg"^1,"m"^-3,"]"))

asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(zprime,aes(x = NOMAD)) +
  geom_abline(slope = 1,intercept = 0,colour="black", na.rm = FALSE, show.legend = FALSE) +
  
  # For varying shapes for SABER and QAA
  geom_point(aes(y = QAA, shape="xx2",fill = "yy2"),
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  geom_point(aes(y = SABER, shape="xx1", fill = "yy1"), 
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  
  scale_shape_manual(name = "Model Type",
                     labels=(c(expression(paste("SABER")),expression(paste("QAAv5")))),
                     values = c(21,23))+
  
  # scale_color_manual(labels=rev(c(expression(paste("SABER")),expression(paste("QAAv5")))),
  #                    values = c("navyblue","goldenrod2"))+
  
  scale_fill_manual(name = "Model Type",
                    labels=c(expression(paste("SABER")),expression(paste("QAAv5"))),
                    values = rev(c("navyblue","goldenrod2")))+
  
  
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
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 0),
        axis.text.y = element_text(size = 15, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
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
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g

ggsave(paste0("./Outputs/nomad_chl_valid_v2.png"), plot = g,
         scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)


#====================
# 2.3.2 Plot [adg443]
#====================
adg.440_saber = data.frame("id" = seq(1, length(fit_param_saber$adg.440),1),
                       "adg.440" = fit_param_saber$adg.440)
adg.440_saber_long = reshape2::melt(data = adg.440_saber, id.vars = "id")

adg.440_qaa = data.frame("id" = seq(1, length(fit_param_qaa$adg.440),1),
                     "adg.440" = fit_param_qaa$adg.440)
adg.440_qaa_long = reshape2::melt(data = adg.440_qaa, id.vars = "id")

adg.440_nomad = data.frame("id" = seq(1, length(iop_param_qc$adg443),1),
                       "adg.440" = iop_param_qc$adg443)
adg.440_nomad_long = reshape2::melt(data = adg.440_nomad, id.vars = "id")

zprime <- merge(adg.440_saber_long,adg.440_qaa_long,by=c("id", "variable"))
zprime <- merge(zprime, adg.440_nomad_long, by=c("id", "variable"))

names(zprime) <- c("id", "param", "SABER", "QAA", "NOMAD")

# Plot the validation
legend_title <- element_blank()
legend_position <- c(0.05, 0.85)

xmin <- 10^-4; xmax <- 10^2; xstp <- xmax/4
xlbl <- expression(paste(italic("a")["dg"](443),italic("NOMAD"), " [", "m"^-1, "]"))
ymin <- xmin; ymax <- xmax ; ystp <- ymax/5
ylbl <-expression(paste(italic("a")["dg"](443),italic("MODEL"), " [", "m"^-1, "]"))

asp_rat <- (xmax-xmin)/(ymax-ymin)

g <- ggplot(zprime,aes(x = NOMAD)) +
  geom_abline(slope = 1,intercept = 0,colour="black", na.rm = FALSE, show.legend = FALSE) +
  
  # For varying shapes for SABER and QAA
  geom_point(aes(y = QAA, shape="xx2",fill = "yy2"),
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  geom_point(aes(y = SABER, shape="xx1", fill = "yy1"), 
             size=3.0, na.rm = FALSE, show.legend = TRUE) +
  
  scale_shape_manual(name = "Model Type",
                     labels=(c(expression(paste("SABER")),expression(paste("QAAv5")))),
                     values = c(21,23))+
  
  # scale_color_manual(labels=rev(c(expression(paste("SABER")),expression(paste("QAAv5")))),
  #                    values = c("navyblue","goldenrod2"))+
  
  scale_fill_manual(name = "Model Type",
                    labels=c(expression(paste("SABER")),expression(paste("QAAv5"))),
                     values = rev(c("navyblue","goldenrod2")))+
  
  
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
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 0),
        axis.text.y = element_text(size = 15, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
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
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g

ggsave(paste0("./Outputs/nomad_adg443_valid_v2.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#----------------------------------------------------
# 2.4 Statistics on the goodness of Fit for inversion 
# retrieved parameters
#----------------------------------------------------




# ########### 2.5 TEST Plot for bbp555
# plot(log10(Fit.optimized.ssobj.nomad$bbp550), 
#      log10(Fit.optimized.qaa.nomad$bbp550), xlim=c(-4,4), ylim=c(-4,4), 
#      pch=19, col="#007500", xlab = "[bbp555]_SABER", ylab = "[bbp555]_QAA", main="[bbp555]")
# abline(0,1, lwd=2)

data(groundbeef)
serving <- groundbeef$serving
fitg <- fitdist(serving, "gamma")
llplot(fitg)

llplot(fitg, expand = 5, pal.col = heat.colors(100), fit.show = TRUE)
llplot(fitg, pal.col = heat.colors(100), fit.show = TRUE)
llplot(fitg, back.col = FALSE, nlev = 25, fit.show = TRUE)
