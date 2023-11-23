#=====================================================================================================
# Saber.sensitivity.analysis.R tests the uncertainty in the forward SABER output (Rrs), that is 
# propagated from the forward model inputs, i.e. [chl], a_dg(443), b_bp(555) and H. The multi-parametric
# uncertainty is obtained from "SOBOL" indice estimation following Monte Carlo Simulations.

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#=====================================================================================================
library(sensitivity)
library(boot)
library(tidymodels)

# Define the forward model whose sensitivity is under scope
model_function <- function(x) {
  new_data <- data.frame(x1 = x[, 1], x2 = x[, 2], x3 = x[, 3], x4 = x[, 4])
  temp = Saber_forward_fast_sensitivity_test(use_true_IOPs = F, use_manual_slope = F, 
                     
                     chl = new_data$x1,
                     a_dg = new_data$x2,
                     bbp.550 = new_data$x3, 
                     
                     slope.parametric = F, 
                     Rrs_input_for_slope = obsdata,
                     
                     z = new_data$x4,
                     rb.fraction = fA.set,
                     
                     
                     wavelength = wavelength,
                     verbose = F)
  return(rowSums(temp))
}

# Create samples of forward model input parameters from truncated Gaussian/Weibull family
set.seed(345)
n = 10000
chl_1 = rnorm.trunc(n = n, mean = 4, sd = 10, min = 0.5, max = 30)
chl_2 = rnorm.trunc(n = n, mean = 4, sd = 10, min = 0.5, max = 30)

adg443_1 = rnorm.trunc(n = n, mean = 1.5, sd = 1, min = 0.25, max = 5)
adg443_2 = rnorm.trunc(n = n, mean = 1.5, sd = 1, min = 0.25, max = 5)

bbp555_1 = rnorm.trunc(n = n, mean = 0.007, sd = 0.002, min = 0.002, max = 0.02)
bbp555_2 = rnorm.trunc(n = n, mean = 0.007, sd = 0.002, min = 0.002, max = 0.02)

H_1 = rnorm.trunc(n = n, mean = 3, sd = 1, min = 1, max = 10)
H_2 = rnorm.trunc(n = n, mean = 3, sd = 1, min = 1, max = 10)


X1 = data.frame(chl_1, adg443_1, bbp555_1, H_1)
#X1 = expand.grid(X1)
names(X1) = c("chl", "adg443", "bbp555", "H")

X2 = data.frame(chl_2, adg443_2, bbp555_2, H_2)
#X2 = expand.grid(X2)
names(X2) = c("chl", "adg443", "bbp555", "H")


#Perform the sensitivity analysis
sobol_results <- sobolSalt(model = model_function, 
                       X1 = X1, 
                       X2 = X2, 
                       nboot = 1000, scheme = "B",conf = 0.95)

#Store the sensitivity analysis outputs in data-frame(s)
sobol_result_df_1st = sobol_results$S*100
sobol_result_df_1st$type = "1st order"
sobol_result_df_1st$id = seq(1,4,1)

sobol_result_df_2nd = sobol_results$S2*100
sobol_result_df_2nd$type = "2nd order"
sobol_result_df_2nd$id = seq(1,6,1)

sobol_result_df_tot = sobol_results$T*100
sobol_result_df_tot$type = "Total"
sobol_result_df_tot$id = seq(1,4,1)


sobol_comb = rbind(sobol_result_df_1st, sobol_result_df_tot)

#Plot the results from sensitivity analysis
chl_lab <- expression(paste("[", italic("chl"), "]"))
adg443_lab <- expression(paste(italic("a")["dg"](443)))
bbp555_lab <- expression(paste(italic("b")["bp"](555)))
H_lab <- expression(paste(italic("H")))


g <- ggplot(sobol_comb, aes(x=id, y=original, fill=type)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()
           ) +
  geom_errorbar(aes(ymin=`min. c.i.`, ymax=`max. c.i.`), width=0.25, size= 1.3, 
                color = "black",
                position=position_dodge(.9))+
  
  scale_x_continuous(name = " ", breaks = seq(1, 4, 1),
              labels = c(chl_lab, adg443_lab, bbp555_lab, H_lab)) +
  ylim(0,100)+
  ylab(expression(paste("Var[",italic("R"),{}[rs],"(",lambda,")] (%)"))) +
  
  scale_shape_manual(name = "legend",values = c(21,23), labels= c("1st order", "Total"))+
  scale_fill_manual(name="legend", values = c("#481567FF", "#20A387FF"), 
                     labels= c("1st order", "Total"))+
  
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
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
g
ggsave(paste0("./outputs/sensitivity_analysis.png"), plot = g,
                scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
