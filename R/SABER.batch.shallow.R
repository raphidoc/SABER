#Pre-FIT of initial values
pre.Fit <- data.frame("C_ph"=seq(1,10,1),
                      "a_cdom.440"=seq(0.5,5,0.5),
                      "a.nap.440"=seq(0.01,0.1,0.01),
                      "z"=seq(0.5,10,1)
                      #,
                      # rb1 = seq(0,1,2*0.05263158),
                      # rb2 = seq(0,1,2*0.05263158),
                      # rb3 = seq(0,1,2*0.05263158),
                      # rb4 = seq(0,1,2*0.05263158),
                      # rb5 = seq(0,1,2*0.05263158)
                      )

pre.Fit.input.LUT <- expand.grid(pre.Fit) #Create pre-Fit parameter space LUT


synthetic_data_shallow = data.frame()
for (i in 1:length(splitdata)) {
  temp = as.data.frame(splitdata[[i]])
  idx = caret::createDataPartition(temp$C_ph, p = 0.05, list = F)
  temp = temp[idx,]
  synthetic_data_shallow = rbind(synthetic_data_shallow, temp)
  
}

synthetic_data_shallow = synthetic_data_shallow[caret::createDataPartition(synthetic_data_shallow$C_ph, 
                                                                           p = 0.25, list = F),]

rrs.forward.SABER <- matrix(nrow = length(synthetic_data_shallow$C_ph),
                                 ncol = length(wavelength),0)

#Create the Progress Bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(synthetic_data_shallow$C_ph), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")


for (i in 1:length(synthetic_data_shallow$C_ph)) { #Create forward SABER LUT
  temp1 <- as.numeric(synthetic_data_shallow[i,])
  temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3], bbp.550 = Fit.input$bbp.550, 
                         verbose = F,z = temp1[4], realdata.exist = FALSE,
                         rb.fraction = fA.set)
  rrs.forward.SABER[i,] <- temp2[[1]]$Rrs
  setTxtProgressBar(pb, i)
  if (i == length(synthetic_data_shallow$C_ph)) {
    cat(paste0("\033[0;32m","###############LUT generation for AM03 FINISHED################","\033[0m","\n"))
  }
  #cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.SABER) - i), " remaining","\033[0m","\n"))
  
}
type_Rrs_below_jacobian = "deep"
Fit.optimized.ssobj.batch <- data.frame(
                                        "zb"=0,
                                        "rb1" = 0,
                                        "rb2" = 0,
                                        "rb3" = 0,
                                        "rb4" = 0,
                                        "rb5" = 0,
                                        "chl.sd"=0, "acdom.440.sd"=0, "anap440.sd"=0,
                                        "zb.sd"=0,
                                        "rb1.sd" = 0,
                                        "rb2.sd" = 0,
                                        "rb3.sd" = 0,
                                        "rb4.sd" = 0,
                                        "rb5.sd" = 0)

#Initial values from pre-fit
if (pop.sd == "unknown" & type_Rrs_below == "shallow") {
  par0 = c(chl = 2, acdom440 = 0.8,
           anap440 = 0.05, z = 2.5, 
           #rb.0 = 0.1,
           rb.1 = 0.3,
           rb.2 = 0.5,
           rb.3 = 0.0,
           rb.4 = 0.2,
           rb.5 = 0.0,
           population.sd = 0.1)
  
  #Autoscale Intital values from pre-Ft
  increament.scale <- 1
  
  lower.bound <- c((par0[1:3] - 0.8*par0[1:3]),z = 0.1,
                   #rb.0 = 0.1,
                   rb.1 = 0,
                   rb.2 = 0,
                   rb.3 = 0,
                   rb.4 = 0,
                   rb.5 = 0,
                   population.sd = 0.0001)
  
  upper.bound <- c((par0[1:3] + 5*par0[1:3]),z = 10,
                   #rb.0 = 1,
                   rb.1 = 1,
                   rb.2 = 1,
                   rb.3 = 1,
                   rb.4 = 1,
                   rb.5 = 1,
                   population.sd = 1)
}
startloc = 1
for (i in startloc:dim(rrs.forward.SABER)[1]) {
  
  ##Do the optimization 
  
  ##BATCH RUN
  inverse_output <- suppressWarnings(solve.objective.inverse.shallow.constrained.batch(
                                            constrain.param = as.numeric(synthetic_data_shallow[i,]),
                                            constrained = T,
                                            initial = as.numeric(par0), 
                                            obsdata = rrs.forward.SABER[i,],
                                            sa.model = "am03", obj.fn =obj.run , 
                                            method.opt = methods.opt[4],
                                            lower.b = lower.bound,
                                            upper.b = upper.bound, 
                                            batch = FALSE, pop.sd = FALSE))
  
  Fit.optimized.ssobj.batch <- rbind(Fit.optimized.ssobj.batch, c(inverse_output[[1]]$estimates,
                                                                  inverse_output[[1]]$`sd(+/-)`))
  cat(paste0("\033[0;32m","###############Calculate analytical Jacobian for each parameter for each wavelength################","\033[0m","\n"))
  
  test = saber.forward.jacobian.analytical(chl = synthetic_data_shallow$C_ph[i],
                                           acdom.440 = synthetic_data_shallow$a_cdom.440[i], 
                                           anap.440 = synthetic_data_shallow$a.nap.440[i],
                                           wavelength = wlrange, verbose = F)
  
  write.csv(test, file = paste0("./Jacobian/jacobian_chl-",synthetic_data_shallow$C_ph[i], "_cdom-", synthetic_data_shallow$a_cdom.440[i], "_nap-",synthetic_data_shallow$a.nap.440[i],".csv"),
              row.names = F, quote = F, sep = ",")
  
  cat(paste0("\033[0;32m","###############analytical Jacobian for each parameter is saved to disc################","\033[0m","\n"))
  
  cat(paste0("\033[0;34m","############### INVERSION FINISHED for ", i, " number spectra ################","\033[0m","\n"))
}

Fit.optimized.ssobj.batch.shallow.final = Fit.optimized.ssobj.batch[-1,]
Fit.optimized.ssobj.batch.shallow.final$zb.actual = synthetic_data_shallow$z

#CHL
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

coplabel <- expression(paste("[",H,"]",italic("synthetic")))
hl.el.label <- expression(paste("[",H,"]",italic("model")))
hl.inel.label <- expression(paste(italic("R")["rs"]("0"^"-", lambda),"HL"[italic("inelastic")]))
albert.label <- expression(paste(italic("R")["rs"]("0"^"-", lambda),"A&M"[italic("2003")]))

xlbl <- coplabel; ylbl <-hl.el.label
ymin <- 0; ymax <- 12 ; ystp <- ymax/4
xmin <- 0; xmax <- 12 ; xstp <- ymax/4
asp_rat <- (xmax-xmin)/(ymax-ymin)


errorb <- qnorm(0.975)*sd(Fit.optimized.ssobj.batch.shallow.final$zb)/sqrt(length(Fit.optimized.ssobj.batch.shallow.final$zb))

idx = is.na(Fit.optimized.ssobj.batch.shallow.final$zb.sd) 
Fit.optimized.ssobj.batch.shallow.final$chl.sd[idx] = errorb

g <- ggplot(Fit.optimized.ssobj.batch.shallow.final,aes(x = zb.actual)) + 
  geom_point(aes(y = zb), shape=21, fill="navyblue", 
             size=5.0, na.rm = T, show.legend = TRUE) +
  geom_abline(slope = 1,linetype="solid", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  geom_errorbar( aes(ymin = zb - (2.5*zb.sd),
                     ymax = zb + (2.5*zb.sd)),
                 color="black",
                 
                 width = 0.1, 
                 size=1.3
                 
  )+
  geom_ribbon(aes(ymin = zb - errorb,
                  ymax = zb + errorb),
              alpha = 0.6,
              fill = "grey70",
              colour="NA"
  )+
  # geom_smooth(size=1.5,level = 0.95,show.legend = F,linetype = "solid", 
  #             color="blue", data = Fit.optimized.ssobj.batch.shallow.final, 
  #             se= F, method = "lm", aes(x=chl.actual, y=chl))+
  # scale_fill_gradientn(colors =cols,breaks=seq(440,710,135),
  #                      limits=c(440,710), labels=paste(seq(440,710,135),"nm")) +
  # scale_colour_gradientn(colors =cols,breaks=seq(400,800,200),
  #                        limits=c(440,710), labels=paste(seq(440,710,135),"nm")) +
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  
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
g
ggsave(paste0("zB.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

