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
splitdata = split.data.frame(pre.Fit.input.LUT, pre.Fit.input.LUT$z)

synthetic_data_shallow = data.frame()
for (i in 1:length(splitdata)) {
  temp = as.data.frame(splitdata[[i]])
  idx = caret::createDataPartition(temp$C_ph, p = 0.05, list = F)
  temp = temp[idx,]
  synthetic_data_shallow = rbind(synthetic_data_shallow, temp)
  
}

synthetic_data_shallow = synthetic_data_shallow[caret::createDataPartition(synthetic_data_shallow$C_ph, 
                                                                           p = 0.25, list = F),]

rownames(synthetic_data_shallow) = seq(1,length(synthetic_data_shallow$C_ph),1)

rrs.forward.SABER <- matrix(nrow = length(synthetic_data_shallow$C_ph),
                            ncol = length(wavelength),0)

#Create the Progress Bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(synthetic_data_shallow$C_ph), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")

type_Rrs_below = "deep"
for (i in 1:length(synthetic_data_shallow$C_ph)) { #Create forward SABER LUT
  temp1 <- as.numeric(synthetic_data_shallow[i,])
  temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3], bbp.550 = Fit.input$bbp.550, 
                         verbose = F,realdata.exist = FALSE,
                         )
  rrs.forward.SABER[i,] <- temp2[[1]]$Rrs
  setTxtProgressBar(pb, i)
  if (i == length(synthetic_data_shallow$C_ph)) {
    cat(paste0("\033[0;32m","###############LUT generation for AM03 FINISHED################","\033[0m","\n"))
  }
  #cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.SABER) - i), " remaining","\033[0m","\n"))
  
}
type_Rrs_below_jacobian = "deep"
Fit.optimized.ssobj.batch <- data.frame("chl"=0, 
                                        "acdom.440"=0,
                                        "anap.440"=0,
                                        "chl.sd"=0, "acdom.440.sd"=0, "anap440.sd"=0)

#Initial values from pre-fit
par0 = c(chl = mean(synthetic_data_shallow$C_ph), acdom440 = mean(synthetic_data_shallow$a_cdom.440), 
         anap440 = mean(synthetic_data_shallow$a.nap.440), population.sd = 0.001)


upper.bound <- par0 + 5*par0
lower.bound <- par0 - 0.9*par0

startloc = 1
for (i in startloc:dim(rrs.forward.SABER)[1]) {
  
  ##Do the optimization 
  
  ##BATCH RUN
  inverse_output <- suppressWarnings(solve.objective.inverse(initial = par0, 
                                                             obsdata = rrs.forward.SABER[i,],
                                                             sa.model = "am03", 
                                                             obj.fn =obj.run , 
                                                             method.opt = methods.opt[4],
                                                             lower.b = lower.bound,
                                                             upper.b = upper.bound, 
                                                             batch = FALSE, pop.sd = FALSE))
  
  Fit.optimized.ssobj.batch <- rbind(Fit.optimized.ssobj.batch, c(inverse_output[[1]]$estimates[1:3],
                                                                  inverse_output[[1]]$`sd(+/-)`[1:3]))
  cat(paste0("\033[0;32m","###############Calculate analytical Jacobian for each parameter for each wavelength################","\033[0m","\n"))
  
  # test = saber.forward.jacobian.analytical(chl = synthetic_data_shallow$C_ph[i],
  #                                          acdom.440 = synthetic_data_shallow$a_cdom.440[i], 
  #                                          anap.440 = synthetic_data_shallow$a.nap.440[i],
  #                                          wavelength = wlrange, verbose = F)
  # 
  # write.csv(test, file = paste0("./Jacobian/deep/jacobian_chl-",synthetic_data_shallow$C_ph[i], "_cdom-", synthetic_data_shallow$a_cdom.440[i], "_nap-",synthetic_data_shallow$a.nap.440[i],".csv"),
  #           row.names = F, quote = F, sep = ",")
  # 
  # cat(paste0("\033[0;32m","###############analytical Jacobian for each parameter is saved to disc################","\033[0m","\n"))
  
  cat(paste0("\033[0;34m","############### INVERSION FINISHED for ", i, " number spectra ################","\033[0m","\n"))
}

Fit.optimized.ssobj.batch.final = Fit.optimized.ssobj.batch[-1,]
Fit.optimized.ssobj.batch.final$chl.actual = synthetic_data_shallow$C_ph
Fit.optimized.ssobj.batch.final$acdom.440.actual = synthetic_data_shallow$a_cdom.440
Fit.optimized.ssobj.batch.final$anap.440.actual = synthetic_data_shallow$a.nap.440

Fit.optimized.ssobj.batch.final <- cbind(Fit.optimized.ssobj.batch.final, 
             synthetic_data_shallow[!names(synthetic_data_shallow) %in% 
                                      names(Fit.optimized.ssobj.batch.final)])

lm_eqn <- function(df){
  m <- lm(anap.440 ~ a.nap.440, df);
  eq <- substitute(italic(m) ==  b %% ","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

#ANAP.440
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

coplabel <- expression(paste(italic("a")["NAP"](lambda=443),italic("synthetic")))
hl.el.label <- expression(paste(italic("a")["NAP"](lambda=443),italic("model")))
hl.inel.label <- expression(paste(italic("R")["rs"]("0"^"-", lambda),"HL"[italic("inelastic")]))
albert.label <- expression(paste(italic("R")["rs"]("0"^"-", lambda),"A&M"[italic("2003")]))

xlbl <- coplabel; ylbl <-hl.el.label
ymin <- 0; ymax <- 0.12 ; ystp <- ymax/4
xmin <- 0; xmax <- 0.12 ; xstp <- ymax/4
asp_rat <- (xmax-xmin)/(ymax-ymin)


errorb <- qnorm(0.975)*sd(Fit.optimized.ssobj.batch.final$anap.440)/sqrt(length(Fit.optimized.ssobj.batch.final$anap.440))

idx = is.na(Fit.optimized.ssobj.batch.final$anap440.sd) 
Fit.optimized.ssobj.batch.final$anap440.sd[idx] = errorb

g <- ggplot(Fit.optimized.ssobj.batch.final,aes(x = a.nap.440)) + 
  geom_point(aes(y = anap.440), shape=21, fill="cyan", 
             size=5.0, na.rm = T, show.legend = TRUE) +
  geom_abline(slope = 1,linetype="solid", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  # geom_text(x = -Inf, y = Inf, hjust = 0, vjust = 1,
  #           label = lm_eqn(Fit.optimized.ssobj.batch.final), parse = TRUE)+
  geom_errorbar( aes(ymin = anap.440 - (2.5*anap440.sd),
                     ymax = anap.440 + (2.5*anap440.sd)),
                color="black",width = 0.0025, size=1.3)+
  geom_ribbon(aes(ymin = a.nap.440 - errorb,
                  ymax = a.nap.440 + errorb),
              alpha = 0.6,
              fill = "grey70",
              colour="NA"
  )+
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
ggsave(paste0("aNAP.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)


#ACDOM440
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

coplabel <- expression(paste(italic("a")["CDOM"](lambda=443),italic("synthetic")))
hl.el.label <- expression(paste(italic("a")["CDOM"](lambda=443),italic("model")))
hl.inel.label <- expression(paste(italic("R")["rs"]("0"^"-", lambda),"HL"[italic("inelastic")]))
albert.label <- expression(paste(italic("R")["rs"]("0"^"-", lambda),"A&M"[italic("2003")]))

xlbl <- coplabel; ylbl <-hl.el.label
ymin <- 0; ymax <- 5 ; ystp <- ymax/5
xmin <- 0; xmax <- 5 ; xstp <- ymax/5
asp_rat <- (xmax-xmin)/(ymax-ymin)


errorb <- qnorm(0.975)*sd(Fit.optimized.ssobj.batch.final$acdom.440)/sqrt(length(Fit.optimized.ssobj.batch.final$acdom.440))

idx = is.na(Fit.optimized.ssobj.batch.final$acdom.440.sd) 
Fit.optimized.ssobj.batch.final$acdom.440.sd[idx] = errorb

g <- ggplot(Fit.optimized.ssobj.batch.final,aes(x = a_cdom.440)) + 
  geom_point(aes(y = acdom.440), shape=21, fill="goldenrod2", 
             size=5.0, na.rm = T, show.legend = TRUE) +
  geom_abline(slope = 1,linetype="solid", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  geom_errorbar( aes(ymin = acdom.440 - (2.5*acdom.440.sd),
                     ymax = acdom.440 + (2.5*acdom.440.sd)),
                 color="black",width = 0.1, size=1.3)+
  geom_ribbon(aes(ymin = a_cdom.440 - errorb,
                  ymax = a_cdom.440 + errorb),
              alpha = 0.6,
              fill = "grey70",
              colour="NA"
  )+

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
ggsave(paste0("aCDOM.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

#CHL
legend_title <- element_blank()
legend_position <- c(0.60, 0.40)

coplabel <- expression(paste("[",italic("chl"),"]",italic("synthetic")))
hl.el.label <- expression(paste("[",italic("chl"),"]",italic("model")))
hl.inel.label <- expression(paste(italic("R")["rs"]("0"^"-", lambda),"HL"[italic("inelastic")]))
albert.label <- expression(paste(italic("R")["rs"]("0"^"-", lambda),"A&M"[italic("2003")]))

xlbl <- coplabel; ylbl <-hl.el.label
ymin <- 0; ymax <- 10 ; ystp <- ymax/5
xmin <- 0; xmax <- 10 ; xstp <- ymax/5
asp_rat <- (xmax-xmin)/(ymax-ymin)


errorb <- qnorm(0.975)*sd(Fit.optimized.ssobj.batch.final$chl)/sqrt(length(Fit.optimized.ssobj.batch.final$chl))

idx = is.na(Fit.optimized.ssobj.batch.final$chl.sd) 
Fit.optimized.ssobj.batch.final$chl.sd[idx] = errorb

g <- ggplot(Fit.optimized.ssobj.batch.final,aes(x = C_ph)) + 
  geom_point(aes(y = chl), shape=21, fill="#568203", 
             size=5.0, na.rm = T, show.legend = TRUE) +
  geom_abline(slope = 1,linetype="solid", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  geom_errorbar( aes(ymin = chl - (2.5*chl.sd),
                     ymax = chl + (2.5*chl.sd)),
                 color="black",
                 
                 width = 0.1, 
                 size=1.3
                 
                 )+
  geom_ribbon(aes(ymin = C_ph - errorb,
                  ymax = C_ph + errorb),
              alpha = 0.6,
              fill = "grey70",
              colour="NA"
  )+

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
ggsave(paste0("chl.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)




