library(brms)

dat <- read.csv("data.csv")
colnames(dat) <- c("d", "v")
#
# mix <- mixture("lognormal", "lognormal")
# mdl_1 <- brm(chl ~ 1, data=HL.deep.iop, family=mix)  # Using the default priors
#
# mdl_ln <- brm(chl ~ 1, data=HL.deep.iop, family="lognormal")
# plot(mdl_ln)
# pp_check(mdl_ln, nsamples = 50)


# dens = density(chl.sample)
#
# q25 <- quantile(chl.sample, .25)
# q75 <- quantile(chl.sample, .75)
#
# dd <- with(dens,data.frame(x,y))
# dd$x[dd$x <0] = 0

# library(ggplot2)
#
# chl.sample = as.data.frame(chl.sample)
#
# legend_title <- element_blank()
# legend_position <- c(0.70, 0.98)
# xmin <- 0; xmax <- signif(max(chl.sample), digits = 1); xstp <- (xmax-xmin)/5
# ymin = 0; ymax = signif((max(dens$y) + 0.2*max(dens$y)), digits = 1); ystp <- (ymax - ymin)/3
# asp_rat = (xmax - xmin)/(ymax - ymin)
# xlbl <- expression(paste(italic("[chl]")))
#
# g <- ggplot(chl.sample, aes(x=chl.sample)) +
#   geom_histogram(aes(y = ..density..), binwidth=density(chl.sample$chl.sample)$bw,
#                  colour="grey", bins=20, alpha =0.7,
#                  size=1, na.rm = FALSE) +
#   geom_density(fill="#568203", alpha = 0.2,
#                na.rm = FALSE) +
#   scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
#                      breaks = seq(xmin, xmax, xstp)) +
#   # scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
#   #                    breaks = seq(ymin, ymax, ystp)) +
#   ylab("Density")+
#   coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
#               ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#   theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         axis.text.x = element_text(size = 20, color = 'black', angle = 0),
#         axis.text.y = element_text(size = 20, color = 'black', angle = 0),
#         axis.title.x = element_text(size = 25),
#         axis.title.y = element_text(size = 25),
#         axis.ticks.length = unit(.25, "cm"),
#         legend.box.just = "right",
#         legend.spacing = unit(-0.5, "cm"),
#         legend.position = legend_position,
#         legend.direction = "vertical",
#         legend.title = element_blank(),
#         legend.text = element_text(colour = "black", size = 25, face = "plain"),
#         legend.background = element_rect(fill = NA, size = 0.5,
#                                          linetype = "solid", colour = 0),
#         legend.key = element_blank(),
#         legend.justification = c("left", "top"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour = "grey",
#                                         size = 0.5, linetype = "dotted"),
#         panel.grid.minor = element_blank(),
#         plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
# g
# ggsave(paste0("./hist_chl.png"), plot = g,
#        scale = 1.5, width = 5.5, height = 4.5, units = "in",dpi = 300)
#
#
#
#   curve(dweibull(x, shape = fit.chl.norm$estimate["shape"],
#                scale = fit.chl.norm$estimate["scale"] ), add = T)


#plot(fit.acdom440.norm) #Infer results of the fit

dens = density(acdom.440.sample)

q25 <- quantile(acdom.440.sample, .25)
q75 <- quantile(acdom.440.sample, .75)

dd <- with(dens,data.frame(x,y))
dd$x[dd$x <0] = 0

library(ggplot2)

acdom.440.sample = as.data.frame(acdom.440.sample)

legend_title <- element_blank()
legend_position <- c(0.70, 0.98)
xmin <- 0; xmax <- signif(max(acdom.440.sample), digits = 1); xstp <- (xmax-xmin)/5
ymin = 0; ymax = signif((max(dens$y) + 0.2*max(dens$y)), digits = 1); ystp <- (ymax - ymin)/3
asp_rat = (xmax - xmin)/(ymax - ymin)
xlbl <- expression(paste(italic("a")["CDOM"](lambda=443),italic("in vivo")))

g <- ggplot(acdom.440.sample, aes(x=acdom.440.sample)) +
  geom_histogram(aes(y = ..density..), binwidth=density(acdom.440.sample$acdom.440.sample)$bw,
                 colour="grey", bins=20, alpha =0.7,
                 size=1, na.rm = FALSE) +
  geom_density(fill="goldenrod2", alpha = 0.2,
               na.rm = FALSE) +
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp)) +
  # scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
  #                    breaks = seq(ymin, ymax, ystp)) +
  ylab("")+
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0),
        axis.text.y = element_text(size = 20, color = 'black', angle = 0),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = legend_position,
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 25, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5,
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey",
                                        size = 0.5, linetype = "dotted"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g
# ggsave(paste0("./hist_acdom440.png"), plot = g,
#        scale = 1.5, width = 5.5, height = 4.5, units = "in",dpi = 300)

# dens = density(anap.440.sample)
#
# q25 <- quantile(anap.440.sample, .25)
# q75 <- quantile(anap.440.sample, .75)
#
# dd <- with(dens,data.frame(x,y))
# dd$x[dd$x <0] = 0
#
# library(ggplot2)
#
# anap.440.sample = as.data.frame(anap.440.sample)
#
# legend_title <- element_blank()
# legend_position <- c(0.70, 0.98)
# xmin <- 0; xmax <- signif(max(anap.440.sample), digits = 1); xstp <- (xmax-xmin)/5
# ymin = 0; ymax = signif((max(dens$y) + 0.2*max(dens$y)), digits = 1); ystp <- (ymax - ymin)/3
# asp_rat = (xmax - xmin)/(ymax - ymin)
# xlbl <- expression(paste(italic("a")["NAP"](lambda=443),italic("in vivo")))
#
# g <- ggplot(anap.440.sample, aes(x=anap.440.sample)) +
#   geom_histogram(aes(y = ..density..), binwidth=density(anap.440.sample$anap.440.sample)$bw,
#                  colour="grey", bins=20, alpha =0.7,
#                  size=1, na.rm = FALSE) +
#   geom_density(fill="cyan", alpha = 0.2,
#                na.rm = FALSE) +
#   scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
#                      breaks = seq(xmin, xmax, xstp)) +
#   # scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
#   #                    breaks = seq(ymin, ymax, ystp)) +
#   ylab("")+
#   coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
#               ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
#   theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         axis.text.x = element_text(size = 20, color = 'black', angle = 0),
#         axis.text.y = element_text(size = 20, color = 'black', angle = 0),
#         axis.title.x = element_text(size = 25),
#         axis.title.y = element_text(size = 25),
#         axis.ticks.length = unit(.25, "cm"),
#         legend.box.just = "right",
#         legend.spacing = unit(-0.5, "cm"),
#         legend.position = legend_position,
#         legend.direction = "vertical",
#         legend.title = element_blank(),
#         legend.text = element_text(colour = "black", size = 25, face = "plain"),
#         legend.background = element_rect(fill = NA, size = 0.5,
#                                          linetype = "solid", colour = 0),
#         legend.key = element_blank(),
#         legend.justification = c("left", "top"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour = "grey",
#                                         size = 0.5, linetype = "dotted"),
#         panel.grid.minor = element_blank(),
#         plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
# g
# ggsave(paste0("./hist_anap440.png"), plot = g,
#        scale = 1.5, width = 5.5, height = 4.5, units = "in",dpi = 300)