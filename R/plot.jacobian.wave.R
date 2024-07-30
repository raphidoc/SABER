library(ggplot2)
jacobian_data = read.csv("./Jacobian/deep/jacobian_chl-1_cdom-0.5_nap-0.04.csv", 
                         header = T)
jacobian_data$abs_CDOM_440 = jacobian_data$abs_CDOM_440*10^-3
jacobian_data$abs_X_440 = jacobian_data$abs_X_440*10^-3
jacobian_data$wl = seq(400,800,10)
jacobian_data$H = rowMeans(jacobian_data[-c(1,4)])
jacobian_data$H = jacobian_data$H*10^-0.5

hl.wl = reshape2::melt(data = jacobian_data, id.vars = "wl")

hl.wl$value = abs(1/hl.wl$value)

hl.wl$value_mod = abs(log10(1/hl.wl$value))
hl.wl$value_mod[hl.wl$variable == "H"] = 1.5*hl.wl$value_mod[hl.wl$variable == "H"]
#hl.wl$value = hl.wl$value/max(hl.wl$value)

g = ggplot(hl.wl, aes(color=variable, y=value_mod, x=wl)) + 
  geom_line(stat="identity", size= 1.3)+
  scale_color_manual(values = rev(c("green4", "yellow3", "blue4", "brown4")), labels =rev(c("[chl]",
                                                                         "adg(443)",
                                                                         "bbp(555)",
                                                                         "H")))+
  xlab("Wavelength [nm]")+
  ylab("Standard error [unscaled]")+
  # scale_fill_manual(labels = c(expression(paste(italic("[chl]"))),
  #                              expression(paste("a"[italic("CDOM")](440))),
  #                              expression(paste("a"[italic("NAP")](440)))),
  #                   values = c(muted("red"),"blue", "green")) +
  
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 0),
        axis.text.y = element_text(size = 15, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        plot.caption = element_text(face="bold", hjust = 0, size = 15, color = 'black'),
        legend.position=c(0.15, 0.95),
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
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g

ggsave(paste0("./outputs/Jacobian_lambda.png"), plot = g,
       scale = 1.5, width = 6, height = 4.5, units = "in",dpi = 300)
