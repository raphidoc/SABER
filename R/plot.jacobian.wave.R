library(ggplot2)
jacobian_data = read.csv("./Jacobian/deep/jacobian_chl-1_cdom-0.5_nap-0.04.csv", 
                         header = T)

jacobian_data$wl = wavelength
#jacobian_data$id = seq(1,41,1)

hl.wl = reshape2::melt(data = jacobian_data, id.vars = "wl")

hl.wl$value = abs(1/hl.wl$value)
hl.wl$value = hl.wl$value/max(hl.wl$value)
g = ggplot(hl.wl, aes(fill=variable, y=abs(log10(value)), x=wl)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("green4", "yellow3", "blue4"), labels = c("[chl]",
                                                                         "aCDOM(443)",
                                                                         "aNAP(443)"))+
  xlab("Wavelength [nm]")+
  ylab("Standard error [%]")+
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
        legend.position=c(0.15, 0.85),
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

ggsave(paste0("./Jacobian_lambda.png"), plot = g,
       scale = 1.5, width = 6, height = 4.5, units = "in",dpi = 300)
