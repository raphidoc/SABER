
#Interpolate to user given lambda interval
wlrange <- seq(400,800,1); wlcount <- length(wlrange)
SAmat <- matrix(nrow = wlcount, ncol = dim(rrs.HL)[1], 0)
for (i in 1:dim(rrs.HL)[1]) {
  SAmat[,i] <- Hmisc::approxExtrap(as.numeric(rrs.HL.wl), rrs.HL[i,],
                                   xout = as.numeric(wlrange), 
                                   na.rm = TRUE, method = "linear")$y
}

rownames(SAmat) <- wlrange

SAmat <- as.data.frame(SAmat)
SAmat$wave <- as.numeric(wlrange)
SAmat[SAmat < 0] <- 0
require(dplyr)
SAmat <- SAmat %>% select(wave, everything())
write.csv(SAmat, file = "./copsdata.csv", sep = " ")

#Assign color to Rrs as per Balasubhranyam et al. 20
rrs.665 <- SAmat[SAmat$wave == "665",]
rrs.560 <- SAmat[SAmat$wave == "560",]
rrs.492 <- SAmat[SAmat$wave == "490",]
rrs.colorset <- rbind(rrs.665, rrs.560, rrs.492)

colorarray <- vector()
for (i in 2:ncol(rrs.colorset)) {
  if ( (rrs.colorset[[i]][rrs.colorset$wave == "665"] >  
        rrs.colorset[[i]][rrs.colorset$wave == "560"]) |  
       (rrs.colorset[[i]][rrs.colorset$wave == "665"] > 0.025)) {
    colorarray[i] <- "brown"
  } else{
    if ( rrs.colorset[[i]][rrs.colorset$wave == "560"] <  rrs.colorset[[i]][rrs.colorset$wave == "490"]) {
      colorarray[i] <- "blue"
    } else{
      colorarray[i] <- "green"
    }
  }
}


#Calculate AVW
rrs.sum <- data.frame("station" = colnames(SAmat)[-1], "rrs.sum" =colSums(SAmat)[-1])
lambda.scale <- vector(); rrs.sum.lambda <- vector()
for (i in 2:ncol(SAmat)) {
  lambda.scale <- SAmat[[i]]/SAmat$wave
  rrs.sum.lambda[i] <- sum(lambda.scale)
}
rrs.sum$rrs.sum.lambda <- rrs.sum.lambda[-1]
rrs.sum$avw <- rrs.sum$rrs.sum/rrs.sum$rrs.sum.lambda

#Calculate NDI (490nm and 665nm)
rrs.490 <- SAmat[SAmat$wave == "490",]
rrs.665 <- SAmat[SAmat$wave == "665",]
ndi.top <- (rrs.665 - rrs.490)[-1]
ndi.down <- (rrs.665 + rrs.490)[-1]
ndi <- ndi.top/ndi.down
ndi <- as.matrix(ndi)
rrs.sum$ndi <- ndi[1,]

#Calculate QWIP
p1 = -8.399885*10^-9 ; p2 = 1.715532*10^-5; p3 = -1.301670*10^-2; p4 = 4.357838*10^0; p5 = -5.449532*10^2
rrs.sum$qwip <- p1*rrs.sum$avw^4 + p2*rrs.sum$avw^3 + p3*rrs.sum$avw^2 + p4*rrs.sum$avw^1 + p5
rrs.sum$qwipscore <- rrs.sum$ndi - rrs.sum$qwip
rrs.sum$abs.qwipscore <- abs(rrs.sum$qwipscore)

rrs.bad <- rrs.sum[rrs.sum$qwipscore > 0.2 | rrs.sum$qwipscore < -0.2,]
rrs.good <- rrs.sum[rrs.sum$qwipscore < 0.2 & rrs.sum$qwipscore > -0.2,]

#assign color
rrs.sum$colorcode <- colorarray[-1]

#add depth
depthdata <- copsdata[4,]
idx <- is.na(depthdata)
depthdata[idx] <- "deep"
depthdata[,41] <- "deep"
depthdata[idx == FALSE] <- "shallow"
rrs.sum$depth <- as.matrix(t(depthdata[1,]))
#errora <- qnorm(0.975)*sd(rrs.sum$qwip)/sqrt(length(rrs.sum$qwip))

xmin <- 510; xmax <- 610; xstp <- 20
xlbl <- "AVW"
ymin <- -1; ymax <- 1 ; ystp <- 0.4
ylbl <- "NDI(490,665)"

g1 <- ggplot(data = rrs.sum)  +
  geom_line(aes(x =avw, y= qwip, color="QWIP"),size=1.3,show.legend = F)+
  geom_ribbon(aes(x =avw, ymax = qwip + 0.2, ymin = qwip - 0.2),
              alpha = 0.0,
              fill = "white",
              colour="purple",
              show.legend = F, linetype="dashed", size=1
  )+
  geom_ribbon(aes(x =avw, ymax = qwip + 0.1, ymin = qwip - 0.1),
              alpha = 0.0,
              fill = "white",
              colour="purple",
              show.legend = F, linetype="dotted", size=1
  )+
  geom_point(aes(x = avw, y=ndi, color=colorcode), size=3, show.legend = T)+ 
  scale_shape_manual(labels=c("deep(>10m)","shallow(<10m)"),
                     values = c(17,19))+
  #scale_fill_discrete("",)
  #labs(x="wavelength [nm]", y=expression(paste(bar(italic("R"))["rs, inelastic"]("0"^"-", lambda)," [sr"^"-1","]")))+
  ggtitle(" ")+
  # geom_text(aes(x=avw,y=ndi,label=as.character(station)),hjust=0, vjust=0, 
  #           angle= 0, check_overlap = T)+
  #scale_color_manual(labels="QWIP", values = c("royalblue"))+
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.caption = element_text(face="bold",size=15, hjust = 0, color = 'black'),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0),
        axis.text.y = element_text(size = 25, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0., 0.95),
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
g1
ggsave(paste0("./QWIP_cops.png"), plot = g1,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)

legend_title <- element_blank()
legend_position <- c(0.70, 0.98)
xmin <- -0.4; xmax <- 0.4; xstp <- 0.2 
xlbl <- "QWIP"
df_means <- plyr::ddply(rrs.sum, dplyr::mean_qwip = mean(qwip))

g <- ggplot(rrs.sum,aes(qwipscore)) + 
  geom_histogram(colour="grey", fill="royalblue",bins=12, alpha =1.0,
                 size=1, na.rm = FALSE, show.legend = TRUE) +
  # geom_density(linetype="solid",
  #              size=1.7, na.rm = FALSE)+
  #annotate("text", x=-0.8, y=4, label= expression(paste(italic(n)," =25")), size=10)+
  #scale_fill_manual(name = "",values = collist)+
  #scale_color_manual(name = "",values = collist)+
  ggtitle(" ")+
  geom_vline(xintercept = mean(rrs.sum$qwip), 
             linetype="dashed", size=2)+
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp)) +
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
        plot.margin = unit(c(0.0,0.9,0.0,0.0), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g
ggsave(paste0("./QWIP_hist.png"), plot = g,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)


png(filename = "./fit1.png", pointsize = 20)
plot(0)
legend(x = "center",          # Position
       legend = c("+/-0.2", "+/-0.1"),  # Legend texts
       lty = c(2, 3),           # Line types
       col = c("purple", "purple"),           # Line colors
       lwd = 2)
dev.off()
ggplot()
require(gridExtra)
require(grid)
legend = cowplot::get_legend(g)
grid.newpage()
grid.draw(legend)
ggsave("legend.png", plot = legend,
       scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
# 