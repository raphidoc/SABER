#======================================================================================
# Functions to create scatterplots for inversion retrieved parameter validation
#======================================================================================


library(ggplot2)
library(ggforce)
library(GGally)
library(ggExtra)


plot_inversion_validation_singlevar_density <- function(input_df, xmin, xmax, xlabel, ylabel, 
                                                        uncertainty = "CI", opacity = 0.3){
  
  xmin = xmin; xmax = xmax; xstp <- xmax/4
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g <- ggplot(input_df,aes(x = actual)) + 
    
    stat_density_2d(aes(y = predicted, fill = ..level..),
                    geom = "polygon", colour="gray", show.legend = F, alpha = 0.5, size=1.1)+
    
    scale_fill_distiller(palette = 2, direction = 1)+
    
    geom_point(aes(y = predicted), shape=21, fill="goldenrod2", 
               size=3.0, na.rm = T, show.legend = F) +
    
    geom_ribbon(aes(ymin = predicted - .data[[uncertainty]],
                    ymax = predicted + .data[[uncertainty]]),
                alpha = opacity,
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


plot_inversion_validation_multvar_encircle <- function(input_df, xmin, xmax, xlabel, ylabel, 
                                                       uncertainty = "CI", opacity = 0.3){
  
  xmin = xmin; xmax = xmax; xstp <- xmax/4
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g <- ggplot(input_df,aes(x = actual)) + 
    # stat_density2d(aes(y = predicted,
    #                    fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
    # scale_fill_continuous(low="green",high="red") +
    
    # stat_density_2d(aes(y = predicted, 
    #                     fill = NA
    #                     ),
    #                 geom = "polygon", colour="black", show.legend = F, alpha = 0.0, size=1.1)+
    # 
    # scale_fill_distiller(palette = 2, direction = 1)+
    
  geom_point(aes(y = predicted, shape = inel_stat, fill = inel_stat), 
             size=3.0, na.rm = T, show.legend = T) +
    
    # stat_ellipse(aes(y = predicted, colour = inel_stat), size=1.1,level = 0.65, segments = 10)+
    # stat_ellipse(aes(y = predicted, colour = inel_stat), size=1.1,level = 0.95, segments = 10)+
    # stat_ellipse(aes(y = predicted, colour = inel_stat), size=1.1,level = 0.99, segments = 10)+
    
    geom_encircle(aes(y = predicted, colour = inel_stat), size=3)+
    
    scale_shape_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                              expression(paste(italic("R")["rs"]^"el")))),
                       values = c(21,23))+
    scale_fill_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                             expression(paste(italic("R")["rs"]^"el")))),
                      values = c("goldenrod2", "navyblue"))+
    scale_colour_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                               expression(paste(italic("R")["rs"]^"el")))),
                        values = c("goldenrod2", "navyblue"))+
    
    geom_ribbon(aes(ymin = predicted - .data[[uncertainty]],
                    ymax = predicted + .data[[uncertainty]],
                    fill = inel_stat),
                alpha = opacity, show.legend = F,
                
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
    annotation_logticks(short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm"))+
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


# library(ks)
# library(tidyverse)
# 
# logticks <- function(datavar,type) {
#   
#   minimum <- 1/10^abs(floor(log10(min(datavar, na.rm=TRUE))))
#   maximum <- 1*10^abs(floor(log10(max(datavar, na.rm=TRUE)))+1)
#   multiple <- floor(log10(maximum/minimum))
#   
#   yourtickvector <- c()
#   
#   if (type=="breaks") {
#     
#     yourtickvector <- c(minimum)
#     
#     for (x in seq(0,multiple)) {
#       
#       andadd <- seq(minimum*10^x,minimum*10^(x+1),minimum*10^x)[-1]
#       
#       yourtickvector <- c(yourtickvector,andadd)
#       
#     }
#     
#   } else if (type=="labels") {
#     
#     for (x in seq(0,multiple)) {
#       
#       andadd <- c(minimum*10^x,rep("",8))
#       
#       yourtickvector <- c(yourtickvector,andadd)
#       
#     }
#     
#     yourtickvector <- c(yourtickvector,minimum*10^multiple)
#     
#   }
#   
#   return(yourtickvector)
#   
# }
# 
# breaks_log10 <- function(x) {
#   low <- floor(log10(min(x)))
#   high <- ceiling(log10(max(x)))
#   
#   10^(seq.int(low, high,by = (high - low)/5))
# }
# 
# breaks_5log10 <- function(x) {
#   low <- floor(log10(min(x)/5))
#   high <- ceiling(log10(max(x)/5))
#   
#   5 * 10^(seq.int(low, high,by = (high - low)/5))
# }
# 
# plot_inversion_validation_multvar_contour <- function(input_df, xmin, xmax, xlabel, ylabel, 
#                                                       uncertainty = "CI", opacity = 0.3){
#   
#   xmin = xmin; xmax = xmax; xstp <- xmax/5
#   
#   asp_rat <- (xmax-xmin)/(ymax-ymin)
#   ## data
#   d <- input_df %>% 
#     as_tibble() 
#   
#   ## density function
#   kd <- ks::kde(d[,1:2], compute.cont=TRUE, h=0.2)
#   
#   ## extract results
#   get_contour <- function(kd_out=kd, prob="5%") {
#     contour_95 <- with(kd_out, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                             z=estimate, levels=cont[prob])[[1]])
#     as_tibble(contour_95) %>% 
#       mutate(prob = prob)
#   }
#   
#   dat_out <- map_dfr(c("10%", "68%", "95%"), ~get_contour(kd, .)) %>% 
#     group_by(prob) %>% 
#     mutate(n_val = 1:n()) %>% 
#     ungroup()
#   
#   ## clean kde output
#   kd_df <- expand_grid(x=kd$eval.points[[1]], y=kd$eval.points[[2]]) %>% 
#     mutate(z = c(kd$estimate %>% t))
#   
#   g <- ggplot(data=kd_df, aes(x, y)) +
#     #geom_contour(aes(z = z), col = "black")+
#     geom_point(data = d, aes(actual, predicted, shape = inel_stat, fill = inel_stat), 
#                alpha = I(0.4), size = I(3)) +
#     geom_path(aes(x, y, group = prob),
#               data=filter(dat_out, !n_val %in% -1:-3), colour = I("black"), size = 0.8) +
#     # 
#     geom_text(aes(label = prob), data =
#                 filter(dat_out, (prob%in% c("10%", "68%") & n_val==15)| 
#                          (prob%in% c("95%") & n_val==30)),
#               colour = I("red"), size =I(5))+
#     # 
#     # geom_encircle(data = d, aes(x = actual, y = predicted, colour = inel_stat), size=3)+
#     
#     scale_shape_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
#                                               expression(paste(italic("R")["rs"]^"el")))),
#                        values = c(21,23))+
#     scale_fill_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
#                                              expression(paste(italic("R")["rs"]^"el")))),
#                       values = c("goldenrod2", "navyblue"))+
#     scale_colour_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
#                                                expression(paste(italic("R")["rs"]^"el")))),
#                         values = c("goldenrod2", "navyblue"))+
#     
#     geom_ribbon(data = d, aes(x = actual, y = predicted, ymin = predicted - .data[[uncertainty]],
#                               ymax = predicted + .data[[uncertainty]],
#                               fill = inel_stat),
#                 alpha = opacity, show.legend = F,
#                 colour="NA"
#     )+
#     geom_abline(slope = 1,linetype="solid", intercept = 0,
#                 colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
#     
#     coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
#                 ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
#     
#     
#     scale_x_log10(
#       #breaks = breaks_log10,
#       breaks = scales::trans_breaks("log10", function(x) 10^x),
#       minor_breaks = breaks_5log10,
#       labels = scales::trans_format("log10", scales::math_format(10^.x))
#     ) +
#     scale_y_log10(
#       #breaks = breaks_log10,
#       breaks = scales::trans_breaks("log10", function(x) 10^x),
#       minor_breaks = breaks_5log10,
#       labels = scales::trans_format("log10", scales::math_format(10^.x))
#     ) +
#     
#     # scale_x_log10(
#     #   breaks = logticks(input_df$actual,"breaks"),
#     #   labels = logticks(input_df$actual,"labels")
#     # ) +
#     # scale_y_log10(
#     #   breaks = logticks(input_df$predicted,"breaks"),
#     #   labels = logticks(input_df$predicted,"labels")
#     # ) +
#     
#     xlab(xlbl)+
#     ylab(ylbl)+
#     #annotation_logticks()+
#     #scale_fill_viridis()+
#     theme_bw() +
#     theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#           axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
#           axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
#           axis.title.x = element_text(size = 25),
#           axis.title.y = element_text(size = 25),
#           axis.ticks.length = unit(.25, "cm"),
#           legend.box.just = "right",
#           legend.spacing = unit(-0.5, "cm"),
#           legend.position = legend_position,
#           #legend.direction = "vertical",
#           legend.title = element_blank(),
#           legend.text = element_text(colour = "black", size = 20, face = "plain"),
#           legend.background = element_rect(fill = NA, size = 0.5, 
#                                            linetype = "solid", colour = 0),
#           legend.key = element_blank(),
#           legend.justification = c("left", "top"),
#           panel.background = element_blank(),
#           panel.grid.major = element_line(colour = "black", 
#                                           size = 0.5, linetype = "dotted"), 
#           panel.grid.minor = element_line(colour = "grey80", 
#                                           linewidth =  0.2, linetype = "solid"),
#           plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
#           legend.direction = "vertical", legend.box = "vertical",
#           legend.text.align = 0,
#           panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
#   return(g)
# }


plot_inversion_validation_multvar_contour <- function(input_df, xmin, xmax, xlabel, ylabel, 
                                                      uncertainty = "CI", opacity = 0.3,
                                                      show_legend){
  
  xmin = xmin; xmax = xmax
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  
  # Define the range and calculate the number of breaks
  min_val <- log10(xmin)
  max_val <- log10(xmax)
  
  num_breaks <- 4
  
  by_value <- (max_val - min_val) / (num_breaks)# - 1)
  
  logrange = seq(-10,10, by = signif(by_value, digits = 0))
  breaks <- 10^(logrange)
  minor_breaks <- rep(1:9, length(breaks))*(10^rep(logrange, each=9))
  
  ## data
  d <- input_df %>% 
    as_tibble() 
  
  ## density function
  kd <- ks::kde(d[,1:2], compute.cont=TRUE, h=0.2)
  
  ## extract results
  get_contour <- function(kd_out=kd, prob="5%") {
    contour_95 <- with(kd_out, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                            z=estimate, levels=cont[prob])[[1]])
    as_tibble(contour_95) %>% 
      mutate(prob = prob)
  }
  
  dat_out <- map_dfr(c("10%", "68%", "95%"), ~get_contour(kd, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup()
  
  ## clean kde output
  kd_df <- expand_grid(x=kd$eval.points[[1]], y=kd$eval.points[[2]]) %>% 
    mutate(z = c(kd$estimate %>% t))
  
  g<-   ggplot(data = d, aes(x = actual, y = predicted,
                             fill = as.factor(inel_stat))) +
    #geom_contour(aes(z = z), col = "black")+
    
    geom_density_2d(data = d, aes(x = actual, y = predicted),na.rm = T,
                    linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
    
    geom_point(data = d, aes(actual, predicted, shape = as.factor(inel_stat),
                             fill = as.factor(inel_stat)), 
               alpha = I(0.4), size = I(3), show.legend = show_legend) +
    # geom_path(aes(x, y, group = prob),
    #           data=filter(dat_out, !n_val %in% -1:-3), 
    #           colour = I("black"), size = 0.8, show.legend = F)+
    # 
    # 
    # geom_text(aes(label = prob), data =
    #             filter(dat_out, (prob%in% c("10%", "68%") & n_val==15)| 
    #                      (prob%in% c("95%") & n_val==30)),
    #           colour = I("red"), size =I(5), show.legend = F)+
    # 
    geom_encircle(data = d, aes(x = actual, y = predicted, 
                                colour = as.factor(inel_stat)), size=3,  
                  linetype = "dotted",  show.legend = F)+
    
    scale_shape_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                              expression(paste(italic("R")["rs"]^"el")))),
                       values = c(21,23))+
    scale_fill_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                             expression(paste(italic("R")["rs"]^"el")))),
                      values = c("goldenrod2", "navyblue"))+
    scale_colour_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                               expression(paste(italic("R")["rs"]^"el")))),
                        values = c("goldenrod2", "navyblue"))+
    
    geom_ribbon(data = d, aes(x = actual, y = predicted, ymin = predicted - .data[[uncertainty]],
                              ymax = predicted + .data[[uncertainty]],
                              fill = as.factor(inel_stat)),
                alpha = opacity, show.legend = F,
                colour="NA"
    )+
    geom_abline(slope = 1,linetype="solid", intercept = 0,
                colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
    
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
    
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(size = 1) +
    
    xlab(xlbl)+
    ylab(ylbl)+
    
    theme_bw() +
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
          panel.grid.minor = element_line(colour = "grey80", 
                                          linewidth =  0.2, linetype = "solid"),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          legend.direction = "vertical", legend.box = "vertical",
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  return(g)
}


# ===================================================================================
# FINAL version for multivariate log scale scatter plot for inversion
# ===================================================================================

plot_inversion_validation_multvar_contour_final <- function(input_df, xmin, xmax, 
                                                            xlabel, ylabel, 
                                                            uncertainty = "CI", opacity = 0.3,
                                                            show_legend){
  ## data
  d <- input_df %>% 
    as_tibble() 
  
  xmin = xmin; xmax = xmax
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  
  # Define the range and calculate the number of breaks
  min_val <- log10(xmin)
  max_val <- log10(xmax)
  
  num_breaks <- 4
  
  by_value <- (max_val - min_val) / (num_breaks)# - 1)
  
  logrange = seq(-10,10, by = signif(by_value, digits = 0))
  breaks <- 10^(logrange)
  minor_breaks <- rep(1:9, length(breaks))*(10^rep(logrange, each=9))
  
  g<-   ggplot( data = d, aes(x = actual, y = predicted, colour = as.factor(inel_stat),
                              fill = as.factor(inel_stat))) +
    #geom_contour(aes(z = z), col = "black")+
    
    geom_density_2d(data = d, aes(x = actual, y = predicted), na.rm = T, bins = 6,
                    linewidth = 0.5,  show.legend = F, size=1.1)+
    
    geom_point(data = d, aes(actual, predicted, shape = as.factor(inel_stat),
                             fill = as.factor(inel_stat)), 
               alpha = I(0.4), size = I(3), show.legend = show_legend) +
    # geom_path(aes(x, y, group = prob),
    #           data=filter(dat_out, !n_val %in% -1:-3), 
    #           colour = I("black"), size = 0.8, show.legend = F)+
    # 
    # 
    # geom_text(aes(label = prob), data =
    #             filter(dat_out, (prob%in% c("10%", "68%") & n_val==15)| 
    #                      (prob%in% c("95%") & n_val==30)),
    #           colour = I("red"), size =I(5), show.legend = F)+
    # 
    # geom_encircle(data = d, aes(x = actual, y = predicted, 
  #                             colour = (inel_stat)), size=3,  
  #               linetype = "dotted",  show.legend = F)+
  
  scale_shape_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                            expression(paste(italic("R")["rs"]^"el")))),
                     values = c(21,23))+
    scale_fill_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                             expression(paste(italic("R")["rs"]^"el")))),
                      values = c("goldenrod2", "navyblue"))+
    scale_colour_manual(name ="", labels=rev(c(expression(paste(italic("R")["rs"]^"obs")),
                                               expression(paste(italic("R")["rs"]^"el")))),
                        values = c("goldenrod2", "navyblue"))+
    
    geom_ribbon(data = d, aes(x = actual, y = predicted, ymin = predicted - .data[[uncertainty]],
                              ymax = predicted + .data[[uncertainty]],
                              fill = as.factor(inel_stat)),
                alpha = opacity, show.legend = F,
                colour="NA"
    )+
    
    geom_abline(slope = 1,linetype="solid", intercept = 0,
                colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
    
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
    
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(size = 1) +
    
    xlab(xlbl)+
    ylab(ylbl)+
    geom_rug(size = 1.1, show.legend = show_legend, alpha = opacity)+
    theme_bw() +
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
          panel.grid.minor = element_line(colour = "grey80", 
                                          linewidth =  0.2, linetype = "solid"),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          legend.direction = "vertical", legend.box = "vertical",
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  g <- ggMarginal(groupFill = T, data = d, type = "densigram", bins = 30,
                  p = g, aes(x = actual, y = predicted))
  
  return(g)
}


# ===================================================================================
# FINAL version for uni-variate log scale scatter plot for inversion (solo cruise)
# ===================================================================================
plot_inversion_validation_singlevar_contour_singlecruise <- function(input_df, xmin, xmax,
                                                                     xlabel, ylabel, 
                                                                     uncertainty = "CI", opacity = 0.3, 
                                                                     plot_col,
                                                                     show_legend){
  
  xmin = xmin; xmax = xmax
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  
  # Define the range and calculate the number of breaks
  min_val <- log10(xmin)
  max_val <- log10(xmax)
  
  num_breaks <- 4
  
  by_value <- (max_val - min_val) / (num_breaks)# - 1)
  
  logrange = seq(-10,10, by = signif(by_value, digits = 0))
  breaks <- 10^(logrange)
  minor_breaks <- rep(1:9, length(breaks))*(10^rep(logrange, each=9))
  
  ## data
  d <- input_df %>% 
    as_tibble() 
  
  ## density function
  kd <- ks::kde(d[,1:2], compute.cont=TRUE, h=0.2)
  
  ## extract results
  get_contour <- function(kd_out=kd, prob="5%") {
    contour_95 <- with(kd_out, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                            z=estimate, levels=cont[prob])[[1]])
    as_tibble(contour_95) %>% 
      mutate(prob = prob)
  }
  
  dat_out <- map_dfr(c("10%", "68%", "95%"), ~get_contour(kd, .)) %>% 
    group_by(prob) %>% 
    mutate(n_val = 1:n()) %>% 
    ungroup()
  
  ## clean kde output
  kd_df <- expand_grid(x=kd$eval.points[[1]], y=kd$eval.points[[2]]) %>% 
    mutate(z = c(kd$estimate %>% t))
  
  
  cols <- as.vector(c("#001bff", "seagreen")) #, "#610000", "#610000"))
  
  cols <- colorRampPalette(cols)(nrow(d))
  
  g<-   ggplot(data=d, aes(x = actual, y = predicted)) +
    #geom_contour(aes(z = z), col = "black")+
    
    
    geom_density_2d(data = d, aes(x = actual, y = predicted),na.rm = T, bins = 6,
                    linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
    
    geom_point(data = d, aes(actual, predicted), shape=21, fill=plot_col,
               alpha = I(0.4), size = I(3), show.legend = show_legend) +
    
    
    geom_ribbon(data = d, aes(x = actual, y = predicted, ymin = predicted - .data[[uncertainty]],
                              ymax = predicted + .data[[uncertainty]]),
                fill = "navyblue",
                alpha = 0.25, show.legend = F,
                colour="NA"
    )+
    geom_rug(size = 1.1, show.legend = T, alpha = opacity, 
             colour = heat.colors(n = length(d$actual)))+
    
    geom_abline(slope = 1,linetype="solid", intercept = 0,
                colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
    
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
    
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(size = 1) +
    
    xlab(xlbl)+
    ylab(ylbl)+
    
    theme_bw() +
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
          panel.grid.minor = element_line(colour = "grey80", 
                                          linewidth =  0.2, linetype = "solid"),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          legend.direction = "vertical", legend.box = "vertical",
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  g <- ggMarginal(fill = plot_col, data = d, type = "densigram", bins = 30, alpha = 0.7,
                  p = g, aes(x = actual, y = predicted))
  return(g)
}


# ===================================================================================
# FINAL version for uni-variate log scale scatter plot for inversion (multi cruise)
# ===================================================================================
plot_inversion_validation_singlevar_contour_multcruise <- function(input_df, plot_insitu=T,
                                                                   xmin, xmax, 
                                                                   xlabel, ylabel, 
                                                                   uncertainty = "CI", opacity = 0.3, 
                                                                   plot_col, 
                                                                   label_list = c("GLORIA","MERMAID","SEABASS", "NOMAD", "AERONET-OC"),
                                                                   shape_list = c(21,22,23,24,25),
                                                                   show_legend, hist_count = 30){
  
  xmin = xmin; xmax = xmax
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  
  # Define the range and calculate the number of breaks
  min_val <- log10(xmin)
  max_val <- log10(xmax)
  
  num_breaks <- 4
  
  by_value <- (max_val - min_val) / (num_breaks)# - 1)
  
  logrange = seq(-10,10, by = signif(by_value, digits = 0))
  breaks <- 10^(logrange)
  minor_breaks <- rep(1:9, length(breaks))*(10^rep(logrange, each=9))
  
  ## data
  d <- input_df %>% 
    as_tibble() 
  
  cols <- as.vector(c("#001bff", "seagreen")) #, "#610000", "#610000"))
  
  cols <- colorRampPalette(cols)(nrow(d))
  
  if (plot_insitu == TRUE) {
    
    g<-   ggplot(data=d, aes(x = actual, y = predicted)) +
      #geom_contour(aes(z = z), col = "black")+
      
      geom_point(data = d, aes(actual, predicted, shape = cruise),  fill=plot_col, 
                 #color = plot_col,
                 alpha = I(opacity), size = I(3), show.legend = show_legend) +
      
      geom_density_2d(data = d, aes(x = actual, y = predicted),na.rm = T, bins = 6,
                      linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
      
      geom_ribbon(data = d, aes(x = actual, y = predicted, ymin = predicted - .data[[uncertainty]],
                                ymax = predicted + .data[[uncertainty]]),
                  fill = "navyblue",
                  alpha = opacity, show.legend = F,
                  colour="NA"
      )+
      scale_shape_manual(name ="", labels=(label_list),
                         values = (shape_list))+
      geom_rug(size = 1.1, show.legend = F, alpha = 0.25,
               color = heat.colors(n = length(d$actual)) 
               
      ) +
      
      geom_abline(slope = 1,linetype="dashed", intercept = 0,
                  colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
      
      geom_smooth(size=1,level = 0.95,show.legend = F,linetype = "solid", color="red4", 
                  se= T, method = "lm")+
      
      coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                  ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
      
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      annotation_logticks(size = 1) +
      
      xlab(xlbl)+
      ylab(ylbl)+
      
      theme_bw() +
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
            panel.grid.minor = element_line(colour = "grey80", 
                                            linewidth =  0.2, linetype = "solid"),
            plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
            legend.direction = "vertical", legend.box = "vertical",
            legend.text.align = 0,
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
    
    g <- ggMarginal(fill = plot_col, data = d, type = "densigram", bins = hist_count, alpha = 0.7,
                    p = g, aes(x = actual, y = predicted))
    
  } else {
    g<-   ggplot(data=d, aes(x = actual, y = predicted)) +
      #geom_contour(aes(z = z), col = "black")+
      
      geom_point(data = d, aes(actual, predicted, shape = cruise),  fill=plot_col,
                 #color = plot_col,
                 alpha = I(opacity), size = I(3), show.legend = show_legend) +
      
      geom_density_2d(data = d, aes(x = actual, y = predicted),na.rm = T, bins = 6,
                      linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
      
      geom_ribbon(data = d, aes(x = actual, y = predicted, ymin = predicted - .data[[uncertainty]],
                                ymax = predicted + .data[[uncertainty]]),
                  fill = "navyblue",
                  alpha = opacity, show.legend = F,
                  colour="NA"
      )+
      scale_shape_manual(name ="", labels=(c("IOCCG2003","Loisel2023")),
                         values = c(21,22))+
      geom_rug(size = 1.1, show.legend = F, alpha = 0.25,
               color = heat.colors(n = length(d$actual)) 
               
      ) +
      
      geom_abline(slope = 1,linetype="solid", intercept = 0,
                  colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
      geom_smooth(size=1,level = 0.95, show.legend = F,linetype = "solid", color="red4", 
                  se= T, method = "lm")+
      
      coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                  ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
      
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      annotation_logticks(size = 1) +
      
      xlab(xlbl)+
      ylab(ylbl)+
      
      theme_bw() +
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
            panel.grid.minor = element_line(colour = "grey80", 
                                            linewidth =  0.2, linetype = "solid"),
            plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
            legend.direction = "vertical", legend.box = "vertical",
            legend.text.align = 0,
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
    
    g <- ggMarginal(fill = plot_col, data = d, type = "densigram", bins = hist_count, alpha = 0.7,
                    p = g, aes(x = actual, y = predicted))
  }
  
  
  return(g)
}

#=======================================================================
#SAME as last function but for OWT type
#=======================================================================
plot_inversion_validation_singlevar_contour_multcruise_OWT <- function(input_df, plot_insitu=T,
                                                                       xmin, xmax, 
                                                                       xlabel, ylabel, 
                                            col_lab = col_list,
                                            uncertainty = "CI", opacity = 0.3, 
                                            plot_col, 
                                            label_list = c(paste0("OWT", seq(1,10,1))),
                                         #shape_list = c(21,22,23,24,25),
                                           shape_list = c(21,22,23,24,25),
                                           cruise_lab = c("GLORIA","MERMAID","SEABASS", "NOMAD",
                                                                                      "AERONET-OC"),
                                           show_legend, hist_count = 30){
  
  xmin = xmin; xmax = xmax
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  
  # Define the range and calculate the number of breaks
  min_val <- log10(xmin)
  max_val <- log10(xmax)
  
  num_breaks <- 4
  
  by_value <- (max_val - min_val) / (num_breaks)# - 1)
  
  logrange = seq(-10,10, by = signif(by_value, digits = 0))
  breaks <- 10^(logrange)
  minor_breaks <- rep(1:9, length(breaks))*(10^rep(logrange, each=9))
  
  ## data
  d <- input_df %>% 
    as_tibble() 
  
  cols <- as.vector(c("#001bff", "seagreen")) #, "#610000", "#610000"))
  
  cols <- colorRampPalette(cols)(nrow(d))
  
  if (plot_insitu == TRUE) {
    
    g<-   ggplot(data=d, aes(x = actual, y = predicted)) +
      #geom_contour(aes(z = z), col = "black")+
      
      geom_point(data = d, aes(actual, predicted, shape = cruise, fill=as.factor(clust_id)), 
                 #color = plot_col,
                 alpha = I(opacity), size = I(3), show.legend = show_legend) +
      
      geom_density_2d(data = d, aes(x = actual, y = predicted),na.rm = T, bins = 6,
                      linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
      
      geom_ribbon(data = d, aes(x = actual, y = predicted, ymin = predicted - .data[[uncertainty]],
                                ymax = predicted + .data[[uncertainty]]),
                  fill = "navyblue",
                  alpha = opacity, show.legend = F,
                  colour="NA"
      )+
      scale_shape_manual(name ="", labels=(cruise_lab),
                         values = (shape_list))+
      scale_fill_manual(name ="", values = col_lab,  labels = label_list)+
      geom_rug(size = 1.1, show.legend = F, alpha = 0.25,
               color = heat.colors(n = length(d$actual)) 
               
      ) +
      
      geom_abline(slope = 1,linetype="dashed", intercept = 0,
                  colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
      
      geom_smooth(size=1,level = 0.95,show.legend = F,linetype = "solid", color="red4", 
                  se= T, method = "lm")+
      
      coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                  ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
      
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      annotation_logticks(size = 1) +
      
      xlab(xlbl)+
      ylab(ylbl)+
      
      theme_bw() +
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
            panel.grid.minor = element_line(colour = "grey80", 
                                            linewidth =  0.2, linetype = "solid"),
            plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
            legend.direction = "vertical", legend.box = "vertical",
            legend.text.align = 0,
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
    g <- g + guides(fill = "none")
    
    g <- ggMarginal(fill = plot_col, data = d, type = "densigram", bins = hist_count, alpha = 0.7,
                    p = g, aes(x = actual, y = predicted))
    
  } else {
    
    g<-   ggplot(data=d, aes(x = actual, y = predicted)) +
      #geom_contour(aes(z = z), col = "black")+
      
      geom_point(data = d, aes(actual, predicted, shape = cruise, fill=as.factor(clust_id)), 
                 #color = plot_col,
                 alpha = I(opacity), size = I(3), show.legend = show_legend) +
      
      geom_density_2d(data = d, aes(x = actual, y = predicted),na.rm = T, bins = 6,
                      linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
      
      geom_ribbon(data = d, aes(x = actual, y = predicted, ymin = predicted - .data[[uncertainty]],
                                ymax = predicted + .data[[uncertainty]]),
                  fill = "navyblue",
                  alpha = opacity, show.legend = F,
                  colour="NA"
      )+
      scale_shape_manual(name ="", labels=(cruise_lab),
                         values = (shape_list))+
      scale_fill_manual(name ="", values = col_lab,  labels = label_list)+
      geom_rug(size = 1.1, show.legend = F, alpha = 0.25,
               color = heat.colors(n = length(d$actual)) 
               
      ) +
      
      geom_abline(slope = 1,linetype="dashed", intercept = 0,
                  colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
      
      geom_smooth(size=1,level = 0.95,show.legend = F,linetype = "solid", color="red4", 
                  se= T, method = "lm")+
      
      coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                  ylim = c(ymin, ymax), expand = FALSE, clip = "on")+
      
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(breaks = breaks, minor_breaks = minor_breaks,
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      annotation_logticks(size = 1) +
      
      xlab(xlbl)+
      ylab(ylbl)+
      
      theme_bw() +
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
            panel.grid.minor = element_line(colour = "grey80", 
                                            linewidth =  0.2, linetype = "solid"),
            plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
            legend.direction = "vertical", legend.box = "vertical",
            legend.text.align = 0,
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
    
    g <- g + guides(fill = "none")
    
    g <- ggMarginal(fill = plot_col, data = d, type = "densigram", bins = hist_count, alpha = 0.7,
                    p = g, aes(x = actual, y = predicted))
    
  }
  
  
  return(g)
}



# ===================================================================================
# FINAL version for uni-variate linear scale scatter plot for inversion
# ===================================================================================
plot_inversion_validation_singlevar_linear_contour <- function(input_df, 
                                                               xmin, xmax, xlabel, 
                                                               ylabel, 
                                                               uncertainty = "CI", opacity = 0.3, 
                                                               plot_col, xstp,
                                                               #label_list = c("GLORIA","MERMAID","SEABASS", "NOMAD", "AERONET-OC"),
                                                               #shape_list = c(21,22,23,24,25),
                                                               show_legend, hist_count = 30){
  
  ## data
  d <- input_df %>% 
    as_tibble() 
  
  ymin = xmin; ymax = xmax
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g<-   ggplot(data=d, aes(x = H_actual, y = H_predicted)) +
    #geom_contour(aes(z = z), col = "black")+
    
    
    geom_density_2d(data = d, aes(x = H_actual, y = H_predicted),na.rm = T, bins = 6,
                    linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
    
    geom_point(data = d, aes(H_actual, H_predicted), shape=21, fill=plot_col,
               alpha = I(opacity), size = I(3), show.legend = show_legend) +
    
    
    geom_ribbon(data = d, aes(x = H_actual, y = H_predicted, 
                              ymin = (H_predicted - .data[[uncertainty]]),
                              ymax = (H_predicted + .data[[uncertainty]])),
                fill = "navyblue",
                alpha = opacity, show.legend = F,
                colour="NA"
    )+
    
    geom_rug(size = 1.1, show.legend = T, alpha = opacity, 
             colour = heat.colors(n = length(d$H_actual)))+
    
    geom_abline(slope = 1,linetype="dashed", intercept = 0,
                colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
    
    geom_smooth(size=1,level = 0.95,show.legend = F,linetype = "solid", color="red4", 
                se= T, method = "lm")+
    
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp)) +
    scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                       breaks = seq(ymin, ymax, ystp)) +
    
    theme_bw() +
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
          panel.grid.minor = element_line(colour = "grey80", 
                                          linewidth =  0.2, linetype = "solid"),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          legend.direction = "vertical", legend.box = "vertical",
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  g <- ggMarginal(fill = plot_col, data = d, type = "densigram", bins = hist_count, alpha = 0.7,
                  p = g, aes(x = H_actual, y = H_predicted))
  return(g)
}