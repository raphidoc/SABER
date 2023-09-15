# Load required libraries
library(ggplot2)
library(ggdist)
library(fitdistrplus)
library(viridis)
library(gridExtra)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(logspline)
library(qqplotr)
library(diptest)
library(mclust)
library(flexmix)
#------------------------------------------------------------
# Functions to calculate densities and quantiles
#------------------------------------------------------------
#Function for mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Function to calculate probability density functions
estimate_prior_density <- function(input_list){
  
  dx = seq(from=min(input_list$data), to=max(input_list$data), length.out=length(input_list$data))
  
  density_data <- data.frame(
    
    data_obs = input_list$data,
    
    x = seq(from=min(input_list$data), to=max(input_list$data), length.out=length(input_list$data)),
    
    weibull_density = dweibull(dx, 
                               shape = input_list$fit_weibull$estimate[1], 
                               scale = input_list$fit_weibull$estimate[2]),
    
    lognormal_density = dlnorm(dx, 
                               meanlog  = input_list$fit_lnorm$estimate[1], 
                               sdlog = input_list$fit_lnorm$estimate[2]),
    
    gamma_density = dgamma(dx,
                           shape = input_list$fit_gamma$estimate[1], 
                           rate = input_list$fit_gamma$estimate[2])
  )
  
  return(density_data)
  
}

#Function to calculate quantiles from drawn distributions
estimate_prior_quantile <- function(input_list){
  
  
  quantile_data <- data.frame(
    
    x = input_list$data,
    
    weibull_qq = qweibull(ppoints(input_list$data), 
                          shape = input_list$fit_weibull$estimate[1], 
                          scale = input_list$fit_weibull$estimate[2]),
    
    lognormal_qq = qlnorm(ppoints(input_list$data), 
                          meanlog  = input_list$fit_lnorm$estimate[1], 
                          sdlog = input_list$fit_lnorm$estimate[2]),
    
    gamma_qq = qgamma(ppoints(input_list$data),
                      shape = input_list$fit_gamma$estimate[1], 
                      rate = input_list$fit_gamma$estimate[2])
  )
  
  return(quantile_data)
  
}

#------------------------------------------------------------
# Functions to Plot priors
#------------------------------------------------------------

#Function to plot histograms using GGPLOT
plot_prior_density <- function(sample_density_df, xlabel, density_multipler, legend){
  
  legend_position <- c(0.55, 0.85)
  
  xmin <- 0; xmax <- max(sample_density_df$x); xstp <- xmax/4
  xmin <- 0; xmax <- signif(max(sample_density_df$x), digits = 1); xstp <- signif(xmax/4, digits = 1)
  xlbl <- xlabel
  ymin <- 0; ymax <- 1 ; ystp <- ymax/5
  ylbl <-"Probability Density"
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  # Histogram with density overlay
  histogram_density_plot <- ggplot(sample_density_df, aes(x=data_obs)) +
    # geom_histogram(aes(y = after_stat(density)),color = "black", fill = "grey", 
    #                bins=2^5, alpha =0.7)+
    
    #ggplot(sample_density_df, aes(x=data_obs)) +
    geom_histogram(stat = "density",color = "black", fill = "grey", 
                   n=2^5, adjust=1, alpha =0.7)+
    
    #                  size=1, na.rm = FALSE ) +
    geom_line(data = sample_density_df, aes(x = x, y = density_multipler[1]*weibull_density, 
                                            color = "xx1"), size = 2, show.legend = legend) +
    geom_line(data = sample_density_df, aes(x = x, y = density_multipler[2]*lognormal_density, 
                                            color = "xx2"), size = 2, show.legend = legend) +
    geom_line(data = sample_density_df, aes(x = x, y = density_multipler[3]*gamma_density, 
                                            color = "xx3"), size = 2, show.legend = legend) +
    
    labs(title = "",
         x = xlbl,
         y = ylbl
    ) +
    scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp))  +
    scale_color_viridis(discrete = T, name = "Distribution", 
                        labels = c("Weibull", "Log-Normal", "Gamma")
    )+
    theme_bw()+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 35, color = 'black', angle = 0),
          axis.text.y = element_text(size = 35, color = 'black', angle = 0),
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.ticks.length = unit(.25, "cm"),
          legend.box.just = "right",
          legend.spacing = unit(-0.5, "cm"),
          legend.position = legend_position,
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 35, face = "plain"),
          legend.background = element_rect(fill = NA, linewidth = 0.5,
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey",
                                          linewidth = 0.5, linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5))+
    guides(linetype = guide_legend(override.aes = list(size = 2)))
  
  return(histogram_density_plot)
  
}

# Function to do the QQ-plot with ggplot
plot_prior_qq <- function(input_list, sample_quantile_df, xlabel, ylabel, density_multipler, legend){
  legend_title <- element_blank()
  legend_position <- c(0.05, 0.85)
  
  xmin <- 0; xmax <- signif(max(sample_quantile_df$x), digits = 1); xstp <- signif(xmax/4, digits = 1)
  xlbl <- xlabel
  ymin <- xmin; ymax <- xmax ; ystp <- xstp
  ylbl <-ylabel
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  qqplot_combined <- ggplot(sample_quantile_df, aes(sample = x)) +
    
    stat_qq(distribution = qweibull, dparams = list(shape = input_list$fit_weibull$estimate[1],
                                                    scale = input_list$fit_weibull$estimate[2]),
            aes(shape = "xx1", fill= "xx1", color="xx1"), show.legend = legend, size=2) +

    stat_qq(distribution = qlnorm, dparams = list(meanlog  = input_list$fit_lnorm$estimate[1],
                                                  sdlog = input_list$fit_lnorm$estimate[2]),
            aes(shape = "xx2", fill = "xx2", color="xx2"), show.legend = legend, size=2) +

    stat_qq(distribution = qgamma, dparams = list(shape = input_list$fit_gamma$estimate[1],
                                                  rate = input_list$fit_gamma$estimate[2]),
            aes(shape = "xx3", fill= "xx3", color="xx3"), show.legend = legend, size=2) +
    
    
    geom_abline(slope = 1,intercept = 0,colour="black", na.rm = FALSE, show.legend = FALSE, size=1.3) +
    
    
    scale_color_viridis(discrete = T, name = "Distribution", 
                        labels = c("Weibull", "Log-Normal", "Gamma")
    )+
    scale_fill_viridis(discrete = T, name = "Distribution", 
                       labels = c("Weibull", "Log-Normal", "Gamma")
    )+
    scale_shape_manual(name = "Distribution"
                       , values = c(21,23,16), 
                       labels = c("Weibull", "Log-Normal", "Gamma")
    )+
    
    scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp))  +
    
    scale_y_continuous(name =ylbl , limits = c(ymin, ymax),
                       breaks = seq(ymin, ymax, ystp))+
    theme_bw()+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 35, color = 'black', angle = 0),
          axis.text.y = element_text(size = 35, color = 'black', angle = 0),
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.ticks.length = unit(.25, "cm"),
          legend.box.just = "right",
          legend.spacing = unit(-0.5, "cm"),
          legend.position = legend_position,
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 35, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5,
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey",
                                          size = 0.5, linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))+
    guides(linetype = guide_legend(override.aes = list(size = 2)))
  
  return(qqplot_combined)
}

#Function to plot CDFs using GGPLOT
plot_prior_cdf <- function(sample_quantile_df, xlabel, density_multipler, legend){
  legend_title <- element_blank()
  legend_position <- c(0.65, 0.25)
  
  xmin <- 0; xmax <- signif(max(sample_quantile_df$x), digits = 1); xstp <- signif(xmax/4, digits = 1)
  xlbl <- xlabel
  ymin <- 0; ymax <- 1 ; ystp <- ymax/4
  ylbl <-"Cumulative Density"
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  cdfplot_combined <- ggplot(sample_quantile_df, aes(x)) +
    stat_ecdf(aes(), color = "black", size = 1.3) +
    
    stat_ecdf(data = data.frame(x = density_multipler[1]*sample_quantile_df$weibull_qq, group = "xx1"), 
              aes(x = x, color = "xx1"), size = 2, linetype = "dashed", show.legend = legend) +
    
    stat_ecdf(data = data.frame(x = density_multipler[2]*sample_quantile_df$lognormal_qq, group = "xx2"), 
              aes(x = x, color = "xx2"), size = 2, linetype = "dashed", show.legend = legend) +
    
    stat_ecdf(data = data.frame(x = density_multipler[3]*sample_quantile_df$gamma_qq, group = "xx3"), 
              aes(x = x, color = "xx3"), size = 2, linetype = "dashed", show.legend = legend) +
    
    scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp))  +
    
    scale_y_continuous(name =ylbl , limits = c(ymin, ymax),
                       breaks = seq(ymin, ymax, ystp))+
    
    scale_color_viridis(discrete = T, name = "Distribution", 
                        labels = c("Weibull", "Log-Normal", "Gamma")
    )+
    theme_bw()+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 35, color = 'black', angle = 0),
          axis.text.y = element_text(size = 35, color = 'black', angle = 0),
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.ticks.length = unit(.25, "cm"),
          legend.box.just = "right",
          legend.spacing = unit(-0.5, "cm"),
          legend.position = legend_position,
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 35, face = "plain"),
          legend.background = element_rect(fill = NA, linewidth = 0.5,
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey",
                                          linewidth = 0.5, linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5))+
    guides(linetype = guide_legend(override.aes = list(size = 2)))
  
  return(cdfplot_combined)
}

#-----------------------------------------------------------------------
# Sample Prior with different distributions with positive support
#-----------------------------------------------------------------------
weibull_fit <- create.prior.data(use.nomad.prior = T, distrib.fit = "weibull", #Weibull
                                 plot.diag = T, truncate_chl = c(1,30),
                                 truncate_adg = c(0.1,5),
                                 truncate_bbp = c(0.002,0.01),
                                 sample_count = 150)

lognormal_fit <- create.prior.data(use.nomad.prior = T, distrib.fit = "lnorm", #Log-Normal
                                   plot.diag = T, truncate_chl = c(1,30),
                                   truncate_adg = c(0.1,5),
                                   truncate_bbp = c(0.002,0.01),
                                   sample_count = 150)

gamma_fit <- create.prior.data(use.nomad.prior = T, distrib.fit = "gamma", #Gamma
                               plot.diag = T, truncate_chl = c(1,30),
                               truncate_adg = c(0.1,5),
                               truncate_bbp = c(0.002,0.01),
                               sample_count = 150)

#-----------------------------------------------------------------------
#Store the goodness of fit for the distributions fitted to prior samples
#-----------------------------------------------------------------------
fit_stat_chl = gofstat(list(weibull_fit$fit_chl, lognormal_fit$fit_chl, gamma_fit$fit_chl))
fit_stat_adg = gofstat(list(weibull_fit$fit_acdm440, lognormal_fit$fit_acdm440, gamma_fit$fit_acdm440))
fit_stat_bbp = gofstat(list(weibull_fit$fit_bbp555, lognormal_fit$fit_bbp555, gamma_fit$fit_bbp555))


#Automatic fitting data to family of distributions
fit <- fitDist(weibull_fit$obs_chl, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)

summary(fit)
plot(fit)
wp(fit, xlim.all = 5, ylim.all = 5)

bcpe_chl = qGG(ppoints(chldata$data), mu = fit$mu, 
      sigma = fit$sigma, nu = fit$nu
      #, nu = fit$nu.fv 
    #  , tau = fit$tau.fv
    
    )

qqplot(bcpe_chl, chldata$data, xlim=c(0,10), ylim=c(0,10))
abline(0,1)

dbcpe_chl = dGG(chldata$data, mu = fit$mu, 
                sigma = fit$sigma, nu = fit$nu
                #, nu = fit$nu.fv 
                #  , tau = fit$tau.fv
                
)

hist(chldata$data, probability  = T)
lines(dbcpe_chl)

qqplot(bcpe_dg, acdm440data$data, xlim=c(0,2), ylim=c(0,2))
abline(0,1)

#Test for non-unimodality and mixture of distributions
dip_adg  = dip.test(acdm440data$data, simulate.p.value = F, B = 2000)
dip_chl  = dip.test(chldata$data, simulate.p.value = F, B = 2000)
dip_bbp  = dip.test(bbp555data$data, simulate.p.value = F, B = 2000)

#Fit a mixture distribution
fit3 <- flexmix(adg~1, data=subset(data.frame("adg"=acdm440data$data)), k=2, 
                model = list(FLXMRglm(family = "Gamma"), FLXMRglm(family = "Gamma"))
)

#Obtain the x for density
xval = seq(from=0, to=max(acdm440data$data), length.out=length(acdm440data$data))

#Extract the fitted parameters of the obtained mixture distribution
param3 <- parameters(fit3)[[1]]
interc <- param3[1,]
shape <- param3[2,]

lam <- table(clusters(fit3))
lambda <- c(lam[1]/sum(lam), lam[2]/sum(lam))

yval <- lambda[[1]]*dgamma(xval, shape=shape[[1]], rate=interc[[1]]*shape[[1]]) + 
  lambda[[2]]*dgamma(xval, shape=shape[[2]], rate=interc[[2]]*shape[[2]])

hist(acdm440data$data, probability = T, col="gray",  breaks=30, xlim=c(0,2))
lines(xval,yval, col="darkred", lwd=2)

qmix_gamma = lambda[[1]]*qgamma(ppoints(acdm440data$data), shape=shape[[1]], rate=interc[[1]]*shape[[1]]) + 
  lambda[[2]]*qgamma(ppoints(acdm440data$data), shape=shape[[2]], rate=interc[[2]]*shape[[2]])

qqplot(x = qmix_gamma, acdm440data$data, xlim=c(0,1.5), ylim=c(0,1.5))
abline(0,1)
#--------------------------------------------------------------------------------
#Create data-frames with prior data and density functions of fitted distributions
#--------------------------------------------------------------------------------
chldata = list("data" = weibull_fit$obs_chl, "fit_weibull" = weibull_fit$fit_chl,
               "fit_lnorm" = lognormal_fit$fit_chl, "fit_gamma" = gamma_fit$fit_chl)

acdm440data = list("data" = weibull_fit$obs_acdm440, "fit_weibull" = weibull_fit$fit_acdm440,
               "fit_lnorm" = lognormal_fit$fit_acdm440, "fit_gamma" = gamma_fit$fit_acdm440)

bbp555data = list("data" = weibull_fit$obs_bbp550, "fit_weibull" = weibull_fit$fit_bbp555,
               "fit_lnorm" = lognormal_fit$fit_bbp555, "fit_gamma" = gamma_fit$fit_bbp555)


#--------------------------------------------------------------------------------
# Plot the fitted probability distribution functions (Hist & PDF)
#--------------------------------------------------------------------------------
#Calculate the density
sample_density_chl = estimate_prior_density(input_list = chldata)
sample_density_dg = estimate_prior_density(input_list = acdm440data)
sample_density_bbp = estimate_prior_density(input_list = bbp555data)

#Create AXIS labels
chl_label = expression(paste("[", italic("chl"), "] [","mg"^1,"m"^-3,"]"))
adg_label = expression(paste(italic("a")["dg"](443),"[m"^"-1","]"))
bbp_label = expression(paste(italic("b")["bp"](555),"[m"^"-1","]"))

#Plot the histograms 
chl_prior_hist = plot_prior_density(sample_density_df = sample_density_chl, xlabel = chl_label,
                                    legend = TRUE, 
                                    density_multipler = c(1,0.82,1))

adg_prior_hist = plot_prior_density(sample_density_df = sample_density_dg, xlabel = adg_label,
                                    legend = FALSE,
                                    density_multipler = c(1,1,1))

bbp_prior_hist = plot_prior_density(sample_density_df = sample_density_bbp, xlabel = bbp_label,
                                    legend = FALSE,
                                    density_multipler = c(1,1,1))


#--------------------------------------------------------------------------------
# Plot the quantiles of fitted probability distributions  (QQ-Plot)
#--------------------------------------------------------------------------------

#Calculate the quantiles
sample_quantile_chl = estimate_prior_quantile(input_list = chldata)
sample_quantile_dg = estimate_prior_quantile(input_list = acdm440data)
sample_quantile_bbp = estimate_prior_quantile(input_list = bbp555data)

#Create AXIS labels
chl_xlabel = expression(paste("Theoretical quantiles of [", italic("chl"), "] [","mg"^1,"m"^-3,"]"))
chl_ylabel = expression(paste("Observed quantiles of [", italic("chl"), "] [","mg"^1,"m"^-3,"]"))

adg_xlabel = expression(paste("Theoretical quantiles of ",italic("a")["dg"](443),"[m"^"-1","]"))
adg_ylabel = expression(paste("Observed quantiles of ",italic("a")["dg"](443),"[m"^"-1","]"))

bbp_xlabel = expression(paste("Theoretical quantiles of ",italic("b")["bp"](555),"[m"^"-1","]"))
bbp_ylabel = expression(paste("Observed quantiles of ",italic("b")["bp"](555),"[m"^"-1","]"))

#Plot the QQ plots
chl_prior_qq = plot_prior_qq(sample_quantile_df = sample_quantile_chl, input_list = chldata,
                             xlabel = chl_xlabel, ylabel = chl_ylabel,
                             legend = FALSE,
                                    density_multipler = c(1,1,1))

adg_prior_qq = plot_prior_qq(sample_quantile_df = sample_quantile_dg, input_list = acdm440data,
                             xlabel = adg_xlabel, ylabel = adg_ylabel,
                             legend = TRUE,
                                    density_multipler = c(1,1,1))

bbp_prior_qq = plot_prior_qq(sample_quantile_df = sample_quantile_bbp, input_list = bbp555data,
                             xlabel = bbp_xlabel, ylabel = bbp_ylabel,
                             legend = FALSE,
                                    density_multipler = c(1,1,1))

#--------------------------------------------------------------------------------
# Plot the Cumulative Distribution functions of fitted distributions  (CDF)
#--------------------------------------------------------------------------------
#Plot the CDF plots
chl_prior_cdf = plot_prior_cdf(sample_quantile_df = sample_quantile_chl,
                             xlabel = chl_label,
                             legend = FALSE,
                             density_multipler = c(1,1,1))

adg_prior_cdf = plot_prior_cdf(sample_quantile_df = sample_quantile_dg,
                             xlabel = adg_label,
                             legend = FALSE,
                             density_multipler = c(1,1,1))

bbp_prior_cdf = plot_prior_cdf(sample_quantile_df = sample_quantile_bbp,
                             xlabel = bbp_label,
                             legend = TRUE,
                             density_multipler = c(1,1,1))

#--------------------------------------------------------------------------------
#Plot PDF, QQ & CDF plots together
#--------------------------------------------------------------------------------
#combine plots
combined_plots <- arrangeGrob(chl_prior_hist, adg_prior_hist, bbp_prior_hist, #density plots
                               chl_prior_qq, adg_prior_qq, bbp_prior_qq, #QQ plots
                              chl_prior_cdf, adg_prior_cdf, bbp_prior_cdf, #CDF plots
                              ncol = 3)
#save to disc
ggsave('./outputs/prior_new.png',combined_plots,scale = 1.7, width = 20, height = 16, 
       units = "in",dpi = 300)


