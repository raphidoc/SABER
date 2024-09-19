#-------------------------------------------------------------------------
#Create a single Bayes Function for the MCMC sampling
#-------------------------------------------------------------------------
library(mcmc)
library(BayesianTools)
library(plot3D)
library(coda)
library(bayesplot)

source("./R/SABER_forward_fast.R")
source("./R/lee_forward_fast.R")
source("./R/solve.objective.inverse_fast.R")

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

inverse_runBayes <- function(obsdata, rrs_type, max_par, min_par, 
                     param_lab = c("chl","adg443","bbp555", 
                                   "H",  rep(paste0("fa", (1:rb_count))), 
                                   "pop_sd"),
                     constrain_config = c("const_bbp" = F, "cost_bgc"=F, "cost_IOP"=F),
                     qaa_slope, manual_slope, 
                     bbp_const_val, bgc_const_val, iop_const_path,
                     manual_slope_vals= c("s_g"=0.014, "s_d"=0.003, "gamma"=0.5),
                     iter_count, sampler_mcmc,
                     wavelngth_sim, sa_model, hybrid_mode, plot_rrs){
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### EMBRACE THE RANDOMNESS #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  
  
  cat(paste0("\033[0;34m**************************************************************************\033[0m","\n"))
  cat(paste0("\033[0;34m MODEL CONSTRAINTS ::::\033[0m","\n"))
  
  cat(paste0("\033[0;32m","Backscatter constrain: ",constrain_config[1],"\033[0m","\n"))
  if (constrain_config[1] == T) {
    cat(paste0("\033[0;32m","bbp(555) constrained value: ",bbp_const_val,"\033[0m","\n"))
  }
  
  cat(paste0("\033[0;32m","Biogeochemical variable constrain: ",constrain_config[2],"\033[0m","\n"))
  if (constrain_config[2] == T) {
    cat(paste0("\033[0;32m","BGC constrained value: ",bgc_const_val,"\033[0m","\n"))
  }
  cat(paste0("\033[0;32m","IOP constrain: ",constrain_config[3],"\033[0m","\n"))
  if (constrain_config[3] == T) {
    cat(paste0("\033[0;32m","Spectral IOP paths: ",paste0(iop_const_path, collapse = ","),"\033[0m","\n"))
  }
  
  cat(paste0("\033[0;34m==========================================================================\033[0m","\n"))
  cat(paste0("\033[0;34m SPECTRAL SLOPE SETTINGS ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Slope for QAA: ",qaa_slope,"\033[0m","\n"))
  cat(paste0("\033[0;32m","Slope from user: ",manual_slope,"\033[0m","\n"))
  if (manual_slope == TRUE) {
    
    cat(paste0("\033[0;32m","Slope values: ",manual_slope_vals,"\033[0m","\n"))
  }
  cat(paste0("\033[0;34m==========================================================================\033[0m","\n"))
  cat(paste0("\033[0;34m ANCILIARY ::::\033[0m","\n"))
  cat(paste0("\033[0;32m","Rrs type: ",rrs_type,"\033[0m","\n"))
  if (sa_model == "am03") {
    cat(paste0("\033[0;32m","Paramterization for forward model: Albert & Mobley (2003) \033[0m","\n"))
    
  } else {
    cat(paste0("\033[0;32m","Paramterization for forward model: Lee et al., (1998) \033[0m","\n"))
    
  }
  cat(paste0("\033[0;32m","MCMC iteration count:",iter_count,"\033[0m","\n"))
  cat(paste0("\033[0;32m","MCMC Sampler:",sampler_mcmc,"\033[0m","\n"))
  cat(paste0("\033[0;32m","HYBRID mode: ",hybrid_mode,"\033[0m","\n"))
  cat(paste0("\033[0;34m**************************************************************************\033[0m","\n"))
  Sys.sleep(1)
  
  initial = max_par
  #initial_rb_length = length(initial[5:(length(initial)-1)])
  initial_rb_length = 3
  
  #bbp constrained Inversion
  if (constrain_config[1] == T) {
    
    LL_unconstr_bayes = function(pars) {
      
      # Values predicted by the forward model
      if (sa_model == "am03") {
        Gpred = Saber_forward_fast(
          use_true_IOPs = F, 
          a_non_water_path = iop_const_path[1],
          bb_non_water_path = iop_const_path[2],
          
          chl = pars[1], 
          a_dg = pars[2],
          bbp.550 = bbp_const_val,
          
          z = pars[4],
          rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
          
          
          Rrs_input_for_slope = obsdata,
          
          slope.parametric = qaa_slope,
          
          
          use_manual_slope =manual_slope,
          manual_slope =  manual_slope_vals,
          
          verbose = F, wavelength = wavelngth_sim
        )
        
      } else {
        Gpred = lee_forward_fast(
          use_true_IOPs = F, 
          
          chl = bgc_const_val[1], 
          a_dg = bgc_const_val[2],
          bbp.550 = bgc_const_val[3],
          
          z = pars[1],
          rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
          
          
          Rrs_input_for_slope = obsdata,
          
          slope.parametric = qaa_slope,
          
          
          use_manual_slope =manual_slope,
          manual_slope =  manual_slope_vals,
          
          verbose = F, wavelength = wavelngth_sim
        )
      }
      
      # Negative log-likelihood
      smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                        log = TRUE), na.rm = T)
      return(smull)
    }
    
  } else {
    
    #IOP constrained Inversion
    if (constrain_config[3] == T) {
      
      LL_unconstr_bayes = function(pars) {
        
        # Values predicted by the forward model
        if (sa_model == "am03") {
          Gpred = Saber_forward_fast(
            use_true_IOPs = T, 
            a_non_water_path = iop_const_path[1],
            bb_non_water_path = iop_const_path[2],
            
            chl = bgc_const_val[1], 
            a_dg = bgc_const_val[2],
            bbp.550 = bgc_const_val[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            Rrs_input_for_slope = obsdata,
            
            slope.parametric = qaa_slope,
            
            
            use_manual_slope =manual_slope,
            manual_slope =  manual_slope_vals,
            
            verbose = F, wavelength = wavelngth_sim
          )
          
        } else {
          Gpred = lee_forward_fast(
            use_true_IOPs = F, 
            
            chl = bgc_const_val[1], 
            a_dg = bgc_const_val[2],
            bbp.550 = bgc_const_val[3],
            
            z = pars[1],
            rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
            
            
            Rrs_input_for_slope = obsdata,
            
            slope.parametric = qaa_slope,
            
            
            use_manual_slope =manual_slope,
            manual_slope =  manual_slope_vals,
            
            verbose = F, wavelength = wavelngth_sim
          )
        }
        
        # Negative log-likelihood
        smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                          log = TRUE), na.rm = T)
        return(smull)
      }
      
    } else {
      
      #Constrained Inversion
      if (constrain_config[2] == T) {
        
        LL_unconstr_bayes = function(pars) {
          
          # Values predicted by the forward model
          if (sa_model == "am03") {
            Gpred = Saber_forward_fast(
              use_true_IOPs = F, 
              
              chl = bgc_const_val[1], 
              a_dg = bgc_const_val[2],
              bbp.550 = bgc_const_val[3],
              
              z = pars[1],
              rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
              
              
              Rrs_input_for_slope = obsdata,
              
              slope.parametric = qaa_slope,
              
              
              use_manual_slope =manual_slope,
              manual_slope =  manual_slope_vals,
              
              verbose = F, wavelength = wavelngth_sim
            )
            
          } else {
            Gpred = lee_forward_fast(
              use_true_IOPs = F, 
              
              chl = bgc_const_val[1], 
              a_dg = bgc_const_val[2],
              bbp.550 = bgc_const_val[3],
              
              z = pars[1],
              rb.fraction = as.numeric(pars[2:(1+initial_rb_length)]),
              
              
              Rrs_input_for_slope = obsdata,
              
              slope.parametric = qaa_slope,
              
              
              use_manual_slope =manual_slope,
              manual_slope =  manual_slope_vals,
              
              verbose = F, wavelength = wavelngth_sim
            )
          }
          
          # Negative log-likelihood
          smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                            log = TRUE), na.rm = T)
          return(smull)
        }
        
      } else {
        
        #Unconstrained Inversion
        if (all(constrain_config) == F) {
          LL_unconstr_bayes = function(pars) {
            
            # Values predicted by the forward model
            if (sa_model == "am03") {
              Gpred = Saber_forward_fast(
                use_true_IOPs = F, 
                
                chl = pars[1], 
                a_dg = pars[2],
                bbp.550 = pars[3],
                
                z = pars[4],
                rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
                
                
                Rrs_input_for_slope = obsdata,
                
                slope.parametric = qaa_slope,
                
                
                use_manual_slope =manual_slope,
                manual_slope =  manual_slope_vals,
                
                verbose = F, wavelength = wavelngth_sim
              )
              
            } else {
              Gpred = lee_forward_fast(
                use_true_IOPs = F, 
                
                chl = pars[1], 
                a_dg = pars[2],
                bbp.550 = pars[3],
                
                z = pars[4],
                rb.fraction = as.numeric(pars[5:(4+initial_rb_length)]),
                
                
                Rrs_input_for_slope = obsdata,
                
                slope.parametric = qaa_slope,
                
                
                use_manual_slope =manual_slope,
                manual_slope =  manual_slope_vals,
                
                verbose = F, wavelength = wavelngth_sim
              )
            }
            
            # Negative log-likelihood
            smull = sum(dnorm(x = 10000*obsdata, mean = 10000*Gpred[[1]]$Rrs, sd = pars[length(pars)], 
                              log = TRUE), na.rm = T)
            return(smull)
          }
        }
        
      }
      
    }
    
  }
  
  
  
  #Instantiate inversion objective function and the optimization scheme 
  #browser()
  inv_bound = create_init_bound(rrs_inv_type = rrs_type, manual_par0 = T, 
                                constrain.bbp = constrain_config[1], 
                                constrain.shallow.bgc = constrain_config[2], 
                                constrain.shallow.iop = constrain_config[3], 
                                pop.sd =  "unknown",
                                
                                init_par = rowMeans(matrix(c(max_par,min_par), 
                                                           ncol=2)),
                                upper_par = max_par,
                                lower_par = min_par
  )
  
  par0_bayes = inv_bound$par0 ; upper.bound_bayes = inv_bound$upper.bound  
  lower.bound_bayes= inv_bound$lower.bound
  
  cat(paste0("\033[0;32m","MINIMA values are: [chl]=",lower.bound_bayes[1], 
             ",adg(440)=",lower.bound_bayes[2], 
             ",bbp(555)=", lower.bound_bayes[3], ", zB=", 
             lower.bound_bayes[4],",RB={", 
             toString(as.numeric(lower.bound_bayes[5:(length(lower.bound_bayes)-1)])),
             "},population_sigma=", lower.bound_bayes[length(lower.bound_bayes)],"\033[0m","\n"))
  
  cat(paste0("\033[0;32m","MAXIMA values are: [chl]=",upper.bound_bayes[1], 
             ",adg(440)=",upper.bound_bayes[2], 
             ",bbp(555)=", upper.bound_bayes[3], ", zB=", 
             upper.bound_bayes[4],",RB={", 
             toString(as.numeric(upper.bound_bayes[5:(length(upper.bound_bayes)-1)])),
             "},population_sigma=", upper.bound_bayes[length(upper.bound_bayes)],"\033[0m","\n"))
  
  
  if (hybrid_mode == TRUE) {
    
    cat(paste0("\033[0;34m EXECUTING HYBRID MODE ::::\033[0m","\n"))
    cat(paste0("\033[0;41m CAUTION: DONT RUN UNLESS UNCONSTRAINED SHALLOW WATER INVERSION ::::\033[0m","\n"))
    
    inverse_output <- quiet(suppressWarnings(solve.objective.inverse.shallow.final.fast(
      
      wave = wavelngth_sim, 
      
      constrain.shallow.iop = constrain_config[3],  
      
      unconstrained = T,
      
      constrained.bgc = constrain_config[2], 
      constrain.bgc.value = NA,
      
      
      initial = as.numeric(par0_bayes), 
      obsdata = as.numeric(obsdata),
      
      auto_spectral_slope = qaa_slope,
      manual_spectral_slope = manual_slope, 
      
      manual_spectral_slope_vals = manual_slope_vals,
      
      sa.model = sa_model, 
      
      obj.fn = "log-LL",
      method.opt = methods.opt[4],
      
      lower.b = lower.bound_bayes,
      upper.b = upper.bound_bayes, 
    )))
    
    cat(paste0("\033[0;32m Expected parameters obtained from gradient-based optimization::::", 
               paste0(signif(inverse_output[[1]]$estimates, digits = 2), collapse = ","),
               "\033[0m","\n"))
    
    bayessetup <- createBayesianSetup(prior = NULL,
                                      likelihood = LL_unconstr_bayes,
                                      lower = lower.bound_bayes,#[c(1:3,length(lower.bound))] 
                                      best = inverse_output[[1]]$estimates, 
                                      upper = upper.bound_bayes#[c(1:3,length(lower.bound))]
                                      ,
                                      names = param_lab,
                                      parallel = F)
    
  } else {
    
    bayessetup <- createBayesianSetup(prior = NULL,
                                      likelihood = LL_unconstr_bayes,
                                      lower = lower.bound_bayes#[c(1:3,length(lower.bound))]
                                      , 
                                      upper = upper.bound_bayes#[c(1:3,length(lower.bound))]
                                      ,
                                      names = param_lab,
                                      parallel = F)
    
  }
  
  cat(paste0("\033[0;34m @@@@@@@@@@ TEST Bayesian Set-UP @@@@@@@@@ \033[0m","\n"))
  print(checkBayesianSetup(bayessetup)) #Test if the setup is initiated for theta pars
  
  max_iter = iter_count ; burn_in = max_iter/4; chain_burn = burn_in/3
  
  settings = list(iterations = max_iter, message = TRUE, nrChains = 1 #Set MCMC config
                  , burnin=burn_in
  )
  
  #Run MCMC
  cat(paste0("\033[0;34m @@@@@@@@@@ MCMC START @@@@@@@@@ \033[0m","\n"))
  out <- runMCMC(bayesianSetup = bayessetup, settings = settings, sampler = sampler_mcmc )
  
  sd_c1 = apply(out[["chain"]][[1]][-1,1:length(upper.bound_bayes)], 2, sd)
  sd_c2 = apply(out[["chain"]][[2]][-1,1:length(upper.bound_bayes)], 2, sd)
  sd_c3 = apply(out[["chain"]][[3]][-1,1:length(upper.bound_bayes)], 2, sd)
  
  sd_bayes = rowMeans(matrix(c(sd_c1,sd_c2,sd_c3), 
                             ncol=3))
  MAP_bayes = MAP(out, start = burn_in
  )
  bayes_output = c(MAP_bayes$parametersMAP, sd_bayes)
  names(bayes_output) = c(param_lab, paste0("sd_", param_lab))
  cat(paste0("\033[0;34m @@@@@@@@@@ MCMC FINISHED @@@@@@@@@ \033[0m","\n"))
  summary(out, start = burn_in)
  
  
  if (plot_rrs == TRUE) {
    
    
    pred_rrs = Saber_forward_fast(
      use_true_IOPs = F, 
      
      chl = bayes_output[1], 
      a_dg = bayes_output[2],
      bbp.550 = bayes_output[3],
      
      z = bayes_output[4],
      rb.fraction = as.numeric(bayes_output[5:(4+initial_rb_length)]),
      
      
      Rrs_input_for_slope = obsdata,
      
      slope.parametric = qaa_slope,
      
      
      use_manual_slope =manual_slope,
      manual_slope =  manual_slope_vals,
      
      verbose = F, wavelength = wavelngth_sim
    )
    
    
    plotframe.rrs <- data.frame("wave"=wavelngth_sim, 
                                "rrs_est"=pred_rrs[[1]]$Rrs,
                                "rrs_obs"=obsdata)
    #Create labels
    map_vals = signif(MAP_bayes$parametersMAP,digits = 2)
    mcmc_label = paste0("theta[MAP]== {",map_vals[1],"*',",
                        map_vals[2],",",
                        map_vals[3],",",
                        map_vals[4],",",
                        map_vals[5],",",
                        map_vals[6],",",
                        map_vals[7],
                        "'}")
  
    #Create AXIS
    xmin = min(wavelngth_sim); xmax= max(wavelngth_sim); xstp=100
    ymin= 0; ymax=max(plotframe.rrs$rrs_obs)+0.20*max(plotframe.rrs$rrs_obs)
    ystp= signif(ymax/5, digits = 1)
    asp_rat <- (xmax-xmin)/(ymax-ymin)
    
    #Create Plot
    g1 <- ggplot()  + 
      geom_line(data = plotframe.rrs,aes(y=rrs_obs,color="xx1",x = wave),
                size=1.3,show.legend = TRUE) +
      geom_line(data = plotframe.rrs,aes(y=rrs_est,x = wave,color="xx2"), 
                linetype = "dashed",
                size=1.3,show.legend = TRUE)+
      
      scale_colour_manual(values = c("red3", "blue3"),
                           labels = c(expression(paste(italic("R")["rs,obs"])),
                                     
                                     expression(paste(italic("R")["rs,MAP"]))) 
                          ) +
      
      scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), 
                         limits = c(xmin, xmax), 
                         breaks = seq(xmin, xmax, xstp))  +
      scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,
                                                ",", 0^"-",")[", sr^-1,"]")) , 
                         limits = c(ymin, ymax),
                         breaks = seq(ymin, ymax, ystp))+ 
      coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                  ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
     
      annotate("text",x=550, y= 0.9*ymax, label = mcmc_label, parse=T, size=3, color="black")+
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
            axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
            axis.title.x = element_text(size = 25),
            axis.title.y = element_text(size = 25),
            axis.ticks.length = unit(.25, "cm"),
            legend.position=c(0.62, 0.905),
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
            #legend.spacing.y = unit(2.0, 'cm'),
            plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
            legend.text.align = 0,
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
    
    dev.new()
    print(g1)
    
  }
  
  MLE <- data.table("param" = param_lab,
                    "estimates" = MAP_bayes$parametersMAP,
                    "95% C.I.(+/-)" = sd_bayes)
  
  print("The inversion retrieved parameters are:")
  prmatrix(MLE)
  
  cat(paste0("\033[0;33m","###################################################################","\033[0m","\n"))
  cat(paste0("\033[0;39m","########### INVERSION IS COMPLETED #######","\033[0m","\n"))
  cat(paste0("\033[0;32m","###################################################################","\033[0m","\n"))
  
  return(bayes_output)
  
}
