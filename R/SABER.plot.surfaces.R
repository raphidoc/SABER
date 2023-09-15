

#==============================================================================
# Plot the objective function Surfaces governed by forward model inputs
#==============================================================================

##Create Objective function for shallow waters

NLL_unconstr = function(pars, data) {
  
  # Values predicted by the forward model for single RUN
  if (sa.model == "am03") {
    Gpred = Saber_forward_final(
      use_true_IOPs = F, 
      a_non_water_path = IOP_files[idx_a],
      bb_non_water_path = IOP_files[idx_bb],
      
      chl = pars[1], 
      acdom440 =NULL, 
      anap440 =NULL , 
      a_dg = pars[2],
      bbp.550 = pars[3],
      
      z = pars[4],
      use_spectral_rb = F, 
      rb.fraction = pars[5:(length(pars)-1)],
      
      
      realdata = data,
      
      slope.parametric = auto_spectral_slope,
      dg_composite = T,
      
      use_manual_slope =manual_spectral_slope,
      manual_slope = manual_spectral_slope_vals ,
      
      use_spectral_shape_chl = F,
      use_spectral_shape_dg = T,
      
      sicf = F, q_phi = 0.02, 
      fDOM = F,
      verbose = F, plot = F
    )
    
  } else {
    Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                        anap440 =pars[3], bbp.550 = bbp.550,
                        z = pars[4], rb.fraction = pars[5:(length(pars)-1)],
                        verbose = F, realdata = data, plot = F)
  }
  
  rrs_est = Gpred[[1]]$Rrs
  
  if (obj.fn == "log-LL") {
    #ON for unknown sigma
    smull = -sum(dnorm(x = 10000*data, mean = 10000*rrs_est, sd = pars[length(pars)], 
                       log = TRUE))
  }
  
  if (obj.fn == "obj_L98") {
    
    #The Spectral error index from Lee 1999
    
    # Define the spectral regions
    region1 <- which(wavelength >= 400 & wavelength <= 675)
    region2 <- which(wavelength >= 750 & wavelength <= 830)
    
    # Calculate the numerator of the error index
    numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                        sum((data[region2] - rrs_est[region2])^2))
    
    # Calculate the denominator of the error index
    denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
    
    # Calculate the error index
    smull =  numerator / denominator
  }
  
  if (obj.fn == "SSR") {
    # Sum-squared residual of error (SSE)
    smull= sum((data - rrs_est)^2)
  }
  
  return(smull)
}


#Create function for calculation of gradient

NLL_unconstr_grad = function(pars, pars_sd, data) {
  
  # Values predicted by the forward model for single RUN
  if (sa.model == "am03") {
    Gpred = Saber_forward_final(
      use_true_IOPs = F, 
      a_non_water_path = IOP_files[idx_a],
      bb_non_water_path = IOP_files[idx_bb],
      
      chl = pars[1], 
      acdom440 =NULL, 
      anap440 =NULL , 
      a_dg = pars[2],
      bbp.550 = pars[3],
      
      z = pars[4],
      use_spectral_rb = F, 
      rb.fraction = pars[5:(length(pars)-1)],
      
      
      realdata = data,
      
      slope.parametric = auto_spectral_slope,
      dg_composite = T,
      
      use_manual_slope =manual_spectral_slope,
      manual_slope = manual_spectral_slope_vals ,
      
      use_spectral_shape_chl = F,
      use_spectral_shape_dg = T,
      
      sicf = F, q_phi = 0.02, 
      fDOM = F,
      verbose = F, plot = F
    )
    
  } else {
    Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                        anap440 =pars[3], bbp.550 = bbp.550,
                        z = pars[4], rb.fraction = pars[5:(length(pars)-1)],
                        verbose = F, realdata = data, plot = F)
  }
  
  rrs_est = Gpred[[1]]$Rrs
  
  if (obj.fn == "log-LL") {
    #ON for unknown sigma
    smull = -sum(dnorm(x = 10000*data, mean = 10000*rrs_est, sd = pars_sd, 
                       log = TRUE))
  }
  
  if (obj.fn == "obj_L98") {
    
    #The Spectral error index from Lee 1999
    
    # Define the spectral regions
    region1 <- which(wavelength >= 400 & wavelength <= 675)
    region2 <- which(wavelength >= 750 & wavelength <= 830)
    
    # Calculate the numerator of the error index
    numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                        sum((data[region2] - rrs_est[region2])^2))
    
    # Calculate the denominator of the error index
    denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
    
    # Calculate the error index
    smull =  numerator / denominator
  }
  
  if (obj.fn == "SSR") {
    # Sum-squared residual of error (SSE)
    smull= sum((data - rrs_est)^2)
  }
  
  return(smull)
}

#------------------------------------------------------------
##Create Objective function for shallow waters
#------------------------------------------------------------
NLL_shallow = function(pars, data) {
  
  # Values predicted by the forward model for single RUN
  if (sa.model == "am03") {
    Gpred = Saber_forward_final(
      use_true_IOPs = F, 
      a_non_water_path = IOP_files[idx_a],
      bb_non_water_path = IOP_files[idx_bb],
      
      chl = pars[1], 
      acdom440 =NULL, 
      anap440 =NULL , 
      a_dg = pars[2],
      bbp.550 = pars[3],
      
      z = pars[4],
      use_spectral_rb = F, 
      rb.fraction = fA.set,
      
      
      realdata = data, realdata_wave = wavelength,
      
      slope.parametric = auto_spectral_slope,
      dg_composite = T,
      
      use_manual_slope =manual_spectral_slope,
      manual_slope =  manual_spectral_slope_vals,
      
      use_spectral_shape_chl = F,
      use_spectral_shape_dg = T,
      
      sicf = F, q_phi = 0.02, 
      fDOM = F,
      verbose = F, plot = F
    )
    
  } else {
    Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                        anap440 =pars[3], bbp.550 = pars[3],
                        verbose = F, realdata = data, plot = F)
  }
  
  rrs_est = Gpred[[1]]$Rrs
  
  if (obj.fn == "log-LL") {
    #ON for unknown sigma
    smull = -sum(dnorm(x = 10000*data, mean = 10000*rrs_est, sd = pars[length(pars)], 
                       log = TRUE))
  }
  
  if (obj.fn == "obj_L98") {
    
    #The Spectral error index from Lee 1999
    
    # Define the spectral regions
    region1 <- which(wavelength >= 400 & wavelength <= 675)
    region2 <- which(wavelength >= 750 & wavelength <= 830)
    
    # Calculate the numerator of the error index
    numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                        sum((data[region2] - rrs_est[region2])^2))
    
    # Calculate the denominator of the error index
    denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
    
    # Calculate the error index
    smull =  numerator / denominator
  }
  
  if (obj.fn == "SSR") {
    # Sum-squared residual of error (SSE)
    smull= sum((data - rrs_est)^2)
  }
  return(smull)
}


NLL_deep = function(pars, data) {
  
  # Values predicted by the forward model for single RUN
  if (sa.model == "am03") {
    Gpred = Saber_forward_final(
      use_true_IOPs = F, 
      a_non_water_path = IOP_files[idx_a],
      bb_non_water_path = IOP_files[idx_bb],
      
      chl = pars[1], 
      acdom440 =NULL, 
      anap440 =NULL , 
      a_dg = pars[2],
      bbp.550 = pars[3],
      
      z = zB,
      use_spectral_rb = F, 
      rb.fraction = fA.set,
      
      
      realdata = data,
      realdata_wave = wavelength,
      
      slope.parametric = auto_spectral_slope,
      dg_composite = T,
      
      use_manual_slope =manual_spectral_slope,
      manual_slope =  manual_spectral_slope_vals,
      
      use_spectral_shape_chl = F,
      use_spectral_shape_dg = T,
      
      sicf = F, q_phi = 0.02, 
      fDOM = F,
      verbose = F, plot = F
    )
    
  } else {
    Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                        anap440 =pars[3], bbp.550 = pars[3],
                        verbose = F, realdata = data, plot = F)
  }
  
  rrs_est = Gpred[[1]]$Rrs
  
  if (obj.fn == "log-LL") {
    #ON for unknown sigma
    smull = -sum(dnorm(x = 10000*data, mean = 10000*rrs_est, sd = pars[length(pars)], 
                       log = TRUE))
  }
  
  if (obj.fn == "obj_L98") {
    
    #The Spectral error index from Lee 1999
    
    # Define the spectral regions
    region1 <- which(wavelength >= 400 & wavelength <= 675)
    region2 <- which(wavelength >= 750 & wavelength <= 830)
    
    # Calculate the numerator of the error index
    numerator <- sqrt(sum((data[region1] - rrs_est[region1])^2) + 
                        sum((data[region2] - rrs_est[region2])^2))
    
    # Calculate the denominator of the error index
    denominator <- sum(rrs_est[region1]) + sum(rrs_est[region2])
    
    # Calculate the error index
    smull =  numerator / denominator
  }
  
  if (obj.fn == "SSR") {
    # Sum-squared residual of error (SSE)
    smull= sum((data - rrs_est)^2)
  }
  return(smull)
}


#==============================================================================
# Create Bounds for parameters
#==============================================================================
below_type =c("deep", "shallow")
type_Rrs_below = below_type[2]

#TEST PURPOSE
auto_spectral_slope = F
manual_spectral_slope = F 

manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.01160, "gamma"=0.46)

obj = c("log-LL", "SSR", "obj_L98")
sa.model = "am03"
obj.fn = obj[1]

forward.op.am.param.conc.dg_comp_sicf_fdom <- Saber_forward_final(
                                use_true_IOPs = F,
                                
                                chl = Fit.input$chl, 
                                acdom440 = NULL, 
                                anap440 = NULL, 
                                a_dg =   Fit.input$acdom.440 + Fit.input$anap.440  ,
                                bbp.550 = Fit.input$bbp.550, 
                                
                                z=zB,
                                rb.fraction = fA.set,
                                use_spectral_rb = F, 
                                spectral_rb_path = Rb_files[idx_rb],
                                
                                realdata.exist = F,
                                realdata = surface_rrs_translate(Rrs = insitu.data),
                                realdata_wave = wavelength,
                                
                                dg_composite = T,
                                
                                slope.parametric = auto_spectral_slope,
                                use_spectral_shape_chl = F,
                                use_spectral_shape_dg = F,
                                
                                use_manual_slope =F,
                                manual_slope = c("s_g"=0.016, "s_d"=0.01160, "gamma"=0.5),
                                
                                sicf = F, q_phi = 0.05, 
                                
                                fDOM = F,
                                sunzen_Ed = -99, 
                                lat_Ed = 49.02487, lon_Ed = -68.37059,
                                date_time_Ed = "2019-08-18 20:59 GMT", 
                                Ed_fDOM_path = "./data/input-spectra/Ed_HL.csv",
                                use_fDOM_rad = F,
                                
                                verbose = T, plot = F)

#Extract AM03 modeled Rrs
rrs.forward.am.param.conc.dg_comp_sicf_fdom <- forward.op.am.param.conc.dg_comp_sicf_fdom[[1]]$Rrs
obsdata = rrs.forward.am.param.conc.dg_comp_sicf_fdom


#==============================================================================
#Plot the likelihood function over parameters
#==============================================================================

# 1.1 plot the likelihood profile of a single parameter among the fit vector
slopevalues = function(x){return(NLL_deep(pars = c(Fit.input$chl, 
                                                   x,
                                                   #Fit.input$acdom.440+Fit.input$anap.440,
                                                   Fit.input$bbp.550,
                                                   0.00001
                                                   #,x
),
data = obsdata))}
param_sequence = seq(0.5, 50, by = 0.5)

slopelikelihoods = lapply(param_sequence, slopevalues )
plot (param_sequence, slopelikelihoods , type="l", xlab = "", 
      ylab = "Log likelihood", col="red", lwd=2)
abline(v=param_sequence[which.min(slopelikelihoods)], col="navyblue", lwd=2)


# 1.2 plot the likelihood profile of two parameter among the fit vector
library(plot3D)
library(lattice)

log_like_graph_applybase <- function(da){
  chl = da[1]
  adg443 = da[2]
  bbp550 = da[3]
  H = da[4]
  sd = 0.00001
  loglik = NLL_shallow(pars = c(chl,adg443,bbp555
                             ,H
                             ,sd),
                    data = obsdata)
  return(loglik)
}

log_like_graph_applybase <- Vectorize(log_like_graph_applybase)

log_like_graph_outer <- function(adg443, H){
  
  loglik = NLL_shallow(pars = c(Fit.input$chl, adg443, Fit.input$bbp.550, H,0.00001),
                    data = obsdata)
  return(loglik)
}
log_like_graph_outer <- Vectorize(log_like_graph_outer)

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#Parameter Vector
da = data.frame("chl"=seq(0.5,50,length.out=100),
                 "adg443" = seq(0.5,5,length.out=100),
                 "bbp555" = seq(0.001,0.05,length.out=100)
                 , "H" = seq(0.5,10,length.out=100),
                  "sd" = seq(0.0001,0.1,length.out=100)
                )

#da = expand.grid(da)

da$LL = apply(da,1,log_like_graph_applybase)
da[which.min(da$LL),]
Fit.input

# Parameter Mesh creation
M <- mesh(da$adg443, da$H); z <- da$bbp555 #user-defined

x_surf <- M$x
y_surf <- M$y
z_surf = rep.col(x = z, n = length(z))

ll_surf = rep.col(x=da$LL, length(da$LL))

H = da$H; adg443 = da$adg443
z = outer(adg443, H,log_like_graph_outer) #This is the correct approach to calculate the matrix of lklhood
LL = z

colnames(LL) = signif(da$adg443, digits = 3)
rownames(LL) = signif(da$H, digits = 3)

#Get the parameter index with highest Posterior
max_index <- which.min(LL)

# Convert the index to row and column number
row_num <- floor((max_index - 1) / nrow(LL)) + 1
col_num <- max_index %% ncol(LL)
cat("Highest Posterior obtained for adg443 = ", da$adg443[col_num], ", H = ", da$H[row_num], " with posterior value = ", LL[row_num, col_num])

#Save the Posterior Matrix
write.table(LL, file = "./outputs/LL_param_space.csv", sep = ",", quote = F, row.names = T,
            col.names = T)

#----------------------------------------------
#Plot Surfaces
#----------------------------------------------
library(viridis)
png(file = paste0("./outputs/param_space.png"),  units = "in", width = 8, height = 6, res=300)
par(mai = c(0.30, 0.1, 0.1, 1))

#Plot with Surf3D
plot3D::surf3D(x =log10(x_surf) ,y = (y_surf) ,
               #z = log(chain$LP[1:(i+1)]),
               z=log10(LL),
               #xlim = c(-1,1), ylim=c(-1, 1),
               xlab = "", ylab="", zlab = "",
               xlim=c(-0.3,0.6), ylim=c(0,10),
               add = F, phi = 15, bty = "g", type = "b", theta =230, #30,
               ticktype = "detailed", pch = 20,
               col = turbo(100),
               #lighting = TRUE,
               clab = c(" "),
               colkey = T, 
               colvar =log10(LL),
               contour= TRUE,
               cex = c(0.5, 1, 1.5)
)
scatter3D(x = log10(da$adg443[col_num]), y = da$H[row_num], z= log10(LL[col_num, row_num]), 
         theta = 230, phi=15, add=T, col = "darkgreen", pch=23, lwd=5)
legend(x = "bottomright",          # Position
       legend = expression(paste("MAP(", phi["par"],")")),  # Legend texts
       lty = c(NA),           # Line types
       pch=23,
       col = c("darkgreen"),           # Line colors
       fill = c("darkgreen"),
       box.col="white",
       box.lwd=0,
       lwd = 2,
       #inset = c(0,1),
       cex=0.9)

dev.off()
text(cex=1.3,x = -0.22, y = -0.5, srt = -45 , labels = expression(paste(italic(a)["dg"](443))))
text(cex=1.3,x = 0.22, y = -0.5, srt = 45 , labels = expression(paste(italic(H))))
text(cex=1.3,x = -0.30, y = -0.03, srt = -90 , labels = expression(paste(italic(P)["posterior"])))



#Plot with PlotLy
library(plotly)
axx <- list(
  range = c(min(da$adg443), max(da$adg443)),
  title = "adg(443)"
)

axy <- list(
  range = c(min(da$H), max(da$H)),
  title = "H"
)

axz <- list(
  title = "LogLL"
)

fig <- plot_ly(z = ~(LL)) %>% add_surface(
  colorscale = list(c(0,1),c(rainbow(2))),
  contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  )
)
fig <- fig %>% layout(
  scene = list(xaxis=axx,yaxis=axy,zaxis=axz,
               camera=list(
                 eye = list(x=1.87, y=0.88, z=-0.64)
               )
  )
)

fig


