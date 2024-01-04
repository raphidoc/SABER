require(dplyr)
require(readxl)
require(stats4)
require(MASS)
require(dglm)
library(fitdistrplus)
require(Riops)
require(Cops)
library(lmodel2)


# Functions to calculate statistical metrics for inversion ----

# Function for summary statistics for different prior parameters
summary_stat <- function(input_var ) {
  summary_temp = summary(input_var)
  desc_stat = data_frame("min" = summary_temp[1], "max" = summary_temp[length(summary_temp)],
                         "sd" = sd(input_var), "mean" = summary_temp[4], 
                         "median" = summary_temp[3], "mode" = Mode(input_var),
                         "skewness" = skewness(input_var), "kurtosis" = kurtosis(input_var))
  
  return(desc_stat)
}


# Function for goodness of fit for inversion retrieved  parameters for QAA and SABER
goodness_of_fit <- function(actual, predicted, take_log = T) {
  
  
  if (take_log == TRUE) {
    
    #BIAS
    bias_e = Metrics::bias(actual = log10(actual), 
                           predicted = log10(predicted))
    
    
    #%-bias
    p_bias_e = Metrics::percent_bias(actual = log10(actual), 
                                     predicted = log10(predicted))
    
    
    #RMSE
    rmse_e = Metrics::rmse(actual = log10(actual), 
                           predicted = log10(predicted))
    
    
    #R^2
    if (var(actual) == 0) {
      print("Actual has zero variance.")
      actual = actual + rnorm(n = length(actual), 0 ,sd = 0.0001)
    } else if (var(predicted) == 0) {
      print("Predicted has zero variance.")
      predicted = predicted + rnorm(n = length(actual), 0 ,sd = 0.0001)
    } else {
      print("Neither actual nor predicted has zero variance.")
    }
    
    yx.lmodel2_e <- lmodel2(log10(predicted) ~ log10(actual)) 
    
    r_square_e = yx.lmodel2_e$rsquare
    slope_e = yx.lmodel2_e$regression.results$Slope[2]
    
    
  } else {
    
    
    #BIAS
    bias_e = Metrics::bias(actual = (actual), 
                           predicted = (predicted))
    
    
    #%-bias
    p_bias_e = Metrics::percent_bias(actual = (actual), 
                                     predicted = (predicted))
    
    
    #RMSE
    rmse_e = Metrics::rmse(actual = (actual), 
                           predicted = (predicted))
    
    
    # #R^2
    # yx.lmodel2_e <- lmodel2((predicted) ~ (actual)) 
    # 
    # r_square_e = yx.lmodel2_e$rsquare
    # slope_e = yx.lmodel2_e$regression.results$Slope[2]
    
    #R^2
    if (var(actual) == 0) {
      print("Actual has zero variance.")
      actual = actual + rnorm(n = length(actual), 0 ,sd = 0.0001)
    } else if (var(predicted) == 0) {
      print("Predicted has zero variance.")
      predicted = predicted + rnorm(n = length(actual), 0 ,sd = 0.0001)
    } else {
      print("Neither actual nor predicted has zero variance.")
    }
    
    yx.lmodel2_e <- lmodel2(log10(predicted) ~ log10(actual)) 
    
    r_square_e = yx.lmodel2_e$rsquare
    slope_e = yx.lmodel2_e$regression.results$Slope[2]
    
    
    
  }
  
  errorstat = c("bias" = bias_e, "p_bias" =p_bias_e ,
                "rmse" = rmse_e, "r_square" = r_square_e, "slope"=slope_e)
  
  print(errorstat)
  
  return(errorstat)
  
}

#Extended validation metric after Seegers (2018) and Pahlevan (2022)
calc_inversion_metrics <- function(actual, predicted){
  
  
  median_ratio = median(log(predicted/actual))
  
  bias = 100*median_ratio*(exp(abs(median_ratio)-1))
  
  uncertainty = 100*(exp(median(abs(log(predicted/actual))))-1)
  
  rmsle = Metrics::rmsle(actual = actual, predicted = predicted)
  
  #R^2
  if (is.na(var(actual)) == TRUE) {
    
    print("NA variance")
    r_square_e = NA
    slope_e = NA
    
  } else {
    
    if (var(actual) == 0) {
      print("Actual has zero variance.")
      actual = actual + rnorm(n = length(actual), 0 ,sd = 0.0001)
    } else if (var(predicted) == 0) {
      print("Predicted has zero variance.")
      predicted = predicted + rnorm(n = length(actual), 0 ,sd = 0.0001)
    } else {
      print("Neither actual nor predicted has zero variance.")
    }
    
    yx.lmodel2_e <- lmodel2(log10(predicted) ~ log10(actual)) 
    
    r_square_e = yx.lmodel2_e$rsquare
    slope_e = yx.lmodel2_e$regression.results$Slope[2]
    
  }
  
  errorstat = c("median_ratio" = median_ratio, "bias" =bias ,
                "uncertainity" = uncertainty, "rmsle" = rmsle, "R_square" = r_square_e, "slope" =  slope_e)
  
  print(errorstat)
  
  return(errorstat)
  
}

# Perform the statistical analysis for validation ----
## For [chl] ----


### in situ ----


#### SABER ----

chl_na_idx = is.nan(inv_chl_insitu_saber$actual)
zeropos = which(inv_chl_insitu_saber$actual[chl_na_idx == F] == 0)



chl_insitu_saber_acc = calc_inversion_metrics(
  actual = inv_chl_insitu_saber$actual[chl_na_idx == F]%>%.[-zeropos],
  predicted = inv_chl_insitu_saber$predicted[chl_na_idx == F]%>% .[-zeropos])

chl_insitu_saber_gof = goodness_of_fit(
  actual = inv_chl_insitu_saber$actual[chl_na_idx == F]%>%.[-zeropos],
  predicted = inv_chl_insitu_saber$predicted[chl_na_idx == F]%>% .[-zeropos], take_log = F)

# Statistics individually for OWT
chl_insitu_saber_acc_owt <- inv_chl_insitu_saber[chl_na_idx == F,]%>%.[-zeropos,] %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


chl_insitu_saber_acc_owt = as.data.frame(do.call(rbind,chl_insitu_saber_acc_owt[["metrics"]]))
chl_insitu_saber_acc_owt$clust_id = rownames(chl_insitu_saber_acc_owt)
chl_insitu_saber_acc_owt$rrs_type = "insitu"
chl_insitu_saber_acc_owt$inv_model = "SABER"
chl_insitu_saber_acc_owt$var = "chl"
chl_insitu_saber_acc_owt <- chl_insitu_saber_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())


#### QAA ----
chl_na_idx = is.nan(inv_chl_insitu_qaa$actual)
zeropos = which(inv_chl_insitu_qaa$actual[chl_na_idx == F] == 0)

chl_insitu_qaa_acc = calc_inversion_metrics(
  actual = inv_chl_insitu_qaa$actual[chl_na_idx == F]%>%.[-zeropos],
  predicted = inv_chl_insitu_qaa$predicted[chl_na_idx == F]%>% .[-zeropos])

chl_insitu_saber_gof = goodness_of_fit(
  actual = inv_chl_insitu_qaa$actual[chl_na_idx == F]%>%.[-zeropos],
  predicted = inv_chl_insitu_qaa$predicted[chl_na_idx == F]%>% .[-zeropos], take_log = F)

# Statistics individually for OWT
chl_insitu_qaa_acc_owt <- inv_chl_insitu_qaa[chl_na_idx == F,]%>%.[-zeropos,] %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


chl_insitu_qaa_acc_owt = as.data.frame(do.call(rbind,chl_insitu_qaa_acc_owt[["metrics"]]))
chl_insitu_qaa_acc_owt$clust_id = rownames(chl_insitu_qaa_acc_owt)
chl_insitu_qaa_acc_owt$rrs_type = "insitu"
chl_insitu_qaa_acc_owt$inv_model = "QAA"
chl_insitu_qaa_acc_owt$var = "chl"
chl_insitu_qaa_acc_owt <- chl_insitu_qaa_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())

### synthetic----

#### SABER----
chl_synth_saber_acc = calc_inversion_metrics(actual = inv_chl_synthetic_saber$actual,
                       predicted = inv_chl_synthetic_saber$predicted)

chl_synth_saber_gof = goodness_of_fit(actual = inv_chl_synthetic_saber$actual,
                                             predicted = inv_chl_synthetic_saber$predicted, 
                                      take_log = F)
# Statistics individually for OWT
chl_synthetic_saber_acc_owt <- inv_chl_synthetic_saber %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


chl_synthetic_saber_acc_owt = as.data.frame(do.call(rbind,chl_synthetic_saber_acc_owt[["metrics"]]))
chl_synthetic_saber_acc_owt$clust_id = rownames(chl_synthetic_saber_acc_owt)
chl_synthetic_saber_acc_owt$rrs_type = "synthetic"
chl_synthetic_saber_acc_owt$inv_model = "SABER"
chl_synthetic_saber_acc_owt$var = "chl"
chl_synthetic_saber_acc_owt <- chl_synthetic_saber_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())

#### QAA----
chl_synth_qaa_acc = calc_inversion_metrics(actual = inv_chl_synthetic_qaa$actual,
                       predicted = inv_chl_synthetic_qaa$predicted)

chl_synth_qaa_gof = goodness_of_fit(actual = inv_chl_synthetic_qaa$actual,
                                           predicted = inv_chl_synthetic_qaa$predicted, take_log = F)

# Statistics individually for OWT
chl_synthetic_qaa_acc_owt <- inv_chl_synthetic_qaa %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


chl_synthetic_qaa_acc_owt = as.data.frame(do.call(rbind,chl_synthetic_qaa_acc_owt[["metrics"]]))
chl_synthetic_qaa_acc_owt$clust_id = rownames(chl_synthetic_qaa_acc_owt)
chl_synthetic_qaa_acc_owt$rrs_type = "synthetic"
chl_synthetic_qaa_acc_owt$inv_model = "QAA"
chl_synthetic_qaa_acc_owt$var = "chl"
chl_synthetic_qaa_acc_owt <- chl_synthetic_qaa_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())


## For adg(443) ----


### in situ ----


#### SABER ----
adg_na_idx = is.na(inv_adg443_insitu_saber$actual)
zeropos = which(inv_adg443_insitu_saber$actual[adg_na_idx == F] == 0)
adg443_insitu_saber_acc = calc_inversion_metrics(actual = 
                                                inv_adg443_insitu_saber$actual[adg_na_idx == F]%>%.[-zeropos],
                       predicted = inv_adg443_insitu_saber$predicted[adg_na_idx == F]%>%.[-zeropos])

adg443_insitu_saber_acc_gof = goodness_of_fit(actual =inv_adg443_insitu_saber$actual[adg_na_idx == F]%>%.[-zeropos],
                      predicted = inv_adg443_insitu_saber$predicted[adg_na_idx == F]%>%.[-zeropos], 
                      take_log = F)

# Statistics individually for OWT
adg443_insitu_saber_acc_owt <- inv_adg443_insitu_saber[adg_na_idx == F,]%>%.[-zeropos,] %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


adg443_insitu_saber_acc_owt = as.data.frame(do.call(rbind,adg443_insitu_saber_acc_owt[["metrics"]]))
adg443_insitu_saber_acc_owt$clust_id = rownames(adg443_insitu_saber_acc_owt)
adg443_insitu_saber_acc_owt$rrs_type = "insitu"
adg443_insitu_saber_acc_owt$inv_model = "SABER"
adg443_insitu_saber_acc_owt$var = "adg443"
adg443_insitu_saber_acc_owt <- adg443_insitu_saber_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())

#### QAA ----
adg_na_idx = is.na(inv_adg443_insitu_qaa$actual)
zeropos = which(inv_adg443_insitu_qaa$actual[adg_na_idx == F] == 0)
adg443_insitu_qaa_acc = calc_inversion_metrics(actual = 
                                                inv_adg443_insitu_qaa$actual[adg_na_idx == F]%>%.[-zeropos],
                       predicted = inv_adg443_insitu_qaa$predicted[adg_na_idx == F]%>%.[-zeropos])

adg443_insitu_qaa_gof = goodness_of_fit(actual = 
                        inv_adg443_insitu_qaa$actual[adg_na_idx == F]%>%.[-zeropos],
                       predicted = inv_adg443_insitu_qaa$predicted[adg_na_idx == F]%>%.[-zeropos], 
                       take_log = F)

# Statistics individually for OWT
adg443_insitu_qaa_acc_owt <- inv_adg443_insitu_qaa[adg_na_idx == F,]%>%.[-zeropos,] %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


adg443_insitu_qaa_acc_owt = as.data.frame(do.call(rbind,adg443_insitu_qaa_acc_owt[["metrics"]]))
adg443_insitu_qaa_acc_owt$clust_id = rownames(adg443_insitu_qaa_acc_owt)
adg443_insitu_qaa_acc_owt$rrs_type = "insitu"
adg443_insitu_qaa_acc_owt$inv_model = "QAA"
adg443_insitu_qaa_acc_owt$var = "adg443"
adg443_insitu_qaa_acc_owt <- adg443_insitu_qaa_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())

### synthetic ----

#### SABER ----
adg443_synth_saber_acc = calc_inversion_metrics(actual = 
                                                   inv_adg443_synthetic_saber$actual,
                       predicted = inv_adg443_synthetic_saber$predicted)

adg443_synth_saber_gof = goodness_of_fit(actual = 
                                                  inv_adg443_synthetic_saber$actual,
                                                predicted = inv_adg443_synthetic_saber$predicted)

# Statistics individually for OWT
adg443_synthetic_saber_acc_owt <- inv_adg443_synthetic_saber %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


adg443_synthetic_saber_acc_owt = as.data.frame(do.call(rbind,adg443_synthetic_saber_acc_owt[["metrics"]]))
adg443_synthetic_saber_acc_owt$clust_id = rownames(adg443_synthetic_saber_acc_owt)
adg443_synthetic_saber_acc_owt$rrs_type = "synthetic"
adg443_synthetic_saber_acc_owt$inv_model = "SABER"
adg443_synthetic_saber_acc_owt$var = "adg443"
adg443_synthetic_saber_acc_owt <- adg443_synthetic_saber_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())

#### QAA ----
adg443_synth_qaa_acc = calc_inversion_metrics(actual = 
                                                  inv_adg443_synthetic_qaa$actual,
                       predicted = inv_adg443_synthetic_qaa$predicted)

adg443_synth_qaa_gof = goodness_of_fit(actual = 
                                                inv_adg443_synthetic_qaa$actual,
                                              predicted = inv_adg443_synthetic_qaa$predicted)

# Statistics individually for OWT
adg443_synthetic_qaa_acc_owt <- inv_adg443_synthetic_qaa %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


adg443_synthetic_qaa_acc_owt = as.data.frame(do.call(rbind,adg443_synthetic_qaa_acc_owt[["metrics"]]))
adg443_synthetic_qaa_acc_owt$clust_id = rownames(adg443_synthetic_qaa_acc_owt)
adg443_synthetic_qaa_acc_owt$rrs_type = "synthetic"
adg443_synthetic_qaa_acc_owt$inv_model = "QAA"
adg443_synthetic_qaa_acc_owt$var = "adg443"
adg443_synthetic_qaa_acc_owt <- adg443_synthetic_qaa_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())


## For bbp(555) ----


### in situ ----

#### SABER ----
bbp_na_idx = is.na(inv_bbp555_insitu_saber$actual)

bbp555_insitu_saber_acc = calc_inversion_metrics(actual = 
                                                   inv_bbp555_insitu_saber$actual[bbp_na_idx == F],
                       predicted = inv_bbp555_insitu_saber$predicted[bbp_na_idx == F])

# Statistics individually for OWT
bbp555_insitu_saber_acc_owt <- inv_bbp555_insitu_saber[bbp_na_idx == F,] %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


bbp555_insitu_saber_acc_owt = as.data.frame(do.call(rbind,bbp555_insitu_saber_acc_owt[["metrics"]]))
bbp555_insitu_saber_acc_owt$clust_id = rownames(bbp555_insitu_saber_acc_owt)
bbp555_insitu_saber_acc_owt$rrs_type = "insitu"
bbp555_insitu_saber_acc_owt$inv_model = "SABER"
bbp555_insitu_saber_acc_owt$var = "bbp555"
bbp555_insitu_saber_acc_owt <- bbp555_insitu_saber_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())

#### QAA ----
bbp555_insitu_qaa_acc = calc_inversion_metrics(actual = 
                                                   inv_bbp555_insitu_qaa$actual[bbp_na_idx == F],
                       predicted = inv_bbp555_insitu_qaa$predicted[bbp_na_idx == F])

# Statistics individually for OWT
bbp555_insitu_qaa_acc_owt <- inv_bbp555_insitu_qaa[bbp_na_idx == F,] %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


bbp555_insitu_qaa_acc_owt = as.data.frame(do.call(rbind,bbp555_insitu_qaa_acc_owt[["metrics"]]))
bbp555_insitu_qaa_acc_owt$clust_id = rownames(bbp555_insitu_qaa_acc_owt)
bbp555_insitu_qaa_acc_owt$rrs_type = "insitu"
bbp555_insitu_qaa_acc_owt$inv_model = "QAA"
bbp555_insitu_qaa_acc_owt$var = "bbp555"
bbp555_insitu_qaa_acc_owt <- bbp555_insitu_qaa_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())


### synthetic ----

#### SABER ----
bbp555_synth_saber_acc = calc_inversion_metrics(actual = inv_bbp555_synthetic_saber$actual,
                       predicted = inv_bbp555_synthetic_saber$predicted)

# Statistics individually for OWT
bbp555_synthetic_saber_acc_owt <- inv_bbp555_synthetic_saber %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


bbp555_synthetic_saber_acc_owt = as.data.frame(do.call(rbind,bbp555_synthetic_saber_acc_owt[["metrics"]]))
bbp555_synthetic_saber_acc_owt$clust_id = rownames(bbp555_synthetic_saber_acc_owt)
bbp555_synthetic_saber_acc_owt$rrs_type = "synthetic"
bbp555_synthetic_saber_acc_owt$inv_model = "SABER"
bbp555_synthetic_saber_acc_owt$var = "bbp555"
bbp555_synthetic_saber_acc_owt <- bbp555_synthetic_saber_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())

#### QAA ----
bbp555_synth_qaa_acc = calc_inversion_metrics(actual = inv_bbp555_synthetic_qaa$actual,
                       predicted = inv_bbp555_synthetic_qaa$predicted)

# Statistics individually for OWT
bbp555_synthetic_qaa_acc_owt <- inv_bbp555_synthetic_qaa %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


bbp555_synthetic_qaa_acc_owt = as.data.frame(do.call(rbind,bbp555_synthetic_qaa_acc_owt[["metrics"]]))
bbp555_synthetic_qaa_acc_owt$clust_id = rownames(bbp555_synthetic_qaa_acc_owt)
bbp555_synthetic_qaa_acc_owt$rrs_type = "synthetic"
bbp555_synthetic_qaa_acc_owt$inv_model = "QAA"
bbp555_synthetic_qaa_acc_owt$var = "bbp555"
bbp555_synthetic_qaa_acc_owt <- bbp555_synthetic_qaa_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())


acc_deep_param = rbind(chl_insitu_saber_acc_owt, chl_insitu_qaa_acc_owt,
                       chl_synthetic_saber_acc_owt, chl_synthetic_qaa_acc_owt,
                       adg443_insitu_saber_acc_owt, adg443_insitu_qaa_acc_owt,
                       adg443_synthetic_saber_acc_owt, adg443_synthetic_qaa_acc_owt,
                       bbp555_insitu_saber_acc_owt, bbp555_insitu_qaa_acc_owt, 
                       bbp555_synthetic_saber_acc_owt, bbp555_synthetic_qaa_acc_owt)

#write.csv(file = "./outputs/accuracy_metrics_OWT.csv", acc_deep_param, quote = F, row.names = F)

## For H ----

#### SABER ----
H_saber_acc = calc_inversion_metrics(actual = 
                                              H_df_synth$H_actual[H_df_synth$sa_model == "BAYES"],
                       predicted = H_df_synth$H_predicted[H_df_synth$sa_model == "BAYES"])

# Statistics individually for OWT
H_saber_acc_owt <- H_df_synth[H_df_synth$sa_model == "BAYES",] %>%
  group_by(clust_id) %>%
  summarize(metrics = list(calc_inversion_metrics(actual, predicted)))


H_saber_acc_owt = as.data.frame(do.call(rbind,H_saber_acc_owt[["metrics"]]))
H_saber_acc_owt$clust_id = rownames(H_saber_acc_owt)
H_saber_acc_owt$rrs_type = "synthetic"
H_saber_acc_owt$inv_model = "SABER"
H_saber_acc_owt$var = "H"
H_saber_acc_owt <- H_saber_acc_owt %>%
  dplyr::select(clust_id,rrs_type, inv_model, var,  everything())

#### HOPE ----
H_hope_acc = calc_inversion_metrics(actual = 
                         H_df_synth$H_actual[H_df_synth$sa_model == "HOPE"],
                       predicted = H_df_synth$H_predicted[H_df_synth$sa_model == "HOPE"])



#Write the accuracy metrics to disc ----
acc_list = rbind(chl_insitu_saber_acc, chl_insitu_qaa_acc,
      chl_synth_saber_acc, chl_synth_qaa_acc,
      
      adg443_insitu_saber_acc, adg443_insitu_qaa_acc,
      adg443_synth_saber_acc, adg443_synth_qaa_acc,
      
      bbp555_insitu_saber_acc, bbp555_insitu_qaa_acc,
      bbp555_synth_saber_acc, bbp555_synth_qaa_acc,
      
      H_saber_acc, H_hope_acc
      )
write.csv(acc_list, "./outputs/accuracy_metrics_v3.csv", quote = F)
