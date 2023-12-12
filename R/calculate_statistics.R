

calc_inversion_metrics <- function(actual, predicted){
  
  
  median_ratio = median(log(predicted/actual))
  
  bias = 100*median_ratio*(exp(abs(median_ratio)-1))
  
  uncertainty = 100*(exp(median(abs(log(predicted/actual))))-1)
  
  rmsle = Metrics::rmsle(actual = actual, predicted = predicted)
  
  
  errorstat = c("median_ratio" = median_ratio, "bias" =bias ,
                "uncertainity" = uncertainty, "rmsle" = rmsle)
  
  return(errorstat)
  
}

#======================================
# Read all the data
#======================================
# #in situ
# inv_chl_insitu_saber = read.csv("./outputs/.csv")
# inv_adg443_insitu_saber
# inv_bbp555_insitu_saber
# 
# 
# inv_chl_insitu_qaa
# inv_adg443_insitu_qaa
# inv_bbp555_insitu_qaa
# 
# #synthetic
# inv_chl_synthetic_saber
# inv_adg443_synthetic_saber
# inv_bbp555_synthetic_saber
# 
# inv_chl_synthetic_qaa
# inv_adg443_synthetic_qaa
# inv_bbp555_synthetic_qaa

#======================================
#For [chl]
#======================================

#---------------
# in situ
#---------------
chl_na_idx = is.na(inv_chl_insitu_saber$actual)

chl_insitu_saber_acc = calc_inversion_metrics(actual = inv_chl_insitu_saber$actual[chl_na_idx == F],
                       predicted = inv_chl_insitu_saber$predicted[chl_na_idx == F])

chl_insitu_qaa_acc = calc_inversion_metrics(actual = inv_chl_insitu_qaa$actual[chl_na_idx == F],
                       predicted = inv_chl_insitu_qaa$predicted[chl_na_idx == F])

#---------------
# synthetic
#---------------
chl_synth_saber_acc = calc_inversion_metrics(actual = inv_chl_synthetic_saber$actual,
                       predicted = inv_chl_synthetic_saber$predicted)

chl_synth_qaa_acc = calc_inversion_metrics(actual = inv_chl_synthetic_qaa$actual,
                       predicted = inv_chl_synthetic_qaa$predicted)

#======================================
#For adg(443)
#======================================

#---------------
# in situ
#---------------
adg_na_idx = is.na(inv_adg443_insitu_saber$actual)

adg443_insitu_saber_acc = calc_inversion_metrics(actual = 
                                                inv_adg443_insitu_saber$actual[adg_na_idx == F],
                       predicted = inv_adg443_insitu_saber$predicted[adg_na_idx == F])

adg443_insitu_qaa_acc = calc_inversion_metrics(actual = 
                                                inv_adg443_insitu_qaa$actual[adg_na_idx == F],
                       predicted = inv_adg443_insitu_qaa$predicted[adg_na_idx == F])

#---------------
# synthetic
#---------------
adg443_synth_saber_acc = calc_inversion_metrics(actual = 
                                                   inv_adg443_synthetic_saber$actual,
                       predicted = inv_adg443_synthetic_saber$predicted)

adg443_synth_qaa_acc = calc_inversion_metrics(actual = 
                                                  inv_adg443_synthetic_qaa$actual,
                       predicted = inv_adg443_synthetic_qaa$predicted)

#======================================
#For bbp(555)
#======================================

#---------------
# in situ
#---------------
bbp_na_idx = is.na(inv_bbp555_insitu_saber$actual)

bbp555_insitu_saber_acc = calc_inversion_metrics(actual = 
                                                   inv_bbp555_insitu_saber$actual[bbp_na_idx == F],
                       predicted = inv_bbp555_insitu_saber$predicted[bbp_na_idx == F])

bbp555_insitu_qaa_acc = calc_inversion_metrics(actual = 
                                                   inv_bbp555_insitu_qaa$actual[bbp_na_idx == F],
                       predicted = inv_bbp555_insitu_qaa$predicted[bbp_na_idx == F])

#---------------
# synthetic
#---------------
bbp555_synth_saber_acc = calc_inversion_metrics(actual = inv_bbp555_synthetic_saber$actual,
                       predicted = inv_bbp555_synthetic_saber$predicted)

bbp555_synth_qaa_acc = calc_inversion_metrics(actual = inv_bbp555_synthetic_qaa$actual,
                       predicted = inv_bbp555_synthetic_qaa$predicted)

#======================================
#For H
#======================================
H_saber_acc = calc_inversion_metrics(actual = 
                                              H_df_synth$H_actual[H_df_synth$sa_model == "BAYES"],
                       predicted = H_df_synth$H_predicted[H_df_synth$sa_model == "BAYES"])

H_hope_acc = calc_inversion_metrics(actual = 
                         H_df_synth$H_actual[H_df_synth$sa_model == "HOPE"],
                       predicted = H_df_synth$H_predicted[H_df_synth$sa_model == "HOPE"])


#======================================
#Write the accuracy metrics to disc
#======================================
acc_list = rbind(chl_insitu_saber_acc, chl_insitu_qaa_acc,
      chl_synth_saber_acc, chl_synth_qaa_acc,
      
      adg443_insitu_saber_acc, adg443_insitu_qaa_acc,
      adg443_synth_saber_acc, adg443_synth_qaa_acc,
      
      bbp555_insitu_saber_acc, bbp555_insitu_qaa_acc,
      bbp555_synth_saber_acc, bbp555_synth_qaa_acc,
      
      H_saber_acc, H_hope_acc
      )
write.csv(acc_list, "./outputs/accuracy_metrics.csv", quote = F)
