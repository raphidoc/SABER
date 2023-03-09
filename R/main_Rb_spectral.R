#==================================================================================================
# main_Rb_spectral.R performs spectral retrieval of bottom reflectance (Rb) based on Radiative 
# transfer in optically complex case-II waters. There are two modes to retrieve Rb:

# 1. Model the Remote-Sensing Reflectance (Rrs) in optically complex waters from user-given [chl],
# CDOM and NAP absorption at 440nm wavelength and depth and subtract the model value from observed 
# Rrs using algebric reparametrization of AM03. 

# 2. Model the Remote-Sensing Reflectance (Rrs) in optically complex waters from user-given spectral
# IOP and depth and subtract the model value using algebric reparametrization of AM03.

#The function to calculate Rb is referred as "Saber_retrieve_rb_wise".
#==================================================================================================

###======================================================================
# 1. Retrieve IOPs, AOPs and BGC variables for the station under process
###======================================================================

#Input the WISE-Man IOP QC sheets for both sets of IOP cages
kildir_iop = read.csv("/home/musk0001/Kildir_IOP.Process_Log.csv", header = T)
saucier_iop = read.csv("/home/musk0001/saucierQC.csv", header = T)

#Input the depth statistics and simulation notes from WISE-Man IOP data
wise_depth_stat = read.csv("/home/musk0001/Depth_wise.csv", header = T)
wise_depth_stat$correct_statname = gsub(wise_depth_stat$Station, pattern = "\\_", replacement = "-") 
wise_depth_stat <- wise_depth_stat %>% dplyr::select(correct_statname, everything())

#Select the QC passed stations for WISE-Man
wise_depth_stat_shallow = wise_depth_stat[wise_depth_stat$Shallow == "T",]
kildir_iop_qc = kildir_iop[kildir_iop$ASPH == "Y" & kildir_iop$HS6 == "Y",]
saucier_iop_qc = saucier_iop[saucier_iop$station.keep == "TRUE",]

qc_statnames_bind = rbind(list(kildir_iop_qc$StationID), 
                                                       list(saucier_iop_qc$Station))

qc_statnames_bind = data.frame("qcstat" = as.character(unlist(qc_statnames_bind)))

shallow_qc_keep_idx = vector()

for (i in 1:length(qc_statnames_bind$qcstat)) {
  tryCatch({
  shallow_qc_keep_idx[i] = grep(x =  wise_depth_stat_shallow$correct_statname, 
                                pattern = qc_statnames_bind$qcstat[i])
  
  if (is.na(shallow_qc_keep_idx[i])) {
    cat(paste0("\033[0;31m","The station ", wise_depth_stat_shallow$correct_statname[i]," is not shallow, coerce to NA","\033[0m","\n"))
    shallow_qc_keep_idx[i] = NA
  } }, 
    error = function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

shallow_qc_stations = qc_statnames_bind$qcstat[!is.na(shallow_qc_keep_idx)]
shallow_qc_stations = shallow_qc_stations[-(3:4)] #To get rid of stations with no Rrs

#=============================================================
## 2. BATCH Rb retrieval for WISE-Man in situ
#=============================================================
Saber_generate_rb_batch(staionlist = shallow_qc_stations)
