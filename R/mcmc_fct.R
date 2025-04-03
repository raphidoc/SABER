#' Create prior distribution
#'
#' @author Soham Mukherjee
#' @param use_prior choose dataset from wich to construct prior
#' @returns input data prior distribution to be used by ?
#' @importFrom readr read_csv
#'

create_prior_distribution <- function(
    use_prior = c("ioccg", "nomad", "wiseman"),
    plot_diag = F,
    distribution_fit = c("weibull", "lognormal", "normal"),
    truncate_chl = c(0.5, 30),
    truncate_adg = c(0.1, 5),
    truncate_bbp = c(0.002, 0.01),
    sample_count = 100
){

  if (use_prior == "ioccg") {

    stop("Not implemented yet")

    # HL_deep_iop = read_IOCCG_data()
    # chl_sample = HL_deep_iop$chl
    # adg_sample = HL_deep_iop$acdom440+HL_deep_iop$anap440
    # bbp_sample = HL_deep_iop$bbp550
    #
    # #Truncate the vector of prior parameters
    # if (!all(is_na(truncate_chl))) {
    #   chl_sample = chl_sample[chl_sample >= truncate_chl[1] & chl_sample <= truncate_chl[2]]
    #
    # }
    #
    # if (!all(is_na(truncate_adg))) {
    #   adg_sample = adg_sample[adg_sample >= truncate_adg[1] & adg_sample <= truncate_adg[2]]
    #
    # }
    #
    # if (!all(is_na(truncate_bbp))) {
    #   bbp_sample = bbp_sample[bbp_sample >= truncate_bbp[1] & bbp_sample <= truncate_bbp[2]]
    #
    # }
    #
    # if (!is_na(sample_count)) {
    #   set_seed(123)
    #   chl_sample = sample(chl_sample, size = sample_count)
    #   adg_sample = sample(adg_sample, size = sample_count)
    #   bbp_sample = sample(bbp_sample, size = sample_count)
    # }
    #
    # #try to fit the user-defined distribution
    # fit_chl_norm <- fitdistrplus::fitdist(chl_sample, distribution_fit)
    #
    # fit_acdm440_norm <- fitdistrplus::fitdist(adg_sample,distribution_fit)
    #
    # fit_bbp550_norm <- fitdistrplus::fitdist(bbp_sample,distribution_fit)

  }

  if (use_prior == "wiseman") {

    stop("Not implemented yet")

    # #Load WISE-Man 2019 in situ Prior data
    # bgc_data <- read_csv("_/data/wise-man/biogeochemistry.csv", header = T)
    # acdom_data <- read_csv("_/data/ag_long.csv", header = T)
    # anap_data <- read_csv("_/data/anap_long.csv", header = T)
    #
    # chl_sample <- bgc_data$chl
    #
    # acdom_440_sample <- acdom_data$ag[acdom_data$wavelength == "440"]
    #
    # anap_440_sample <- anap_data$ad[anap_data$wavelength == "440"]
    #
    # adg_sample = acdom_440_sample + anap_440_sample
    #
    # #Truncate the vector of prior parameters
    # if (!all(is_na(truncate_chl))) {
    #   chl_sample = chl_sample[chl_sample >= truncate_chl[1] & chl_sample <= truncate_chl[2]]
    #
    # }
    #
    # if (!all(is_na(truncate_adg))) {
    #   adg_sample = adg_sample[adg_sample >= truncate_adg[1] & adg_sample <= truncate_adg[2]]
    #
    # }
    #
    # if (!is_na(sample_count)) {
    #   set_seed(123)
    #   chl_sample = sample(chl_sample, size = sample_count)
    #   adg_sample = sample(adg_sample, size = sample_count)
    # }
    #
    # #try to fit the user-defined distribution
    # fit_chl_norm <- fitdistrplus::fitdist(chl_sample, distribution_fit)
    #
    # fit_acdm440_norm <- fitdistrplus::fitdist(adg_sample,distribution_fit)

  }

  browser()

  if (use_prior == "nomad") {

    #Load NOMAD in situ Prior data
    bgc_data <- read_csv("_/data/nomad/dataset_simplified_csv",
                         header = T, sep = ",")


    chl_sample <- bgc_data$chl
    chl_sample = chl_sample[!chl_sample %in% -999]


    acdom_440_sample <- bgc_data$ag443
    acdom_440_sample = acdom_440_sample[!acdom_440_sample %in% -999]


    anap_440_sample <- bgc_data$ad443
    anap_440_sample = anap_440_sample[!anap_440_sample %in% -999]

    adg_sample = acdom_440_sample + anap_440_sample

    bbp_sample = bgc_data$bb555
    bbp_sample = bbp_sample[!bbp_sample %in% -999]

    #Truncate the vector of prior parameters
    if (!is.null(truncate_chl)) {
      chl_sample = chl_sample[chl_sample >= truncate_chl[1] & chl_sample <= truncate_chl[2]]
    }

    if (all(!is_na(truncate_adg))) {
      adg_sample = adg_sample[adg_sample >= truncate_adg[1] & adg_sample <= truncate_adg[2]]
    }

    if (all(!is_na(truncate_bbp))) {
      bbp_sample = bbp_sample[bbp_sample >= truncate_bbp[1] & bbp_sample <= truncate_bbp[2]]
    }

    if (!is_na(sample_count)) {
      set_seed(123)
      chl_sample = sample(chl_sample, size = sample_count)
      adg_sample = sample(adg_sample, size = sample_count)
      bbp_sample = sample(bbp_sample, size = sample_count)
    }

    #try to fit the user-defined distribution
    fit_chl_norm <- fitdistrplus::fitdist(chl_sample, distribution_fit)

    fit_acdm440_norm <- fitdistrplus::fitdist(adg_sample,distribution_fit)

    fit_bbp550_norm <- fitdistrplus::fitdist(bbp_sample,distribution_fit)
  }

  if (plot_diag == TRUE && use_wise_prior == TRUE) {

    plot(fit_chl_norm)
    plot(fit_acdm440_norm)

  } else {
    plot(fit_chl_norm)
    plot(fit_acdm440_norm)
    plot(fit_bbp550_norm)
  }

  if (plot_diag == TRUE && use_wise_prior == TRUE) {

    return(list("fit_chl"=fit_chl_norm, "obs_chl"= chl_sample,
                "fit_acdm440"=fit_acdm440_norm, "obs_acdm440"= adg_sample,
                #"fit_bbp555"=fit_bbp550_norm, "obs_bbp550" = bbp_sample,
                "distribution_fitted" = distribution_fit
    ))

  } else {
    return(list("fit_chl"=fit_chl_norm, "obs_chl"= chl_sample,
                "fit_acdm440"=fit_acdm440_norm, "obs_acdm440"= adg_sample,
                "fit_bbp555"=fit_bbp550_norm, "obs_bbp550" = bbp_sample,
                "distribution_fitted" = distribution_fit
    ))
  }

}

#-------------------------------------------------------------------------
#Create Prior density function
#-------------------------------------------------------------------------
prior_bayes <- function(param, pop_sd = pop_sd, verbose = F){

  chl = param[1]
  acdom_440 = param[2]
  anap_440 = param[3]
  x_sd = param[4]
  chl_prior = dweibull(x=chl,shape = fit_chl_norm$estimate[1],
                       scale =fit_chl_norm$estimate[2] , log = T)

  acdom440_prior = dweibull(x=acdom_440,shape = fit_acdom440_norm$estimate[1],
                            scale =fit_acdom440_norm$estimate[2] , log = T)

  anap440_prior = dweibull(x=anap_440,shape = fit_anap440_norm$estimate[1],
                           scale =fit_anap440_norm$estimate[2] , log = T)

  lklhood_prior = dunif(x = x_sd, min = 0.000001, max  = 0.001, log = T)

  if (verbose == T) {
    if (pop_sd == TRUE) {
      print("Prior of forward model noise is not fitted")
      return(chl_prior+acdom440_prior+anap440_prior)

    } else {
      print("Prior of forward model noise is fitted")
      return(chl_prior+acdom440_prior+anap440_prior+lklhood_prior)
    }
  }
}
