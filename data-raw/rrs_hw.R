## code to prepare `rrs_hw` dataset goes here

library(readr)

rrs_hw <- read_csv(fs::path_package("SABER", "extdata", "2022-07-06_caa_jet-ski.csv"))

usethis::use_data(rrs_hw, overwrite = TRUE)
