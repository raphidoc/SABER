## code to prepare `rrs_hw_gamache` dataset goes here

library(readr)

rrs_hw_gamache <- read_csv(fs::path_package("SABER", "extdata", "2022-07-05_gamache_bay_jet-ski.csv"))

usethis::use_data(rrs_hw_gamache, overwrite = TRUE)
