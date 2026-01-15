## code to prepare `r_b_class_egsl` dataset goes here
library(readr)

r_b_gamache <- read_csv(fs::path_package("SABER", "extdata", "r_b_gamache.csv"))

usethis::use_data(r_b_gamache, overwrite = TRUE)
