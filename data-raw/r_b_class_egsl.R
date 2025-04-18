## code to prepare `r_b_class_egsl` dataset goes here
r_b_class_egsl <- read_csv(fs::path_package("SABER", "extdata", "r_b_class_egsl.csv"))

usethis::use_data(r_b_class_egsl, overwrite = TRUE)
