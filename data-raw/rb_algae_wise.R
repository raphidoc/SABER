## code to prepare `rb_algae_wise` dataset goes here

rb_algae_wise <- read_csv(fs::path_package("SABER", "extdata", "rb_endmembers_algae-wise.csv"))

usethis::use_data(rb_algae_wise, overwrite = TRUE)
