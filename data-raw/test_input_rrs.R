## code to prepare `test_input_rrs` dataset goes here

library(readr)

test_input_rrs <- read_csv(
  "data-raw/test_input_rrs.csv"
)

usethis::use_data(test_input_rrs, overwrite = TRUE)
