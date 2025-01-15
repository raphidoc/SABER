## code to prepare `absorption_O3` dataset goes here

library(readr)

a_O3 <- read_tsv(
  "data-raw/a_O3.txt",
  col_names = c("wavelength", "a_O3")
)

usethis::use_data(a_O3, overwrite = TRUE)
