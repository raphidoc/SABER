## code to prepare `absorption_O2` dataset goes here

library(readr)

a_O2 <- read_tsv(
  "data-raw/a_O2.txt",
  col_names = c("wavelength", "a_O2")
)

usethis::use_data(a_O2, overwrite = TRUE)
