## code to prepare `extraterrestrial_Ed` dataset goes here

library(readr)

F0 <- read_tsv(
  "data-raw/F0.txt",
  col_names = c("wavelength", "F0")
)

usethis::use_data(F0, overwrite = TRUE)
