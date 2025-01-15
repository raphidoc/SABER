## code to prepare `absorption_WV` dataset goes here

library(readr)

a_WV <- read_tsv(
  "data-raw/a_WV.txt",
  col_names = c("wavelength", "a_WV")
  )

usethis::use_data(a_WV, overwrite = TRUE)
