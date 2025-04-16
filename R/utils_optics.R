# Function to translate Rrs0+ to Rrs0-
surface_rrs_translate <- function(Rrs) {
  rrs <- Rrs / (0.52 + 1.7 * Rrs)
  return(rrs)
}

snell_law <- function(theta_view, theta_sun) { # Function to convert above water to under water geometry

  # Index of refrations (real)
  n_air <- 1 # air index of refration (real part)
  n_w <- 1.33 # water index of refration (real part)

  # Angles from the water

  # from deg to rad
  theta_view <- theta_view * (180 / pi) # rad
  theta_sun <- theta_sun * (180 / pi) # rad

  # angles inside the water in rad
  view_w <- asin((n_air / n_w) * sin(theta_view)) # rad
  sun_w <- asin((n_air / n_w) * sin(theta_sun)) # rad

  # Fresnel Law

  rho_L <- (1 / 2) * abs(((sin(theta_view - view_w)^2) / (sin(theta_view + view_w)^2)) + ((tan(theta_view - view_w)^2) / (tan(theta_view + view_w)^2)))
  return(data.frame("view_w" = view_w, "sun_w" = sun_w, "rho_L" = rho_L))
}
