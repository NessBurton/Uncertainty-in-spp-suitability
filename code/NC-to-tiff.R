# kelvin to celsius
to_celsius <- function (x) {
  x - 273.15
}

# read with raster and convert to stars
read_nc <- function(x, path, crs = 27700) {
  out <- raster::brick(file.path(path, x))
  out <- sf::st_set_crs(stars::st_as_stars(out), crs)
  out
}

# read raster, convert to stars and write to tif
nc_to_tif <- function(x, path, path_out, crs = 27700) {
  
  out <- read_nc(x, path, crs)
  if (grepl("min_tasmin|max_tasmax", x)) {
    out <- out - 273.15
  }
  stars::write_stars(out, file.path(path_out, sub("\\.nc$", "\\.tif", x)))
  
  if (grepl("aws", x)) {
    aws_depth <- list()
    aws_depth <- within(aws_depth, {
      thick <- round(st_get_dimension_values(out, "band"), 2)
      to <- cumsum(thick)
      from <- to - thick
      labels <- paste0(from, "-", to, " m")
      units <- "m^3/m^3"
    })
    saveRDS(aws_depth, file.path(path_out, sub("\\.nc$", "\\.RDS", x)))
  }
}