
# date: 16-04-24
# author: VB
# desc: script to process climate data (CHESS baseline 1991-2011 from CEH), convert nc files to tif,
# calculate CMD for ESC

### load libs ----
library(dplyr)
library(stars)
library(raster)
library(ncmeta)

# load in functions from #CEH script
source(file.path(wd,"code","NC-to-tiff.R"))

### working dirs ----

wd <- "~/Documents/Woodland-Trust/Data-Analysis/SoWT2_Spp-uncertainty" # Mac path
dirData <-  paste0(wd,"/data-raw/")
dirScratch <-  paste0(wd,"/data-scratch/")
dirOut <- paste0(wd,"/data-out/")

# list .nc files
path_baseline <- file.path(dirData,"chess_baseline")
base_files <- list.files(path_baseline, pattern = "\\.nc$")

# create metadata
description <- c(
  gdd = "growing degree days and number of days",
  precip = "precipitation",
  pet = "Penman-Monteith potential evapotranspiration, for a well-watered grass surface"
)

metadata <- tibble::tibble(file = base_files) %>%
  tidyr::extract(
    file,
    c("variable", "summary", "temp_resolution", "from", "to"),
    "chess-[[:alpha:]]+_(.+)_[[:alpha:]]{2}_1km_20yr-([[:alpha:]]+)-([[:alpha:]]+)_([[:digit:]]+)-([[:digit:]]+)",
    remove = FALSE
  ) %>%
  mutate(
    from_year = substr(from, 1, 4),
    from_month = substr(from, 5, 6),
    to_year = substr(to, 1, 4),
    to_month = substr(to, 5, 6),
    description = description[variable]
  )

write.csv(metadata, file.path(path_baseline, "metadata.csv"), row.names = FALSE)

### process nc files ----

# gdd (growing degree days, sometimes called accumulated temperature AT)
gdd1991_2011 <- stars::read_ncdf(file.path(path_baseline,metadata$file[4]), 27700)
plot(gdd1991_2011)
stars::write_stars(gdd1991_2011, dsn = paste0(dirScratch,"/chess_gdd_1991_2011_annual.tif"))

# precipitation
prec1991_2011 <- stars::read_ncdf(file.path(path_baseline,metadata$file[8]), 27700)
plot(prec1991_2011)

# potential evapotranspiration
pet1991_2011 <- stars::read_ncdf(file.path(path_baseline,metadata$file[12]), path_baseline, 27700)
plot(pet1991_2011)

# total monthly potential evapotranspiration
# the data provided is mean monthly potential evapotranspiration (per day), it need to be
# transformed to total. Need to multiply by the number of days in a month.

# convert mean monthly potential evapotranspiration to total
pet_to_total <- function (x) {
  month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  dims <- dim(x)
  total <- x * array(rep(month_days, each = prod(dims[1:2])), dims)
  total <- st_set_dimensions(total, "time", tolower(month.abb))
  total
}

pet <- pet_to_total(pet1991_2011)

plot(pet, key.pos = 4, key.width = lcm(1), box_col = "white", col = hcl.colors(10))

stars::write_stars(pet, dsn = paste0(dirScratch,"chess_pet_1991_2011_monthly.tif"))

# total monthly precipitation
# the data provided is mean monthly precipitation (per second)
# needs to be transformed to total - multiply by the number seconds in a month

# convert mean monthly precipitation per second to total
precip_to_total <- function (x) {
  month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  month_seconds <- month_days * 3600 * 24
  dims <- dim(x)
  total <- x * array(rep(month_seconds, each = prod(dims[1:2])), dims)
  total <- st_set_dimensions(total, "time", tolower(month.abb))
  total
}

prec <- precip_to_total(prec1991_2011)

plot(prec, key.pos = 4, key.width = lcm(1), box_col = "white", col = hcl.colors(10))
plot(prec)

stars::write_stars(prec, paste0(dirScratch,"chess_pr_1991_2011_monthly.tif"))
