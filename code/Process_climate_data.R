
# date: 16-04-24
# author: VB
# desc: script to process climate data (CHESS baseline 1991-2011 from CEH), convert nc files to tif,
# calculate CMD for ESC

### working dirs ----

#wd <- "~/Documents/Woodland-Trust/Data-Analysis/SoWT2_Spp-uncertainty" # Mac path
wd <- "C:/Users/vbu/OneDrive - the Woodland Trust/Data-analysis/SoWT2_Spp-uncertainty" # WT path
dirData <-  paste0(wd,"/data-raw/")
dirScratch <-  paste0(wd,"/data-scratch/")
dirOut <- paste0(wd,"/data-out/")

### load libs ----

library(dplyr)
library(stars)
library(raster)
library(terra)
library(ncmeta)

# load in functions from #CEH script
source(file.path(wd,"code","NC-to-tiff.R"))

### loop to process all .nc files, per baseline/scenario ----

# list folders within data-raw
lstFolders <- c("chess_baseline","speed_future_rcp26","speed_future_rcp45","speed_future_rcp85")

for (folder in lstFolders){
  
  folder <- lstFolders[1]
  
  lstFiles <- list.files(paste0(dirData,"/",folder))
  
  # create metadata
  description <- c(
    gdd = "growing degree days and number of days",
    precip = "precipitation",
    pet = "Penman-Monteith potential evapotranspiration, for a well-watered grass surface")
  
  metadata <- tibble::tibble(file = lstFiles) %>%
    tidyr::extract(
    file,
    c("variable", "summary", "temp_resolution", "from", "to"),
    "chess-[[:alpha:]]+_(.+)_[[:alpha:]]{2}_1km_20yr-([[:alpha:]]+)-([[:alpha:]]+)_([[:digit:]]+)-([[:digit:]]+)",
    remove = FALSE) %>%
    mutate(from_year = substr(from, 1, 4),
           from_month = substr(from, 5, 6),
           to_year = substr(to, 1, 4),
           to_month = substr(to, 5, 6),
           description = description[variable])
  
  write.csv(metadata, file.path(dirData,"/",folder, "metadata.csv"), row.names = FALSE)
  
  # process nc files
  
  # gdd (growing degree days, sometimes called accumulated temperature AT)
  gdd1991_2011 <- stars::read_ncdf(file.path(dirData,folder,"/",metadata$file[4]), 27700)
  plot(gdd1991_2011)
  stars::write_stars(gdd1991_2011, dsn = paste0(dirScratch,"/chess_gdd_1991_2011_annual.tif"))
  
  # precipitation
  prec1991_2011 <- stars::read_ncdf(file.path(dirData,folder,"/",metadata$file[8]), 27700)
  plot(prec1991_2011)
  
  # calculate total monthly precipitation
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
  
  dev.off()
  plot(prec, key.pos = 4, key.width = lcm(1), box_col = "white", col = hcl.colors(10))

  stars::write_stars(prec, paste0(dirScratch,"chess_pr_1991_2011_monthly.tif"))
  
  # potential evapotranspiration
  pet1991_2011 <- stars::read_ncdf(file.path(dirData,folder,"/",metadata$file[12]), path_baseline, 27700)
  plot(pet1991_2011)
  
  # calculate total monthly potential evapotranspiration
  # the data provided is mean monthly potential evapotranspiration (per day)
  # it needs to be transformed to total - multiply by the number of days in a month
  
  # convert mean monthly potential evapotranspiration to total
  pet_to_total <- function (x) {
    month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    dims <- dim(x)
    total <- x * array(rep(month_days, each = prod(dims[1:2])), dims)
    total <- st_set_dimensions(total, "time", tolower(month.abb))
    total
  }
  
  pet <- pet_to_total(pet1991_2011)
  
  dev.off()
  plot(pet, key.pos = 4, key.width = lcm(1), box_col = "white", col = hcl.colors(10))
  
  stars::write_stars(pet, dsn = paste0(dirScratch,"chess_pet_1991_2011_monthly.tif"))
  
  # now calculate CMD = climatic moisture deficit, required for ESC
  # use monthly values for pet and prec
  
  # calculate monthly moisture deficit (mMD)
  # mMD = pet - precip
  # positive values = deficit, negative values = surplus
  mMD <- pet - prec
  
  dev.off()
  plot(mMD, key.pos = 4, key.width = lcm(1), box_col = "white", col = hcl.colors(10))
  
  stars::write_stars(mMD, dsn = paste0(dirScratch,"chess_mMD_1991_2011_monthly.tif"), NA_value = NA)
  
  # issue here. can't read in tiff as brick, which then affects function below
  #EtoPrDiff <- brick(paste0(dirScratch,"chess_mMD_1991_2011_monthly.tif"))
  
  # work around - subset each month from stars object and convert to raster
  jan <- raster(mMD[[1]][,,1], crs = 27700)
  feb <- raster(mMD[[1]][,,2], crs = 27700)
  mar <- raster(mMD[[1]][,,3], crs = 27700)
  apr <- raster(mMD[[1]][,,4], crs = 27700)
  may <- raster(mMD[[1]][,,5], crs = 27700)
  jun <- raster(mMD[[1]][,,6], crs = 27700)
  jul <- raster(mMD[[1]][,,7], crs = 27700)
  aug <- raster(mMD[[1]][,,8], crs = 27700)
  sep <- raster(mMD[[1]][,,9], crs = 27700)
  oct <- raster(mMD[[1]][,,10], crs = 27700)
  nov <- raster(mMD[[1]][,,11], crs = 27700)
  dec <- raster(mMD[[1]][,,12], crs = 27700)
  
  EtoPrDiff <- brick(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,dec)
  
  calcCMD <- function(EtoPrDiff) {
    EtoPrDiff <- EtoPrDiff[!is.na(EtoPrDiff)]
    if (length(EtoPrDiff) > 0) {
      CMD <- 0
      accumulator <- 0
      for (a in EtoPrDiff) {
        if (accumulator >= 0) {
          accumulator <- max(accumulator + a, 0)
          if (accumulator > CMD) CMD <- accumulator
        } else {
          accumulator <- 0
        }
      }
    } else {
      CMD <- NA
    }
    
    return(CMD)
  }
  
  rasterCMD <- raster::calc(EtoPrDiff, fun = calcCMD)
  #rasterCMD <- terra::app(EtoPrDiff, fun = calcCMD)
  
  dev.off()
  plot(rasterCMD, col = hcl.colors(10))
  
  # apply adjustment
  diff <- 0.0011*rasterCMD^2 - 0.076*rasterCMD + 0.08465
  
  rasterCMD_adj <- rasterCMD - diff
  
  dev.off()
  plot(rasterCMD_adj, col = hcl.colors(10))
  
  writeRaster(rasterCMD_adj, paste0(dirScratch,"/chess_CMD_adj_1991_2011_baseline.tif") , overwrite=TRUE)
  
}

